subroutine writed(n)
    use params
    use globals
    use interface_mod, only : bqtoq

    implicit none
    integer, intent(in) :: n
    integer i, j, index
    integer, parameter :: sep = 50
    real(8) x_write, y_write, temp_q(4)
    character filename*128

    open(10, file = 'data_q.csv')
    j = int(jmax / 2)
    do i = imin, imax
        !q(cell center)を表示
        temp_q(:) = bqtoq(bq(i, j, :))
        !glid(cell vertex)をq(i,j)の位置(cell center)にするため、cell center周りの4隅を平均する
        x_write = (x(i-1, j-1) + x(i, j-1) + x(i, j) + x(i-1, j)) / 4.0d0
        y_write = (y(i-1, j-1) + y(i, j-1) + y(i, j) + y(i-1, j)) / 4.0d0
        write(10, *) x_write, ',', temp_q(1), ',', temp_q(2), ',', temp_q(3), ',', temp_q(4)
        !& , ',', ',', e(i, j, 1), ',', e(i, j, 2), ',', e(i, j, 3), ',', e(i, j, 4) &
        !& , ',', ',', f(j, i, 1), ',', f(j, i, 2), ',', f(j, i, 3), ',', f(j, i, 4) &
        !& , ',', ',', rec(i, 1), ',', rec(i, 2), ',', rec(i, 3), ',', rec(i, 4), ',', rec(i, 5)
    enddo
    close(10)
    ! x4   x3     q(i,j) = cell_center
    !    q        x1(i-1, j-1), x2(i, j-1)
    ! x1   x2     x3(i, j), x4(i-1, j)

    !data_res.csvフォーマット、書き込み
    if ( n == 1 ) then
        open(11, file = 'data_res.csv')
        write(11, '(a9, a1, a9, a1, a9, a1, a9, a1, a9, a1, a9)') &
        & "time", ',', "n_time =", ',', "res_rho", ',', "res_u", ',', "res_v", ',', "res_p"
        close(11)
    end if
    open(11, file = 'data_res.csv', position = 'append')
    write(11, *) time, ',', n, ',', res(1), ',', res(2), ',', res(3), ',', res(4)
    close(11)


    !!以下make_vtkファイル用
    !imin, imax, jmin, jmax書き出し
    open(12, file = 'GridNum.txt')
    write(12, *) imin, ',', imax, ',', jmin, ',', jmax
    close(12)

    index = int((n-1) / sep) + 1 !n=1~100でindex=1, n=101~200でindex=2
    write(filename,'("Qascii_",i5.5,".dat")') index !index(integer)をfilename(char)に代入し、文字+整数を文字に変換
    open(120+index,file = filename)
    !Qの内容を読み込みます
    rewind(120+index)
    write(120+index,*) 'meshfile.txt'
    do j= jmin-2, jmax+2
        do i= imin-2, imax+2
            write(120+index,*) bq(i,j,1), bq(i,j,2), bq(i,j,3), bq(i,j,4)
        enddo
    enddo
    close(120+index)

    !ntime, time書き出し
    !time.txtをフォーマット
    if ( n==1 ) then
        open(50, file = 'time.txt')
        write(50, *) 'filename,     n,      time'
        write(50, *)
        write(50, *)
        close(50)
    end if

    open(50, file = 'time.txt', position = 'append')
    if (index /= int((n) / sep) + 1) then !index_n /= index_n+1 (cf.n=100,200...)
        backspace(50) !前のループで書かれたnに空白を上書き
        backspace(50) !前のループで書かれたfilenameとかに空白を上書き
        write(50, '(a15, a1, i5, a1, f15.10)') filename, ',', n, ',', time
        write(50, *) n, ',', sep !make_vtkでnmaxを読み込むため
        write(50, *) !filenameとかが消されないよう空白行を追加
    else
        backspace(50)
        backspace(50)
        write(50, '(a15, a1, i5, a1, f15.10)') filename, ',', n, ',', time
        write(50, *) n
    endif
    close(50)

    !n=2,4,8を見るため改変
    !if ( n==100 ) then
        !index = 1
    !else if (n==130) then
        !index = 2
    !else if (n==200) then
        !index = 3
    !endif
    !if (n==100 .or. n==130 .or. n==200) then
        !write(filename,'("Qascii_",i1.1,".dat")') index !index(integer)をfilename(char)に代入し、文字+整数を文字に変換
        !open(20+index,file = filename)
        !!Qの内容を読み込みます
        !rewind(20+index)
        !write(20+index,*) 'meshfile.txt'
        !do j= jmin-2, jmax+2
            !do i= imin-2, imax+2
                !write(20+index,*) bq(i,j,1), bq(i,j,2), bq(i,j,3), bq(i,j,4)
            !enddo
        !enddo
        !close(20+index)
    !endif
    !

end subroutine writed
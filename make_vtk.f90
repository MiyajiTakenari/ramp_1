program visualize
    implicit none
    integer :: fo=20, i=0, j=0, GridNum=18, Grid=19, Qfi=100, Timefi=30, &
                & index=0, max_index, imin, imax, jmin, jmax, ntime, sep, iost
    character filename*128, meshfile*64
    real(8),allocatable :: x(:, :), y(:, :), bq(:, :, :)
    real(8) :: xcenter, ycenter, rho, u, v, p, t
    real(8),parameter :: GAMMA = 1.4d0


    !マニュアル著者の環境では，格子数を"100 100"のように記載したファイルを用意して，
    !格子数を前処理，流れ場の計算，後処理，のどのプログラムからも読み出せるようにしています．
    !そのため，格子数を読みだすコードを書いています．!!格子数ファイルがGridNum.txt
    open(GridNum,file = 'GridNum.txt')
    read(GridNum,*) imin, imax, jmin, jmax
    close(GridNum)
    write(*,*) 'imin=', imin, 'imax=', imax, 'jmin=', jmin, 'jmax=', jmax

    !ntime, time読み込み, templeteに追加
    open(Timefi,file = 'time.txt')
    do !ファイル終了条件(iost<0)が検出されたらexit、それまでntimeに行を読み込み
        read(Timefi, *, iostat = iost) ntime, sep
        if ( iost < 0 ) exit
    end do
    close(Timefi)
    max_index = int((ntime-1) / sep) + 1 !n=1~100でindex=1, n=101~200でindex=2
    write(*,*) 'ntime=', ntime, 'max_index=', max_index, 'n_seprate=', sep

    allocate (x(imin-3:imax+2, jmin-3:jmax+2), y(imin-3:imax+2, jmin-3:jmax+2))
    allocate (bq(imin-2:imax+2, jmin-2:jmax+2, 4))


    !メッシュの情報を読み込みます．
    !HOGEHOGEをメッシュの座標を記載しているファイル名に置き換えて使用してください
    !※マニュアル著者は，メッシュの座標を記載したファイルと計算結果(基本量Qの内容)が対応し，後処理の際に
    !読み込むメッシュの情報に間違いがないように，Qの結果を出力するファイルの1行目に使用したメッシュファイル名を出力しています
    !
    !# Qを出力するファイルの内容
    !   meshfile.txt (←メッシュの座標を記載したファイル名)
    !   1.0000         2.0000        3.0000       4.0000 (←以下がQの内容)
    !   1.0000         2.0000        3.0000       4.0000
    !   1.0000         2.0000        3.0000       4.0000
    !                           ・
    !                           ・
    !                           ・
    !これに従って読み込むので，40行目はQを出力したファイルの1行目，すなわちメッシュの座標を記載したファイル名を読み込む，ということになります
    !以上の内容はマニュアル著者が利便性のために行った細工なので，ご自身の環境に合わせてファイル読み込みを行ってください
    write(filename,'("Qascii_",i5.5,".dat")') index !index(integer)をfilename(char)に代入し、文字+整数を文字に変換
    open(Qfi,file = filename)
    read(Qfi,*) meshfile
    write(*,*) meshfile
    close(Qfi) !templeteに追加

    !読み込んだファイル名を用いてファイルオープン，座標を読み込んでいきます
    open(Grid,file = meshfile)
    do j= jmin-3, jmax+2
        do i= imin-3, imax+2
            read(Grid,*) x(i,j), y(i,j)
        enddo
    enddo
    close(Grid)

    !Qの内容を読み込みます.繰り返し部分追加
    do index = 0, max_index !!Qasciiの個数で変える
        write(filename,'("Qascii_",i5.5,".dat")') index !index(integer)をfilename(char)に代入し、文字+整数を文字に変換
        open(Qfi,file = filename)
        rewind(Qfi) !templeteに追加
        read(Qfi,*) meshfile !templeteに追加, meshfile.txtの名前分ずらす
        do j= jmin-2, jmax+2
            do i= imin-2, imax+2
                read(Qfi,*) bq(i,j,1), bq(i,j,2), bq(i,j,3), bq(i,j,4)
            enddo
        enddo
        close(Qfi)
        Qfi = Qfi + 1 !追加

        !paraviewに読み込ませるVTKファイルの名前です．適宜変更してください
        write(filename,'("ramp_",i5.5,".vtk")') index
        open(fo,file = filename)

        !VTKフォーマットに従い，69行目までは指定されています．ただし66行目のFUGAは任意のファイル名なので適宜変更してください
        write(fo,"('# vtk DataFile Version 3.0')")
        write(fo,"('CELLCENTER')") !FUGA
        write(fo,"('ASCII')")
        write(fo,"('DATASET STRUCTURED_GRID')") !非構造格子を用いる場合は 'DATASET UNSTRUCTURED GRID'　となります
        write(fo,"('DIMENSIONS',3(1x,i3))") imax-imin+1, jmax-jmin+1,1 !マニュアル著者は2次元の計算を行っているのでz座標の格子数は1です
        write(fo,"('POINTS ',i9,' float')") (imax-imin+1) * (jmax-jmin+1)    !3次元計算を行う方は Nx*Ny*Nz となります
        do j= jmin, jmax
            do i= imin, imax
                xcenter  = 0.25d0 *(x(i,j) + x(i-1,j) + x(i,j-1) + x(i-1,j-1))
                ycenter  = 0.25d0 *(y(i,j) + y(i-1,j) + y(i,j-1) + y(i-1,j-1))
                write(fo,'(3(f9.4,1x))') xcenter, ycenter, 0.0d0
            enddo
        enddo


        !81~83行目まではvtkフォーマット既定のものです．そのあとは出力した座標に対応した値を出力してください
        write(fo,"('POINT_DATA',i9)") (imax-imin+1) * (jmax-jmin+1)
        write(fo,"('SCALARS rho float')")
        write(fo,"('LOOKUP_TABLE default')")
        do j= jmin, jmax
            do i= imin, imax
                write(fo,"(f10.7)") bq(i,j,1)
            enddo
        enddo

        !速度のようなベクトル場に対しては　VECTORS　という形式を指定することもできるようです
        write(fo,"('SCALARS U float')")
        write(fo,"('LOOKUP_TABLE default')")
        do j= jmin, jmax
            do i= imin, imax
                rho = bq(i,j,1)
                u   = bq(i,j,2) / bq(i,j,1)
                v   = bq(i,j,3) / bq(i,j,1)
                p   = (GAMMA-1.0d0) * (bq(i,j,4) - 0.5d0 * rho * (u**2.0d0 + v**2.0d0))
                t   = p / rho
                write(fo,"(f10.7)") u
            enddo
        enddo

        write(fo,"('SCALARS V float')")
        write(fo,"('LOOKUP_TABLE default')")
        do j= jmin, jmax
            do i= imin, imax
                rho = bq(i,j,1)
                u   = bq(i,j,2) / bq(i,j,1)
                v   = bq(i,j,3) / bq(i,j,1)
                p   = (GAMMA-1.0d0) * (bq(i,j,4) - 0.5d0 * rho * (u**2.0d0 + v**2.0d0))
                t   = p / rho
                write(fo,"(f10.7)") v
            enddo
        enddo

        write(fo,"('SCALARS p float')")
        write(fo,"('LOOKUP_TABLE default')")
        do j= jmin, jmax
            do i= imin, imax
                rho = bq(i,j,1)
                u   = bq(i,j,2) / bq(i,j,1)
                v   = bq(i,j,3) / bq(i,j,1)
                p   = (GAMMA-1.0d0) * (bq(i,j,4) - 0.5d0 * rho * (u**2.0d0 + v**2.0d0))
                t   = p / rho
                write(fo,"(f10.7)")p
            enddo
        enddo
        close(fo)
        fo = fo + 1 !追加
    end do


    deallocate(x,y,bq)

end program visualize
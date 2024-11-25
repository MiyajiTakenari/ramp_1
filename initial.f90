subroutine init
    use params
    use globals
    use interface_mod, only : qtobq

    implicit none
    integer i, j

    !初期条件Q
    !call glid
    do i = imin-2, imax+2
        do j = jmin-2, jmax+2
                !-0.5d0 + dx * x(i, j) <= 0.0d0
                bq(i, j, 1:4) = qtobq(1.0d0, 2.0d0, 0.0d0, 1.0d0/gamma)
        enddo
    enddo

    open(20,file = 'Qascii_00000.dat')
    !Qの内容を読み込みます
    rewind(20)
    write(20, *) 'meshfile.txt'
    do j= jmin-2, jmax+2
        do i= imin-2, imax+2
            write(20, *) bq(i,j,1), bq(i,j,2), bq(i,j,3), bq(i,j,4)
        enddo
    enddo
    close(20)

end subroutine init
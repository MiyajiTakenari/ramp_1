subroutine alloc(al_flag)
    use params
    use globals

    implicit none
    integer, intent(in) :: al_flag
    if (al_flag == 1) then
        !表示するqつまりQはmin~max
        allocate (x(imin-3:imax+2, jmin-3:jmax+2), y(imin-3:imax+2, jmin-3:jmax+2))
        allocate (bq(imin-2:imax+2, jmin-2:jmax+2, 4))
        allocate (e(imin-1:imax, jmin-1:jmax, 4))
        allocate (f(imin-1:imax, jmin-1:jmax, 4))
        allocate (mx(imin-3:imax+2, jmin-2:jmax+2), my(imin-3:imax+2, jmin-2:jmax+2))
        allocate (nx(imin-2:imax+2, jmin-3:jmax+2), ny(imin-2:imax+2, jmin-3:jmax+2))
        allocate (s_j(imin:imax, jmin:jmax))
        allocate (bq_n(imin:imax, jmin:jmax, 4))
        !B.Cのfuncion ave_n, ave_mでmetricsの配列の下限値使ってる。下限値変更時注意
    else if (al_flag == 0) then
        deallocate (x, y, bq, e, f, mx, my, nx, ny, s_j, bq_n)
    end if

end subroutine alloc
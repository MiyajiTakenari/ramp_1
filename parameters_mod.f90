module params
    implicit none
    integer, parameter :: imin=1, imax=100, jmin=1, jmax=60, nmax=3000
    real(8), parameter :: ex_time = 0.2d0, &
    & cfl = 0.15d0, gamma = 1.4d0, beta = 2.5d0, phi = 1.0d0 / 3.0d0, &
    & pi = acos(-1.0d0), theta = (5.0d0 / 180.0d0) * pi, res_min = 10.0d0 ** (-4.0d0)
    !dx,dy注意：range_x/(Qの表示数)
    !imin,jminを変更した場合,interface_modのave_m, ave_nの仮引数の下限と
    !boundaryのave_m, ave_nの仮引数の下限を変更すること
    !range_yは坂道無しの大きい方
    !range_x = 1.0d0, & !dx = range_x / dble(imax-imin+1), &
    !& range_y = 0.6d0, & !dy = range_y / dble(jmax-jmin+1), &
end module params
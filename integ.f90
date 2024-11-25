subroutine integ
    use params
    use globals
    use interface_mod, only : xflux, yflux, bound

    implicit none
    integer i, j

    !2段階2次精度スキーム
    bq_n(imin:imax, jmin:jmax, 1:4) = bq(imin:imax, jmin:jmax, 1:4)
    call xflux
    call yflux
    do i = imin, imax
        do j = jmin, jmax
            bq(i, j, 1:4) = bq(i, j, 1:4) - (1.0d0 / s_j(i, j)) * (dt / 2.0d0) &
            & * ((e(i, j, 1:4) - e(i-1, j, 1:4)) + (f(i, j, 1:4) - f(i, j-1, 1:4)))
        end do
    end do

    call bound
    call xflux
    call yflux
    do i = imin, imax
        do j = jmin, jmax
            bq(i, j, 1:4) = bq_n(i, j, 1:4) - (1.0d0 / s_j(i, j)) * dt &
            & * ((e(i, j, 1:4) - e(i-1, j, 1:4)) + (f(i, j, 1:4) - f(i, j-1, 1:4)))
        end do
    end do

end subroutine integ
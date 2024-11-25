subroutine calc_res
    use params
    use globals

    implicit none
    integer i, j

    !res_xはx方向の残差,それがy本(jmin:jmax)ある
    !res_yはy方向の残差,それがx本(imin:imax)ある
    res(:) = 0.0d0
    do j = jmin, jmax
        do i = imin, imax
            res(:) = res(:) + abs(bq(i, j, :) - bq_n(i, j, :)) &
            & / (dble((imax - imin + 1)*(jmax - jmin + 1)) * dt)
        end do
    end do

end subroutine calc_res
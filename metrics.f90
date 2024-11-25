subroutine metrics
    use globals
    use params

    implicit none
    integer i, j

    do i = imin-3, imax+2
        do j = jmin-2, jmax+2
            mx(i, j) = (y(i, j) - y(i, j-1))
            my(i, j) = -(x(i, j) - x(i, j-1))
        end do
    end do

    do j = jmin-3, jmax+2
        do i = imin-2, imax+2
            nx(i, j) = -(y(i, j) - y(i-1, j))
            ny(i, j) = (x(i, j) - x(i-1, j))
        end do
    end do

    do i = imin, imax
        do j = jmin, jmax
            s_j(i, j) = ny(i, j) * mx(i, j) - nx(i, j) * my(i, j)
        end do
    end do
end subroutine metrics
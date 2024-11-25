subroutine glid
    use params
    use globals

    implicit none
    integer i, j
    real(8) dx
    real(8), parameter :: range_x = 1.0d0, range_y = 0.6d0, ramp_x = 0.4d0

    dx = range_x / dble(imax-imin+1)

    open(13,file = 'meshfile.txt')
    do j = jmin-3, jmax+2
        do i = imin-3, imax+2
            !物理空間初期分布
            x(i, j) = dx * dble(i)
            if ( x(i, j)<= ramp_x ) then
                y(i, j) = (range_y / dble(jmax-jmin+1)) * dble(j)
            else
                y(i, j) = (range_y / dble(jmax-jmin+1)) * dble(j) &
                & + tan(theta) * (1.0d0 - dble(j) / dble(jmax-jmin+1)) * (x(i, j) - ramp_x)
            endif
            write(13,*) x(i,j), ',',  y(i,j)
        end do
    end do
    close(13)

end subroutine glid
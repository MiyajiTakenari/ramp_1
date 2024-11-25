subroutine glid
    use params
    use globals

    implicit none
    integer i, j
    real(8) dx, range_y, h
    real(8), parameter :: range_x = 1.5d0, ramp_x = 0.2d0, ramp_leng = 1.0d0, swa = (26.619d0 / 180.0d0) * pi

    dx = range_x / dble(imax-imin+1)
    h = ((tan(swa) - tan(theta))) ** 2.0d0 / (2.0d0 * tan(swa) + tan(theta) * (tan(swa) ** 2.0d0 - 1.0d0))
    range_y = h + ramp_leng * tan(theta)

    open(13,file = 'meshfile.txt')
    do j = jmin-3, jmax+2
        do i = imin-3, imax+2
            !物理空間初期分布
            x(i, j) = dx * dble(i)
            if ( x(i, j)<= ramp_x ) then
                y(i, j) = (range_y / dble(jmax-jmin+1)) * dble(j)
            elseif( x(i, j) >= ramp_x + ramp_leng ) then
                y(i, j) = (h / dble(jmax-jmin+1)) * dble(j) + ramp_leng * tan(theta)
            else
                y(i, j) = (range_y / dble(jmax-jmin+1)) * dble(j) &
                & + tan(theta) * (1.0d0 - dble(j) / dble(jmax-jmin+1)) * (x(i, j) - ramp_x)
            endif
            write(13,*) x(i,j), ',',  y(i,j)
        end do
    end do
    close(13)

end subroutine glid
function minmod(a, b) result(z)
    implicit none
    real(8), intent(in) :: a(4), b(4)
    real(8) z(4)
    integer i
    do i = 1, 4
        z(i) = sign(1.0d0, a(i)) * max(0.0d0, min(abs(a(i)), b(i) * sign(1.0d0, a(i))))
    end do
end function minmod
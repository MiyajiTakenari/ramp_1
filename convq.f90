function qtobq(rho, u, v, p) result(big_q)
    use params
    use globals

    implicit none
    real(8), intent(in) :: rho, u, v, p
    real(8) big_q(4)
    big_q(1) = rho
    big_q(2) = rho * u
    big_q(3) = rho * v
    big_q(4) = (p/(gamma-1.0d0)) + 0.5d0 * (rho * u * u + rho * v * v)
end function qtobq

function bqtoq(big_q) result(q)
    use params
    use globals

    implicit none
    real(8), intent(in) :: big_q(4)
    real(8) q(4)
    q(1) = big_q(1)
    q(2) = big_q(2)/big_q(1)
    q(3) = big_q(3)/big_q(1)
    q(4) = (gamma-1.0d0) * (big_q(4) - 0.5d0 * big_q(1) * ((big_q(2)/big_q(1)) ** 2.0d0 + (big_q(3)/big_q(1)) ** 2.0d0))
end function bqtoq
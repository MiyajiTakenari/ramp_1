module interface_mod
    interface

        function qtobq(rho, u, v, p) result(big_q)
            real(8), intent(in) :: rho, u, v, p
            real(8) big_q(4)
        end function qtobq

        function bqtoq(big_q) result(q)
            real(8), intent(in) :: big_q(4)
            real(8) q(4)
        end function bqtoq

        subroutine alloc(al_flag)
            integer, intent(in) :: al_flag
        end subroutine alloc

        subroutine glid
        end subroutine glid

        subroutine init
        end subroutine init

        subroutine metrics
        endsubroutine metrics

        function ave_m(met, i, j) result(ave)
            use params
            integer, intent(in) :: i, j
            real(8), intent(in) :: met(imin-3:, jmin-2:)
            real(8) ave
        end function ave_m

        function ave_n(met, i, j) result(ave)
            use params
            integer, intent(in) :: i, j
            real(8), intent(in) :: met(imin-2:, jmin-3:)
            real(8) ave
        end function ave_n

        subroutine bound
        end subroutine bound

        subroutine cflc
        end subroutine cflc

        function minmod(a, b) result(z)
            real(8), intent(in) :: a(4), b(4)
            real(8) z(4)
        end function minmod

        subroutine delbar(qim, qi, qip, bar_plus, bar_minus)
            real(8), intent(in) :: qim(4), qi(4), qip(4)
            real(8), intent(out) :: bar_plus(4), bar_minus(4)
        end subroutine delbar

        subroutine xflux
        end subroutine xflux

        subroutine yflux
        end subroutine yflux

        subroutine integ
        end subroutine integ

        subroutine calc_res
        end subroutine calc_res

        subroutine writed(n)
            integer, intent(in) :: n
        end subroutine writed

    end interface
end module interface_mod
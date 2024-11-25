subroutine delbar(qim, qi, qip, bar_plus, bar_minus)
    use params
    use interface_mod, only : minmod

    implicit none
    real(8), intent(in) :: qim(4), qi(4), qip(4)
    real(8), intent(out) :: bar_plus(4), bar_minus(4)
    real(8) plus(4), minus(4)
    plus(:) = qip(:) - qi(:)
    minus(:) = qi(:) - qim(:)
    bar_plus(:) = minmod(plus(:), beta * minus(:))
    bar_minus(:) = minmod(minus(:), beta * plus(:))
end subroutine delbar
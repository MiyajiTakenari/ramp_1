program euler_2D_B
    use params
    use globals
    use interface_mod

    implicit none
    integer :: n = 1

    !計算開始
    call alloc(1)
    call glid
    call init
    call metrics
    !時間進める(n>nmaxかres(1:4)<res_minでexit)
    do while (n <= nmax .and. (res(1) >= res_min .or. res(2) >= res_min &
        & .or. res(3) >= res_min .or. res(4) >= res_min))
        call bound
        call cflc
        call integ
        call calc_res
        call writed(n)
        !write (*, *) 'dt = ', dt
        n = n + 1
    end do
    write(*,'(a9, i8)') 'ntime = ', n - 1

    call alloc(0)

end program euler_2D_B
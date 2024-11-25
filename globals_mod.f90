module globals
    real(8), save :: dt, time = 0.0d0, res(1:4) = 1.0d0!res_init>res_exitであればよい
    !real(8), save :: rec(-1:100, 50)
    !(x, y, 成分)
    real(8), allocatable, save :: x(:, :), y(:, :), bq(:, :, :), e(:, :, :), f(:, :, :), &
    & mx(:, :), my(:, :), nx(:, :), ny(:, :), s_j(:, :), &
    & bq_n(:, :, :)
    !integer, save :: exit_flag = 0
end module globals
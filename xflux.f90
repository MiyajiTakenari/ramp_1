subroutine xflux
    use params
    use globals
    use interface_mod, only : bqtoq, delbar

    implicit none
    integer i, j
    real(8) temp_q(4), bar_m(4), bar_p(4), & !MUSCLで使う
    & rho_l, su_l, sv_l, p_l, & !MUSCLの答え_left
    & rho_r, su_r, sv_r, p_r, & !MUSCLの答え_right
    & u_l, u_r, v_l, v_r, mx_h, my_h, & !large_u, large_v
    & a_l, a_r, cm, ul_p, ur_m, pl_p, pr_m, & !2,4,3で使う
    & rus_r, rus_d, mdot, s, & !5で使う
    & e_l, e_r, h_l, h_r, eb(4) !Flux評価で使う

    !位置進めるE_bar(-1:imax)
    !fluxの計算
    do j = jmin-1, jmax
        do i = imin-1, imax
            !MUSCL approach
            !qlを計算
            call delbar(bqtoq(bq(i-1,j,:)), bqtoq(bq(i,j,:)), bqtoq(bq(i+1,j,:)), bar_p, bar_m)
            temp_q(:) = bqtoq(bq(i,j,:)) + 0.25d0 * (1.0d0 - phi) * bar_m(:) + 0.25d0 * (1.0d0 + phi) * bar_p(:)
            rho_l = temp_q(1)
            su_l = temp_q(2)
            sv_l = temp_q(3)
            p_l = temp_q(4)
            !qrを計算
            call delbar(bqtoq(bq(i,j,:)), bqtoq(bq(i+1,j,:)), bqtoq(bq(i+2,j,:)), bar_p, bar_m)
            temp_q(:) = bqtoq(bq(i+1,j,:)) - 0.25d0 * (1.0d0 - phi) * bar_p(:) - 0.25d0 * (1.0d0 + phi) * bar_m(:)
            rho_r= temp_q(1)
            su_r = temp_q(2)
            sv_r = temp_q(3)
            p_r = temp_q(4)

            !large_u, laege_vのl,r定義
            mx_h = mx(i, j) / sqrt(mx(i, j) ** 2.0d0 + my(i, j) ** 2.0d0)
            my_h = my(i, j) / sqrt(mx(i, j) ** 2.0d0 + my(i, j) ** 2.0d0)
            u_l = mx_h * su_l + my_h * sv_l
            u_r = mx_h * su_r + my_h * sv_r
            v_l = -my_h * su_l + mx_h * sv_l
            v_r = -my_h * su_r + mx_h * sv_r

            !2,4,3を計算
            a_l = (2.0d0 * (p_l / rho_l)) / ((p_l/rho_l) + (p_r/rho_r))
            a_r = (2.0d0 * (p_r / rho_r)) / ((p_l/rho_l) + (p_r/rho_r))
            cm = max(sqrt(gamma * p_l / rho_l), sqrt(gamma * p_r / rho_r))

            !!rec
            !rec(i, 1) = su_l
            !rec(i, 2) = sv_l
            !rec(i, 3) = mx_h
            !rec(i, 4) = my_h
            !rec(i, 5) = u_l
            !rec(i, 6) = u_r
            !rec(i, 7) = v_l
            !rec(i, 8) = v_r
            !rec(i, 9) = a_l
            !rec(i, 10) = a_r
            !rec(i, 11) = cm
            !rec(i, 12) = rho_l
            !rec(i, 13) = p_l

            if (abs(u_l) <= cm) then
                ul_p = ((a_l * (u_l + cm) ** 2.0d0) / (4.0d0 * cm)) + ((1.0d0 - a_l) * (u_l + abs(u_l)) / 2.0d0)
                pl_p = (p_l * (2.0d0 - (u_l/cm)) * ((u_l/cm) + 1.0d0) ** 2.0d0) / 4.0d0
            else
                ul_p = (u_l + abs(u_l)) / 2.0d0
                pl_p = p_l * (u_l + abs(u_l)) / (2.0d0 * u_l)
            end if

            if (abs(u_r) <= cm) then
                ur_m = ((-a_r * (u_r - cm) ** 2.0d0) / (4.0d0 * cm)) + ((1.0d0 - a_r) * (u_r - abs(u_r)) / 2.0d0)
                pr_m = (p_r * (2.0d0 + (u_r/cm)) * ((u_r/cm) - 1.0d0) ** 2.0d0) / 4.0d0
            else
                ur_m = (u_r - abs(u_r)) / 2.0d0
                pr_m = p_r * (u_r - abs(u_r)) / (2.0d0 * u_r)
            end if

            !5,1を計算
            rus_r = ul_p * (rho_l * u_l) + ur_m * (rho_r * u_r)
            mdot = ul_p * rho_l + ur_m * rho_r
            rus_d = 0.5d0 * (mdot * (u_l + u_r) - abs(mdot) * (u_r - u_l))
            s = 0.5d0 * (1.0d0 + min(1.0d0, 10.0d0 * abs(p_r - p_l) / min(p_l, p_r)))

            !E_barを計算
            eb(1) = mdot
            eb(2) = s * rus_r + (1.0d0 - s) * rus_d + pl_p + pr_m
            eb(3) = 0.5d0 * (mdot * (v_l + v_r) - abs(mdot) * (v_r - v_l))
            e_l = (p_l / (gamma - 1.0d0)) + 0.5d0 * rho_l * (su_l ** 2.0d0 + sv_l ** 2.0d0)
            e_r = (p_r / (gamma - 1.0d0)) + 0.5d0 * rho_r * (su_r ** 2.0d0 + sv_r ** 2.0d0)
            h_l = (e_l + p_l) / rho_l
            h_r = (e_r + p_r) / rho_r
            eb(4) = 0.5d0 * (mdot * (h_l + h_r) - abs(mdot) * (h_r - h_l))

            !e_tilde = e を計算
            e(i, j, 1) = sqrt(mx(i, j) ** 2.0d0 + my(i, j) ** 2.0d0) * eb(1)
            e(i, j, 2) = sqrt(mx(i, j) ** 2.0d0 + my(i, j) ** 2.0d0) * (mx_h * eb(2) - my_h * eb(3))
            e(i, j, 3) = sqrt(mx(i, j) ** 2.0d0 + my(i, j) ** 2.0d0) * (my_h * eb(2) + mx_h * eb(3))
            e(i, j, 4) = sqrt(mx(i, j) ** 2.0d0 + my(i, j) ** 2.0d0) * eb(4)

        enddo
    enddo
end subroutine xflux
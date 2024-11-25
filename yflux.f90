subroutine yflux
    use params
    use globals
    use interface_mod, only : bqtoq, delbar

    implicit none
    integer i, j
    real(8) temp_q(4), bar_m(4), bar_p(4), & !MUSCLで使う
    & rho_l, su_l, sv_l, p_l, & !MUSCLの答え_left
    & rho_r, su_r, sv_r, p_r, & !MUSCLの答え_right
    & u_l, u_r, v_l, v_r, nx_h, ny_h, & !large_u, large_v
    & a_l, a_r, cm, vl_p, vr_m, pl_p, pr_m, & !2,4,3で使う
    & rvs_r, rvs_d, mdot, s, & !5で使う
    & e_l, e_r, h_l, h_r, fb(4) !Flux評価で使う

    !位置進めるF_bar(-1:jmax)
    !fluxの計算
    do i = imin-1, imax
        do j = jmin-1, jmax
            !MUSCL approach
            !qlを計算
            call delbar(bqtoq(bq(i,j-1,:)), bqtoq(bq(i,j,:)), bqtoq(bq(i,j+1,:)), bar_p, bar_m)
            temp_q(:) = bqtoq(bq(i,j,:)) + 0.25d0 * (1.0d0 - phi) * bar_m(:) + 0.25d0 * (1.0d0 + phi) * bar_p(:)
            rho_l = temp_q(1)
            su_l = temp_q(2)
            sv_l = temp_q(3)
            p_l = temp_q(4)
            !qrを計算
            call delbar(bqtoq(bq(i,j,:)), bqtoq(bq(i,j+1,:)), bqtoq(bq(i,j+2,:)), bar_p, bar_m)
            temp_q(:) = bqtoq(bq(i,j+1,:)) - 0.25d0 * (1.0d0 - phi) * bar_p(:) - 0.25d0 * (1.0d0 + phi) * bar_m(:)
            rho_r = temp_q(1)
            su_r = temp_q(2)
            sv_r = temp_q(3)
            p_r = temp_q(4)

            !large_u, laege_vのl,r定義
            nx_h = nx(i, j) / sqrt(nx(i, j) ** 2.0d0 + ny(i, j) ** 2.0d0)
            ny_h = ny(i, j) / sqrt(nx(i, j) ** 2.0d0 + ny(i, j) ** 2.0d0)
            u_l = ny_h * su_l - nx_h * sv_l
            u_r = ny_h * su_r - nx_h * sv_r
            v_l = nx_h * su_l + ny_h * sv_l
            v_r = nx_h * su_r + ny_h * sv_r

            !2,4,3を計算
            a_l = (2.0d0 * (p_l / rho_l)) / ((p_l/rho_l) + (p_r/rho_r))
            a_r = (2.0d0 * (p_r / rho_r)) / ((p_l/rho_l) + (p_r/rho_r))
            cm = max(sqrt(gamma * p_l / rho_l), sqrt(gamma * p_r / rho_r))

            !!rec
            !rec(j, 21) = su_l
            !rec(j, 22) = sv_l
            !rec(j, 23) = nx_h
            !rec(j, 24) = ny_h
            !rec(j, 25) = u_l
            !rec(j, 26) = u_r
            !rec(j, 27) = v_l
            !rec(j, 28) = v_r
            !rec(j, 29) = a_l
            !rec(j, 30) = a_r
            !rec(j, 31) = cm
            !rec(j, 32) = rho_l
            !rec(j, 33) = p_l

            if (abs(v_l) <= cm) then
                vl_p = ((a_l * (v_l + cm) ** 2.0d0) / (4.0d0 * cm)) + ((1.0d0 - a_l) * (v_l + abs(v_l)) / 2.0d0)
                pl_p = (p_l * (2.0d0 - (v_l/cm)) * ((v_l/cm) + 1.0d0) ** 2.0d0) / 4.0d0
            else
                vl_p = (v_l + abs(v_l)) / 2.0d0
                pl_p = p_l * (v_l + abs(v_l)) / (2.0d0 * v_l)
            endif

            if (abs(v_r) <= cm) then
                vr_m = ((-a_r * (v_r - cm) ** 2.0d0) / (4.0d0 * cm)) + ((1.0d0 - a_r) * (v_r - abs(v_r)) / 2.0d0)
                pr_m = (p_r * (2.0d0 + (v_r/cm)) * ((v_r/cm) - 1.0d0) ** 2.0d0) / 4.0d0
            else
                vr_m = (v_r - abs(v_r)) / 2.0d0
                pr_m = p_r * (v_r - abs(v_r)) / (2.0d0 * v_r)
            endif

            !5,1を計算
            rvs_r = vl_p * (rho_l * v_l) + vr_m * (rho_r * v_r)
            mdot = vl_p * rho_l + vr_m * rho_r
            rvs_d = 0.5d0 * (mdot * (v_l + v_r) - abs(mdot) * (v_r - v_l))
            s = 0.5d0 * (1.0d0 + min(1.0d0, 10.0d0 * abs(p_r - p_l) / min(p_l, p_r)))

            !F_barを計算
            fb(1) = mdot
            fb(2) = 0.5d0 * (mdot * (u_l + u_r) - abs(mdot) * (u_r - u_l))
            fb(3) = s * rvs_r + (1.0d0 - s) * rvs_d + pl_p + pr_m
            e_l = (p_l / (gamma - 1.0d0)) + 0.5d0 * rho_l * (su_l ** 2.0d0 + sv_l ** 2.0d0)
            e_r = (p_r / (gamma - 1.0d0)) + 0.5d0 * rho_r * (su_r ** 2.0d0 + sv_r ** 2.0d0)
            h_l = (e_l + p_l) / rho_l
            h_r = (e_r + p_r) / rho_r
            fb(4) = 0.5d0 * (mdot * (h_l + h_r) - abs(mdot) * (h_r - h_l))

            !f_tilde = f を計算
            f(i, j, 1) = sqrt(nx(i, j) ** 2.0d0 + ny(i, j) ** 2.0d0) * fb(1)
            f(i, j, 2) = sqrt(nx(i, j) ** 2.0d0 + ny(i, j) ** 2.0d0) * (ny_h * fb(2) + nx_h * fb(3))
            f(i, j, 3) = sqrt(nx(i, j) ** 2.0d0 + ny(i, j) ** 2.0d0) * (-nx_h * fb(2) + ny_h * fb(3))
            f(i, j, 4) = sqrt(nx(i, j) ** 2.0d0 + ny(i, j) ** 2.0d0) * fb(4)

        enddo
    enddo
end subroutine yflux
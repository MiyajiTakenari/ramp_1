function ave_m(met, i, j) result(ave)
    use params
    implicit none
    integer, intent(in) :: i, j
    real(8), intent(in) :: met(imin-3:, jmin-2:)
    real(8) ave
    ave = (met(i, j) + met(i-1, j)) / 2.0d0
end function ave_m

function ave_n(met, i, j) result(ave)
    use params
    implicit none
    integer, intent(in) :: i, j
    real(8), intent(in) :: met(imin-2:, jmin-3:)
    real(8) ave
    ave = (met(i, j) + met(i, j-1)) / 2.0d0
end function ave_n

subroutine bound
    use globals
    use params
    use interface_mod, only : qtobq, bqtoq, ave_m, ave_n

    implicit none
    real(8) temp_q(4), q_bc(4), bu, bv, u, v
    integer i, j
    !境界条件: qからq_bcを出す、その後 q_bc→bq
        !BD4
    !BD1    !BD2
        !BD3

    !BD3
    !slip condition
    !j=j+1, u_j = -u_j+1, rho_j = rho_j+1, e_j = e_j+1
    do i = imin-2, imax+2
        !rho_0 = rho_1, u_0 = -u_1, e_0 = e_1
        !temp_q(2) = u_1 = u, temp_q(3) = v_1 = v
        temp_q(:) = bqtoq(bq(i, jmin, :))
        q_bc(1) = temp_q(1)
        u = temp_q(2)
        v = temp_q(3)
        ! bu=U_1, bv=UU_1を求め、q_bc(2)=u_0, q_bc(3)=v_0を求める
        bu = ave_n(nx, i, jmin) * u + ave_n(ny, i, jmin) * v
        bv = -ave_n(ny, i, jmin) * u + ave_n(nx, i, jmin) * v
        q_bc(2) = (-ave_n(nx, i, jmin-1) * bu - ave_n(ny, i, jmin-1) * bv)&
                & / (ave_n(nx, i, jmin-1) ** 2.0d0 + ave_n(ny, i, jmin-1) ** 2.0d0)
        q_bc(3) = (-ave_n(ny, i, jmin-1) * bu + ave_n(nx, i, jmin-1) * bv)&
                & / (ave_n(nx, i, jmin-1) ** 2.0d0 + ave_n(ny, i, jmin-1) ** 2.0d0)
        !pは適当、e=bq(i, 0, 4)はe_0 = e_1
        bq(i, jmin-1, :) = qtobq(q_bc(1), q_bc(2), q_bc(3), 0.0d0)
        bq(i, jmin-1, 4) = bq(i, jmin, 4)

        !rho_-1 = rho_2, u_-1 = -u_2, e_-1 = e_2
        !temp_q(2) = u_2 = u, temp_q(3) = v_2 = v
        temp_q(:) = bqtoq(bq(i, jmin+1, :))
        q_bc(1) = temp_q(1)
        u = temp_q(2)
        v = temp_q(3)
        ! bu=U_2, bv=UU_2を求め、q_bc(2)=u_-1, q_bc(3)=v_-1を求める
        bu = ave_n(nx, i, jmin+1) * u + ave_n(ny, i, jmin+1) * v
        bv = -ave_n(ny, i, jmin+1) * u + ave_n(nx, i, jmin+1) * v
        q_bc(2) = (-ave_n(nx, i, jmin-2) * bu - ave_n(ny, i, jmin-2) * bv)&
                & / (ave_n(nx, i, jmin-2) ** 2.0d0 + ave_n(ny, i, jmin-2) ** 2.0d0)
        q_bc(3) = (-ave_n(ny, i, jmin-2) * bu + ave_n(nx, i, jmin-2) * bv)&
                & / (ave_n(nx, i, jmin-2) ** 2.0d0 + ave_n(ny, i, jmin-2) ** 2.0d0)
        !pは適当、e=bq(i, -1, 4)はe_-1 = e_2
        bq(i, jmin-2, :) = qtobq(q_bc(1), q_bc(2), q_bc(3), 0.0d0)
        bq(i, jmin-2, 4) = bq(i, jmin+1, 4)
    end do


    !BD4
    !slip condition
    !j+1=j, u_j+1 = -u_j, rho_j+1 = rho_j, e_j+1 = e_j
    do i = imin-2, imax+2
        !101=100, u_101 = u_100, rho_101 = rho_100, e_101 = e_100
        !temp_q(2) = u_100 = u, temp_q(3) = v_100 = v
        temp_q(:) = bqtoq(bq(i, jmax, :))
        q_bc(1) = temp_q(1)
        q_bc(2) = temp_q(2)
        q_bc(3) = temp_q(3)
        q_bc(4) = temp_q(4)
        bq(i, jmax+1, :) = qtobq(q_bc(1), q_bc(2), q_bc(3), q_bc(4))

        !102=99, u_102 = u_99, rho_102 = rho_99, e_102 = e_99
        !temp_q(2) = u_99 = u, temp_q(3) = v_99 = v
        temp_q(:) = bqtoq(bq(i, jmax-1, :))
        q_bc(1) = temp_q(1)
        q_bc(2) = temp_q(2)
        q_bc(3) = temp_q(3)
        q_bc(4) = temp_q(4)
        bq(i, jmax+2, :) = qtobq(q_bc(1), q_bc(2), q_bc(3), q_bc(4))
    end do


    !BD1
    !slip condition
    !i=i+1, u_i = u_i+1, uu_i = uu_i+1, rho_i = rho_i+1, e_i = e_i+1
    do j = jmin-2, jmax+2
        !0=1, u_0 = -u_1, uu_0 = uu_1, rho_0 = rho_1, e_0 = e_1
        !temp_q(2) = u_1 = u, temp_q(3) = v_1 = v
        temp_q(:) = bqtoq(bq(imin, j, :))
        q_bc(1) = temp_q(1)
        q_bc(2) = temp_q(2)
        q_bc(3) = temp_q(3)
        q_bc(4) = temp_q(4)
        bq(imin-1, j, :) = qtobq(q_bc(1), q_bc(2), q_bc(3), q_bc(4))

        !-1=2, u_-1 = u_2, uu_-1 = uu_2, rho_-1 = rho_2, e_-1 = e_2
        !rho_-1 = rho_2
        !temp_q(2) = u_2 = u, temp_q(3) = v_2 = v
        temp_q(:) = bqtoq(bq(imin+1, j, :))
        q_bc(1) = temp_q(1)
        q_bc(2) = temp_q(2)
        q_bc(3) = temp_q(3)
        q_bc(4) = temp_q(4)
        bq(imin-2, j, :) = qtobq(q_bc(1), q_bc(2), q_bc(3), q_bc(4))
    end do


    !BD2
    !slip condition
    !i+1=i, u_i+1 = -u_i, uu_i+1 = uu_i, rho_i+1 = rho_i, e_i+1 = e_i
    do j = jmin-2, jmax+2
        !101=100, u_101 = u_100, uu_101 = uu_100, rho_101 = rho_100, e_101 = e_100
        !temp_q(2) = u_100 = u, temp_q(3) = v_100 = v
        temp_q(:) = bqtoq(bq(imax, j, :))
        q_bc(1) = temp_q(1)
        q_bc(2) = temp_q(2)
        q_bc(3) = temp_q(3)
        q_bc(4) = temp_q(4)
        bq(imax+1, j, :) = qtobq(q_bc(1), q_bc(2), q_bc(3), q_bc(4))

        !102=99, u_102 = u_99, uu_102 = uu_99, rho_102 = rho_99, e_102 = e_99
        !temp_q(2) = u_99 = u, temp_q(3) = v_99 = v
        temp_q(:) = bqtoq(bq(imax-1, j, :))
        q_bc(1) = temp_q(1)
        q_bc(2) = temp_q(2)
        q_bc(3) = temp_q(3)
        q_bc(4) = temp_q(4)
        bq(imax+2, j, :) = qtobq(q_bc(1), q_bc(2), q_bc(3), q_bc(4))
    end do

end subroutine bound

!BD3ひな型
!!slip condition
!!j=j+1, u_j = -u_j+1, rho_j = rho_j+1, e_j = e_j+1
!do i = -2, 102
    !!rho_j = rho_j+1
    !!temp_q(2) = u_j+1 = u, temp_q(3) = v_j+1 = v
    !temp_q(:) = bqtoq(bq(i, j+1, :))
    !q_bc(1) = temp_q(1)
    !u = temp_q(2)
    !v = temp_q(3)
    !! bu=U_j+1, bv=UU_j+1を求め、q_bc(2)=u_j, q_bc(3)=v_jを求める
    !bu = ave_n(nx, i, j+1) * u + ave_n(ny, i, j+1) * v
    !bv = -ave_n(ny, i, j+1) * u + ave_n(nx, i, j+1) * v
    !q_bc(2) = (-ave_n(nx, i, j) * bu - ave_n(ny, i, j) * bv) / (ave_n(nx, i, j) ** 2.0d0 + ave_n(ny, i, j) ** 2.0d0)
    !q_bc(3) = (-ave_n(ny, i, j) * bu + ave_n(nx, i, j) * bv) / (ave_n(nx, i, j) ** 2.0d0 + ave_n(ny, i, j) ** 2.0d0)

    !!pは適当、e=bq(i, j, 4)はe_j = e_j+1
    !bq(i, j, :) = qtobq(q_bc(1), q_bc(2), q_bc(3), 0.0d0)
    !bq(i, j, 4) = bq(i, j+1, 4)
!end do

!BD4ひな型
!!slip condition
!!j+1=j, u_j+1 = -u_j, rho_j+1 = rho_j, e_j+1 = e_j
!do i = -2, 102
    !!rho_j+1 = rho_j
    !!temp_q(2) = u_j = u, temp_q(3) = v_j = v
    !temp_q(:) = bqtoq(bq(i, j, :))
    !q_bc(1) = temp_q(1)
    !u = temp_q(2)
    !v = temp_q(3)
    !! bu=U_j, bv=UU_jを求め、q_bc(2)=u_j+1, q_bc(3)=v_j+1を求める
    !bu = ave_n(nx, i, j) * u + ave_n(ny, i, j) * v
    !bv = -ave_n(ny, i, j) * u + ave_n(nx, i, j) * v
    !q_bc(2) = (-ave_n(nx, i, j+1) * bu - ave_n(ny, i, j+1) * bv) / (ave_n(nx, i, j+1) ** 2.0d0 + ave_n(ny, i, j+1) ** 2.0d0)
    !q_bc(3) = (-ave_n(ny, i, j+1) * bu + ave_n(nx, i, j+1) * bv) / (ave_n(nx, i, j+1) ** 2.0d0 + ave_n(ny, i, j+1) ** 2.0d0)

    !!pは適当、e=bq(i, j+1, 4)はe_j+1 = e_j
    !bq(i, j+1, :) = qtobq(q_bc(1), q_bc(2), q_bc(3), 0.0d0)
    !bq(i, j+1, 4) = bq(i, j, 4)
!end do

!BD1ひな型
!!slip condition
!!i=i+1, u_i = -u_i+1, uu_i = uu_i+1, rho_i = rho_i+1, e_i = e_i+1
!do j = -2, 102
    !!rho_i = rho_i+1
    !!temp_q(2) = u_i+1 = u, temp_q(3) = v_i+1 = v
    !temp_q(:) = bqtoq(bq(i+1, j, :))
    !q_bc(1) = temp_q(1)
    !u = temp_q(2)
    !v = temp_q(3)
    !! bu=U_i+1, bv=UU_i+1を求め、q_bc(2)=u_i, q_bc(3)=v_iを求める
    !bu = ave_m(mx, i+1, j) * u + ave_m(my, i+1, j) * v
    !bv = -ave_m(my, i+1, j) * u + ave_m(mx, i+1, j) * v
    !q_bc(2) = (-ave_m(mx, i, j) * bu - ave_m(my, i, j) * bv) / (ave_m(mx, i, j) ** 2.0d0 + ave_m(my, i, j) ** 2.0d0)
    !q_bc(3) = (-ave_m(my, i, j) * bu + ave_m(mx, i, j) * bv) / (ave_m(mx, i, j) ** 2.0d0 + ave_m(my, i, j) ** 2.0d0)

    !!pは適当、e=bq(i, j, 4)はe_i = e_i+1
    !bq(i, j, :) = qtobq(q_bc(1), q_bc(2), q_bc(3), 0.0d0)
    !bq(i, j, 4) = bq(i+1, j, 4)
!end do

!BD2ひな型
!!slip condition
!!i+1=i, u_i+1 = -u_i, uu_i+1 = uu_i, rho_i+1 = rho_i, e_i+1 = e_i
!do j = -2, 102
    !!rho_i+1 = rho_i
    !!temp_q(2) = u_i = u, temp_q(3) = v_i = v
    !temp_q(:) = bqtoq(bq(i, j, :))
    !q_bc(1) = temp_q(1)
    !u = temp_q(2)
    !v = temp_q(3)
    !! bu=U_i, bv=UU_iを求め、q_bc(2)=u_i+1, q_bc(3)=v_i+1を求める
    !bu = ave_m(mx, i, j) * u + ave_m(my, i, j) * v
    !bv = -ave_m(my, i, j) * u + ave_m(mx, i, j) * v
    !q_bc(2) = (-ave_m(mx, i+1, j) * bu - ave_m(my, i+1, j) * bv) / (ave_m(mx, i+1, j) ** 2.0d0 + ave_m(my, i+1, j) ** 2.0d0)
    !q_bc(3) = (-ave_m(my, i+1, j) * bu + ave_m(mx, i+1, j) * bv) / (ave_m(mx, i+1, j) ** 2.0d0 + ave_m(my, i+1, j) ** 2.0d0)

    !!pは適当、e=bq(i+1, j, 4)はe_i+1 = e_i
    !bq(i+1, j, :) = qtobq(q_bc(1), q_bc(2), q_bc(3), 0.0d0)
    !bq(i+1, j, 4) = bq(i, j, 4)
!end do
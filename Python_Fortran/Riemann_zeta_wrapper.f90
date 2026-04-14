!-----------------------------------------------------------------------
! Thin wrapper subroutines around Riemann_zeta_module for f2py.
! These accept double precision complex arrays, convert to/from quad
! precision internally, and call the module functions element-wise.
!-----------------------------------------------------------------------

subroutine Riemann_zeta_wrap(s, f, n)
    use Riemann_zeta_module
    implicit none
    integer, intent(in) :: n
    complex(8), intent(in) :: s(n)
    complex(8), intent(out) :: f(n)
!f2py intent(in) s
!f2py intent(out) f
!f2py intent(hide), depend(s) :: n = shape(s, 0)

    integer, parameter :: qp = selected_real_kind(33, 4931)
    complex(kind=qp) :: s_qp, f_qp
    integer :: k

    do k = 1, n
        s_qp = cmplx(real(s(k), kind=qp), aimag(s(k)), kind=qp)
        f_qp = Riemann_zeta(s_qp)
        f(k) = cmplx(real(f_qp, kind=8), real(aimag(f_qp), kind=8), kind=8)
    end do
end subroutine Riemann_zeta_wrap

subroutine Riemann_zeta_prime_wrap(s, f1, f2, n)
    use Riemann_zeta_module
    implicit none
    integer, intent(in) :: n
    complex(8), intent(in) :: s(n)
    complex(8), intent(out) :: f1(n), f2(n)
!f2py intent(in) s
!f2py intent(out) f1, f2
!f2py intent(hide), depend(s) :: n = shape(s, 0)

    integer, parameter :: qp = selected_real_kind(33, 4931)
    complex(kind=qp) :: s_qp, g(1:2)
    integer :: k

    do k = 1, n
        s_qp = cmplx(real(s(k), kind=qp), aimag(s(k)), kind=qp)
        g = Riemann_zeta_prime(s_qp)
        f1(k) = cmplx(real(g(1), kind=8), real(aimag(g(1)), kind=8), kind=8)
        f2(k) = cmplx(real(g(2), kind=8), real(aimag(g(2)), kind=8), kind=8)
    end do
end subroutine Riemann_zeta_prime_wrap

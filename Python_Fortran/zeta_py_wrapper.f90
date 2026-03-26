module zeta_py_wrapper
  use Riemann_zeta_module, only: Riemann_zeta
  implicit none
contains

  subroutine riemann_zeta_scalar(z_in, z_out)
    complex(kind=8), intent(in)  :: z_in
    complex(kind=8), intent(out) :: z_out
    complex(kind=16) :: s, f

    s = cmplx(real(z_in, kind=16), aimag(z_in), kind=16)
    f = Riemann_zeta(s)
    z_out = cmplx(real(f, kind=8), aimag(f), kind=8)
  end subroutine riemann_zeta_scalar

  subroutine riemann_zeta_vec(z_in, z_out)
    complex(kind=8), intent(in)  :: z_in(:)
    complex(kind=8), intent(out) :: z_out(size(z_in))
    integer :: k
    complex(kind=16) :: s, f

    do k = 1, size(z_in)
      s = cmplx(real(z_in(k), kind=16), aimag(z_in(k)), kind=16)
      f = Riemann_zeta(s)
      z_out(k) = cmplx(real(f, kind=8), aimag(f), kind=8)
    end do
  end subroutine riemann_zeta_vec

end module zeta_py_wrapper

module sphum_mod

use constants_mod,      only: rdgas, rvgas
use sat_vapor_pres_mod, only: escomp
use fms_mod,            only: WARNING
use land_debug_mod,     only: check_temp_range

implicit none
private

public :: qscomp

! ==== module constants ======================================================
real, parameter :: d622 = rdgas/rvgas
real, parameter :: d378 = 1.0-d622
real, parameter :: del_temp = 0.1 ! temperature increment for q_sat derivative calc.

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine qscomp(T, p, qsat, DqsatDT )
  real, intent(in) :: T    ! temperature
  real, intent(in) :: p    ! pressure
  real, intent(out):: qsat ! saturated specific humidity
  real, intent(out), optional :: DqsatDT ! deriv of specific humidity w.r.t. T

  real :: esat ! sat. water vapor pressure

  call check_temp_range(T,'qscomp','temperature')

  ! calculate saturated specific humidity
  call escomp(T,esat)
  qsat = d622*esat /(p-d378*esat )

  ! if requested, calculate the derivative of qsat w.r.t. temperature
  if (present(DqsatDT)) then
     call escomp(T+del_temp,esat)
     DqsatDT = (d622*esat/(p-d378*esat)-qsat)/del_temp
  endif
end subroutine qscomp

end module sphum_mod

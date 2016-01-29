! -*- mode: f90; indent-tabs-mode: nil; f90-do-indent:3; f90-if-indent:3; f90-type-indent:3; f90-program-indent:2; f90-associate-indent:0; f90-continuation-indent:5  -*-

module marbl_oxygen

  use kinds_mod, only : log_kind
  use kinds_mod, only : int_kind
  use kinds_mod, only : r8

  use constants, only : c0

  implicit none

  private

  public :: &
       schmidt_o2, &
       o2sat, &
       o2sat_scalar

contains

  !*****************************************************************************

  function SCHMIDT_O2(nx, ny, SST, LAND_MASK)

    ! !DESCRIPTION:
    !  Compute Schmidt number of O2 in seawater as function of SST
    !  where LAND_MASK is true. Give zero where LAND_MASK is false.
    !
    !  ref : Keeling et al, Global Biogeochem. Cycles, Vol. 12,
    !        No. 1, pp. 141-163, March 1998
    !
    ! !REVISION HISTORY:
    !  same as module

    ! !INPUT PARAMETERS:
    integer(int_kind), intent(in) :: nx, ny
    real (r8), intent(in) :: SST(nx, ny)

    logical (log_kind), intent(in) :: LAND_MASK(nx, ny)

    ! !OUTPUT PARAMETERS:

    real (r8) :: SCHMIDT_O2(nx, ny)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    real (r8), parameter :: &
         a = 1638.0_r8, &
         b = 81.83_r8, &
         c = 1.483_r8, &
         d = 0.008004_r8

    where (LAND_MASK)
       SCHMIDT_O2 = a + SST * (-b + SST * (c + SST * (-d)))
    elsewhere
       SCHMIDT_O2 = c0
    end where

  end function SCHMIDT_O2

  !*****************************************************************************

  function O2SAT(nx, ny, SST, SSS, LAND_MASK)

    !  Computes oxygen saturation concentration at 1 atm total pressure
    !  in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
    !  in permil) where LAND_MASK is true. Give zero where LAND_MASK is false.
    !
    !  FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
    !  THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
    !
    !  *** NOTE: THE "A_3*TS^2" TERM (IN THE PAPER) IS INCORRECT. ***
    !  *** IT SHOULD NOT BE THERE.                                ***
    !
    !  O2SAT IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
    !  0 permil <= S <= 42 permil
    !  CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil,
    !  O2SAT = 282.015 mmol/m^3
    !
    ! !INPUT PARAMETERS:

    integer(int_kind), intent(in) :: nx, ny
    real (r8), intent(in) :: SST(nx, ny) ! sea surface temperature (C)
    real (r8), intent(in) :: SSS(nx, ny) ! sea surface salinity (psu)

    logical (log_kind), intent(in) :: LAND_MASK(nx, ny)

    ! !OUTPUT PARAMETERS:

    real (r8) :: O2SAT(nx, ny)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer(int_kind) :: i, j

    do j = 1, ny
       do i = 1, nx
          if (LAND_MASK(i, j)) then
             O2SAT(i, j) = O2SAT_scalar(SST(i, j), SSS(i, j))
          else
             O2SAT(i, j) = c0
          end if
       end do
    end do

  end function O2SAT

  !-----------------------------------------------------------------------
  function O2SAT_scalar(SST, SSS)

    ! !DESCRIPTION:
    !
    !  Computes oxygen saturation concentration at 1 atm total pressure
    !  in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
    !  in permil) where LAND_MASK is true. Give zero where LAND_MASK is false.
    !
    !  FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
    !  THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
    !
    !  *** NOTE: THE "A_3*TS^2" TERM (IN THE PAPER) IS INCORRECT. ***
    !  *** IT SHOULD NOT BE THERE.                                ***
    !
    !  O2SAT IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
    !  0 permil <= S <= 42 permil
    !  CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil,
    !  O2SAT = 282.015 mmol/m^3
    !
    ! !REVISION HISTORY:
    !  same as module
    use constants, only : T0_Kelvin

    ! !INPUT PARAMETERS:

    real (r8), intent(in) :: SST    ! sea surface temperature (C)
    real (r8), intent(in) :: SSS    ! sea surface salinity (psu)

    ! !OUTPUT PARAMETERS:

    real (r8) :: O2SAT_scalar
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    real (r8) :: TS
    real (r8) :: O2SAT

    !  coefficients in expansion
    real (r8), parameter :: &
         a_0 = 2.00907_r8, &
         a_1 = 3.22014_r8, &
         a_2 = 4.05010_r8, &
         a_3 = 4.94457_r8, &
         a_4 = -2.56847E-1_r8, &
         a_5 = 3.88767_r8, &
         b_0 = -6.24523E-3_r8, &
         b_1 = -7.37614E-3_r8, &
         b_2 = -1.03410E-2_r8, &
         b_3 = -8.17083E-3_r8, &
         c_0 = -4.88682E-7_r8

    TS = log( ((T0_Kelvin + 25.0_r8) - SST) / (T0_Kelvin + SST) )

    O2SAT = exp(a_0+TS*(a_1+TS*(a_2+TS*(a_3+TS*(a_4+TS*a_5)))) + &
         SSS*( (b_0+TS*(b_1+TS*(b_2+TS*b_3))) + SSS*c_0 ))


    !  Convert from ml/l to mmol/m^3
    O2SAT_scalar = O2SAT / 0.0223916_r8

  end function O2SAT_scalar

end module marbl_oxygen

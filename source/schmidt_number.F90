module schmidt_number

  implicit none

  private

  public :: &
       schmidt_co2,  &
       schmidt_co2_surf

contains
  
  !*****************************************************************************

  function schmidt_co2(SST, LAND_MASK)
    !
    !  Compute Schmidt number of CO2 in seawater as function of SST
    !  where LAND_MASK is true. Give zero where LAND_MASK is false.
    !
    !  ref : Wanninkhof, J. Geophys. Res, Vol. 97, No. C5,
    !  pp. 7373-7382, May 15, 1992
    !
    
    use kinds_mod, only: r8, log_kind
    use constants, only: c0
    use blocks, only: nx_block, ny_block

    ! INPUT PARAMETERS:
    real (r8), dimension(nx_block, ny_block), intent(in) :: SST
    logical (log_kind), dimension(nx_block, ny_block), intent(in) :: LAND_MASK
   
    ! OUTPUT PARAMETERS:
    real (r8), dimension(nx_block, ny_block) :: SCHMIDT_CO2

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    real (r8), parameter :: &
         a = 2073.1_r8, &
         b = 125.62_r8, &
         c = 3.6276_r8, &
         d = 0.043219_r8

    where (LAND_MASK)
       SCHMIDT_CO2 = a + SST * (-b + SST * (c + SST * (-d)))
    elsewhere
       SCHMIDT_CO2 = c0
    end where

  end function schmidt_co2

  !*****************************************************************************

  function schmidt_co2_surf(n, SST, LAND_MASK)
    !
    !  Compute Schmidt number of CO2 in seawater as function of SST
    !  where LAND_MASK is true. Give zero where LAND_MASK is false.
    !
    !  ref : Wanninkhof, J. Geophys. Res, Vol. 97, No. C5,
    !  pp. 7373-7382, May 15, 1992
    !
    
    use kinds_mod, only: r8, log_kind, int_kind
    use constants, only: c0

    ! INPUT PARAMETERS:
    integer(int_kind), intent(in) :: n
    real (r8), dimension(n), intent(in) :: SST
    logical (log_kind), dimension(n), intent(in) :: LAND_MASK
   
    ! OUTPUT PARAMETERS:
    real (r8), dimension(n) :: SCHMIDT_CO2_SURF

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    real (r8), parameter :: &
         a = 2073.1_r8, &
         b = 125.62_r8, &
         c = 3.6276_r8, &
         d = 0.043219_r8

    where (LAND_MASK)
       SCHMIDT_CO2_SURF = a + SST * (-b + SST * (c + SST * (-d)))
    elsewhere
       SCHMIDT_CO2_SURF = c0
    end where

  end function schmidt_co2_surf

  !*****************************************************************************
  
end module schmidt_number

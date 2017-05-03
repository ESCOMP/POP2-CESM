module schmidt_number

  implicit none

  private

  public :: schmidt_co2

contains
  
  !*****************************************************************************

  function schmidt_co2(SST_IN, LAND_MASK)
    !
    !  Compute Schmidt number of CO2 in seawater as function of SST
    !  where LAND_MASK is true. Give zero where LAND_MASK is false.
    !
    !  range of validity of fit is -2:40
    !
    !  Ref : Wanninkhof 2014, Relationship between wind speed
    !        and gas exchange over the ocean revisited,
    !        Limnol. Oceanogr.: Methods, 12,
    !        doi:10.4319/lom.2014.12.351
    !
    
    use kinds_mod, only: r8, log_kind, int_kind
    use constants, only: c0
    use blocks, only: nx_block, ny_block

    ! INPUT PARAMETERS:
    real (r8)         , intent(in) :: SST_IN(nx_block,ny_block)
    logical (log_kind), intent(in) :: LAND_MASK(nx_block,ny_block)
   
    ! OUTPUT PARAMETERS:
    real (r8) :: SCHMIDT_CO2(nx_block, ny_block)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer(int_kind)    :: i, j
    real (r8)            :: SST(nx_block,ny_block)

    real (r8), parameter :: a = 2116.8_r8
    real (r8), parameter :: b = -136.25_r8
    real (r8), parameter :: c =    4.7353_r8
    real (r8), parameter :: d =   -0.092307_r8
    real (r8), parameter :: e =    0.0007555_r8

    !-----------------------------------------------------------------------

    do j = 1, ny_block
       do i = 1, nx_block
          if (LAND_MASK(i,j)) then
             SST(i,j) = max(-2.0_r8, min(40.0_r8, SST_IN(i,j)))
             SCHMIDT_CO2(i,j) = a + SST(i,j) * (b + SST(i,j) * (c + SST(i,j) * (d + SST(i,j) * e)))
          else
             SCHMIDT_CO2(i,j) = c0
          endif
       end do
    end do

  end function schmidt_co2

  !*****************************************************************************
  
end module schmidt_number

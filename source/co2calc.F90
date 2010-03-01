MODULE co2calc

  !-----------------------------------------------------------------------------
  !   based upon OCMIP2 co2calc
  !
  !   CVS:$Id: co2calc.F90 941 2006-05-12 21:36:48Z klindsay $
  !   CVS:$Name$
  !-----------------------------------------------------------------------------

  USE blocks,      ONLY : nx_block
  USE domain_size, ONLY : max_blocks_clinic
  USE kinds_mod

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  !   public/private declarations
  !-----------------------------------------------------------------------------

  PRIVATE
  PUBLIC :: co2calc_row

  !-----------------------------------------------------------------------------
  !   module parameters
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !   The current setting of xacc, a tolerance critera, will result in co2star
  !   being accurate to 3 significant figures (xx.y). Making xacc bigger will
  !   result in faster convergence also, but this is not recommended (xacc of
  !   10**-9 drops precision to 2 significant figures).
  !-----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: xacc = 1e-10_r8
  INTEGER(KIND=int_kind), PARAMETER :: maxit = 100

  !-----------------------------------------------------------------------------
  !   declarations for function coefficients & species concentrations
  !-----------------------------------------------------------------------------

  REAL(KIND=r8), DIMENSION(nx_block,max_blocks_clinic) :: &
       k0, k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, ff, &
       bt, st, ft, dic, ta, pt, sit

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  SUBROUTINE co2calc_row(iblock, mask, t, s, dic_in, ta_in, pt_in, sit_in, &
       phlo, phhi, ph, xco2_in, atmpres, co2star, dco2star, pCO2surf, dpco2)

    !---------------------------------------------------------------------------
    !   SUBROUTINE CO2CALC
    !
    !   PURPOSE : Calculate delta co2*, etc. from total alkalinity, total CO2,
    !             temp, salinity (s), etc.
    !---------------------------------------------------------------------------

    USE constants, ONLY : c0, c1, c1p5, c10, c1000, rho_sw, T0_Kelvin

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: iblock
    LOGICAL(KIND=log_kind), DIMENSION(nx_block), INTENT(IN) :: mask
    REAL(KIND=r8), DIMENSION(nx_block), INTENT(IN) :: &
         t,        & ! temperature (degrees C)
         s,        & ! salinity (PSU)
         dic_in,   & ! total inorganic carbon (nmol/cm^3)
         ta_in,    & ! total alkalinity (neq/cm^3)
         pt_in,    & ! inorganic phosphate (nmol/cm^3)
         sit_in,   & ! inorganic silicate (nmol/cm^3)
         phlo,     & ! lower limit of pH range
         phhi,     & ! upper limit of pH range
         xco2_in,  & ! atmospheric mole fraction CO2 in dry air (ppmv)
         atmpres     ! atmospheric pressure (atmosphere)

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=r8), DIMENSION(nx_block), INTENT(OUT) :: &
         ph,       & ! computed ph values, for initial guess on next time step
         co2star,  & ! CO2*water (nmol/cm^3)
         dco2star, & ! delta CO2 (nmol/cm^3)
         pco2surf, & ! oceanic pCO2 (ppmv)
         dpco2       ! Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind) :: i

    REAL(KIND=r8) :: &
         mass_to_vol,  & ! (mol/kg) -> (mmol/m^3)
         vol_to_mass,  & ! (mmol/m^3) -> (mol/kg)
         tk,           & ! temperature (K)
         is,           & ! ionic strength
         scl,          & ! chlorinity
         co2starair,   & ! co2star saturation
         tk100, tk1002, invtk, dlogtk, is2, sqrtis, &
         s2, sqrts, s15, htotal2

    REAL(KIND=r8), DIMENSION(nx_block) :: &
         xco2,         & ! atmospheric CO2 (atm)
         htotal,       & ! free concentration of H ion
         x1, x2          ! bounds on htotal for solver

    !---------------------------------------------------------------------------
    !   check for existence of ocean points
    !---------------------------------------------------------------------------

    IF (COUNT(mask) == 0) THEN
       ph          = c0
       co2star     = c0
       dco2star    = c0
       pCO2surf    = c0
       dpCO2       = c0
       RETURN
    END IF

    !---------------------------------------------------------------------------
    !   set unit conversion factors
    !---------------------------------------------------------------------------

    mass_to_vol = 1e6_r8 * rho_sw
    vol_to_mass = c1 / mass_to_vol

    !---------------------------------------------------------------------------
    !   convert tracer units to per mass & xco2 from uatm to atm
    !---------------------------------------------------------------------------

    DO i = 1,nx_block
       IF (mask(i)) THEN
          dic(i,iblock)  = dic_in(i)  * vol_to_mass
          ta(i,iblock)   = ta_in(i)   * vol_to_mass
          pt(i,iblock)   = max(pt_in(i),c0)   * vol_to_mass
          sit(i,iblock)  = max(sit_in(i),c0)  * vol_to_mass
          xco2(i) = xco2_in(i) * 1e-6_r8

          !---------------------------------------------------------------------
          !   Calculate all constants needed to convert between various
          !   measured carbon species. References for each equation are
          !   noted in the code.  Once calculated, the constants are stored
          !   and passed in the common block "const". The original version
          !   of this code was based on the code by Dickson in Version 2 of
          !   "Handbook of Methods for the Analysis of the Various Parameters
          !   of the Carbon Dioxide System in Seawater", DOE, 1994 (SOP No. 3,
          !   p25-26).
          !   Derive simple terms used more than once
          !---------------------------------------------------------------------

          tk       = T0_Kelvin + t(i)
          tk100    = tk * 1e-2_r8
          tk1002   = tk100 * tk100
          invtk    = c1 / tk
          dlogtk   = LOG(tk)

          is       = 19.924_r8 * s(i) / (c1000 - 1.005_r8 * s(i))
          is2      = is * is
          sqrtis   = SQRT(is)
          sqrts    = SQRT(s(i))
          s15      = s(i) ** c1p5
          s2       = s(i) ** 2
          scl      = s(i) / 1.80655_r8

          !---------------------------------------------------------------------
          !   f = k0(1-pH2O)*correction term for non-ideality
          !   Weiss & Price (1980, Mar. Chem., 8, 347-359;
          !                 Eq 13 with table 6 values)
          !---------------------------------------------------------------------

          ff(i,iblock) = EXP(-162.8301_r8 + 218.2968_r8 / tk100 + &
               90.9241_r8 * LOG(tk100) - 1.47696_r8 * tk1002 + &
               s(i) * (.025695_r8 - .025225_r8 * tk100 + &
               0.0049867_r8 * tk1002))

          !---------------------------------------------------------------------
          !   K0 from Weiss 1974
          !---------------------------------------------------------------------

          k0(i,iblock) = EXP(93.4517_r8 / tk100 - 60.2409_r8 + &
               23.3585_r8 * LOG(tk100) + s(i) * (.023517_r8 - &
               0.023656_r8 * tk100 + 0.0047036_r8 * tk1002))

          !---------------------------------------------------------------------
          !   k1 = [H][HCO3]/[H2CO3]
          !   k2 = [H][CO3]/[HCO3]
          !   Millero p.664 (1995) using Mehrbach et al. data on seawater scale
          !---------------------------------------------------------------------

          k1(i,iblock) = c10 ** ((-c1) * (3670.7_r8 * invtk - 62.008_r8 + &
               9.7944_r8 * dlogtk - 0.0118_r8 * s(i) + &
               0.000116_r8 * s2))

          k2(i,iblock) = c10 ** ((-c1) * (1394.7_r8 * invtk + 4.777_r8 - &
               0.0184_r8 * s(i) + 0.000118_r8 * s2))

          !---------------------------------------------------------------------
          !   kb = [H][BO2]/[HBO2]
          !   Millero p.669 (1995) using data from Dickson (1990)
          !---------------------------------------------------------------------

          kb(i,iblock) = EXP((-8966.90_r8 - 2890.53_r8 * sqrts - &
               77.942_r8 * s(i) + 1.728_r8 * s15 - &
               0.0996_r8 * s2) * invtk + (148.0248_r8 + &
               137.1942_r8 * sqrts + 1.62142_r8 * s(i)) + &
               (-24.4344_r8 - 25.085_r8 * sqrts - &
               0.2474_r8 * s(i)) * dlogtk + &
               0.053105_r8 * sqrts * tk)

          !---------------------------------------------------------------------
          !   k1p = [H][H2PO4]/[H3PO4]
          !   DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
          !---------------------------------------------------------------------

          k1p(i,iblock) = EXP(-4576.752_r8 * invtk + 115.525_r8 - &
               18.453_r8 * dlogtk + &
               (-106.736_r8 * invtk + 0.69171_r8) * sqrts + &
               (-0.65643_r8 * invtk - 0.01844_r8) * s(i))

          !---------------------------------------------------------------------
          !   k2p = [H][HPO4]/[H2PO4]
          !   DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
          !---------------------------------------------------------------------

          k2p(i,iblock) = EXP(-8814.715_r8 * invtk + 172.0883_r8 - &
               27.927_r8 * dlogtk + &
               (-160.340_r8 * invtk + 1.3566_r8) * sqrts + &
               (0.37335_r8 * invtk - 0.05778_r8) * s(i))

          !---------------------------------------------------------------------
          !   k3p = [H][PO4]/[HPO4]
          !   DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
          !---------------------------------------------------------------------

          k3p(i,iblock) = EXP(-3070.75_r8 * invtk - 18.141_r8 + &
               (17.27039_r8 * invtk + 2.81197_r8) * sqrts + &
               (-44.99486_r8 * invtk - 0.09984_r8) * s(i))

          !---------------------------------------------------------------------
          !   ksi = [H][SiO(OH)3]/[Si(OH)4]
          !   Millero p.671 (1995) using data from Yao and Millero (1995)
          !---------------------------------------------------------------------

          ksi(i,iblock) = EXP(-8904.2_r8 * invtk + 117.385_r8 - &
               19.334_r8 * dlogtk + &
               (-458.79_r8 * invtk + 3.5913_r8) * sqrtis + &
               (188.74_r8 * invtk - 1.5998_r8) * is + &
               (-12.1652_r8 * invtk + 0.07871_r8) * is2 + &
               LOG(c1 - 0.001005_r8 * s(i)))

          !---------------------------------------------------------------------
          !   kw = [H][OH]
          !   Millero p.670 (1995) using composite data
          !---------------------------------------------------------------------

          kw(i,iblock) = EXP(-13847.26_r8 * invtk + 148.9652_r8 - &
               23.6521_r8 * dlogtk + (118.67_r8 * invtk - &
               5.977_r8 + 1.0495_r8 * dlogtk) * sqrts - &
               0.01615_r8 * s(i))

          !---------------------------------------------------------------------
          !   ks = [H][SO4]/[HSO4]
          !   Dickson (1990, J. chem. Thermodynamics 22, 113)
          !---------------------------------------------------------------------

          ks(i,iblock) = EXP(-4276.1_r8 * invtk + 141.328_r8 - &
               23.093_r8 * dlogtk + (-13856.0_r8 * invtk + &
               324.57_r8 - 47.986_r8 * dlogtk) * sqrtis + &
               (35474.0_r8 * invtk - 771.54_r8 + &
               114.723_r8 * dlogtk) * is - &
               2698.0_r8 * invtk * is ** c1p5 + &
               1776.0_r8 * invtk * is2 + &
               LOG(c1 - 0.001005_r8 * s(i)))

          !---------------------------------------------------------------------
          !   kf = [H][F]/[HF]
          !   Dickson and Riley (1979) -- change pH scale to total
          !---------------------------------------------------------------------

          kf(i,iblock) = EXP(1590.2_r8 * invtk - 12.641_r8 + &
               1.525_r8 * sqrtis + &
               LOG(c1 - 0.001005_r8 * s(i)) +  &
               LOG(c1 + (0.1400_r8 / 96.062_r8) * (scl) / ks(i,iblock)))

          !---------------------------------------------------------------------
          !   Calculate concentrations for borate, sulfate, and fluoride
          !   bt : Uppstrom (1974)
          !   st : Morris & Riley (1966)
          !   ft : Riley (1965)
          !---------------------------------------------------------------------

          bt(i,iblock) = 0.000232_r8 * scl / 10.811_r8
          st(i,iblock) = 0.14_r8 * scl / 96.062_r8
          ft(i,iblock) = 0.000067_r8 * scl / 18.9984_r8

          x1(i) = c10 ** (-phhi(i))
          x2(i) = c10 ** (-phlo(i))

       END IF ! if mask
    END DO ! i loop

    !---------------------------------------------------------------------------
    !   If DIC and TA are known then either a root finding or iterative
    !   method must be used to calculate htotal. In this case we use
    !   the Newton-Raphson "safe" method taken from "Numerical Recipes"
    !   (function "rtsafe.f" with error trapping removed).
    !
    !   As currently set, this procedure iterates about 12 times. The
    !   x1 and x2 values set below will accomodate ANY oceanographic
    !   values. If an initial guess of the pH is known, then the
    !   number of iterations can be reduced to about 5 by narrowing
    !   the gap between x1 and x2. It is recommended that the first
    !   few time steps be run with x1 and x2 set as below. After that,
    !   set x1 and x2 to the previous value of the pH +/- ~0.5.
    !---------------------------------------------------------------------------

    CALL drtsafe_row(iblock, mask, x1, x2, xacc, htotal)

    !---------------------------------------------------------------------------
    !   Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
    !   ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
    !---------------------------------------------------------------------------

    DO i = 1,nx_block
       IF (mask(i)) THEN

          htotal2 = htotal(i) ** 2
          co2star(i) = dic(i,iblock) * htotal2 / &
               (htotal2 + k1(i,iblock) * htotal(i) + k1(i,iblock) * k2(i,iblock))
          co2starair = xco2(i) * ff(i,iblock) * atmpres(i)
          dco2star(i) = co2starair - co2star(i)
          ph(i) = -LOG10(htotal(i))

          !---------------------------------------------------------------------
          !   Add two output arguments for storing pCO2surf
          !   Should we be using K0 or ff for the solubility here?
          !---------------------------------------------------------------------

          pCO2surf(i) = co2star(i) / ff(i,iblock)
          dpCO2(i)    = pCO2surf(i) - xco2(i) * atmpres(i)

          !---------------------------------------------------------------------
          !   Convert units of output arguments
          !   Note: pCO2surf and dpCO2 are calculated in atm above.
          !---------------------------------------------------------------------

          co2star(i)  = co2star(i) * mass_to_vol
          dco2star(i) = dco2star(i) * mass_to_vol

          pCO2surf(i) = pCO2surf(i) * 1e6_r8
          dpCO2(i)    = dpCO2(i) * 1e6_r8

       ELSE ! if mask

          ph(i)       = c0
          co2star(i)  = c0
          dco2star(i) = c0
          pCO2surf(i) = c0
          dpCO2(i)    = c0

       END IF ! if mask
    END DO ! i loop

  END SUBROUTINE co2calc_row

  !*****************************************************************************

  SUBROUTINE talk_row(iblock, mask, x, fn, df)

    !---------------------------------------------------------------------------
    !   This routine computes TA as a function of DIC, htotal and constants.
    !   It also calculates the derivative of this function with respect to
    !   htotal. It is used in the iterative solution for htotal. In the call
    !   "x" is the input value for htotal, "fn" is the calculated value for
    !   TA and "df" is the value for dTA/dhtotal.
    !---------------------------------------------------------------------------

    USE constants, ONLY : c1, c2, c3

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: iblock
    LOGICAL(KIND=log_kind), DIMENSION(nx_block), INTENT(IN) :: mask
    REAL(KIND=r8), DIMENSION(nx_block), INTENT(IN) :: x

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=r8), DIMENSION(nx_block), INTENT(OUT) :: fn, df

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind) :: i

    REAL(KIND=r8) :: &
         x1, x2, x3, k12, k12p, k123p, a, a2, da, b, b2, db, c

    !---------------------------------------------------------------------------

    DO i = 1,nx_block
       IF (mask(i)) THEN
          x1 = x(i)
          x2 = x1 * x1
          x3 = x2 * x1
          k12 = k1(i,iblock) * k2(i,iblock)
          k12p = k1p(i,iblock) * k2p(i,iblock)
          k123p = k12p * k3p(i,iblock)
          a = x3 + k1p(i,iblock) * x2 + k12p * x1 + k123p
          a2 = a * a
          da = c3 * x2 + c2 * k1p(i,iblock) * x1 + k12p
          b = x2 + k1(i,iblock) * x1 + k12
          b2 = b * b
          db = c2 * x1 + k1(i,iblock)
          c = c1 + st(i,iblock) / ks(i,iblock)

          !---------------------------------------------------------------------
          !   fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
          !---------------------------------------------------------------------

          fn(i) = k1(i,iblock) * x1 * dic(i,iblock) / b + &
               c2 * dic(i,iblock) * k12 / b + &
               bt(i,iblock) / (c1 + x1 / kb(i,iblock)) + &
               kw(i,iblock) / x1 + &
               pt(i,iblock) * k12p * x1 / a + &
               c2 * pt(i,iblock) * k123p / a + &
               sit(i,iblock) / (c1 + x1 / ksi(i,iblock)) - &
               x1 / c - &
               st(i,iblock) / (c1 + ks(i,iblock) / x1 / c) - &
               ft(i,iblock) / (c1 + kf(i,iblock) / x1) - &
               pt(i,iblock) * x3 / a - &
               ta(i,iblock)

          !---------------------------------------------------------------------
          !   df = d(fn)/dx
          !---------------------------------------------------------------------

          df(i) = ((k1(i,iblock) * dic(i,iblock) * b) - k1(i,iblock) * x1 * dic(i,iblock) * db) / b2 - &
               c2 * dic(i,iblock) * k12 * db / b2 - &
               bt(i,iblock) / kb(i,iblock) / (c1 + x1 / kb(i,iblock)) ** 2 - &
               kw(i,iblock) / x2 + &
               (pt(i,iblock) * k12p * (a - x1 * da)) / a2 - &
               c2 * pt(i,iblock) * k123p * da / a2 - &
               sit(i,iblock) / ksi(i,iblock) / (c1 + x1 / ksi(i,iblock)) ** 2 - &
               c1 / c + &
               st(i,iblock) * (c1 + ks(i,iblock) / x1 / c) ** (-2) * (ks(i,iblock) / c / x2) + &
               ft(i,iblock) * (c1 + kf(i,iblock) / x1) ** (-2) * kf(i,iblock) / x2 - &
               pt(i,iblock) * x2 * (c3 * a - x1 * da) / a2

       END IF ! if mask
    END DO ! i loop

  END SUBROUTINE talk_row

  !*****************************************************************************

  SUBROUTINE drtsafe_row(iblock, mask_in, x1, x2, xacc, soln)

    !---------------------------------------------------------------------------
    !   Vectorized version of drtsafe, which was a modified version of
    !      Numerical Recipes algorithm.
    !   Keith Lindsay, Oct 1999
    !
    !   Algorithm comment :
    !      Iteration from Newtons method is used unless it leaves
    !      bracketing interval or the dx is > 0.5 the previous dx.
    !      In that case, bisection method is used.
    !---------------------------------------------------------------------------

    USE constants, ONLY : c0, c2, p5
    !USE shr_sys_mod, ONLY : shr_sys_abort

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: iblock
    LOGICAL(KIND=log_kind), DIMENSION(nx_block), INTENT(IN) :: mask_in
    REAL(KIND=r8), DIMENSION(nx_block), INTENT(IN) :: x1, x2
    REAL(KIND=r8), INTENT(IN) :: xacc

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=r8), DIMENSION(nx_block), INTENT(OUT) :: soln

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    LOGICAL(KIND=log_kind) :: leave_bracket, dx_decrease
    LOGICAL(KIND=log_kind), DIMENSION(nx_block) :: mask
    INTEGER(KIND=int_kind) ::  i, it
    REAL(KIND=r8) :: temp
    REAL(KIND=r8), DIMENSION(nx_block) :: xlo, xhi, flo, fhi, f, df, dxold, dx

    !---------------------------------------------------------------------------
    !   bracket root at each location and set up first iteration
    !---------------------------------------------------------------------------

    mask = mask_in

    CALL talk_row(iblock, mask, x1, flo, df)
    CALL talk_row(iblock, mask, x2, fhi, df)

    DO i = 1,nx_block
       IF (mask(i)) THEN
          IF (flo(i) .LT. c0) THEN
             xlo(i) = x1(i)
             xhi(i) = x2(i)
          ELSE
             xlo(i) = x2(i)
             xhi(i) = x1(i)
             temp = flo(i)
             flo(i) = fhi(i)
             fhi(i) = temp
          END IF
          soln(i) = p5 * (xlo(i) + xhi(i))
          dxold(i) = ABS(xlo(i) - xhi(i))
          dx(i) = dxold(i)
       END IF
    END DO

    CALL talk_row(iblock, mask, soln, f, df)

    !---------------------------------------------------------------------------
    !   perform iterations, zeroing mask when a location has converged
    !---------------------------------------------------------------------------

    DO it = 1,maxit
       DO i = 1,nx_block
          IF (mask(i)) THEN
             leave_bracket = ((soln(i) - xhi(i)) * df(i) - f(i)) * &
                  ((soln(i) - xlo(i)) * df(i) - f(i)) .GE. 0
             dx_decrease = ABS(c2 * f(i)) .LE. ABS(dxold(i) * df(i))
             IF (leave_bracket .OR. .NOT. dx_decrease) THEN
                dxold(i) = dx(i)
                dx(i) = p5 * (xhi(i) - xlo(i))
                soln(i) = xlo(i) + dx(i)
                IF (xlo(i) .EQ. soln(i)) mask(i) = .FALSE.
             ELSE
                dxold(i) = dx(i)
                dx(i) = -f(i) / df(i)
                temp = soln(i)
                soln(i) = soln(i) + dx(i)
                IF (temp .EQ. soln(i)) mask(i) = .FALSE.
             END IF
             IF (ABS(dx(i)) .LT. xacc) mask(i) = .FALSE.
          END IF
       END DO

       IF (.NOT. ANY(mask)) RETURN

       CALL talk_row(iblock, mask, soln, f, df)

       DO i = 1,nx_block
          IF (mask(i)) THEN
             IF (f(i) .LT. c0) THEN
                xlo(i) = soln(i)
                flo(i) = f(i)
             ELSE
                xhi(i) = soln(i)
                fhi(i) = f(i)
             END IF
          END IF
       END DO

    END DO ! iteration loop

    !CALL shr_sys_abort('lack of convergence in drtsafe_row')

  END SUBROUTINE drtsafe_row

  !*****************************************************************************

END MODULE co2calc

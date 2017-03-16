!*********************************************************************************
! Part of the MEDINA: MECCA - KPP Fortran to CUDA source-to-source pre-processor
!*********************************************************************************
!
! Extracted Minimal Definitions of machine precision and physical constants 
! as Fortran PARAMETERs for MESSy
!
! Author: Michail Alvanos
!
! Original Authors:
!  - Rolf Sander     
!  - Patrick Joeckel 
!

MODULE messy_main_constants_mem

  ! PHYSICAL CONSTANTS (CODATA Recommended Values, 2010, 
  ! http://physics.nist.gov/cuu/Constants/)
  REAL(dp), PARAMETER :: pi       = 3.14159265358979323846_dp
  REAL(dp), PARAMETER :: R_gas    = 8.3144621_dp     ! R [J/K/mol]
  REAL(dp), PARAMETER :: h_Planck = 6.62606957E34_dp ! Planck constant [Js]
  REAL(dp), PARAMETER :: c_light  = 2.99792458E8_dp  ! speed of light [m/s]
  REAL(dp), PARAMETER :: stbo     = 5.670373E-8_dp   ! Stephan-Boltzmann constant [W/m2/K4]
  REAL(dp), PARAMETER :: N_A      = 6.02214129E23_dp ! Avogadro constant [1/mol]
  REAL(dp), PARAMETER :: N_A_kmol = 6.02214129E26_dp ! Avogadro constant [1/kmol] 
#ifndef MESSYIDTC
  REAL(dp), PARAMETER :: g        = 9.80665_dp       ! gravity acceleration [m/s2]
#else
  REAL(dp), PARAMETER :: g        = 9.80616_dp       ! gravity acceleration [m/s2]
#endif
  REAL(dp), PARAMETER :: T0       = 298.15_dp        ! standard temperature [K]
  REAL(dp), PARAMETER :: T0_INV   = 1._DP / T0       ! 1/T0 [1/K]
#ifndef MESSYIDTC
  REAL(dp), PARAMETER :: atm2Pa   = 101325._dp       ! conversion from [atm] to [Pa]
#else
  REAL(dp), PARAMETER :: atm2Pa   = 100000._dp       ! conversion from [atm] to [Pa]
#endif

END MODULE messy_main_constants_mem


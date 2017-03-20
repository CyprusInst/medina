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

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)  
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
  REAL(DP), PARAMETER :: TINY_DP = TINY(0._dp) ! mz_rs_20060114
  REAL(DP), PARAMETER :: HUGE_DP = HUGE(0._dp) ! mz_rs_20100409


  ! PHYSICAL CONSTANTS (CODATA Recommended Values, 2010, 
  ! http://physics.nist.gov/cuu/Constants/)
  REAL(dp), PARAMETER :: pi       = 3.14159265358979323846_dp
  REAL(dp), PARAMETER :: R_gas    = 8.3144621_dp     ! R [J/K/mol]
  REAL(dp), PARAMETER :: h_Planck = 6.62606957E34_dp ! Planck constant [Js]
  REAL(dp), PARAMETER :: c_light  = 2.99792458E8_dp  ! speed of light [m/s]
  REAL(dp), PARAMETER :: stbo     = 5.670373E-8_dp   ! Stephan-Boltzmann constant [W/m2/K4]
  REAL(dp), PARAMETER :: N_A      = 6.02214129E23_dp ! Avogadro constant [1/mol]
  REAL(dp), PARAMETER :: N_A_kmol = 6.02214129E26_dp ! Avogadro constant [1/kmol] 
  REAL(dp), PARAMETER :: g        = 9.80665_dp       ! gravity acceleration [m/s2]
  REAL(dp), PARAMETER :: T0       = 298.15_dp        ! standard temperature [K]
  REAL(dp), PARAMETER :: T0_INV   = 1._DP / T0       ! 1/T0 [1/K]
  REAL(dp), PARAMETER :: atm2Pa   = 101325._dp       ! conversion from [atm] to [Pa]

END MODULE messy_main_constants_mem


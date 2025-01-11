!*********************************************************************************
! Part of the MEDINA: MECCA - KPP Fortran to CUDA source-to-source pre-processor
!*********************************************************************************
!
! Definitions that all photolysis submodels have in common. Used for extraction
! of definitions in accelerated code.
!
! Author: Michail Alvanos
!
! Original Authors:
!  - Rolf Sander     
!

MODULE messy_cmn_photol_mem

  IMPLICIT NONE

  ! ip_* = index of photolysis
  INTEGER, PUBLIC, PARAMETER :: &
    ip_O2       =  1, ip_O3P      =  2, ip_O1D      =  3, ip_H2O2     =  4, &
    ip_NO2      =  5, ip_NO2O     =  6, ip_NOO2     =  7, ip_N2O5     =  8, &
    ip_HNO3     =  9, ip_HNO4     = 10, ip_PAN      = 11, ip_HONO     = 12, &
    ip_CH3OOH   = 13, ip_COH2     = 14, ip_CHOH     = 15, ip_CH3CO3H  = 16, &
    ip_CH3CHO   = 17, ip_CH3COCH3 = 18, ip_MGLYOX   = 19, ip_HOCl     = 20, &
    ip_OClO     = 21, ip_Cl2O2    = 22, ip_ClNO3    = 23, ip_ClNO2    = 24, &
    ip_Cl2      = 25, ip_BrO      = 26, ip_HOBr     = 27, ip_BrCl     = 28, &
    ip_BrNO3    = 29, ip_BrNO2    = 30, ip_Br2      = 31, ip_CCl4     = 32, &
    ip_CH3Cl    = 33, ip_CH3CCl3  = 34, ip_CFCl3    = 35, ip_CF2Cl2   = 36, &
    ip_CH3Br    = 37, ip_CF2ClBr  = 38, ip_CF3Br    = 39, ip_CH3I     = 40, &
    ip_C3H7I    = 41, ip_CH2ClI   = 42, ip_CH2I2    = 43, ip_IO       = 44, &
    ip_HOI      = 45, ip_I2       = 46, ip_ICl      = 47, ip_IBr      = 48, &
    ip_INO2     = 49, ip_INO3     = 50, ip_SO2      = 51, ip_SO3      = 52, &
    ip_OCS      = 53, ip_CS2      = 54, ip_H2O      = 55, ip_N2O      = 56, &
    ip_NO       = 57, ip_CO2      = 58, ip_HCl      = 59, ip_CHCl2Br  = 60, &
    ip_CHClBr2  = 61, ip_CH2ClBr  = 62, ip_CH2Br2   = 63, ip_CHBr3    = 64, &
    ip_SF6      = 65, ip_NO3NOO   = 66, ip_ClONO2   = 67, ip_MACR     = 68, &
    ip_MVK      = 69, ip_GLYOX    = 70, ip_HOCH2CHO = 71, ip_CH4      = 72, &
    ip_O2_b1b2  = 73, ip_O2_b1    = 74, ip_O2_b2    = 75, ip_O3PO1D   = 76, &
    ip_O3Pp     = 77, ip_H2O1D    = 78, ip_N2       = 79, ip_N2_b1    = 80, &
    ip_N2_b2    = 81, ip_N2_b3    = 82, ip_NN2D     = 83, ip_NOp      = 84, &
    ip_Op_em    = 85, ip_O2p_em   = 86, ip_Op_O_em  = 87, ip_N2p_em   = 88, &
    ip_Np_N_em  = 89, ip_Np_N2D_em= 90, ip_N_N2D_em = 91, ip_Op_em_b  = 92, &
    ip_se_O2_b1 = 93, ip_se_O2_b2 = 94, ip_se_N2_b1 = 95, ip_se_N2_b2 = 96, &
    ip_se_N2_b3 = 97, ip_se_N2_b4 = 98, ip_se_Op_em = 99, ip_O2_aurq  =100, &
    ip_N2_aurq  =101, ip_H2SO4    =102, ip_C3O2     =103, ip_CH3NO3   =104, &
    ip_CH3O2NO2 =105, ip_CH3ONO   =106, ip_CH3O2    =107, ip_HCOOH    =108, &
    ! ju_jg_20140124+
    ip_HO2NO2   =109, ip_OHNO3    =110, ip_BrONO2   =111, ip_CH3OCl   =112, &
    ip_MEO2NO2  =113, ip_CHF2Cl   =114, ip_F113    = 115
    ! ju_jg_20140124-

  ! IP_MAX must be set to the highest ip_* value from the definitions above:
  INTEGER, PUBLIC, PARAMETER :: IP_MAX = 115



END MODULE messy_cmn_photol_mem

!*****************************************************************************

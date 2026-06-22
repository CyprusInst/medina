/*************************************************************
 *
 *    kpp_integrate_cuda_prototype.cu
 *    Prototype file for kpp CUDA kernel
 *
 *    Copyright 2016 The Cyprus Institute
 *
 *    Developers: Michail Alvanos - m.alvanos@cyi.ac.cy
 *                Giannis Ashiotis
 *                Theodoros Christoudias - christoudias@cyi.ac.cy
 *
 ********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include "cuda.h"

#define NSPEC 74
#define NVAR 73
#define NFIX 1
#define NREACT 77
#define LU_NONZERO 455
#define NBSIZE 140
#define BLOCKSIZE 64

//#define MAX_VL_GLO 12288 /* elements that will pass in each call */

#define REDUCTION_SIZE_1 64
#define REDUCTION_SIZE_2 32

#define R_gas 8.3144621
#define N_A 6.02214129e+23
#define atm2Pa 101325.0

#define ip_O2 0
#define ip_O3P 1
#define ip_O1D 2
#define ip_H2O2 3
#define ip_NO2 4
#define ip_NO2O 5
#define ip_NOO2 6
#define ip_N2O5 7
#define ip_HNO3 8
#define ip_HNO4 9
#define ip_PAN 10
#define ip_HONO 11
#define ip_CH3OOH 12
#define ip_COH2 13
#define ip_CHOH 14
#define ip_CH3CO3H 15
#define ip_CH3CHO 16
#define ip_CH3COCH3 17
#define ip_MGLYOX 18
#define ip_HOCl 19
#define ip_OClO 20
#define ip_Cl2O2 21
#define ip_ClNO3 22
#define ip_ClNO2 23
#define ip_Cl2 24
#define ip_BrO 25
#define ip_HOBr 26
#define ip_BrCl 27
#define ip_BrNO3 28
#define ip_BrNO2 29
#define ip_Br2 30
#define ip_CCl4 31
#define ip_CH3Cl 32
#define ip_CH3CCl3 33
#define ip_CFCl3 34
#define ip_CF2Cl2 35
#define ip_CH3Br 36
#define ip_CF2ClBr 37
#define ip_CF3Br 38
#define ip_CH3I 39
#define ip_C3H7I 40
#define ip_CH2ClI 41
#define ip_CH2I2 42
#define ip_IO 43
#define ip_HOI 44
#define ip_I2 45
#define ip_ICl 46
#define ip_IBr 47
#define ip_INO2 48
#define ip_INO3 49
#define ip_SO2 50
#define ip_SO3 51
#define ip_OCS 52
#define ip_CS2 53
#define ip_H2O 54
#define ip_N2O 55
#define ip_NO 56
#define ip_CO2 57
#define ip_HCl 58
#define ip_CHCl2Br 59
#define ip_CHClBr2 60
#define ip_CH2ClBr 61
#define ip_CH2Br2 62
#define ip_CHBr3 63
#define ip_SF6 64
#define ip_NO3NOO 65
#define ip_ClONO2 66
#define ip_MACR 67
#define ip_MVK 68
#define ip_GLYOX 69
#define ip_HOCH2CHO 70
#define ip_CH4 71
#define ip_O2_b1b2 72
#define ip_O2_b1 73
#define ip_O2_b2 74
#define ip_O3PO1D 75
#define ip_O3Pp 76
#define ip_H2O1D 77
#define ip_N2 78
#define ip_N2_b1 79
#define ip_N2_b2 80
#define ip_N2_b3 81
#define ip_NN2D 82
#define ip_NOp 83
#define ip_Op_em 84
#define ip_O2p_em 85
#define ip_Op_O_em 86
#define ip_N2p_em 87
#define ip_Np_N_em 88
#define ip_Np_N2D_em 89
#define ip_N_N2D_em 90
#define ip_Op_em_b 91
#define ip_se_O2_b1 92
#define ip_se_O2_b2 93
#define ip_se_N2_b1 94
#define ip_se_N2_b2 95
#define ip_se_N2_b3 96
#define ip_se_N2_b4 97
#define ip_se_Op_em 98
#define ip_O2_aurq 99
#define ip_N2_aurq 100
#define ip_H2SO4 101
#define ip_C3O2 102
#define ip_CH3NO3 103
#define ip_CH3O2NO2 104
#define ip_CH3ONO 105
#define ip_CH3O2 106
#define ip_HCOOH 107
#define ip_HO2NO2 108
#define ip_OHNO3 109
#define ip_BrONO2 110
#define ip_CH3OCl 111
#define ip_MEO2NO2 112
#define ip_CHF2Cl 113
#define ip_F113 114

#define ind_N2O5 0
#define ind_HNO3 1
#define ind_CO 2
#define ind_CH3CO3 3
#define ind_HOCl 4
#define ind_OClO 5
#define ind_HOBr 6
#define ind_H2SO4 7
#define ind_CH3SO3H 8
#define ind_NO3m_cs 9
#define ind_Hp_cs 10
#define ind_LHOC3H6O2 11
#define ind_LossG2104 12
#define ind_LossG2107 13
#define ind_LossG2110 14
#define ind_LossHO2 15
#define ind_LossO1D 16
#define ind_LossO3 17
#define ind_LossO3Br 18
#define ind_LossO3Cl 19
#define ind_LossO3H 20
#define ind_LossO3Nn 21
#define ind_LossO3O 22
#define ind_LossOH 23
#define ind_ProdLBr 24
#define ind_ProdLCl 25
#define ind_ProdO3 26
#define ind_ProdSBr 27
#define ind_O3P 28
#define ind_H2OH2O 29
#define ind_Cl2O2 30
#define ind_H 31
#define ind_CH3Br 32
#define ind_CHCl2Br 33
#define ind_CHClBr2 34
#define ind_CH2ClBr 35
#define ind_CH2Br2 36
#define ind_CHBr3 37
#define ind_CH3SO3 38
#define ind_CH3CCl3 39
#define ind_H2 40
#define ind_O1D 41
#define ind_CH4 42
#define ind_C2H6 43
#define ind_SO2 44
#define ind_C2H4 45
#define ind_CH3CHO 46
#define ind_H2O 47
#define ind_NO 48
#define ind_C2H2 49
#define ind_DMSO 50
#define ind_CH3OOH 51
#define ind_H2O2 52
#define ind_NO3 53
#define ind_Cl2 54
#define ind_NO2 55
#define ind_BrCl 56
#define ind_HBr 57
#define ind_HCl 58
#define ind_Br2 59
#define ind_CH3SO2 60
#define ind_BrNO3 61
#define ind_HCHO 62
#define ind_ClNO3 63
#define ind_DMS 64
#define ind_OH 65
#define ind_O3 66
#define ind_Cl 67
#define ind_Br 68
#define ind_ClO 69
#define ind_CH3O2 70
#define ind_HO2 71
#define ind_BrO 72
#define ind_O2 73
#define ind_N -1
#define ind_N2D -1
#define ind_N2 -1
#define ind_NH3 -1
#define ind_N2O -1
#define ind_HONO -1
#define ind_HOONO -1
#define ind_HNO4 -1
#define ind_NH2 -1
#define ind_HNO -1
#define ind_NHOH -1
#define ind_NH2O -1
#define ind_NH2OH -1
#define ind_CO2 -1
#define ind_HCOOH -1
#define ind_CH2OO -1
#define ind_CH3 -1
#define ind_CH3O -1
#define ind_HOCH2O2 -1
#define ind_CH3OH -1
#define ind_HOCH2OOH -1
#define ind_HOCH2OH -1
#define ind_CH3ONO -1
#define ind_CH3NO3 -1
#define ind_CH3O2NO2 -1
#define ind_HOCH2O2NO2 -1
#define ind_LCARBON -1
#define ind_HCOCO3 -1
#define ind_HCOCO3A -1
#define ind_GLYOX -1
#define ind_HCOCO2H -1
#define ind_CHOOCHO -1
#define ind_HCOCO3H -1
#define ind_HCOCH2O2 -1
#define ind_HOCH2CO3 -1
#define ind_HOOCH2CO3 -1
#define ind_CH3CO2H -1
#define ind_HOCH2CHO -1
#define ind_CH3CO3H -1
#define ind_HOCH2CO2H -1
#define ind_HOCH2CO3H -1
#define ind_C2H5O2 -1
#define ind_HOCH2CH2O -1
#define ind_HOCH2CH2O2 -1
#define ind_C2H5OOH -1
#define ind_ETHGLY -1
#define ind_HYETHO2H -1
#define ind_PAN -1
#define ind_PHAN -1
#define ind_ETHOHNO3 -1
#define ind_C33CO -1
#define ind_CHOCOCH2O2 -1
#define ind_HCOCH2CO3 -1
#define ind_ALCOCH2OOH -1
#define ind_MGLYOX -1
#define ind_HOCH2COCHO -1
#define ind_HCOCH2CHO -1
#define ind_HOCH2COCO2H -1
#define ind_HCOCH2CO2H -1
#define ind_HCOCH2CO3H -1
#define ind_CH3COCH2O2 -1
#define ind_HOC2H4CO3 -1
#define ind_C3H6 -1
#define ind_CH3COCH3 -1
#define ind_ACETOL -1
#define ind_HYPERACET -1
#define ind_HOC2H4CO2H -1
#define ind_HOC2H4CO3H -1
#define ind_IC3H7O2 -1
#define ind_HYPROPO2 -1
#define ind_C3H8 -1
#define ind_IC3H7OOH -1
#define ind_HYPROPO2H -1
#define ind_C3PAN2 -1
#define ind_NOA -1
#define ind_C3PAN1 -1
#define ind_PRONO3BO2 -1
#define ind_IC3H7NO3 -1
#define ind_PR2O2HNO3 -1
#define ind_C312COCO3 -1
#define ind_C4CODIAL -1
#define ind_CO23C3CHO -1
#define ind_C312COCO3H -1
#define ind_CO2H3CHO -1
#define ind_MACO3 -1
#define ind_BIACETO2 -1
#define ind_CHOC3COO2 -1
#define ind_C44O2 -1
#define ind_CO2H3CO3 -1
#define ind_MACR -1
#define ind_MVK -1
#define ind_BIACET -1
#define ind_MACO2H -1
#define ind_MVKOH -1
#define ind_MACO3H -1
#define ind_BIACETOH -1
#define ind_C413COOOH -1
#define ind_BIACETOOH -1
#define ind_C44OOH -1
#define ind_CO2H3CO3H -1
#define ind_MACRO2 -1
#define ind_MEK -1
#define ind_HO12CO3C4 -1
#define ind_MACROH -1
#define ind_MACROOH -1
#define ind_NC4H10 -1
#define ind_C312COPAN -1
#define ind_MPAN -1
#define ind_LMEKO2 -1
#define ind_LHMVKABO2 -1
#define ind_LMVKOHABO2 -1
#define ind_LMEKOOH -1
#define ind_LHMVKABOOH -1
#define ind_LMVKOHABOOH -1
#define ind_LC4H9O2 -1
#define ind_LC4H9OOH -1
#define ind_LC4H9NO3 -1
#define ind_CHOC3COCO3 -1
#define ind_CO23C4CO3 -1
#define ind_CO13C4CHO -1
#define ind_CO23C4CHO -1
#define ind_C513CO -1
#define ind_CHOC3COOOH -1
#define ind_CO23C4CO3H -1
#define ind_C511O2 -1
#define ind_C512O2 -1
#define ind_C513O2 -1
#define ind_C5H8 -1
#define ind_HCOC5 -1
#define ind_C511OOH -1
#define ind_C512OOH -1
#define ind_C513OOH -1
#define ind_ISOPBO2 -1
#define ind_ISOPDO2 -1
#define ind_C59O2 -1
#define ind_ISOPAOH -1
#define ind_ISOPBOH -1
#define ind_ISOPDOH -1
#define ind_ISOPBOOH -1
#define ind_ISOPDOOH -1
#define ind_C59OOH -1
#define ind_C514OOH -1
#define ind_C514O2 -1
#define ind_CHOC3COPAN -1
#define ind_C5PAN9 -1
#define ind_NC4CHO -1
#define ind_NISOPO2 -1
#define ind_ISOPBNO3 -1
#define ind_ISOPDNO3 -1
#define ind_NISOPOOH -1
#define ind_C514NO3 -1
#define ind_LHC4ACCO3 -1
#define ind_LHC4ACCHO -1
#define ind_LHC4ACCO2H -1
#define ind_LHC4ACCO3H -1
#define ind_LISOPACO2 -1
#define ind_LC578O2 -1
#define ind_LISOPACOOH -1
#define ind_LC578OOH -1
#define ind_LC5PAN1719 -1
#define ind_LISOPACNO3 -1
#define ind_LNISO3 -1
#define ind_LNISOOH -1
#define ind_CO235C5CHO -1
#define ind_CO235C6O2 -1
#define ind_C614CO -1
#define ind_CO235C6OOH -1
#define ind_C614O2 -1
#define ind_C614OOH -1
#define ind_C614NO3 -1
#define ind_CO235C6CO3 -1
#define ind_CO235C6CHO -1
#define ind_C235C6CO3H -1
#define ind_C716O2 -1
#define ind_C716OOH -1
#define ind_ROO6R4P -1
#define ind_ROO6R5P -1
#define ind_ROO6R3O -1
#define ind_C721O2 -1
#define ind_C722O2 -1
#define ind_ROO6R3O2 -1
#define ind_ROO6R5O2 -1
#define ind_C721OOH -1
#define ind_C722OOH -1
#define ind_ROO6R3OOH -1
#define ind_C7PAN3 -1
#define ind_ROO6R3NO3 -1
#define ind_C8BCO2 -1
#define ind_C721CO3 -1
#define ind_C8BCCO -1
#define ind_C8BCOOH -1
#define ind_C721CHO -1
#define ind_NORPINIC -1
#define ind_C721CO3H -1
#define ind_C85O2 -1
#define ind_C89O2 -1
#define ind_C811O2 -1
#define ind_C86O2 -1
#define ind_C812O2 -1
#define ind_C813O2 -1
#define ind_C8BC -1
#define ind_C85OOH -1
#define ind_C811OOH -1
#define ind_C86OOH -1
#define ind_C812OOH -1
#define ind_C813OOH -1
#define ind_C89OOH -1
#define ind_C810OOH -1
#define ind_C810O2 -1
#define ind_C8BCNO3 -1
#define ind_C721PAN -1
#define ind_C89NO3 -1
#define ind_C810NO3 -1
#define ind_C85CO3 -1
#define ind_NOPINDCO -1
#define ind_C85CO3H -1
#define ind_NOPINDO2 -1
#define ind_C89CO3 -1
#define ind_C811CO3 -1
#define ind_NOPINONE -1
#define ind_NOPINOO -1
#define ind_NORPINAL -1
#define ind_C89CO2H -1
#define ind_NOPINDOOH -1
#define ind_RO6R3P -1
#define ind_C89CO3H -1
#define ind_PINIC -1
#define ind_C811CO3H -1
#define ind_C96O2 -1
#define ind_C97O2 -1
#define ind_C98O2 -1
#define ind_C96OOH -1
#define ind_C97OOH -1
#define ind_C98OOH -1
#define ind_C89PAN -1
#define ind_C9PAN2 -1
#define ind_C811PAN -1
#define ind_C96NO3 -1
#define ind_C98NO3 -1
#define ind_C109CO -1
#define ind_PINALO2 -1
#define ind_PINALOOH -1
#define ind_C109O2 -1
#define ind_C96CO3 -1
#define ind_C106O2 -1
#define ind_APINENE -1
#define ind_BPINENE -1
#define ind_PINAL -1
#define ind_APINAOO -1
#define ind_APINBOO -1
#define ind_MENTHEN6ONE -1
#define ind_PINONIC -1
#define ind_C109OOH -1
#define ind_PERPINONIC -1
#define ind_C106OOH -1
#define ind_BPINAO2 -1
#define ind_OH2MENTHEN6ONE -1
#define ind_RO6R1O2 -1
#define ind_ROO6R1O2 -1
#define ind_RO6R3O2 -1
#define ind_OHMENTHEN6ONEO2 -1
#define ind_BPINAOOH -1
#define ind_RO6R1OOH -1
#define ind_RO6R3OOH -1
#define ind_ROO6R1OOH -1
#define ind_PINALNO3 -1
#define ind_C10PAN2 -1
#define ind_C106NO3 -1
#define ind_BPINANO3 -1
#define ind_RO6R1NO3 -1
#define ind_RO6R3NO3 -1
#define ind_ROO6R1NO3 -1
#define ind_LAPINABO2 -1
#define ind_LAPINABOOH -1
#define ind_LNAPINABO2 -1
#define ind_LNBPINABO2 -1
#define ind_LAPINABNO3 -1
#define ind_LNAPINABOOH -1
#define ind_LNBPINABOOH -1
#define ind_ClNO2 -1
#define ind_CCl4 -1
#define ind_CH3Cl -1
#define ind_CF2Cl2 -1
#define ind_CFCl3 -1
#define ind_BrNO2 -1
#define ind_CF3Br -1
#define ind_CF2ClBr -1
#define ind_I -1
#define ind_I2 -1
#define ind_IO -1
#define ind_OIO -1
#define ind_I2O2 -1
#define ind_HI -1
#define ind_HOI -1
#define ind_HIO3 -1
#define ind_INO2 -1
#define ind_INO3 -1
#define ind_CH3I -1
#define ind_CH2I2 -1
#define ind_C3H7I -1
#define ind_ICl -1
#define ind_CH2ClI -1
#define ind_IBr -1
#define ind_S -1
#define ind_SO -1
#define ind_SO3 -1
#define ind_SH -1
#define ind_OCS -1
#define ind_SF6 -1
#define ind_Hg -1
#define ind_HgO -1
#define ind_HgCl -1
#define ind_HgCl2 -1
#define ind_HgBr -1
#define ind_HgBr2 -1
#define ind_ClHgBr -1
#define ind_BrHgOBr -1
#define ind_ClHgOBr -1
#define ind_RGM_cs -1
#define ind_IPART -1
#define ind_Dummy -1
#define ind_O3s -1
#define ind_LO3s -1
#define ind_LHOC3H6OOH -1
#define ind_ISO2 -1
#define ind_ISON -1
#define ind_ISOOH -1
#define ind_MVKO2 -1
#define ind_MVKOOH -1
#define ind_NACA -1
#define ind_Op -1
#define ind_O2p -1
#define ind_Np -1
#define ind_N2p -1
#define ind_NOp -1
#define ind_Hp -1
#define ind_em -1
#define ind_kJmol -1
#define ind_RH2O -1
#define ind_RNOy -1
#define ind_RCly -1
#define ind_RBr -1
#define ind_CFCl3_c -1
#define ind_CF2Cl2_c -1
#define ind_N2O_c -1
#define ind_CH3CCl3_c -1
#define ind_CF2ClBr_c -1
#define ind_CF3Br_c -1
#define ind_LTERP -1
#define ind_LALK4 -1
#define ind_LALK5 -1
#define ind_LARO1 -1
#define ind_LARO2 -1
#define ind_LOLE1 -1
#define ind_LOLE2 -1
#define ind_LfPOG01 -1
#define ind_LfPOG02 -1
#define ind_LfPOG03 -1
#define ind_LfPOG04 -1
#define ind_LbbPOG01 -1
#define ind_LbbPOG02 -1
#define ind_LbbPOG03 -1
#define ind_LbbPOG04 -1
#define ind_LfSOGsv01 -1
#define ind_LbbSOGsv01 -1
#define ind_LfSOGiv01 -1
#define ind_LfSOGiv02 -1
#define ind_LfSOGiv03 -1
#define ind_LbbSOGiv01 -1
#define ind_LbbSOGiv02 -1
#define ind_LbbSOGiv03 -1
#define ind_LbSOGv01 -1
#define ind_LbSOGv02 -1
#define ind_LbSOGv03 -1
#define ind_LbSOGv04 -1
#define ind_LbOSOGv01 -1
#define ind_LbOSOGv02 -1
#define ind_LbOSOGv03 -1
#define ind_LaSOGv01 -1
#define ind_LaSOGv02 -1
#define ind_LaSOGv03 -1
#define ind_LaSOGv04 -1
#define ind_LaOSOGv01 -1
#define ind_LaOSOGv02 -1
#define ind_LaOSOGv03 -1
#define ind_O2_a01 -1
#define ind_O3_a01 -1
#define ind_OH_a01 -1
#define ind_HO2_a01 -1
#define ind_H2O_a01 -1
#define ind_H2O2_a01 -1
#define ind_NH3_a01 -1
#define ind_NO_a01 -1
#define ind_NO2_a01 -1
#define ind_NO3_a01 -1
#define ind_HONO_a01 -1
#define ind_HNO3_a01 -1
#define ind_HNO4_a01 -1
#define ind_N2O5_a01 -1
#define ind_CH3OH_a01 -1
#define ind_HCOOH_a01 -1
#define ind_HCHO_a01 -1
#define ind_CH3O2_a01 -1
#define ind_CH3OOH_a01 -1
#define ind_CO2_a01 -1
#define ind_CH3CO2H_a01 -1
#define ind_PAN_a01 -1
#define ind_C2H5O2_a01 -1
#define ind_CH3CHO_a01 -1
#define ind_CH3COCH3_a01 -1
#define ind_Cl_a01 -1
#define ind_Cl2_a01 -1
#define ind_HCl_a01 -1
#define ind_HOCl_a01 -1
#define ind_Br_a01 -1
#define ind_Br2_a01 -1
#define ind_HBr_a01 -1
#define ind_HOBr_a01 -1
#define ind_BrCl_a01 -1
#define ind_I2_a01 -1
#define ind_IO_a01 -1
#define ind_HI_a01 -1
#define ind_HOI_a01 -1
#define ind_ICl_a01 -1
#define ind_IBr_a01 -1
#define ind_HIO3_a01 -1
#define ind_SO2_a01 -1
#define ind_H2SO4_a01 -1
#define ind_DMS_a01 -1
#define ind_DMSO_a01 -1
#define ind_Hg_a01 -1
#define ind_HgO_a01 -1
#define ind_HgOH_a01 -1
#define ind_HgOHOH_a01 -1
#define ind_HgOHCl_a01 -1
#define ind_HgCl2_a01 -1
#define ind_HgBr2_a01 -1
#define ind_HgSO3_a01 -1
#define ind_ClHgBr_a01 -1
#define ind_BrHgOBr_a01 -1
#define ind_ClHgOBr_a01 -1
#define ind_O2m_a01 -1
#define ind_OHm_a01 -1
#define ind_Hp_a01 -1
#define ind_NH4p_a01 -1
#define ind_NO2m_a01 -1
#define ind_NO3m_a01 -1
#define ind_NO4m_a01 -1
#define ind_CO3m_a01 -1
#define ind_HCOOm_a01 -1
#define ind_HCO3m_a01 -1
#define ind_CH3COOm_a01 -1
#define ind_Clm_a01 -1
#define ind_Cl2m_a01 -1
#define ind_ClOm_a01 -1
#define ind_ClOHm_a01 -1
#define ind_Brm_a01 -1
#define ind_Br2m_a01 -1
#define ind_BrOm_a01 -1
#define ind_BrOHm_a01 -1
#define ind_BrCl2m_a01 -1
#define ind_Br2Clm_a01 -1
#define ind_Im_a01 -1
#define ind_IO2m_a01 -1
#define ind_IO3m_a01 -1
#define ind_ICl2m_a01 -1
#define ind_IClBrm_a01 -1
#define ind_IBr2m_a01 -1
#define ind_SO3m_a01 -1
#define ind_SO3mm_a01 -1
#define ind_SO4m_a01 -1
#define ind_SO4mm_a01 -1
#define ind_SO5m_a01 -1
#define ind_HSO3m_a01 -1
#define ind_HSO4m_a01 -1
#define ind_HSO5m_a01 -1
#define ind_CH3SO3m_a01 -1
#define ind_CH2OHSO3m_a01 -1
#define ind_Hgp_a01 -1
#define ind_Hgpp_a01 -1
#define ind_HgOHp_a01 -1
#define ind_HgClp_a01 -1
#define ind_HgCl3m_a01 -1
#define ind_HgCl4mm_a01 -1
#define ind_HgBrp_a01 -1
#define ind_HgBr3m_a01 -1
#define ind_HgBr4mm_a01 -1
#define ind_HgSO32mm_a01 -1
#define ind_D1O_a01 -1
#define ind_D2O_a01 -1
#define ind_DAHp_a01 -1
#define ind_DA_a01 -1
#define ind_DAm_a01 -1
#define ind_DGtAi_a01 -1
#define ind_DGtAs_a01 -1
#define ind_PROD1_a01 -1
#define ind_PROD2_a01 -1
#define ind_Nap_a01 -1
#define ind_LossG2106 -1
#define ind_LossG3103 -1
#define ind_LossG3105 -1
#define ind_LossG3106 -1
#define ind_LossG3201 -1
#define ind_LossG3201a -1
#define ind_LossG3201aKG -1
#define ind_LossG3201b -1
#define ind_LossG3201bKG -1
#define ind_LossG3202 -1
#define ind_LossG4110 -1
#define ind_LossJ3101 -1
#define ind_LossJ3103a -1
#define ind_LossO3Cln -1
#define ind_LossO3Hn -1
#define ind_LossO3N -1
#define ind_LossO3N2 -1
#define ind_LossO3R -1
#define ind_ProdHO2 -1
#define ind_ProdMeO2 -1
#define ind_ProdRO2 -1
#define ind_ProdSCl -1

#define ihs_N2O5_H2O 0
#define ihs_HOCl_HCl 1
#define ihs_ClNO3_HCl 2
#define ihs_ClNO3_H2O 3
#define ihs_N2O5_HCl 4
#define ihs_ClNO3_HBr 5
#define ihs_BrNO3_HCl 6
#define ihs_HOCl_HBr 7
#define ihs_HOBr_HCl 8
#define ihs_HOBr_HBr 9
#define ihs_BrNO3_H2O 10
#define ihs_Hg 11
#define ihs_RGM 12

#define iht_N2O5 0
#define iht_HNO3 1
#define iht_Hg 2
#define iht_RGM 3

#define k_s   (8.42E-13 )
#define k_t   (1.75E-12)
#define k_p   (1.24E-13 )
#define k_rohro   (1.6E-13 )
#define k_adp   (0.45E-11 )
#define k_ads   (3.0E-11)
#define k_adt   (5.5E-11 )
#define k_adsecprim   (3.0E-11 )
#define k_adtertprim   (5.7E-11 )
#define f_soh   (3.44   )
#define f_toh   (2.68 )
#define f_sooh   (7. )
#define f_tooh   (7. )
#define f_ono2   (0.04  )
#define f_ch2ono2   (0.2 )
#define f_cpan  (.25 )
#define f_allyl   (3.6 )
#define f_alk  (1.23 )
#define f_cho   (0.55 )
#define f_co2h   (1.67 )
#define f_co   (0.73 )
#define f_o   (8.15 )
#define f_pch2oh   (1.29)
#define f_tch2oh   (0.53)
#define a_pan   (0.56       )
#define a_cho   (0.31    )
#define a_coch3   (0.76 )
#define a_ch2ono2   (0.47   )
#define a_ch2oh   (1.7   )
#define a_ch2ooh   (0.21 )
#define a_coh   (2.2        )
#define a_cooh   (2.2    )
#define a_co2h   (0.25)

#define ifun 0
#define ijac 1
#define istp 2
#define iacc 3
#define irej 4
#define idec 5
#define isol 6
#define isng 7
#define itexit 0
#define ihexit 1

#define ZERO 0.0
#define ONE 1.0
#define HALF 0.5


/*
 * Fortran to C macros 
 * GPU-friendly array deffinition 
 * i:VL_GLO, j:NVAR 
 *
 */
#define conc(i,j)    conc[(j)*VL_GLO+(i)]
#define khet_st(i,j) khet_st[(j)*VL_GLO+(i)]
#define khet_tr(i,j) khet_tr[(j)*VL_GLO+(i)]
#define jx(i,j)      jx[j*VL_GLO+i]
#define istatus(i,j) istatus[(j)*(VL_GLO)+(i)]
#define rstatus(i,j) rstatus[(j)*(VL_GLO)+(i)]


#define ROUND128(X)  (X + (128 - 1)) & ~(128 - 1)

#define rconst(i,j)  rconst[(j)]


/* Temporary arrays allocated in stack */
#define var(i,j)     var[(j)]
#define fix(i,j)     fix[(j)]
#define jcb(i,j)     jcb[(j)]
#define varDot(i,j)  varDot[j]
#define varNew(i,j) varNew[(j)]
#define Fcn0(i,j)   Fcn0[(j)]
#define Fcn(i,j)    Fcn[(j)]
#define Fcn(i,j)    Fcn[(j)]
#define dFdT(i,j)   dFdT[(j)]
#define varErr(i,j) varErr[(j)]
#define K(i,j,k) K[(j)*(NVAR)+(k)]
#define jac0(i,j)    jac0[(j)]
#define Ghimj(i,j)   Ghimj[(j)]


/* Enable debug flags for GPU */
//#define DEBUG

#ifdef DEBUG
#define GPU_DEBUG()\
    gpuErrchk( cudaPeekAtLastError()   ); \
    gpuErrchk( cudaDeviceSynchronize() ); 

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

static inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#else 
/* If debug flags are disabled */
#define GPU_DEBUG()
#define gpuErrchk(ans) ans
#endif

/** prefetches into L1 cache */
__device__ inline void prefetch_gl1(const void *p) {
#if __CUDA_ARCH__ <= 300
        asm("prefetch.global.L1 [%0];": :"l"(p));
#endif
}
__device__ inline void prefetch_ll1(const void *p) {
#if __CUDA_ARCH__ <= 300
        asm("prefetch.local.L1 [%0];": :"l"(p));
#endif
}

/** prefetches into L2 cache */
__device__ inline void prefetch_gl2(const void *p) {
#if __CUDA_ARCH__ <= 300
        asm("prefetch.global.L2 [%0];": :"l"(p));
#endif
}
__device__ inline void prefetch_ll2(const void *p) {
#if __CUDA_ARCH__ <= 300
        asm("prefetch.local.L2 [%0];": :"l"(p));
#endif
}

#if defined(__SINGLEPREC)
typedef float REAL;
#else
typedef double REAL;
#endif


// __device__ void  update_rconst(const REAL * __restrict__ var,
//                                const REAL * __restrict__ khet_st, const REAL * __restrict__ khet_tr,
//                                const REAL * __restrict__ jx, REAL * __restrict__ rconst,
//                                const REAL * __restrict__ temp_gpu,
//                                const REAL * __restrict__ press_gpu,
//                                const REAL * __restrict__ cair_gpu,
//                                const int VL_GLO);

/* This runs on CPU */
REAL machine_eps_flt()
{
    REAL machEps = 1.0f;

    do
    {
        machEps /= 2.0f;
        // If next epsilon yields 1, then break, because current
        // epsilon is the machine epsilon.
    }
    while ((REAL)(1.0 + (machEps/2.0)) != 1.0);

    return machEps;
}

/* This runs on GPU */
__device__ REAL machine_eps_flt_cuda()
{
    typedef union
    {
        long  i64;
        REAL f64;
    } flt_64;

    flt_64 s;

    s.f64 = 1.;
    s.i64++;
    return (s.f64 - 1.);
}

__device__  static REAL alpha_AN(const int n, const int ro2type, const REAL temp, const REAL cair){
    REAL alpha=2.E-22, beta=1.0, Yinf_298K=0.43,  F=0.41, m0=0., minf=8.0;
    REAL Y0_298K, Y0_298K_tp, Yinf_298K_t, zeta, k_ratio, alpha_a;
    /*  IF (ro2type = 1) THEN   m = 0.4                !   primary RO2
        ELSE IF (ro2type = 2) THEN  m = 1.                 !   secondary RO2
        ELSE IF (ro2type = 3) THEN  m = 0.3                !   tertiary RO2
        ELSE  m = 1.
  */
    REAL m = 1.;
    Y0_298K     = alpha*exp(beta*n);
    Y0_298K_tp  = Y0_298K *cair *pow((temp/298.),(- m0));
    Yinf_298K_t = Yinf_298K * pow((temp/298.),(- minf));
    zeta        = 1/(1+ pow(log10(Y0_298K_tp/Yinf_298K_t),2));
    k_ratio     = (Y0_298K_tp/(1+ Y0_298K_tp/Yinf_298K_t))*pow(F,zeta);
    alpha_a    = k_ratio/(1+ k_ratio) *m;
    return alpha_a;
}
__device__  static REAL alpha_AN(const int n, const int ro2type, const int bcarb, const int gcarb, const int abic, const REAL temp, const REAL cair){
    REAL alpha=2.E-22, beta=1.0, Yinf_298K=0.43,  F=0.41, m0=0., minf=8.0;
    REAL Y0_298K, Y0_298K_tp, Yinf_298K_t, zeta, k_ratio, alpha_a;
    REAL bcf=1., gcf=1., abf=1.;
    REAL m = 1.; //According to Teng, ref3189

if (bcarb == 1) { bcf = 0.19; }// derived from Praske, ref3190: alpha_AN = 0.03 for the secondary HMKO2 relative to alpha_AN for 6C RO2 (0.16)
if (gcarb == 1) {gcf = 0.44; }// derived from Praske, ref3190: alpha_AN = 0.07 for the primary HMKO2 relative to alpha_AN for 6C RO2 (0.16)
if (abic == 1) { abf = 0.24; }// derived from the ratio of AN- yield for toluene from Elrod et al. (ref3180), 5.5 0x1.9206e69676542p+ 229t & 
                              // 200 torr, and this SAR for linear alkyl RO2 with 9 heavy atoms, 23.3%

    Y0_298K     = alpha*exp(beta*n);
    Y0_298K_tp  = Y0_298K *cair *pow((temp/298.),(- m0));
    Yinf_298K_t = Yinf_298K * pow((temp/298.),(- minf));
    zeta        = 1/(1+ pow(log10(Y0_298K_tp/Yinf_298K_t),2));
    k_ratio     = (Y0_298K_tp/(1+ Y0_298K_tp/Yinf_298K_t))*pow(F,zeta);
    alpha_a    = k_ratio/(1+ k_ratio) *m*bcf*gcf*abf;
    return alpha_a;
}

__device__  static REAL k_RO2_HO2(const REAL temp, const int nC){
    return 2.91e-13*exp(1300./temp)*(1.-exp(-0.245*nC)); // ref1630
}
__device__ REAL ros_ErrorNorm(REAL * __restrict__ var, REAL * __restrict__ varNew, REAL * __restrict__ varErr,
                                const REAL * __restrict__ absTol, const REAL * __restrict__ relTol,
                                const int vectorTol )
{
    REAL err, scale, varMax;


    err = ZERO;

    if (vectorTol){
        for (int i=0;i<NVAR - 16;i+=16){
            prefetch_ll1(&varErr[i]);
            prefetch_ll1(&absTol[i]);
            prefetch_ll1(&relTol[i]);
            prefetch_ll1(&var[i]);
            prefetch_ll1(&varNew[i]);
        }

        for (int i=0; i<NVAR; i++)
        {
            varMax = fmax(fabs(var[i]),fabs(varNew[i]));
            scale = absTol[i]+ relTol[i]*varMax;

            err += pow((REAL)varErr[i]/scale,2.0);
        }
        err  = sqrt((REAL) err/NVAR);
    }else{
        for (int i=0;i<NVAR - 16;i+=16){
            prefetch_ll1(&varErr[i]);
            prefetch_ll1(&var[i]);
            prefetch_ll1(&varNew[i]);
        }

        for (int i=0; i<NVAR; i++)
        {
            varMax = fmax(fabs(var[i]),fabs(varNew[i]));

            scale = absTol[0]+ relTol[0]*varMax;
            err += pow((REAL)varErr[i]/scale,2.0);
        }
        err  = sqrt((REAL) err/NVAR);
    }

    return err;

}

__device__ void kppSolve(const REAL * __restrict__ Ghimj, REAL * __restrict__ K, 
                         const int istage, const int ros_S ){
    int index = blockIdx.x*blockDim.x+threadIdx.x;

       K = &K[istage*NVAR];

        K[9] = K[9]- Ghimj[32]*K[0];
        K[10] = K[10]- Ghimj[34]*K[0];
        K[17] = K[17]- Ghimj[54]*K[0];
        K[47] = K[47]- Ghimj[155]*K[29]- Ghimj[156]*K[32]- Ghimj[157]*K[33]- Ghimj[158]*K[34]- Ghimj[159]*K[35]- Ghimj[160]*K[36]- Ghimj[161]*K[37]  - Ghimj[162]*K[39]- Ghimj[163]*K[40]- Ghimj[164]*K[41];
        K[55] = K[55]- Ghimj[200]*K[48];
        K[56] = K[56]- Ghimj[206]*K[54];
        K[57] = K[57]- Ghimj[214]*K[45]- Ghimj[215]*K[46]- Ghimj[216]*K[49]- Ghimj[217]*K[51];
        K[58] = K[58]- Ghimj[225]*K[42]- Ghimj[226]*K[43]- Ghimj[227]*K[45]- Ghimj[228]*K[46]- Ghimj[229]*K[49]- Ghimj[230]*K[51]- Ghimj[231]*K[52];
        K[59] = K[59]- Ghimj[239]*K[56];
        K[60] = K[60]- Ghimj[248]*K[53];
        K[61] = K[61]- Ghimj[257]*K[55];
        K[62] = K[62]- Ghimj[263]*K[50]- Ghimj[264]*K[51]- Ghimj[265]*K[53]- Ghimj[266]*K[61];
        K[63] = K[63]- Ghimj[276]*K[55]- Ghimj[277]*K[61];
        K[64] = K[64]- Ghimj[283]*K[53]- Ghimj[284]*K[61]- Ghimj[285]*K[63];
        K[65] = K[65]- Ghimj[292]*K[32]- Ghimj[293]*K[33]- Ghimj[294]*K[34]- Ghimj[295]*K[35]- Ghimj[296]*K[36]- Ghimj[297]*K[37]- Ghimj[298]*K[39]  - Ghimj[299]*K[40]- Ghimj[300]*K[41]- Ghimj[301]*K[44]- Ghimj[302]*K[47]- Ghimj[303]*K[50]- Ghimj[304]*K[51]- Ghimj[305]*K[52] - Ghimj[306]*K[57]- Ghimj[307]*K[58]- Ghimj[308]*K[59]- Ghimj[309]*K[60]- Ghimj[310]*K[61]- Ghimj[311]*K[62]- Ghimj[312]*K[63] - Ghimj[313]*K[64];
        K[66] = K[66]- Ghimj[322]*K[28]- Ghimj[323]*K[41]- Ghimj[324]*K[47]- Ghimj[325]*K[52]- Ghimj[326]*K[57]- Ghimj[327]*K[58]- Ghimj[328]*K[60]  - Ghimj[329]*K[61]- Ghimj[330]*K[62]- Ghimj[331]*K[63]- Ghimj[332]*K[64]- Ghimj[333]*K[65];
        K[67] = K[67]- Ghimj[341]*K[39]- Ghimj[342]*K[42]- Ghimj[343]*K[43]- Ghimj[344]*K[45]- Ghimj[345]*K[46]- Ghimj[346]*K[48]- Ghimj[347]*K[49]  - Ghimj[348]*K[51]- Ghimj[349]*K[52]- Ghimj[350]*K[54]- Ghimj[351]*K[56]- Ghimj[352]*K[58]- Ghimj[353]*K[59]- Ghimj[354]*K[61] - Ghimj[355]*K[62]- Ghimj[356]*K[63]- Ghimj[357]*K[64]- Ghimj[358]*K[65]- Ghimj[359]*K[66];
        K[68] = K[68]- Ghimj[366]*K[32]- Ghimj[367]*K[33]- Ghimj[368]*K[34]- Ghimj[369]*K[35]- Ghimj[370]*K[36]- Ghimj[371]*K[37]- Ghimj[372]*K[45]  - Ghimj[373]*K[46]- Ghimj[374]*K[48]- Ghimj[375]*K[49]- Ghimj[376]*K[51]- Ghimj[377]*K[54]- Ghimj[378]*K[56]- Ghimj[379]*K[57] - Ghimj[380]*K[59]- Ghimj[381]*K[61]- Ghimj[382]*K[62]- Ghimj[383]*K[63]- Ghimj[384]*K[64]- Ghimj[385]*K[65]- Ghimj[386]*K[66] - Ghimj[387]*K[67];
        K[69] = K[69]- Ghimj[393]*K[30]- Ghimj[394]*K[48]- Ghimj[395]*K[55]- Ghimj[396]*K[61]- Ghimj[397]*K[63]- Ghimj[398]*K[66]- Ghimj[399]*K[67]  - Ghimj[400]*K[68];
        K[70] = K[70]- Ghimj[405]*K[42]- Ghimj[406]*K[43]- Ghimj[407]*K[49]- Ghimj[408]*K[50]- Ghimj[409]*K[51]- Ghimj[410]*K[60]- Ghimj[411]*K[61]  - Ghimj[412]*K[63]- Ghimj[413]*K[64]- Ghimj[414]*K[65]- Ghimj[415]*K[66]- Ghimj[416]*K[67]- Ghimj[417]*K[68]- Ghimj[418]*K[69];
        K[71] = K[71]- Ghimj[422]*K[31]- Ghimj[423]*K[38]- Ghimj[424]*K[40]- Ghimj[425]*K[44]- Ghimj[426]*K[50]- Ghimj[427]*K[52]- Ghimj[428]*K[60]  - Ghimj[429]*K[61]- Ghimj[430]*K[62]- Ghimj[431]*K[63]- Ghimj[432]*K[64]- Ghimj[433]*K[65]- Ghimj[434]*K[66]- Ghimj[435]*K[67] - Ghimj[436]*K[68]- Ghimj[437]*K[69]- Ghimj[438]*K[70];
        K[72] = K[72]- Ghimj[441]*K[48]- Ghimj[442]*K[55]- Ghimj[443]*K[61]- Ghimj[444]*K[63]- Ghimj[445]*K[64]- Ghimj[446]*K[65]- Ghimj[447]*K[66]  - Ghimj[448]*K[67]- Ghimj[449]*K[68]- Ghimj[450]*K[69]- Ghimj[451]*K[70]- Ghimj[452]*K[71];
        K[72] = K[72]/ Ghimj[453];
        K[71] = (K[71]- Ghimj[440]*K[72])/(Ghimj[439]);
        K[70] = (K[70]- Ghimj[420]*K[71]- Ghimj[421]*K[72])/(Ghimj[419]);
        K[69] = (K[69]- Ghimj[402]*K[70]- Ghimj[403]*K[71]- Ghimj[404]*K[72])/(Ghimj[401]);
        K[68] = (K[68]- Ghimj[389]*K[69]- Ghimj[390]*K[70]- Ghimj[391]*K[71]- Ghimj[392]*K[72])/(Ghimj[388]);
        K[67] = (K[67]- Ghimj[361]*K[68]- Ghimj[362]*K[69]- Ghimj[363]*K[70]- Ghimj[364]*K[71]- Ghimj[365]*K[72])/(Ghimj[360]);
        K[66] = (K[66]- Ghimj[335]*K[67]- Ghimj[336]*K[68]- Ghimj[337]*K[69]- Ghimj[338]*K[70]- Ghimj[339]*K[71]- Ghimj[340]*K[72])/(Ghimj[334]);
        K[65] = (K[65]- Ghimj[315]*K[66]- Ghimj[316]*K[67]- Ghimj[317]*K[68]- Ghimj[318]*K[69]- Ghimj[319]*K[70]- Ghimj[320]*K[71]- Ghimj[321]  *K[72])/(Ghimj[314]);
        K[64] = (K[64]- Ghimj[287]*K[65]- Ghimj[288]*K[67]- Ghimj[289]*K[68]- Ghimj[290]*K[69]- Ghimj[291]*K[72])/(Ghimj[286]);
        K[63] = (K[63]- Ghimj[279]*K[67]- Ghimj[280]*K[68]- Ghimj[281]*K[69]- Ghimj[282]*K[72])/(Ghimj[278]);
        K[62] = (K[62]- Ghimj[268]*K[63]- Ghimj[269]*K[64]- Ghimj[270]*K[65]- Ghimj[271]*K[67]- Ghimj[272]*K[68]- Ghimj[273]*K[69]- Ghimj[274]  *K[70]- Ghimj[275]*K[72])/(Ghimj[267]);
        K[61] = (K[61]- Ghimj[259]*K[63]- Ghimj[260]*K[68]- Ghimj[261]*K[69]- Ghimj[262]*K[72])/(Ghimj[258]);
        K[60] = (K[60]- Ghimj[250]*K[61]- Ghimj[251]*K[63]- Ghimj[252]*K[64]- Ghimj[253]*K[65]- Ghimj[254]*K[66]- Ghimj[255]*K[67]- Ghimj[256]  *K[68])/(Ghimj[249]);
        K[59] = (K[59]- Ghimj[241]*K[61]- Ghimj[242]*K[63]- Ghimj[243]*K[65]- Ghimj[244]*K[67]- Ghimj[245]*K[68]- Ghimj[246]*K[69]- Ghimj[247]  *K[72])/(Ghimj[240]);
        K[58] = (K[58]- Ghimj[233]*K[62]- Ghimj[234]*K[64]- Ghimj[235]*K[65]- Ghimj[236]*K[67]- Ghimj[237]*K[68]- Ghimj[238]*K[71])/(Ghimj[232]);
        K[57] = (K[57]- Ghimj[219]*K[62]- Ghimj[220]*K[64]- Ghimj[221]*K[65]- Ghimj[222]*K[67]- Ghimj[223]*K[68]- Ghimj[224]*K[71])/(Ghimj[218]);
        K[56] = (K[56]- Ghimj[208]*K[59]- Ghimj[209]*K[63]- Ghimj[210]*K[67]- Ghimj[211]*K[68]- Ghimj[212]*K[69]- Ghimj[213]*K[72])/(Ghimj[207]);
        K[55] = (K[55]- Ghimj[202]*K[61]- Ghimj[203]*K[63]- Ghimj[204]*K[69]- Ghimj[205]*K[72])/(Ghimj[201]);
        K[54] = (K[54]- Ghimj[195]*K[56]- Ghimj[196]*K[63]- Ghimj[197]*K[67]- Ghimj[198]*K[68]- Ghimj[199]*K[69])/(Ghimj[194]);
        K[53] = (K[53]- Ghimj[189]*K[61]- Ghimj[190]*K[63]- Ghimj[191]*K[64]- Ghimj[192]*K[67]- Ghimj[193]*K[68])/(Ghimj[188]);
        K[52] = (K[52]- Ghimj[185]*K[65]- Ghimj[186]*K[67]- Ghimj[187]*K[71])/(Ghimj[184]);
        K[51] = (K[51]- Ghimj[182]*K[67]- Ghimj[183]*K[68])/(Ghimj[181]);
        K[50] = (K[50]- Ghimj[178]*K[64]- Ghimj[179]*K[65]- Ghimj[180]*K[72])/(Ghimj[177]);
        K[49] = (K[49]- Ghimj[175]*K[67]- Ghimj[176]*K[68])/(Ghimj[174]);
        K[48] = (K[48]- Ghimj[172]*K[69]- Ghimj[173]*K[72])/(Ghimj[171]);
        K[47] = (K[47]- Ghimj[166]*K[52]- Ghimj[167]*K[57]- Ghimj[168]*K[58]- Ghimj[169]*K[65]- Ghimj[170]*K[71])/(Ghimj[165]);
        K[46] = (K[46]- Ghimj[153]*K[67]- Ghimj[154]*K[68])/(Ghimj[152]);
        K[45] = (K[45]- Ghimj[150]*K[67]- Ghimj[151]*K[68])/(Ghimj[149]);
        K[44] = (K[44]- Ghimj[146]*K[50]- Ghimj[147]*K[60]- Ghimj[148]*K[65])/(Ghimj[145]);
        K[43] = (K[43]- Ghimj[144]*K[67])/(Ghimj[143]);
        K[42] = (K[42]- Ghimj[142]*K[67])/(Ghimj[141]);
        K[41] = (K[41]- Ghimj[140]*K[47])/(Ghimj[139]);
        K[40] = (K[40]- Ghimj[138]*K[65])/(Ghimj[137]);
        K[39] = (K[39]- Ghimj[136]*K[65])/(Ghimj[135]);
        K[38] = (K[38]- Ghimj[132]*K[60]- Ghimj[133]*K[66]- Ghimj[134]*K[71])/(Ghimj[131]);
        K[37] = (K[37]- Ghimj[130]*K[65])/(Ghimj[129]);
        K[36] = (K[36]- Ghimj[128]*K[65])/(Ghimj[127]);
        K[35] = (K[35]- Ghimj[126]*K[65])/(Ghimj[125]);
        K[34] = (K[34]- Ghimj[124]*K[65])/(Ghimj[123]);
        K[33] = (K[33]- Ghimj[122]*K[65])/(Ghimj[121]);
        K[32] = (K[32]- Ghimj[120]*K[65])/(Ghimj[119]);
        K[31] = (K[31]- Ghimj[117]*K[40]- Ghimj[118]*K[65])/(Ghimj[116]);
        K[30] = (K[30]- Ghimj[115]*K[69])/(Ghimj[114]);
        K[29] = (K[29]- Ghimj[113]*K[47])/(Ghimj[112]);
        K[28] = (K[28]- Ghimj[111]*K[41])/(Ghimj[110]);
        K[27] = (K[27]- Ghimj[104]*K[33]- Ghimj[105]*K[34]- Ghimj[106]*K[35]- Ghimj[107]*K[36]- Ghimj[108]*K[37]- Ghimj[109]*K[65])/(Ghimj[103]);
        K[26] = (K[26]- Ghimj[99]*K[50]- Ghimj[100]*K[59]- Ghimj[101]*K[60]- Ghimj[102]*K[65])/(Ghimj[98]);
        K[25] = (K[25]- Ghimj[96]*K[39]- Ghimj[97]*K[65])/(Ghimj[95]);
        K[24] = (K[24]- Ghimj[93]*K[32]- Ghimj[94]*K[65])/(Ghimj[92]);
        K[23] = (K[23]- Ghimj[90]*K[65]- Ghimj[91]*K[66])/(Ghimj[89]);
        K[22] = (K[22]- Ghimj[87]*K[41]- Ghimj[88]*K[47])/(Ghimj[86]);
        K[21] = (K[21]- Ghimj[83]*K[48]- Ghimj[84]*K[69]- Ghimj[85]*K[72])/(Ghimj[82]);
        K[20] = (K[20]- Ghimj[79]*K[65]- Ghimj[80]*K[66]- Ghimj[81]*K[71])/(Ghimj[78]);
        K[19] = (K[19]- Ghimj[74]*K[48]- Ghimj[75]*K[69]- Ghimj[76]*K[70]- Ghimj[77]*K[72])/(Ghimj[73]);
        K[18] = (K[18]- Ghimj[69]*K[48]- Ghimj[70]*K[69]- Ghimj[71]*K[70]- Ghimj[72]*K[72])/(Ghimj[68]);
        K[17] = (K[17]- Ghimj[56]*K[41]- Ghimj[57]*K[44]- Ghimj[58]*K[47]- Ghimj[59]*K[53]- Ghimj[60]*K[60]- Ghimj[61]*K[64]- Ghimj[62]*K[65]  - Ghimj[63]*K[66]- Ghimj[64]*K[69]- Ghimj[65]*K[70]- Ghimj[66]*K[71]- Ghimj[67]*K[72])/(Ghimj[55]);
        K[16] = (K[16]- Ghimj[52]*K[41]- Ghimj[53]*K[47])/(Ghimj[51]);
        K[15] = (K[15]- Ghimj[49]*K[66]- Ghimj[50]*K[71])/(Ghimj[48]);
        K[14] = (K[14]- Ghimj[47]*K[71])/(Ghimj[46]);
        K[13] = (K[13]- Ghimj[44]*K[66]- Ghimj[45]*K[71])/(Ghimj[43]);
        K[12] = (K[12]- Ghimj[41]*K[65]- Ghimj[42]*K[66])/(Ghimj[40]);
        K[11] = (K[11]- Ghimj[37]*K[45]- Ghimj[38]*K[67]- Ghimj[39]*K[68])/(Ghimj[36]);
        K[10] = K[10]/ Ghimj[35];
        K[9] = K[9]/ Ghimj[33];
        K[8] = (K[8]- Ghimj[28]*K[38]- Ghimj[29]*K[50]- Ghimj[30]*K[65]- Ghimj[31]*K[71])/(Ghimj[27]);
        K[7] = (K[7]- Ghimj[25]*K[44]- Ghimj[26]*K[65])/(Ghimj[24]);
        K[6] = (K[6]- Ghimj[19]*K[59]- Ghimj[20]*K[65]- Ghimj[21]*K[70]- Ghimj[22]*K[71]- Ghimj[23]*K[72])/(Ghimj[18]);
        K[5] = (K[5]- Ghimj[16]*K[69]- Ghimj[17]*K[72])/(Ghimj[15]);
        K[4] = (K[4]- Ghimj[13]*K[69]- Ghimj[14]*K[71])/(Ghimj[12]);
        K[3] = (K[3]- Ghimj[9]*K[46]- Ghimj[10]*K[67]- Ghimj[11]*K[68])/(Ghimj[8]);
        K[2] = (K[2]- Ghimj[5]*K[62]- Ghimj[6]*K[67]- Ghimj[7]*K[68])/(Ghimj[4]);
        K[1] = (K[1]- Ghimj[2]*K[53]- Ghimj[3]*K[64])/(Ghimj[1]);
        K[0] = K[0]/ Ghimj[0];
}

__device__ void ros_Solve(REAL * __restrict__ Ghimj, REAL * __restrict__ K, int &Nsol, const int istage, const int ros_S)
{

    #pragma unroll 4 
    for (int i=0;i<LU_NONZERO-16;i+=16){
        prefetch_ll1(&Ghimj[i]);
    }

    kppSolve(Ghimj, K, istage, ros_S);
    Nsol++;
}

__device__ void kppDecomp(REAL *Ghimj, int VL_GLO)
{
    REAL a=0.0;

 REAL dummy, W_0, W_1, W_2, W_3, W_4, W_5, W_6, W_7, W_8, W_9, W_10, W_11, W_12, W_13, W_14, W_15, W_16, W_17, W_18, W_19, W_20, W_21, W_22, W_23, W_24, W_25, W_26, W_27, W_28, W_29, W_30, W_31, W_32, W_33, W_34, W_35, W_36, W_37, W_38, W_39, W_40, W_41, W_42, W_43, W_44, W_45, W_46, W_47, W_48, W_49, W_50, W_51, W_52, W_53, W_54, W_55, W_56, W_57, W_58, W_59, W_60, W_61, W_62, W_63, W_64, W_65, W_66, W_67, W_68, W_69, W_70, W_71, W_72, W_73;

        W_0 = Ghimj[32];
        W_9 = Ghimj[33];
        a = - W_0/ Ghimj[0];
        W_0 = -a;
        Ghimj[32] = W_0;
        Ghimj[33] = W_9;
        W_0 = Ghimj[34];
        W_10 = Ghimj[35];
        a = - W_0/ Ghimj[0];
        W_0 = -a;
        Ghimj[34] = W_0;
        Ghimj[35] = W_10;
        W_0 = Ghimj[54];
        W_17 = Ghimj[55];
        W_41 = Ghimj[56];
        W_44 = Ghimj[57];
        W_47 = Ghimj[58];
        W_53 = Ghimj[59];
        W_60 = Ghimj[60];
        W_64 = Ghimj[61];
        W_65 = Ghimj[62];
        W_66 = Ghimj[63];
        W_69 = Ghimj[64];
        W_70 = Ghimj[65];
        W_71 = Ghimj[66];
        W_72 = Ghimj[67];
        a = - W_0/ Ghimj[0];
        W_0 = -a;
        Ghimj[54] = W_0;
        Ghimj[55] = W_17;
        Ghimj[56] = W_41;
        Ghimj[57] = W_44;
        Ghimj[58] = W_47;
        Ghimj[59] = W_53;
        Ghimj[60] = W_60;
        Ghimj[61] = W_64;
        Ghimj[62] = W_65;
        Ghimj[63] = W_66;
        Ghimj[64] = W_69;
        Ghimj[65] = W_70;
        Ghimj[66] = W_71;
        Ghimj[67] = W_72;
        W_29 = Ghimj[155];
        W_32 = Ghimj[156];
        W_33 = Ghimj[157];
        W_34 = Ghimj[158];
        W_35 = Ghimj[159];
        W_36 = Ghimj[160];
        W_37 = Ghimj[161];
        W_39 = Ghimj[162];
        W_40 = Ghimj[163];
        W_41 = Ghimj[164];
        W_47 = Ghimj[165];
        W_52 = Ghimj[166];
        W_57 = Ghimj[167];
        W_58 = Ghimj[168];
        W_65 = Ghimj[169];
        W_71 = Ghimj[170];
        a = - W_29/ Ghimj[112];
        W_29 = -a;
        W_47 = W_47+ a *Ghimj[113];
        a = - W_32/ Ghimj[119];
        W_32 = -a;
        W_65 = W_65+ a *Ghimj[120];
        a = - W_33/ Ghimj[121];
        W_33 = -a;
        W_65 = W_65+ a *Ghimj[122];
        a = - W_34/ Ghimj[123];
        W_34 = -a;
        W_65 = W_65+ a *Ghimj[124];
        a = - W_35/ Ghimj[125];
        W_35 = -a;
        W_65 = W_65+ a *Ghimj[126];
        a = - W_36/ Ghimj[127];
        W_36 = -a;
        W_65 = W_65+ a *Ghimj[128];
        a = - W_37/ Ghimj[129];
        W_37 = -a;
        W_65 = W_65+ a *Ghimj[130];
        a = - W_39/ Ghimj[135];
        W_39 = -a;
        W_65 = W_65+ a *Ghimj[136];
        a = - W_40/ Ghimj[137];
        W_40 = -a;
        W_65 = W_65+ a *Ghimj[138];
        a = - W_41/ Ghimj[139];
        W_41 = -a;
        W_47 = W_47+ a *Ghimj[140];
        Ghimj[155] = W_29;
        Ghimj[156] = W_32;
        Ghimj[157] = W_33;
        Ghimj[158] = W_34;
        Ghimj[159] = W_35;
        Ghimj[160] = W_36;
        Ghimj[161] = W_37;
        Ghimj[162] = W_39;
        Ghimj[163] = W_40;
        Ghimj[164] = W_41;
        Ghimj[165] = W_47;
        Ghimj[166] = W_52;
        Ghimj[167] = W_57;
        Ghimj[168] = W_58;
        Ghimj[169] = W_65;
        Ghimj[170] = W_71;
        W_48 = Ghimj[200];
        W_55 = Ghimj[201];
        W_61 = Ghimj[202];
        W_63 = Ghimj[203];
        W_69 = Ghimj[204];
        W_72 = Ghimj[205];
        a = - W_48/ Ghimj[171];
        W_48 = -a;
        W_69 = W_69+ a *Ghimj[172];
        W_72 = W_72+ a *Ghimj[173];
        Ghimj[200] = W_48;
        Ghimj[201] = W_55;
        Ghimj[202] = W_61;
        Ghimj[203] = W_63;
        Ghimj[204] = W_69;
        Ghimj[205] = W_72;
        W_54 = Ghimj[206];
        W_56 = Ghimj[207];
        W_59 = Ghimj[208];
        W_63 = Ghimj[209];
        W_67 = Ghimj[210];
        W_68 = Ghimj[211];
        W_69 = Ghimj[212];
        W_72 = Ghimj[213];
        a = - W_54/ Ghimj[194];
        W_54 = -a;
        W_56 = W_56+ a *Ghimj[195];
        W_63 = W_63+ a *Ghimj[196];
        W_67 = W_67+ a *Ghimj[197];
        W_68 = W_68+ a *Ghimj[198];
        W_69 = W_69+ a *Ghimj[199];
        Ghimj[206] = W_54;
        Ghimj[207] = W_56;
        Ghimj[208] = W_59;
        Ghimj[209] = W_63;
        Ghimj[210] = W_67;
        Ghimj[211] = W_68;
        Ghimj[212] = W_69;
        Ghimj[213] = W_72;
        W_45 = Ghimj[214];
        W_46 = Ghimj[215];
        W_49 = Ghimj[216];
        W_51 = Ghimj[217];
        W_57 = Ghimj[218];
        W_62 = Ghimj[219];
        W_64 = Ghimj[220];
        W_65 = Ghimj[221];
        W_67 = Ghimj[222];
        W_68 = Ghimj[223];
        W_71 = Ghimj[224];
        a = - W_45/ Ghimj[149];
        W_45 = -a;
        W_67 = W_67+ a *Ghimj[150];
        W_68 = W_68+ a *Ghimj[151];
        a = - W_46/ Ghimj[152];
        W_46 = -a;
        W_67 = W_67+ a *Ghimj[153];
        W_68 = W_68+ a *Ghimj[154];
        a = - W_49/ Ghimj[174];
        W_49 = -a;
        W_67 = W_67+ a *Ghimj[175];
        W_68 = W_68+ a *Ghimj[176];
        a = - W_51/ Ghimj[181];
        W_51 = -a;
        W_67 = W_67+ a *Ghimj[182];
        W_68 = W_68+ a *Ghimj[183];
        Ghimj[214] = W_45;
        Ghimj[215] = W_46;
        Ghimj[216] = W_49;
        Ghimj[217] = W_51;
        Ghimj[218] = W_57;
        Ghimj[219] = W_62;
        Ghimj[220] = W_64;
        Ghimj[221] = W_65;
        Ghimj[222] = W_67;
        Ghimj[223] = W_68;
        Ghimj[224] = W_71;
        W_42 = Ghimj[225];
        W_43 = Ghimj[226];
        W_45 = Ghimj[227];
        W_46 = Ghimj[228];
        W_49 = Ghimj[229];
        W_51 = Ghimj[230];
        W_52 = Ghimj[231];
        W_58 = Ghimj[232];
        W_62 = Ghimj[233];
        W_64 = Ghimj[234];
        W_65 = Ghimj[235];
        W_67 = Ghimj[236];
        W_68 = Ghimj[237];
        W_71 = Ghimj[238];
        a = - W_42/ Ghimj[141];
        W_42 = -a;
        W_67 = W_67+ a *Ghimj[142];
        a = - W_43/ Ghimj[143];
        W_43 = -a;
        W_67 = W_67+ a *Ghimj[144];
        a = - W_45/ Ghimj[149];
        W_45 = -a;
        W_67 = W_67+ a *Ghimj[150];
        W_68 = W_68+ a *Ghimj[151];
        a = - W_46/ Ghimj[152];
        W_46 = -a;
        W_67 = W_67+ a *Ghimj[153];
        W_68 = W_68+ a *Ghimj[154];
        a = - W_49/ Ghimj[174];
        W_49 = -a;
        W_67 = W_67+ a *Ghimj[175];
        W_68 = W_68+ a *Ghimj[176];
        a = - W_51/ Ghimj[181];
        W_51 = -a;
        W_67 = W_67+ a *Ghimj[182];
        W_68 = W_68+ a *Ghimj[183];
        a = - W_52/ Ghimj[184];
        W_52 = -a;
        W_65 = W_65+ a *Ghimj[185];
        W_67 = W_67+ a *Ghimj[186];
        W_71 = W_71+ a *Ghimj[187];
        Ghimj[225] = W_42;
        Ghimj[226] = W_43;
        Ghimj[227] = W_45;
        Ghimj[228] = W_46;
        Ghimj[229] = W_49;
        Ghimj[230] = W_51;
        Ghimj[231] = W_52;
        Ghimj[232] = W_58;
        Ghimj[233] = W_62;
        Ghimj[234] = W_64;
        Ghimj[235] = W_65;
        Ghimj[236] = W_67;
        Ghimj[237] = W_68;
        Ghimj[238] = W_71;
        W_56 = Ghimj[239];
        W_59 = Ghimj[240];
        W_61 = Ghimj[241];
        W_63 = Ghimj[242];
        W_65 = Ghimj[243];
        W_67 = Ghimj[244];
        W_68 = Ghimj[245];
        W_69 = Ghimj[246];
        W_72 = Ghimj[247];
        a = - W_56/ Ghimj[207];
        W_56 = -a;
        W_59 = W_59+ a *Ghimj[208];
        W_63 = W_63+ a *Ghimj[209];
        W_67 = W_67+ a *Ghimj[210];
        W_68 = W_68+ a *Ghimj[211];
        W_69 = W_69+ a *Ghimj[212];
        W_72 = W_72+ a *Ghimj[213];
        Ghimj[239] = W_56;
        Ghimj[240] = W_59;
        Ghimj[241] = W_61;
        Ghimj[242] = W_63;
        Ghimj[243] = W_65;
        Ghimj[244] = W_67;
        Ghimj[245] = W_68;
        Ghimj[246] = W_69;
        Ghimj[247] = W_72;
        W_53 = Ghimj[248];
        W_60 = Ghimj[249];
        W_61 = Ghimj[250];
        W_63 = Ghimj[251];
        W_64 = Ghimj[252];
        W_65 = Ghimj[253];
        W_66 = Ghimj[254];
        W_67 = Ghimj[255];
        W_68 = Ghimj[256];
        a = - W_53/ Ghimj[188];
        W_53 = -a;
        W_61 = W_61+ a *Ghimj[189];
        W_63 = W_63+ a *Ghimj[190];
        W_64 = W_64+ a *Ghimj[191];
        W_67 = W_67+ a *Ghimj[192];
        W_68 = W_68+ a *Ghimj[193];
        Ghimj[248] = W_53;
        Ghimj[249] = W_60;
        Ghimj[250] = W_61;
        Ghimj[251] = W_63;
        Ghimj[252] = W_64;
        Ghimj[253] = W_65;
        Ghimj[254] = W_66;
        Ghimj[255] = W_67;
        Ghimj[256] = W_68;
        W_55 = Ghimj[257];
        W_61 = Ghimj[258];
        W_63 = Ghimj[259];
        W_68 = Ghimj[260];
        W_69 = Ghimj[261];
        W_72 = Ghimj[262];
        a = - W_55/ Ghimj[201];
        W_55 = -a;
        W_61 = W_61+ a *Ghimj[202];
        W_63 = W_63+ a *Ghimj[203];
        W_69 = W_69+ a *Ghimj[204];
        W_72 = W_72+ a *Ghimj[205];
        Ghimj[257] = W_55;
        Ghimj[258] = W_61;
        Ghimj[259] = W_63;
        Ghimj[260] = W_68;
        Ghimj[261] = W_69;
        Ghimj[262] = W_72;
        W_50 = Ghimj[263];
        W_51 = Ghimj[264];
        W_53 = Ghimj[265];
        W_61 = Ghimj[266];
        W_62 = Ghimj[267];
        W_63 = Ghimj[268];
        W_64 = Ghimj[269];
        W_65 = Ghimj[270];
        W_67 = Ghimj[271];
        W_68 = Ghimj[272];
        W_69 = Ghimj[273];
        W_70 = Ghimj[274];
        W_72 = Ghimj[275];
        a = - W_50/ Ghimj[177];
        W_50 = -a;
        W_64 = W_64+ a *Ghimj[178];
        W_65 = W_65+ a *Ghimj[179];
        W_72 = W_72+ a *Ghimj[180];
        a = - W_51/ Ghimj[181];
        W_51 = -a;
        W_67 = W_67+ a *Ghimj[182];
        W_68 = W_68+ a *Ghimj[183];
        a = - W_53/ Ghimj[188];
        W_53 = -a;
        W_61 = W_61+ a *Ghimj[189];
        W_63 = W_63+ a *Ghimj[190];
        W_64 = W_64+ a *Ghimj[191];
        W_67 = W_67+ a *Ghimj[192];
        W_68 = W_68+ a *Ghimj[193];
        a = - W_61/ Ghimj[258];
        W_61 = -a;
        W_63 = W_63+ a *Ghimj[259];
        W_68 = W_68+ a *Ghimj[260];
        W_69 = W_69+ a *Ghimj[261];
        W_72 = W_72+ a *Ghimj[262];
        Ghimj[263] = W_50;
        Ghimj[264] = W_51;
        Ghimj[265] = W_53;
        Ghimj[266] = W_61;
        Ghimj[267] = W_62;
        Ghimj[268] = W_63;
        Ghimj[269] = W_64;
        Ghimj[270] = W_65;
        Ghimj[271] = W_67;
        Ghimj[272] = W_68;
        Ghimj[273] = W_69;
        Ghimj[274] = W_70;
        Ghimj[275] = W_72;
        W_55 = Ghimj[276];
        W_61 = Ghimj[277];
        W_63 = Ghimj[278];
        W_67 = Ghimj[279];
        W_68 = Ghimj[280];
        W_69 = Ghimj[281];
        W_72 = Ghimj[282];
        a = - W_55/ Ghimj[201];
        W_55 = -a;
        W_61 = W_61+ a *Ghimj[202];
        W_63 = W_63+ a *Ghimj[203];
        W_69 = W_69+ a *Ghimj[204];
        W_72 = W_72+ a *Ghimj[205];
        a = - W_61/ Ghimj[258];
        W_61 = -a;
        W_63 = W_63+ a *Ghimj[259];
        W_68 = W_68+ a *Ghimj[260];
        W_69 = W_69+ a *Ghimj[261];
        W_72 = W_72+ a *Ghimj[262];
        Ghimj[276] = W_55;
        Ghimj[277] = W_61;
        Ghimj[278] = W_63;
        Ghimj[279] = W_67;
        Ghimj[280] = W_68;
        Ghimj[281] = W_69;
        Ghimj[282] = W_72;
        W_53 = Ghimj[283];
        W_61 = Ghimj[284];
        W_63 = Ghimj[285];
        W_64 = Ghimj[286];
        W_65 = Ghimj[287];
        W_67 = Ghimj[288];
        W_68 = Ghimj[289];
        W_69 = Ghimj[290];
        W_72 = Ghimj[291];
        a = - W_53/ Ghimj[188];
        W_53 = -a;
        W_61 = W_61+ a *Ghimj[189];
        W_63 = W_63+ a *Ghimj[190];
        W_64 = W_64+ a *Ghimj[191];
        W_67 = W_67+ a *Ghimj[192];
        W_68 = W_68+ a *Ghimj[193];
        a = - W_61/ Ghimj[258];
        W_61 = -a;
        W_63 = W_63+ a *Ghimj[259];
        W_68 = W_68+ a *Ghimj[260];
        W_69 = W_69+ a *Ghimj[261];
        W_72 = W_72+ a *Ghimj[262];
        a = - W_63/ Ghimj[278];
        W_63 = -a;
        W_67 = W_67+ a *Ghimj[279];
        W_68 = W_68+ a *Ghimj[280];
        W_69 = W_69+ a *Ghimj[281];
        W_72 = W_72+ a *Ghimj[282];
        Ghimj[283] = W_53;
        Ghimj[284] = W_61;
        Ghimj[285] = W_63;
        Ghimj[286] = W_64;
        Ghimj[287] = W_65;
        Ghimj[288] = W_67;
        Ghimj[289] = W_68;
        Ghimj[290] = W_69;
        Ghimj[291] = W_72;
        W_32 = Ghimj[292];
        W_33 = Ghimj[293];
        W_34 = Ghimj[294];
        W_35 = Ghimj[295];
        W_36 = Ghimj[296];
        W_37 = Ghimj[297];
        W_39 = Ghimj[298];
        W_40 = Ghimj[299];
        W_41 = Ghimj[300];
        W_44 = Ghimj[301];
        W_47 = Ghimj[302];
        W_50 = Ghimj[303];
        W_51 = Ghimj[304];
        W_52 = Ghimj[305];
        W_57 = Ghimj[306];
        W_58 = Ghimj[307];
        W_59 = Ghimj[308];
        W_60 = Ghimj[309];
        W_61 = Ghimj[310];
        W_62 = Ghimj[311];
        W_63 = Ghimj[312];
        W_64 = Ghimj[313];
        W_65 = Ghimj[314];
        W_66 = Ghimj[315];
        W_67 = Ghimj[316];
        W_68 = Ghimj[317];
        W_69 = Ghimj[318];
        W_70 = Ghimj[319];
        W_71 = Ghimj[320];
        W_72 = Ghimj[321];
        a = - W_32/ Ghimj[119];
        W_32 = -a;
        W_65 = W_65+ a *Ghimj[120];
        a = - W_33/ Ghimj[121];
        W_33 = -a;
        W_65 = W_65+ a *Ghimj[122];
        a = - W_34/ Ghimj[123];
        W_34 = -a;
        W_65 = W_65+ a *Ghimj[124];
        a = - W_35/ Ghimj[125];
        W_35 = -a;
        W_65 = W_65+ a *Ghimj[126];
        a = - W_36/ Ghimj[127];
        W_36 = -a;
        W_65 = W_65+ a *Ghimj[128];
        a = - W_37/ Ghimj[129];
        W_37 = -a;
        W_65 = W_65+ a *Ghimj[130];
        a = - W_39/ Ghimj[135];
        W_39 = -a;
        W_65 = W_65+ a *Ghimj[136];
        a = - W_40/ Ghimj[137];
        W_40 = -a;
        W_65 = W_65+ a *Ghimj[138];
        a = - W_41/ Ghimj[139];
        W_41 = -a;
        W_47 = W_47+ a *Ghimj[140];
        a = - W_44/ Ghimj[145];
        W_44 = -a;
        W_50 = W_50+ a *Ghimj[146];
        W_60 = W_60+ a *Ghimj[147];
        W_65 = W_65+ a *Ghimj[148];
        a = - W_47/ Ghimj[165];
        W_47 = -a;
        W_52 = W_52+ a *Ghimj[166];
        W_57 = W_57+ a *Ghimj[167];
        W_58 = W_58+ a *Ghimj[168];
        W_65 = W_65+ a *Ghimj[169];
        W_71 = W_71+ a *Ghimj[170];
        a = - W_50/ Ghimj[177];
        W_50 = -a;
        W_64 = W_64+ a *Ghimj[178];
        W_65 = W_65+ a *Ghimj[179];
        W_72 = W_72+ a *Ghimj[180];
        a = - W_51/ Ghimj[181];
        W_51 = -a;
        W_67 = W_67+ a *Ghimj[182];
        W_68 = W_68+ a *Ghimj[183];
        a = - W_52/ Ghimj[184];
        W_52 = -a;
        W_65 = W_65+ a *Ghimj[185];
        W_67 = W_67+ a *Ghimj[186];
        W_71 = W_71+ a *Ghimj[187];
        a = - W_57/ Ghimj[218];
        W_57 = -a;
        W_62 = W_62+ a *Ghimj[219];
        W_64 = W_64+ a *Ghimj[220];
        W_65 = W_65+ a *Ghimj[221];
        W_67 = W_67+ a *Ghimj[222];
        W_68 = W_68+ a *Ghimj[223];
        W_71 = W_71+ a *Ghimj[224];
        a = - W_58/ Ghimj[232];
        W_58 = -a;
        W_62 = W_62+ a *Ghimj[233];
        W_64 = W_64+ a *Ghimj[234];
        W_65 = W_65+ a *Ghimj[235];
        W_67 = W_67+ a *Ghimj[236];
        W_68 = W_68+ a *Ghimj[237];
        W_71 = W_71+ a *Ghimj[238];
        a = - W_59/ Ghimj[240];
        W_59 = -a;
        W_61 = W_61+ a *Ghimj[241];
        W_63 = W_63+ a *Ghimj[242];
        W_65 = W_65+ a *Ghimj[243];
        W_67 = W_67+ a *Ghimj[244];
        W_68 = W_68+ a *Ghimj[245];
        W_69 = W_69+ a *Ghimj[246];
        W_72 = W_72+ a *Ghimj[247];
        a = - W_60/ Ghimj[249];
        W_60 = -a;
        W_61 = W_61+ a *Ghimj[250];
        W_63 = W_63+ a *Ghimj[251];
        W_64 = W_64+ a *Ghimj[252];
        W_65 = W_65+ a *Ghimj[253];
        W_66 = W_66+ a *Ghimj[254];
        W_67 = W_67+ a *Ghimj[255];
        W_68 = W_68+ a *Ghimj[256];
        a = - W_61/ Ghimj[258];
        W_61 = -a;
        W_63 = W_63+ a *Ghimj[259];
        W_68 = W_68+ a *Ghimj[260];
        W_69 = W_69+ a *Ghimj[261];
        W_72 = W_72+ a *Ghimj[262];
        a = - W_62/ Ghimj[267];
        W_62 = -a;
        W_63 = W_63+ a *Ghimj[268];
        W_64 = W_64+ a *Ghimj[269];
        W_65 = W_65+ a *Ghimj[270];
        W_67 = W_67+ a *Ghimj[271];
        W_68 = W_68+ a *Ghimj[272];
        W_69 = W_69+ a *Ghimj[273];
        W_70 = W_70+ a *Ghimj[274];
        W_72 = W_72+ a *Ghimj[275];
        a = - W_63/ Ghimj[278];
        W_63 = -a;
        W_67 = W_67+ a *Ghimj[279];
        W_68 = W_68+ a *Ghimj[280];
        W_69 = W_69+ a *Ghimj[281];
        W_72 = W_72+ a *Ghimj[282];
        a = - W_64/ Ghimj[286];
        W_64 = -a;
        W_65 = W_65+ a *Ghimj[287];
        W_67 = W_67+ a *Ghimj[288];
        W_68 = W_68+ a *Ghimj[289];
        W_69 = W_69+ a *Ghimj[290];
        W_72 = W_72+ a *Ghimj[291];
        Ghimj[292] = W_32;
        Ghimj[293] = W_33;
        Ghimj[294] = W_34;
        Ghimj[295] = W_35;
        Ghimj[296] = W_36;
        Ghimj[297] = W_37;
        Ghimj[298] = W_39;
        Ghimj[299] = W_40;
        Ghimj[300] = W_41;
        Ghimj[301] = W_44;
        Ghimj[302] = W_47;
        Ghimj[303] = W_50;
        Ghimj[304] = W_51;
        Ghimj[305] = W_52;
        Ghimj[306] = W_57;
        Ghimj[307] = W_58;
        Ghimj[308] = W_59;
        Ghimj[309] = W_60;
        Ghimj[310] = W_61;
        Ghimj[311] = W_62;
        Ghimj[312] = W_63;
        Ghimj[313] = W_64;
        Ghimj[314] = W_65;
        Ghimj[315] = W_66;
        Ghimj[316] = W_67;
        Ghimj[317] = W_68;
        Ghimj[318] = W_69;
        Ghimj[319] = W_70;
        Ghimj[320] = W_71;
        Ghimj[321] = W_72;
        W_28 = Ghimj[322];
        W_41 = Ghimj[323];
        W_47 = Ghimj[324];
        W_52 = Ghimj[325];
        W_57 = Ghimj[326];
        W_58 = Ghimj[327];
        W_60 = Ghimj[328];
        W_61 = Ghimj[329];
        W_62 = Ghimj[330];
        W_63 = Ghimj[331];
        W_64 = Ghimj[332];
        W_65 = Ghimj[333];
        W_66 = Ghimj[334];
        W_67 = Ghimj[335];
        W_68 = Ghimj[336];
        W_69 = Ghimj[337];
        W_70 = Ghimj[338];
        W_71 = Ghimj[339];
        W_72 = Ghimj[340];
        a = - W_28/ Ghimj[110];
        W_28 = -a;
        W_41 = W_41+ a *Ghimj[111];
        a = - W_41/ Ghimj[139];
        W_41 = -a;
        W_47 = W_47+ a *Ghimj[140];
        a = - W_47/ Ghimj[165];
        W_47 = -a;
        W_52 = W_52+ a *Ghimj[166];
        W_57 = W_57+ a *Ghimj[167];
        W_58 = W_58+ a *Ghimj[168];
        W_65 = W_65+ a *Ghimj[169];
        W_71 = W_71+ a *Ghimj[170];
        a = - W_52/ Ghimj[184];
        W_52 = -a;
        W_65 = W_65+ a *Ghimj[185];
        W_67 = W_67+ a *Ghimj[186];
        W_71 = W_71+ a *Ghimj[187];
        a = - W_57/ Ghimj[218];
        W_57 = -a;
        W_62 = W_62+ a *Ghimj[219];
        W_64 = W_64+ a *Ghimj[220];
        W_65 = W_65+ a *Ghimj[221];
        W_67 = W_67+ a *Ghimj[222];
        W_68 = W_68+ a *Ghimj[223];
        W_71 = W_71+ a *Ghimj[224];
        a = - W_58/ Ghimj[232];
        W_58 = -a;
        W_62 = W_62+ a *Ghimj[233];
        W_64 = W_64+ a *Ghimj[234];
        W_65 = W_65+ a *Ghimj[235];
        W_67 = W_67+ a *Ghimj[236];
        W_68 = W_68+ a *Ghimj[237];
        W_71 = W_71+ a *Ghimj[238];
        a = - W_60/ Ghimj[249];
        W_60 = -a;
        W_61 = W_61+ a *Ghimj[250];
        W_63 = W_63+ a *Ghimj[251];
        W_64 = W_64+ a *Ghimj[252];
        W_65 = W_65+ a *Ghimj[253];
        W_66 = W_66+ a *Ghimj[254];
        W_67 = W_67+ a *Ghimj[255];
        W_68 = W_68+ a *Ghimj[256];
        a = - W_61/ Ghimj[258];
        W_61 = -a;
        W_63 = W_63+ a *Ghimj[259];
        W_68 = W_68+ a *Ghimj[260];
        W_69 = W_69+ a *Ghimj[261];
        W_72 = W_72+ a *Ghimj[262];
        a = - W_62/ Ghimj[267];
        W_62 = -a;
        W_63 = W_63+ a *Ghimj[268];
        W_64 = W_64+ a *Ghimj[269];
        W_65 = W_65+ a *Ghimj[270];
        W_67 = W_67+ a *Ghimj[271];
        W_68 = W_68+ a *Ghimj[272];
        W_69 = W_69+ a *Ghimj[273];
        W_70 = W_70+ a *Ghimj[274];
        W_72 = W_72+ a *Ghimj[275];
        a = - W_63/ Ghimj[278];
        W_63 = -a;
        W_67 = W_67+ a *Ghimj[279];
        W_68 = W_68+ a *Ghimj[280];
        W_69 = W_69+ a *Ghimj[281];
        W_72 = W_72+ a *Ghimj[282];
        a = - W_64/ Ghimj[286];
        W_64 = -a;
        W_65 = W_65+ a *Ghimj[287];
        W_67 = W_67+ a *Ghimj[288];
        W_68 = W_68+ a *Ghimj[289];
        W_69 = W_69+ a *Ghimj[290];
        W_72 = W_72+ a *Ghimj[291];
        a = - W_65/ Ghimj[314];
        W_65 = -a;
        W_66 = W_66+ a *Ghimj[315];
        W_67 = W_67+ a *Ghimj[316];
        W_68 = W_68+ a *Ghimj[317];
        W_69 = W_69+ a *Ghimj[318];
        W_70 = W_70+ a *Ghimj[319];
        W_71 = W_71+ a *Ghimj[320];
        W_72 = W_72+ a *Ghimj[321];
        Ghimj[322] = W_28;
        Ghimj[323] = W_41;
        Ghimj[324] = W_47;
        Ghimj[325] = W_52;
        Ghimj[326] = W_57;
        Ghimj[327] = W_58;
        Ghimj[328] = W_60;
        Ghimj[329] = W_61;
        Ghimj[330] = W_62;
        Ghimj[331] = W_63;
        Ghimj[332] = W_64;
        Ghimj[333] = W_65;
        Ghimj[334] = W_66;
        Ghimj[335] = W_67;
        Ghimj[336] = W_68;
        Ghimj[337] = W_69;
        Ghimj[338] = W_70;
        Ghimj[339] = W_71;
        Ghimj[340] = W_72;
        W_39 = Ghimj[341];
        W_42 = Ghimj[342];
        W_43 = Ghimj[343];
        W_45 = Ghimj[344];
        W_46 = Ghimj[345];
        W_48 = Ghimj[346];
        W_49 = Ghimj[347];
        W_51 = Ghimj[348];
        W_52 = Ghimj[349];
        W_54 = Ghimj[350];
        W_56 = Ghimj[351];
        W_58 = Ghimj[352];
        W_59 = Ghimj[353];
        W_61 = Ghimj[354];
        W_62 = Ghimj[355];
        W_63 = Ghimj[356];
        W_64 = Ghimj[357];
        W_65 = Ghimj[358];
        W_66 = Ghimj[359];
        W_67 = Ghimj[360];
        W_68 = Ghimj[361];
        W_69 = Ghimj[362];
        W_70 = Ghimj[363];
        W_71 = Ghimj[364];
        W_72 = Ghimj[365];
        a = - W_39/ Ghimj[135];
        W_39 = -a;
        W_65 = W_65+ a *Ghimj[136];
        a = - W_42/ Ghimj[141];
        W_42 = -a;
        W_67 = W_67+ a *Ghimj[142];
        a = - W_43/ Ghimj[143];
        W_43 = -a;
        W_67 = W_67+ a *Ghimj[144];
        a = - W_45/ Ghimj[149];
        W_45 = -a;
        W_67 = W_67+ a *Ghimj[150];
        W_68 = W_68+ a *Ghimj[151];
        a = - W_46/ Ghimj[152];
        W_46 = -a;
        W_67 = W_67+ a *Ghimj[153];
        W_68 = W_68+ a *Ghimj[154];
        a = - W_48/ Ghimj[171];
        W_48 = -a;
        W_69 = W_69+ a *Ghimj[172];
        W_72 = W_72+ a *Ghimj[173];
        a = - W_49/ Ghimj[174];
        W_49 = -a;
        W_67 = W_67+ a *Ghimj[175];
        W_68 = W_68+ a *Ghimj[176];
        a = - W_51/ Ghimj[181];
        W_51 = -a;
        W_67 = W_67+ a *Ghimj[182];
        W_68 = W_68+ a *Ghimj[183];
        a = - W_52/ Ghimj[184];
        W_52 = -a;
        W_65 = W_65+ a *Ghimj[185];
        W_67 = W_67+ a *Ghimj[186];
        W_71 = W_71+ a *Ghimj[187];
        a = - W_54/ Ghimj[194];
        W_54 = -a;
        W_56 = W_56+ a *Ghimj[195];
        W_63 = W_63+ a *Ghimj[196];
        W_67 = W_67+ a *Ghimj[197];
        W_68 = W_68+ a *Ghimj[198];
        W_69 = W_69+ a *Ghimj[199];
        a = - W_56/ Ghimj[207];
        W_56 = -a;
        W_59 = W_59+ a *Ghimj[208];
        W_63 = W_63+ a *Ghimj[209];
        W_67 = W_67+ a *Ghimj[210];
        W_68 = W_68+ a *Ghimj[211];
        W_69 = W_69+ a *Ghimj[212];
        W_72 = W_72+ a *Ghimj[213];
        a = - W_58/ Ghimj[232];
        W_58 = -a;
        W_62 = W_62+ a *Ghimj[233];
        W_64 = W_64+ a *Ghimj[234];
        W_65 = W_65+ a *Ghimj[235];
        W_67 = W_67+ a *Ghimj[236];
        W_68 = W_68+ a *Ghimj[237];
        W_71 = W_71+ a *Ghimj[238];
        a = - W_59/ Ghimj[240];
        W_59 = -a;
        W_61 = W_61+ a *Ghimj[241];
        W_63 = W_63+ a *Ghimj[242];
        W_65 = W_65+ a *Ghimj[243];
        W_67 = W_67+ a *Ghimj[244];
        W_68 = W_68+ a *Ghimj[245];
        W_69 = W_69+ a *Ghimj[246];
        W_72 = W_72+ a *Ghimj[247];
        a = - W_61/ Ghimj[258];
        W_61 = -a;
        W_63 = W_63+ a *Ghimj[259];
        W_68 = W_68+ a *Ghimj[260];
        W_69 = W_69+ a *Ghimj[261];
        W_72 = W_72+ a *Ghimj[262];
        a = - W_62/ Ghimj[267];
        W_62 = -a;
        W_63 = W_63+ a *Ghimj[268];
        W_64 = W_64+ a *Ghimj[269];
        W_65 = W_65+ a *Ghimj[270];
        W_67 = W_67+ a *Ghimj[271];
        W_68 = W_68+ a *Ghimj[272];
        W_69 = W_69+ a *Ghimj[273];
        W_70 = W_70+ a *Ghimj[274];
        W_72 = W_72+ a *Ghimj[275];
        a = - W_63/ Ghimj[278];
        W_63 = -a;
        W_67 = W_67+ a *Ghimj[279];
        W_68 = W_68+ a *Ghimj[280];
        W_69 = W_69+ a *Ghimj[281];
        W_72 = W_72+ a *Ghimj[282];
        a = - W_64/ Ghimj[286];
        W_64 = -a;
        W_65 = W_65+ a *Ghimj[287];
        W_67 = W_67+ a *Ghimj[288];
        W_68 = W_68+ a *Ghimj[289];
        W_69 = W_69+ a *Ghimj[290];
        W_72 = W_72+ a *Ghimj[291];
        a = - W_65/ Ghimj[314];
        W_65 = -a;
        W_66 = W_66+ a *Ghimj[315];
        W_67 = W_67+ a *Ghimj[316];
        W_68 = W_68+ a *Ghimj[317];
        W_69 = W_69+ a *Ghimj[318];
        W_70 = W_70+ a *Ghimj[319];
        W_71 = W_71+ a *Ghimj[320];
        W_72 = W_72+ a *Ghimj[321];
        a = - W_66/ Ghimj[334];
        W_66 = -a;
        W_67 = W_67+ a *Ghimj[335];
        W_68 = W_68+ a *Ghimj[336];
        W_69 = W_69+ a *Ghimj[337];
        W_70 = W_70+ a *Ghimj[338];
        W_71 = W_71+ a *Ghimj[339];
        W_72 = W_72+ a *Ghimj[340];
        Ghimj[341] = W_39;
        Ghimj[342] = W_42;
        Ghimj[343] = W_43;
        Ghimj[344] = W_45;
        Ghimj[345] = W_46;
        Ghimj[346] = W_48;
        Ghimj[347] = W_49;
        Ghimj[348] = W_51;
        Ghimj[349] = W_52;
        Ghimj[350] = W_54;
        Ghimj[351] = W_56;
        Ghimj[352] = W_58;
        Ghimj[353] = W_59;
        Ghimj[354] = W_61;
        Ghimj[355] = W_62;
        Ghimj[356] = W_63;
        Ghimj[357] = W_64;
        Ghimj[358] = W_65;
        Ghimj[359] = W_66;
        Ghimj[360] = W_67;
        Ghimj[361] = W_68;
        Ghimj[362] = W_69;
        Ghimj[363] = W_70;
        Ghimj[364] = W_71;
        Ghimj[365] = W_72;
        W_32 = Ghimj[366];
        W_33 = Ghimj[367];
        W_34 = Ghimj[368];
        W_35 = Ghimj[369];
        W_36 = Ghimj[370];
        W_37 = Ghimj[371];
        W_45 = Ghimj[372];
        W_46 = Ghimj[373];
        W_48 = Ghimj[374];
        W_49 = Ghimj[375];
        W_51 = Ghimj[376];
        W_54 = Ghimj[377];
        W_56 = Ghimj[378];
        W_57 = Ghimj[379];
        W_59 = Ghimj[380];
        W_61 = Ghimj[381];
        W_62 = Ghimj[382];
        W_63 = Ghimj[383];
        W_64 = Ghimj[384];
        W_65 = Ghimj[385];
        W_66 = Ghimj[386];
        W_67 = Ghimj[387];
        W_68 = Ghimj[388];
        W_69 = Ghimj[389];
        W_70 = Ghimj[390];
        W_71 = Ghimj[391];
        W_72 = Ghimj[392];
        a = - W_32/ Ghimj[119];
        W_32 = -a;
        W_65 = W_65+ a *Ghimj[120];
        a = - W_33/ Ghimj[121];
        W_33 = -a;
        W_65 = W_65+ a *Ghimj[122];
        a = - W_34/ Ghimj[123];
        W_34 = -a;
        W_65 = W_65+ a *Ghimj[124];
        a = - W_35/ Ghimj[125];
        W_35 = -a;
        W_65 = W_65+ a *Ghimj[126];
        a = - W_36/ Ghimj[127];
        W_36 = -a;
        W_65 = W_65+ a *Ghimj[128];
        a = - W_37/ Ghimj[129];
        W_37 = -a;
        W_65 = W_65+ a *Ghimj[130];
        a = - W_45/ Ghimj[149];
        W_45 = -a;
        W_67 = W_67+ a *Ghimj[150];
        W_68 = W_68+ a *Ghimj[151];
        a = - W_46/ Ghimj[152];
        W_46 = -a;
        W_67 = W_67+ a *Ghimj[153];
        W_68 = W_68+ a *Ghimj[154];
        a = - W_48/ Ghimj[171];
        W_48 = -a;
        W_69 = W_69+ a *Ghimj[172];
        W_72 = W_72+ a *Ghimj[173];
        a = - W_49/ Ghimj[174];
        W_49 = -a;
        W_67 = W_67+ a *Ghimj[175];
        W_68 = W_68+ a *Ghimj[176];
        a = - W_51/ Ghimj[181];
        W_51 = -a;
        W_67 = W_67+ a *Ghimj[182];
        W_68 = W_68+ a *Ghimj[183];
        a = - W_54/ Ghimj[194];
        W_54 = -a;
        W_56 = W_56+ a *Ghimj[195];
        W_63 = W_63+ a *Ghimj[196];
        W_67 = W_67+ a *Ghimj[197];
        W_68 = W_68+ a *Ghimj[198];
        W_69 = W_69+ a *Ghimj[199];
        a = - W_56/ Ghimj[207];
        W_56 = -a;
        W_59 = W_59+ a *Ghimj[208];
        W_63 = W_63+ a *Ghimj[209];
        W_67 = W_67+ a *Ghimj[210];
        W_68 = W_68+ a *Ghimj[211];
        W_69 = W_69+ a *Ghimj[212];
        W_72 = W_72+ a *Ghimj[213];
        a = - W_57/ Ghimj[218];
        W_57 = -a;
        W_62 = W_62+ a *Ghimj[219];
        W_64 = W_64+ a *Ghimj[220];
        W_65 = W_65+ a *Ghimj[221];
        W_67 = W_67+ a *Ghimj[222];
        W_68 = W_68+ a *Ghimj[223];
        W_71 = W_71+ a *Ghimj[224];
        a = - W_59/ Ghimj[240];
        W_59 = -a;
        W_61 = W_61+ a *Ghimj[241];
        W_63 = W_63+ a *Ghimj[242];
        W_65 = W_65+ a *Ghimj[243];
        W_67 = W_67+ a *Ghimj[244];
        W_68 = W_68+ a *Ghimj[245];
        W_69 = W_69+ a *Ghimj[246];
        W_72 = W_72+ a *Ghimj[247];
        a = - W_61/ Ghimj[258];
        W_61 = -a;
        W_63 = W_63+ a *Ghimj[259];
        W_68 = W_68+ a *Ghimj[260];
        W_69 = W_69+ a *Ghimj[261];
        W_72 = W_72+ a *Ghimj[262];
        a = - W_62/ Ghimj[267];
        W_62 = -a;
        W_63 = W_63+ a *Ghimj[268];
        W_64 = W_64+ a *Ghimj[269];
        W_65 = W_65+ a *Ghimj[270];
        W_67 = W_67+ a *Ghimj[271];
        W_68 = W_68+ a *Ghimj[272];
        W_69 = W_69+ a *Ghimj[273];
        W_70 = W_70+ a *Ghimj[274];
        W_72 = W_72+ a *Ghimj[275];
        a = - W_63/ Ghimj[278];
        W_63 = -a;
        W_67 = W_67+ a *Ghimj[279];
        W_68 = W_68+ a *Ghimj[280];
        W_69 = W_69+ a *Ghimj[281];
        W_72 = W_72+ a *Ghimj[282];
        a = - W_64/ Ghimj[286];
        W_64 = -a;
        W_65 = W_65+ a *Ghimj[287];
        W_67 = W_67+ a *Ghimj[288];
        W_68 = W_68+ a *Ghimj[289];
        W_69 = W_69+ a *Ghimj[290];
        W_72 = W_72+ a *Ghimj[291];
        a = - W_65/ Ghimj[314];
        W_65 = -a;
        W_66 = W_66+ a *Ghimj[315];
        W_67 = W_67+ a *Ghimj[316];
        W_68 = W_68+ a *Ghimj[317];
        W_69 = W_69+ a *Ghimj[318];
        W_70 = W_70+ a *Ghimj[319];
        W_71 = W_71+ a *Ghimj[320];
        W_72 = W_72+ a *Ghimj[321];
        a = - W_66/ Ghimj[334];
        W_66 = -a;
        W_67 = W_67+ a *Ghimj[335];
        W_68 = W_68+ a *Ghimj[336];
        W_69 = W_69+ a *Ghimj[337];
        W_70 = W_70+ a *Ghimj[338];
        W_71 = W_71+ a *Ghimj[339];
        W_72 = W_72+ a *Ghimj[340];
        a = - W_67/ Ghimj[360];
        W_67 = -a;
        W_68 = W_68+ a *Ghimj[361];
        W_69 = W_69+ a *Ghimj[362];
        W_70 = W_70+ a *Ghimj[363];
        W_71 = W_71+ a *Ghimj[364];
        W_72 = W_72+ a *Ghimj[365];
        Ghimj[366] = W_32;
        Ghimj[367] = W_33;
        Ghimj[368] = W_34;
        Ghimj[369] = W_35;
        Ghimj[370] = W_36;
        Ghimj[371] = W_37;
        Ghimj[372] = W_45;
        Ghimj[373] = W_46;
        Ghimj[374] = W_48;
        Ghimj[375] = W_49;
        Ghimj[376] = W_51;
        Ghimj[377] = W_54;
        Ghimj[378] = W_56;
        Ghimj[379] = W_57;
        Ghimj[380] = W_59;
        Ghimj[381] = W_61;
        Ghimj[382] = W_62;
        Ghimj[383] = W_63;
        Ghimj[384] = W_64;
        Ghimj[385] = W_65;
        Ghimj[386] = W_66;
        Ghimj[387] = W_67;
        Ghimj[388] = W_68;
        Ghimj[389] = W_69;
        Ghimj[390] = W_70;
        Ghimj[391] = W_71;
        Ghimj[392] = W_72;
        W_30 = Ghimj[393];
        W_48 = Ghimj[394];
        W_55 = Ghimj[395];
        W_61 = Ghimj[396];
        W_63 = Ghimj[397];
        W_66 = Ghimj[398];
        W_67 = Ghimj[399];
        W_68 = Ghimj[400];
        W_69 = Ghimj[401];
        W_70 = Ghimj[402];
        W_71 = Ghimj[403];
        W_72 = Ghimj[404];
        a = - W_30/ Ghimj[114];
        W_30 = -a;
        W_69 = W_69+ a *Ghimj[115];
        a = - W_48/ Ghimj[171];
        W_48 = -a;
        W_69 = W_69+ a *Ghimj[172];
        W_72 = W_72+ a *Ghimj[173];
        a = - W_55/ Ghimj[201];
        W_55 = -a;
        W_61 = W_61+ a *Ghimj[202];
        W_63 = W_63+ a *Ghimj[203];
        W_69 = W_69+ a *Ghimj[204];
        W_72 = W_72+ a *Ghimj[205];
        a = - W_61/ Ghimj[258];
        W_61 = -a;
        W_63 = W_63+ a *Ghimj[259];
        W_68 = W_68+ a *Ghimj[260];
        W_69 = W_69+ a *Ghimj[261];
        W_72 = W_72+ a *Ghimj[262];
        a = - W_63/ Ghimj[278];
        W_63 = -a;
        W_67 = W_67+ a *Ghimj[279];
        W_68 = W_68+ a *Ghimj[280];
        W_69 = W_69+ a *Ghimj[281];
        W_72 = W_72+ a *Ghimj[282];
        a = - W_66/ Ghimj[334];
        W_66 = -a;
        W_67 = W_67+ a *Ghimj[335];
        W_68 = W_68+ a *Ghimj[336];
        W_69 = W_69+ a *Ghimj[337];
        W_70 = W_70+ a *Ghimj[338];
        W_71 = W_71+ a *Ghimj[339];
        W_72 = W_72+ a *Ghimj[340];
        a = - W_67/ Ghimj[360];
        W_67 = -a;
        W_68 = W_68+ a *Ghimj[361];
        W_69 = W_69+ a *Ghimj[362];
        W_70 = W_70+ a *Ghimj[363];
        W_71 = W_71+ a *Ghimj[364];
        W_72 = W_72+ a *Ghimj[365];
        a = - W_68/ Ghimj[388];
        W_68 = -a;
        W_69 = W_69+ a *Ghimj[389];
        W_70 = W_70+ a *Ghimj[390];
        W_71 = W_71+ a *Ghimj[391];
        W_72 = W_72+ a *Ghimj[392];
        Ghimj[393] = W_30;
        Ghimj[394] = W_48;
        Ghimj[395] = W_55;
        Ghimj[396] = W_61;
        Ghimj[397] = W_63;
        Ghimj[398] = W_66;
        Ghimj[399] = W_67;
        Ghimj[400] = W_68;
        Ghimj[401] = W_69;
        Ghimj[402] = W_70;
        Ghimj[403] = W_71;
        Ghimj[404] = W_72;
        W_42 = Ghimj[405];
        W_43 = Ghimj[406];
        W_49 = Ghimj[407];
        W_50 = Ghimj[408];
        W_51 = Ghimj[409];
        W_60 = Ghimj[410];
        W_61 = Ghimj[411];
        W_63 = Ghimj[412];
        W_64 = Ghimj[413];
        W_65 = Ghimj[414];
        W_66 = Ghimj[415];
        W_67 = Ghimj[416];
        W_68 = Ghimj[417];
        W_69 = Ghimj[418];
        W_70 = Ghimj[419];
        W_71 = Ghimj[420];
        W_72 = Ghimj[421];
        a = - W_42/ Ghimj[141];
        W_42 = -a;
        W_67 = W_67+ a *Ghimj[142];
        a = - W_43/ Ghimj[143];
        W_43 = -a;
        W_67 = W_67+ a *Ghimj[144];
        a = - W_49/ Ghimj[174];
        W_49 = -a;
        W_67 = W_67+ a *Ghimj[175];
        W_68 = W_68+ a *Ghimj[176];
        a = - W_50/ Ghimj[177];
        W_50 = -a;
        W_64 = W_64+ a *Ghimj[178];
        W_65 = W_65+ a *Ghimj[179];
        W_72 = W_72+ a *Ghimj[180];
        a = - W_51/ Ghimj[181];
        W_51 = -a;
        W_67 = W_67+ a *Ghimj[182];
        W_68 = W_68+ a *Ghimj[183];
        a = - W_60/ Ghimj[249];
        W_60 = -a;
        W_61 = W_61+ a *Ghimj[250];
        W_63 = W_63+ a *Ghimj[251];
        W_64 = W_64+ a *Ghimj[252];
        W_65 = W_65+ a *Ghimj[253];
        W_66 = W_66+ a *Ghimj[254];
        W_67 = W_67+ a *Ghimj[255];
        W_68 = W_68+ a *Ghimj[256];
        a = - W_61/ Ghimj[258];
        W_61 = -a;
        W_63 = W_63+ a *Ghimj[259];
        W_68 = W_68+ a *Ghimj[260];
        W_69 = W_69+ a *Ghimj[261];
        W_72 = W_72+ a *Ghimj[262];
        a = - W_63/ Ghimj[278];
        W_63 = -a;
        W_67 = W_67+ a *Ghimj[279];
        W_68 = W_68+ a *Ghimj[280];
        W_69 = W_69+ a *Ghimj[281];
        W_72 = W_72+ a *Ghimj[282];
        a = - W_64/ Ghimj[286];
        W_64 = -a;
        W_65 = W_65+ a *Ghimj[287];
        W_67 = W_67+ a *Ghimj[288];
        W_68 = W_68+ a *Ghimj[289];
        W_69 = W_69+ a *Ghimj[290];
        W_72 = W_72+ a *Ghimj[291];
        a = - W_65/ Ghimj[314];
        W_65 = -a;
        W_66 = W_66+ a *Ghimj[315];
        W_67 = W_67+ a *Ghimj[316];
        W_68 = W_68+ a *Ghimj[317];
        W_69 = W_69+ a *Ghimj[318];
        W_70 = W_70+ a *Ghimj[319];
        W_71 = W_71+ a *Ghimj[320];
        W_72 = W_72+ a *Ghimj[321];
        a = - W_66/ Ghimj[334];
        W_66 = -a;
        W_67 = W_67+ a *Ghimj[335];
        W_68 = W_68+ a *Ghimj[336];
        W_69 = W_69+ a *Ghimj[337];
        W_70 = W_70+ a *Ghimj[338];
        W_71 = W_71+ a *Ghimj[339];
        W_72 = W_72+ a *Ghimj[340];
        a = - W_67/ Ghimj[360];
        W_67 = -a;
        W_68 = W_68+ a *Ghimj[361];
        W_69 = W_69+ a *Ghimj[362];
        W_70 = W_70+ a *Ghimj[363];
        W_71 = W_71+ a *Ghimj[364];
        W_72 = W_72+ a *Ghimj[365];
        a = - W_68/ Ghimj[388];
        W_68 = -a;
        W_69 = W_69+ a *Ghimj[389];
        W_70 = W_70+ a *Ghimj[390];
        W_71 = W_71+ a *Ghimj[391];
        W_72 = W_72+ a *Ghimj[392];
        a = - W_69/ Ghimj[401];
        W_69 = -a;
        W_70 = W_70+ a *Ghimj[402];
        W_71 = W_71+ a *Ghimj[403];
        W_72 = W_72+ a *Ghimj[404];
        Ghimj[405] = W_42;
        Ghimj[406] = W_43;
        Ghimj[407] = W_49;
        Ghimj[408] = W_50;
        Ghimj[409] = W_51;
        Ghimj[410] = W_60;
        Ghimj[411] = W_61;
        Ghimj[412] = W_63;
        Ghimj[413] = W_64;
        Ghimj[414] = W_65;
        Ghimj[415] = W_66;
        Ghimj[416] = W_67;
        Ghimj[417] = W_68;
        Ghimj[418] = W_69;
        Ghimj[419] = W_70;
        Ghimj[420] = W_71;
        Ghimj[421] = W_72;
        W_31 = Ghimj[422];
        W_38 = Ghimj[423];
        W_40 = Ghimj[424];
        W_44 = Ghimj[425];
        W_50 = Ghimj[426];
        W_52 = Ghimj[427];
        W_60 = Ghimj[428];
        W_61 = Ghimj[429];
        W_62 = Ghimj[430];
        W_63 = Ghimj[431];
        W_64 = Ghimj[432];
        W_65 = Ghimj[433];
        W_66 = Ghimj[434];
        W_67 = Ghimj[435];
        W_68 = Ghimj[436];
        W_69 = Ghimj[437];
        W_70 = Ghimj[438];
        W_71 = Ghimj[439];
        W_72 = Ghimj[440];
        a = - W_31/ Ghimj[116];
        W_31 = -a;
        W_40 = W_40+ a *Ghimj[117];
        W_65 = W_65+ a *Ghimj[118];
        a = - W_38/ Ghimj[131];
        W_38 = -a;
        W_60 = W_60+ a *Ghimj[132];
        W_66 = W_66+ a *Ghimj[133];
        W_71 = W_71+ a *Ghimj[134];
        a = - W_40/ Ghimj[137];
        W_40 = -a;
        W_65 = W_65+ a *Ghimj[138];
        a = - W_44/ Ghimj[145];
        W_44 = -a;
        W_50 = W_50+ a *Ghimj[146];
        W_60 = W_60+ a *Ghimj[147];
        W_65 = W_65+ a *Ghimj[148];
        a = - W_50/ Ghimj[177];
        W_50 = -a;
        W_64 = W_64+ a *Ghimj[178];
        W_65 = W_65+ a *Ghimj[179];
        W_72 = W_72+ a *Ghimj[180];
        a = - W_52/ Ghimj[184];
        W_52 = -a;
        W_65 = W_65+ a *Ghimj[185];
        W_67 = W_67+ a *Ghimj[186];
        W_71 = W_71+ a *Ghimj[187];
        a = - W_60/ Ghimj[249];
        W_60 = -a;
        W_61 = W_61+ a *Ghimj[250];
        W_63 = W_63+ a *Ghimj[251];
        W_64 = W_64+ a *Ghimj[252];
        W_65 = W_65+ a *Ghimj[253];
        W_66 = W_66+ a *Ghimj[254];
        W_67 = W_67+ a *Ghimj[255];
        W_68 = W_68+ a *Ghimj[256];
        a = - W_61/ Ghimj[258];
        W_61 = -a;
        W_63 = W_63+ a *Ghimj[259];
        W_68 = W_68+ a *Ghimj[260];
        W_69 = W_69+ a *Ghimj[261];
        W_72 = W_72+ a *Ghimj[262];
        a = - W_62/ Ghimj[267];
        W_62 = -a;
        W_63 = W_63+ a *Ghimj[268];
        W_64 = W_64+ a *Ghimj[269];
        W_65 = W_65+ a *Ghimj[270];
        W_67 = W_67+ a *Ghimj[271];
        W_68 = W_68+ a *Ghimj[272];
        W_69 = W_69+ a *Ghimj[273];
        W_70 = W_70+ a *Ghimj[274];
        W_72 = W_72+ a *Ghimj[275];
        a = - W_63/ Ghimj[278];
        W_63 = -a;
        W_67 = W_67+ a *Ghimj[279];
        W_68 = W_68+ a *Ghimj[280];
        W_69 = W_69+ a *Ghimj[281];
        W_72 = W_72+ a *Ghimj[282];
        a = - W_64/ Ghimj[286];
        W_64 = -a;
        W_65 = W_65+ a *Ghimj[287];
        W_67 = W_67+ a *Ghimj[288];
        W_68 = W_68+ a *Ghimj[289];
        W_69 = W_69+ a *Ghimj[290];
        W_72 = W_72+ a *Ghimj[291];
        a = - W_65/ Ghimj[314];
        W_65 = -a;
        W_66 = W_66+ a *Ghimj[315];
        W_67 = W_67+ a *Ghimj[316];
        W_68 = W_68+ a *Ghimj[317];
        W_69 = W_69+ a *Ghimj[318];
        W_70 = W_70+ a *Ghimj[319];
        W_71 = W_71+ a *Ghimj[320];
        W_72 = W_72+ a *Ghimj[321];
        a = - W_66/ Ghimj[334];
        W_66 = -a;
        W_67 = W_67+ a *Ghimj[335];
        W_68 = W_68+ a *Ghimj[336];
        W_69 = W_69+ a *Ghimj[337];
        W_70 = W_70+ a *Ghimj[338];
        W_71 = W_71+ a *Ghimj[339];
        W_72 = W_72+ a *Ghimj[340];
        a = - W_67/ Ghimj[360];
        W_67 = -a;
        W_68 = W_68+ a *Ghimj[361];
        W_69 = W_69+ a *Ghimj[362];
        W_70 = W_70+ a *Ghimj[363];
        W_71 = W_71+ a *Ghimj[364];
        W_72 = W_72+ a *Ghimj[365];
        a = - W_68/ Ghimj[388];
        W_68 = -a;
        W_69 = W_69+ a *Ghimj[389];
        W_70 = W_70+ a *Ghimj[390];
        W_71 = W_71+ a *Ghimj[391];
        W_72 = W_72+ a *Ghimj[392];
        a = - W_69/ Ghimj[401];
        W_69 = -a;
        W_70 = W_70+ a *Ghimj[402];
        W_71 = W_71+ a *Ghimj[403];
        W_72 = W_72+ a *Ghimj[404];
        a = - W_70/ Ghimj[419];
        W_70 = -a;
        W_71 = W_71+ a *Ghimj[420];
        W_72 = W_72+ a *Ghimj[421];
        Ghimj[422] = W_31;
        Ghimj[423] = W_38;
        Ghimj[424] = W_40;
        Ghimj[425] = W_44;
        Ghimj[426] = W_50;
        Ghimj[427] = W_52;
        Ghimj[428] = W_60;
        Ghimj[429] = W_61;
        Ghimj[430] = W_62;
        Ghimj[431] = W_63;
        Ghimj[432] = W_64;
        Ghimj[433] = W_65;
        Ghimj[434] = W_66;
        Ghimj[435] = W_67;
        Ghimj[436] = W_68;
        Ghimj[437] = W_69;
        Ghimj[438] = W_70;
        Ghimj[439] = W_71;
        Ghimj[440] = W_72;
        W_48 = Ghimj[441];
        W_55 = Ghimj[442];
        W_61 = Ghimj[443];
        W_63 = Ghimj[444];
        W_64 = Ghimj[445];
        W_65 = Ghimj[446];
        W_66 = Ghimj[447];
        W_67 = Ghimj[448];
        W_68 = Ghimj[449];
        W_69 = Ghimj[450];
        W_70 = Ghimj[451];
        W_71 = Ghimj[452];
        W_72 = Ghimj[453];
        a = - W_48/ Ghimj[171];
        W_48 = -a;
        W_69 = W_69+ a *Ghimj[172];
        W_72 = W_72+ a *Ghimj[173];
        a = - W_55/ Ghimj[201];
        W_55 = -a;
        W_61 = W_61+ a *Ghimj[202];
        W_63 = W_63+ a *Ghimj[203];
        W_69 = W_69+ a *Ghimj[204];
        W_72 = W_72+ a *Ghimj[205];
        a = - W_61/ Ghimj[258];
        W_61 = -a;
        W_63 = W_63+ a *Ghimj[259];
        W_68 = W_68+ a *Ghimj[260];
        W_69 = W_69+ a *Ghimj[261];
        W_72 = W_72+ a *Ghimj[262];
        a = - W_63/ Ghimj[278];
        W_63 = -a;
        W_67 = W_67+ a *Ghimj[279];
        W_68 = W_68+ a *Ghimj[280];
        W_69 = W_69+ a *Ghimj[281];
        W_72 = W_72+ a *Ghimj[282];
        a = - W_64/ Ghimj[286];
        W_64 = -a;
        W_65 = W_65+ a *Ghimj[287];
        W_67 = W_67+ a *Ghimj[288];
        W_68 = W_68+ a *Ghimj[289];
        W_69 = W_69+ a *Ghimj[290];
        W_72 = W_72+ a *Ghimj[291];
        a = - W_65/ Ghimj[314];
        W_65 = -a;
        W_66 = W_66+ a *Ghimj[315];
        W_67 = W_67+ a *Ghimj[316];
        W_68 = W_68+ a *Ghimj[317];
        W_69 = W_69+ a *Ghimj[318];
        W_70 = W_70+ a *Ghimj[319];
        W_71 = W_71+ a *Ghimj[320];
        W_72 = W_72+ a *Ghimj[321];
        a = - W_66/ Ghimj[334];
        W_66 = -a;
        W_67 = W_67+ a *Ghimj[335];
        W_68 = W_68+ a *Ghimj[336];
        W_69 = W_69+ a *Ghimj[337];
        W_70 = W_70+ a *Ghimj[338];
        W_71 = W_71+ a *Ghimj[339];
        W_72 = W_72+ a *Ghimj[340];
        a = - W_67/ Ghimj[360];
        W_67 = -a;
        W_68 = W_68+ a *Ghimj[361];
        W_69 = W_69+ a *Ghimj[362];
        W_70 = W_70+ a *Ghimj[363];
        W_71 = W_71+ a *Ghimj[364];
        W_72 = W_72+ a *Ghimj[365];
        a = - W_68/ Ghimj[388];
        W_68 = -a;
        W_69 = W_69+ a *Ghimj[389];
        W_70 = W_70+ a *Ghimj[390];
        W_71 = W_71+ a *Ghimj[391];
        W_72 = W_72+ a *Ghimj[392];
        a = - W_69/ Ghimj[401];
        W_69 = -a;
        W_70 = W_70+ a *Ghimj[402];
        W_71 = W_71+ a *Ghimj[403];
        W_72 = W_72+ a *Ghimj[404];
        a = - W_70/ Ghimj[419];
        W_70 = -a;
        W_71 = W_71+ a *Ghimj[420];
        W_72 = W_72+ a *Ghimj[421];
        a = - W_71/ Ghimj[439];
        W_71 = -a;
        W_72 = W_72+ a *Ghimj[440];
        Ghimj[441] = W_48;
        Ghimj[442] = W_55;
        Ghimj[443] = W_61;
        Ghimj[444] = W_63;
        Ghimj[445] = W_64;
        Ghimj[446] = W_65;
        Ghimj[447] = W_66;
        Ghimj[448] = W_67;
        Ghimj[449] = W_68;
        Ghimj[450] = W_69;
        Ghimj[451] = W_70;
        Ghimj[452] = W_71;
        Ghimj[453] = W_72;
}

__device__ void ros_Decomp(REAL * __restrict__ Ghimj, int &Ndec, int VL_GLO)
{
    kppDecomp(Ghimj, VL_GLO);
    Ndec++;
}

__device__ void ros_PrepareMatrix(REAL &H, int direction, REAL gam, REAL *jac0, REAL *Ghimj,  int &Nsng, int &Ndec, int VL_GLO)
{
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    int ising, nConsecutive;
    REAL ghinv;
    
        ghinv = ONE/(direction*H*gam);
        for (int i=0; i<LU_NONZERO; i++)
            Ghimj[i] = -jac0[i];

        Ghimj[0] += ghinv;
        Ghimj[1] += ghinv;
        Ghimj[4] += ghinv;
        Ghimj[8] += ghinv;
        Ghimj[12] += ghinv;
        Ghimj[15] += ghinv;
        Ghimj[18] += ghinv;
        Ghimj[24] += ghinv;
        Ghimj[27] += ghinv;
        Ghimj[33] += ghinv;
        Ghimj[35] += ghinv;
        Ghimj[36] += ghinv;
        Ghimj[40] += ghinv;
        Ghimj[43] += ghinv;
        Ghimj[46] += ghinv;
        Ghimj[48] += ghinv;
        Ghimj[51] += ghinv;
        Ghimj[55] += ghinv;
        Ghimj[68] += ghinv;
        Ghimj[73] += ghinv;
        Ghimj[78] += ghinv;
        Ghimj[82] += ghinv;
        Ghimj[86] += ghinv;
        Ghimj[89] += ghinv;
        Ghimj[92] += ghinv;
        Ghimj[95] += ghinv;
        Ghimj[98] += ghinv;
        Ghimj[103] += ghinv;
        Ghimj[110] += ghinv;
        Ghimj[112] += ghinv;
        Ghimj[114] += ghinv;
        Ghimj[116] += ghinv;
        Ghimj[119] += ghinv;
        Ghimj[121] += ghinv;
        Ghimj[123] += ghinv;
        Ghimj[125] += ghinv;
        Ghimj[127] += ghinv;
        Ghimj[129] += ghinv;
        Ghimj[131] += ghinv;
        Ghimj[135] += ghinv;
        Ghimj[137] += ghinv;
        Ghimj[139] += ghinv;
        Ghimj[141] += ghinv;
        Ghimj[143] += ghinv;
        Ghimj[145] += ghinv;
        Ghimj[149] += ghinv;
        Ghimj[152] += ghinv;
        Ghimj[165] += ghinv;
        Ghimj[171] += ghinv;
        Ghimj[174] += ghinv;
        Ghimj[177] += ghinv;
        Ghimj[181] += ghinv;
        Ghimj[184] += ghinv;
        Ghimj[188] += ghinv;
        Ghimj[194] += ghinv;
        Ghimj[201] += ghinv;
        Ghimj[207] += ghinv;
        Ghimj[218] += ghinv;
        Ghimj[232] += ghinv;
        Ghimj[240] += ghinv;
        Ghimj[249] += ghinv;
        Ghimj[258] += ghinv;
        Ghimj[267] += ghinv;
        Ghimj[278] += ghinv;
        Ghimj[286] += ghinv;
        Ghimj[314] += ghinv;
        Ghimj[334] += ghinv;
        Ghimj[360] += ghinv;
        Ghimj[388] += ghinv;
        Ghimj[401] += ghinv;
        Ghimj[419] += ghinv;
        Ghimj[439] += ghinv;
        Ghimj[453] += ghinv;
        Ghimj[454] += ghinv;
        ros_Decomp(Ghimj, Ndec, VL_GLO);
}

__device__ void Jac_sp(const REAL * __restrict__ var, const REAL * __restrict__ fix,
                 const REAL * __restrict__ rconst, REAL * __restrict__ jcb, int &Njac, const int VL_GLO)
{
    int index = blockIdx.x*blockDim.x+threadIdx.x;

 REAL dummy, B_0, B_1, B_2, B_3, B_4, B_5, B_6, B_7, B_8, B_9, B_10, B_11, B_12, B_13, B_14, B_15, B_16, B_17, B_18, B_19, B_20, B_21, B_22, B_23, B_24, B_25, B_26, B_27, B_28, B_29, B_30, B_31, B_32, B_33, B_34, B_35, B_36, B_37, B_38, B_39, B_40, B_41, B_42, B_43, B_44, B_45, B_46, B_47, B_48, B_49, B_50, B_51, B_52, B_53, B_54, B_55, B_56, B_57, B_58, B_59, B_60, B_61, B_62, B_63, B_64, B_65, B_66, B_67, B_68, B_69, B_70, B_71, B_72, B_73, B_74, B_75, B_76, B_77, B_78, B_79, B_80, B_81, B_82, B_83, B_84, B_85, B_86, B_87, B_88, B_89, B_90, B_91, B_92, B_93, B_94, B_95, B_96, B_97, B_98, B_99, B_100, B_101, B_102, B_103, B_104, B_105, B_106, B_107, B_108, B_109, B_110, B_111, B_112, B_113, B_114, B_115, B_116, B_117, B_118, B_119, B_120, B_121, B_122, B_123, B_124, B_125, B_126, B_127, B_128, B_129, B_130, B_131, B_132, B_133, B_134, B_135, B_136, B_137, B_138, B_139;


    Njac++;

        B_0 = rconst(index,0)*fix[0];
        B_2 = rconst(index,1)*fix[0];
        B_4 = rconst(index,2)*fix[0];
        B_6 = rconst(index,3)*var[66];
        B_7 = rconst(index,3)*var[65];
        B_8 = rconst(index,4)*var[65];
        B_9 = rconst(index,4)*var[40];
        B_10 = rconst(index,5)*var[71];
        B_11 = rconst(index,5)*var[66];
        B_12 = rconst(index,6)*var[71];
        B_13 = rconst(index,6)*var[65];
        B_14 = rconst(index,7)*2*var[71];
        B_15 = rconst(index,8)*var[47];
        B_16 = rconst(index,8)*var[41];
        B_17 = 1.8e-12*var[65];
        B_18 = 1.8e-12*var[52];
        B_19 = rconst(index,10)*2*var[47];
        B_20 = 1e+06;
        B_21 = rconst(index,12)*var[67];
        B_22 = rconst(index,12)*var[66];
        B_23 = rconst(index,13)*2*var[69];
        B_24 = rconst(index,14)*2*var[69];
        B_25 = rconst(index,15)*2*var[69];
        B_26 = rconst(index,16)*2*var[69];
        B_27 = rconst(index,17);
        B_28 = rconst(index,18)*var[67];
        B_29 = rconst(index,18)*var[52];
        B_30 = rconst(index,19)*var[71];
        B_31 = rconst(index,19)*var[69];
        B_32 = rconst(index,20)*var[65];
        B_33 = rconst(index,20)*var[58];
        B_34 = rconst(index,21)*var[69];
        B_35 = rconst(index,21)*var[48];
        B_36 = rconst(index,22)*var[69];
        B_37 = rconst(index,22)*var[55];
        B_38 = rconst(index,23);
        B_39 = rconst(index,24)*var[67];
        B_40 = rconst(index,24)*var[63];
        B_41 = rconst(index,25)*var[67];
        B_42 = rconst(index,25)*var[42];
        B_43 = rconst(index,26)*var[67];
        B_44 = rconst(index,26)*var[62];
        B_45 = 5.9e-11*var[67];
        B_46 = 5.9e-11*var[51];
        B_47 = rconst(index,28)*var[70];
        B_48 = rconst(index,28)*var[69];
        B_49 = rconst(index,29)*var[65];
        B_50 = rconst(index,29)*var[39];
        B_51 = rconst(index,30)*var[67];
        B_52 = rconst(index,30)*var[45];
        B_53 = 8e-11*var[67];
        B_54 = 8e-11*var[46];
        B_55 = rconst(index,32)*var[67];
        B_56 = rconst(index,32)*var[49];
        B_57 = rconst(index,33)*var[67];
        B_58 = rconst(index,33)*var[43];
        B_59 = rconst(index,34)*var[68];
        B_60 = rconst(index,34)*var[66];
        B_61 = 2.7e-12*2*var[72];
        B_62 = rconst(index,36)*2*var[72];
        B_63 = rconst(index,37)*var[71];
        B_64 = rconst(index,37)*var[68];
        B_65 = rconst(index,38)*var[72];
        B_66 = rconst(index,38)*var[71];
        B_67 = rconst(index,39)*var[65];
        B_68 = rconst(index,39)*var[57];
        B_69 = rconst(index,40)*var[65];
        B_70 = rconst(index,40)*var[59];
        B_71 = 4.9e-11*var[68];
        B_72 = 4.9e-11*var[61];
        B_73 = rconst(index,42)*var[72];
        B_74 = rconst(index,42)*var[48];
        B_75 = rconst(index,43)*var[72];
        B_76 = rconst(index,43)*var[55];
        B_77 = rconst(index,44);
        B_78 = rconst(index,45)*var[68];
        B_79 = rconst(index,45)*var[62];
        B_80 = rconst(index,46)*var[68];
        B_81 = rconst(index,46)*var[51];
        B_82 = rconst(index,47)*var[72];
        B_83 = rconst(index,47)*var[70];
        B_84 = rconst(index,48)*var[72];
        B_85 = rconst(index,48)*var[70];
        B_86 = rconst(index,49)*var[65];
        B_87 = rconst(index,49)*var[32];
        B_88 = rconst(index,50)*var[68];
        B_89 = rconst(index,50)*var[45];
        B_90 = rconst(index,51)*var[68];
        B_91 = rconst(index,51)*var[46];
        B_92 = rconst(index,52)*var[68];
        B_93 = rconst(index,52)*var[49];
        B_94 = rconst(index,53)*var[65];
        B_95 = rconst(index,53)*var[37];
        B_96 = rconst(index,54)*var[65];
        B_97 = rconst(index,54)*var[36];
        B_98 = 3.32e-15*var[68];
        B_99 = 3.32e-15*var[56];
        B_100 = 1.1e-15*var[68];
        B_101 = 1.1e-15*var[54];
        B_102 = rconst(index,57)*var[67];
        B_103 = rconst(index,57)*var[59];
        B_104 = rconst(index,58)*var[72];
        B_105 = rconst(index,58)*var[69];
        B_106 = rconst(index,59)*var[72];
        B_107 = rconst(index,59)*var[69];
        B_108 = rconst(index,60)*var[72];
        B_109 = rconst(index,60)*var[69];
        B_110 = 1.45e-11*var[67];
        B_111 = 1.45e-11*var[56];
        B_112 = rconst(index,62)*var[65];
        B_113 = rconst(index,62)*var[33];
        B_114 = rconst(index,63)*var[65];
        B_115 = rconst(index,63)*var[34];
        B_116 = rconst(index,64)*var[65];
        B_117 = rconst(index,64)*var[35];
        B_118 = rconst(index,65)*var[65];
        B_119 = rconst(index,65)*var[44];
        B_120 = rconst(index,66)*var[65];
        B_121 = rconst(index,66)*var[64];
        B_122 = rconst(index,67)*var[65];
        B_123 = rconst(index,67)*var[64];
        B_124 = rconst(index,68)*var[64];
        B_125 = rconst(index,68)*var[53];
        B_126 = 1e-10*var[65];
        B_127 = 1e-10*var[50];
        B_128 = rconst(index,70);
        B_129 = 3e-13*var[66];
        B_130 = 3e-13*var[60];
        B_131 = 5e-11*var[71];
        B_132 = 5e-11*var[38];
        B_133 = 3.3e-10*var[67];
        B_134 = 3.3e-10*var[64];
        B_135 = rconst(index,74)*var[68];
        B_136 = rconst(index,74)*var[64];
        B_137 = 4.4e-13*var[72];
        B_138 = 4.4e-13*var[64];
        B_139 = rconst(index,76);
        jcb[0] = - B_139;
        jcb[1] = 0;
        jcb[2] = B_124;
        jcb[3] = B_125;
        jcb[4] = 0;
        jcb[5] = B_43+ B_78;
        jcb[6] = B_44;
        jcb[7] = B_79;
        jcb[8] = 0;
        jcb[9] = B_53+ B_90;
        jcb[10] = B_54;
        jcb[11] = B_91;
        jcb[12] = 0;
        jcb[13] = B_30;
        jcb[14] = B_31;
        jcb[15] = 0;
        jcb[16] = B_25+ B_104;
        jcb[17] = B_105;
        jcb[18] = 0;
        jcb[19] = B_69;
        jcb[20] = B_70;
        jcb[21] = B_82;
        jcb[22] = B_65;
        jcb[23] = B_66+ B_83;
        jcb[24] = 0;
        jcb[25] = B_118;
        jcb[26] = B_119;
        jcb[27] = 0;
        jcb[28] = B_131;
        jcb[29] = 0.4*B_126;
        jcb[30] = 0.4*B_127;
        jcb[31] = B_132;
        jcb[32] = 2*B_139;
        jcb[33] = 0;
        jcb[34] = 2*B_139;
        jcb[35] = 0;
        jcb[36] = 0;
        jcb[37] = 0.666667*B_51+ 0.666667*B_88;
        jcb[38] = 0.666667*B_52;
        jcb[39] = 0.666667*B_89;
        jcb[40] = 0;
        jcb[41] = B_6;
        jcb[42] = B_7;
        jcb[43] = 0;
        jcb[44] = B_10;
        jcb[45] = B_11;
        jcb[46] = 0;
        jcb[47] = B_14;
        jcb[48] = 0;
        jcb[49] = B_10;
        jcb[50] = B_11;
        jcb[51] = 0;
        jcb[52] = B_15;
        jcb[53] = B_16;
        jcb[54] = 3*B_139;
        jcb[55] = 0;
        jcb[56] = B_15;
        jcb[57] = B_118;
        jcb[58] = B_16;
        jcb[59] = B_124;
        jcb[60] = B_129;
        jcb[61] = B_125+ B_137;
        jcb[62] = B_6+ B_119;
        jcb[63] = B_7+ B_10+ B_130;
        jcb[64] = 2*B_23+ 2*B_24+ B_25+ B_47+ B_104+ 2*B_106+ 2*B_108;
        jcb[65] = B_48+ B_84;
        jcb[66] = B_11;
        jcb[67] = 2*B_61+ 2*B_62+ B_85+ B_105+ 2*B_107+ 2*B_109+ B_138;
        jcb[68] = 0;
        jcb[69] = B_73;
        jcb[70] = B_106+ B_108;
        jcb[71] = B_84;
        jcb[72] = 2*B_61+ 2*B_62+ B_74+ B_85+ B_107+ B_109;
        jcb[73] = 0;
        jcb[74] = B_34;
        jcb[75] = B_35+ B_47+ B_106+ B_108;
        jcb[76] = B_48;
        jcb[77] = B_107+ B_109;
        jcb[78] = 0;
        jcb[79] = B_6;
        jcb[80] = B_7+ B_10;
        jcb[81] = B_11;
        jcb[82] = 0;
        jcb[83] = B_34+ B_73;
        jcb[84] = B_35;
        jcb[85] = B_74;
        jcb[86] = 0;
        jcb[87] = B_15;
        jcb[88] = B_16;
        jcb[89] = 0;
        jcb[90] = B_6;
        jcb[91] = B_7;
        jcb[92] = 0;
        jcb[93] = B_86;
        jcb[94] = B_87;
        jcb[95] = 0;
        jcb[96] = 3*B_49;
        jcb[97] = 3*B_50;
        jcb[98] = 0;
        jcb[99] = 0.6*B_126;
        jcb[100] = B_69;
        jcb[101] = B_128;
        jcb[102] = B_70+ 0.6*B_127;
        jcb[103] = 0;
        jcb[104] = B_112;
        jcb[105] = 2*B_114;
        jcb[106] = B_116;
        jcb[107] = 2*B_96;
        jcb[108] = 3*B_94;
        jcb[109] = 3*B_95+ 2*B_97+ B_113+ 2*B_115+ B_117;
        jcb[110] = - B_2;
        jcb[111] = B_0;
        jcb[112] = - B_20;
        jcb[113] = B_19;
        jcb[114] = - B_27;
        jcb[115] = B_26;
        jcb[116] = - B_4;
        jcb[117] = B_8;
        jcb[118] = B_9;
        jcb[119] = - B_86;
        jcb[120] = - B_87;
        jcb[121] = - B_112;
        jcb[122] = - B_113;
        jcb[123] = - B_114;
        jcb[124] = - B_115;
        jcb[125] = - B_116;
        jcb[126] = - B_117;
        jcb[127] = - B_96;
        jcb[128] = - B_97;
        jcb[129] = - B_94;
        jcb[130] = - B_95;
        jcb[131] = - B_131;
        jcb[132] = B_129;
        jcb[133] = B_130;
        jcb[134] = - B_132;
        jcb[135] = - B_49;
        jcb[136] = - B_50;
        jcb[137] = - B_8;
        jcb[138] = - B_9;
        jcb[139] = - B_0- B_15;
        jcb[140] = - B_16;
        jcb[141] = - B_41;
        jcb[142] = - B_42;
        jcb[143] = - B_57;
        jcb[144] = - B_58;
        jcb[145] = - B_118;
        jcb[146] = 0.6*B_126;
        jcb[147] = B_128;
        jcb[148] = - B_119+ 0.6*B_127;
        jcb[149] = - B_51- B_88;
        jcb[150] = - B_52;
        jcb[151] = - B_89;
        jcb[152] = - B_53- B_90;
        jcb[153] = - B_54;
        jcb[154] = - B_91;
        jcb[155] = 2*B_20;
        jcb[156] = B_86;
        jcb[157] = B_112;
        jcb[158] = B_114;
        jcb[159] = B_116;
        jcb[160] = B_96;
        jcb[161] = B_94;
        jcb[162] = B_49;
        jcb[163] = B_8;
        jcb[164] = - B_15;
        jcb[165] = - B_16- 2*B_19;
        jcb[166] = B_17;
        jcb[167] = B_67;
        jcb[168] = B_32;
        jcb[169] = B_9+ B_12+ B_18+ B_33+ B_50+ B_68+ B_87+ B_95+ B_97+ B_113+ B_115+ B_117;
        jcb[170] = B_13;
        jcb[171] = - B_34- B_73;
        jcb[172] = - B_35;
        jcb[173] = - B_74;
        jcb[174] = - B_55- B_92;
        jcb[175] = - B_56;
        jcb[176] = - B_93;
        jcb[177] = - B_126;
        jcb[178] = B_122+ B_137;
        jcb[179] = B_123- B_127;
        jcb[180] = B_138;
        jcb[181] = - B_45- B_80;
        jcb[182] = - B_46;
        jcb[183] = - B_81;
        jcb[184] = - B_17- B_28;
        jcb[185] = - B_18;
        jcb[186] = - B_29;
        jcb[187] = B_14;
        jcb[188] = - B_124;
        jcb[189] = B_71;
        jcb[190] = B_39;
        jcb[191] = - B_125;
        jcb[192] = B_40;
        jcb[193] = B_72;
        jcb[194] = - B_100;
        jcb[195] = B_110;
        jcb[196] = B_39;
        jcb[197] = B_40+ B_111;
        jcb[198] = - B_101;
        jcb[199] = B_23;
        jcb[200] = B_34+ B_73;
        jcb[201] = - B_36- B_75;
        jcb[202] = B_77;
        jcb[203] = B_38;
        jcb[204] = B_35- B_37;
        jcb[205] = B_74- B_76;
        jcb[206] = B_100;
        jcb[207] = - B_98- B_110;
        jcb[208] = B_102;
        jcb[209] = 0;
        jcb[210] = B_103- B_111;
        jcb[211] = - B_99+ B_101;
        jcb[212] = B_108;
        jcb[213] = B_109;
        jcb[214] = B_88;
        jcb[215] = B_90;
        jcb[216] = B_92;
        jcb[217] = B_80;
        jcb[218] = - B_67;
        jcb[219] = B_78;
        jcb[220] = B_135;
        jcb[221] = - B_68;
        jcb[222] = 0;
        jcb[223] = B_63+ B_79+ B_81+ B_89+ B_91+ B_93+ B_136;
        jcb[224] = B_64;
        jcb[225] = B_41;
        jcb[226] = B_57;
        jcb[227] = B_51;
        jcb[228] = B_53;
        jcb[229] = B_55;
        jcb[230] = B_45;
        jcb[231] = B_28;
        jcb[232] = - B_32;
        jcb[233] = B_43;
        jcb[234] = B_133;
        jcb[235] = - B_33;
        jcb[236] = B_29+ B_42+ B_44+ B_46+ B_52+ B_54+ B_56+ B_58+ B_134;
        jcb[237] = 0;
        jcb[238] = 0;
        jcb[239] = B_98;
        jcb[240] = - B_69- B_102;
        jcb[241] = B_71;
        jcb[242] = 0;
        jcb[243] = - B_70;
        jcb[244] = - B_103;
        jcb[245] = B_72+ B_99;
        jcb[246] = 0;
        jcb[247] = B_62;
        jcb[248] = B_124;
        jcb[249] = - B_128- B_129;
        jcb[250] = 0;
        jcb[251] = 0;
        jcb[252] = B_120+ B_125+ B_133+ B_135;
        jcb[253] = B_121;
        jcb[254] = - B_130;
        jcb[255] = B_134;
        jcb[256] = B_136;
        jcb[257] = B_75;
        jcb[258] = - B_71- B_77;
        jcb[259] = 0;
        jcb[260] = - B_72;
        jcb[261] = 0;
        jcb[262] = B_76;
        jcb[263] = B_126;
        jcb[264] = B_45;
        jcb[265] = B_124;
        jcb[266] = 0;
        jcb[267] = - B_43- B_78;
        jcb[268] = 0;
        jcb[269] = B_120+ B_125+ B_133+ B_135;
        jcb[270] = B_121+ B_127;
        jcb[271] = - B_44+ B_46+ B_134;
        jcb[272] = - B_79+ B_136;
        jcb[273] = B_47;
        jcb[274] = B_48+ B_82+ B_84;
        jcb[275] = B_83+ B_85;
        jcb[276] = B_36;
        jcb[277] = 0;
        jcb[278] = - B_38- B_39;
        jcb[279] = - B_40;
        jcb[280] = 0;
        jcb[281] = B_37;
        jcb[282] = 0;
        jcb[283] = - B_124;
        jcb[284] = 0;
        jcb[285] = 0;
        jcb[286] = - B_120- B_122- B_125- B_133- B_135- B_137;
        jcb[287] = - B_121- B_123;
        jcb[288] = - B_134;
        jcb[289] = - B_136;
        jcb[290] = 0;
        jcb[291] = - B_138;
        jcb[292] = - B_86;
        jcb[293] = - B_112;
        jcb[294] = - B_114;
        jcb[295] = - B_116;
        jcb[296] = - B_96;
        jcb[297] = - B_94;
        jcb[298] = - B_49;
        jcb[299] = - B_8;
        jcb[300] = 2*B_15;
        jcb[301] = - B_118;
        jcb[302] = 2*B_16;
        jcb[303] = - B_126;
        jcb[304] = B_45;
        jcb[305] = - B_17;
        jcb[306] = - B_67;
        jcb[307] = - B_32;
        jcb[308] = - B_69;
        jcb[309] = 0;
        jcb[310] = 0;
        jcb[311] = 0;
        jcb[312] = 0;
        jcb[313] = - B_120- B_122;
        jcb[314] = - B_6- B_9- B_12- B_18- B_33- B_50- B_68- B_70- B_87- B_95- B_97- B_113- B_115- B_117- B_119- B_121- B_123 - B_127;
        jcb[315] = - B_7+ B_10;
        jcb[316] = B_46;
        jcb[317] = 0;
        jcb[318] = 0;
        jcb[319] = 0;
        jcb[320] = B_11- B_13;
        jcb[321] = 0;
        jcb[322] = B_2;
        jcb[323] = 0;
        jcb[324] = 0;
        jcb[325] = 0;
        jcb[326] = 0;
        jcb[327] = 0;
        jcb[328] = - B_129;
        jcb[329] = 0;
        jcb[330] = 0;
        jcb[331] = 0;
        jcb[332] = 0;
        jcb[333] = - B_6;
        jcb[334] = - B_7- B_10- B_21- B_59- B_130;
        jcb[335] = - B_22;
        jcb[336] = - B_60;
        jcb[337] = 0;
        jcb[338] = 0;
        jcb[339] = - B_11;
        jcb[340] = 0;
        jcb[341] = 3*B_49;
        jcb[342] = - B_41;
        jcb[343] = - B_57;
        jcb[344] = - B_51;
        jcb[345] = - B_53;
        jcb[346] = B_34;
        jcb[347] = - B_55;
        jcb[348] = - B_45;
        jcb[349] = - B_28;
        jcb[350] = B_100;
        jcb[351] = B_98- B_110;
        jcb[352] = B_32;
        jcb[353] = - B_102;
        jcb[354] = 0;
        jcb[355] = - B_43;
        jcb[356] = - B_39;
        jcb[357] = - B_133;
        jcb[358] = B_33+ 3*B_50;
        jcb[359] = - B_21;
        jcb[360] = - B_22- B_29- B_40- B_42- B_44- B_46- B_52- B_54- B_56- B_58- B_103- B_111- B_134;
        jcb[361] = B_99+ B_101;
        jcb[362] = 2*B_24+ B_25+ B_35+ B_47+ B_106;
        jcb[363] = B_48;
        jcb[364] = 0;
        jcb[365] = B_107;
        jcb[366] = B_86;
        jcb[367] = B_112;
        jcb[368] = 2*B_114;
        jcb[369] = B_116;
        jcb[370] = 2*B_96;
        jcb[371] = 3*B_94;
        jcb[372] = - B_88;
        jcb[373] = - B_90;
        jcb[374] = B_73;
        jcb[375] = - B_92;
        jcb[376] = - B_80;
        jcb[377] = - B_100;
        jcb[378] = - B_98+ B_110;
        jcb[379] = B_67;
        jcb[380] = B_69+ B_102;
        jcb[381] = - B_71;
        jcb[382] = - B_78;
        jcb[383] = 0;
        jcb[384] = - B_135+ B_137;
        jcb[385] = B_68+ B_70+ B_87+ 3*B_95+ 2*B_97+ B_113+ 2*B_115+ B_117;
        jcb[386] = - B_59;
        jcb[387] = B_103+ B_111;
        jcb[388] = - B_60- B_63- B_72- B_79- B_81- B_89- B_91- B_93- B_99- B_101- B_136;
        jcb[389] = B_104+ B_106;
        jcb[390] = B_84;
        jcb[391] = - B_64;
        jcb[392] = 2*B_61+ B_74+ B_85+ B_105+ B_107+ B_138;
        jcb[393] = 2*B_27;
        jcb[394] = - B_34;
        jcb[395] = - B_36;
        jcb[396] = 0;
        jcb[397] = B_38;
        jcb[398] = B_21;
        jcb[399] = B_22;
        jcb[400] = 0;
        jcb[401] = - 2*B_23- 2*B_24- 2*B_25- 2*B_26- B_30- B_35- B_37- B_47- B_104- B_106- B_108;
        jcb[402] = - B_48;
        jcb[403] = - B_31;
        jcb[404] = - B_105- B_107- B_109;
        jcb[405] = B_41;
        jcb[406] = B_57;
        jcb[407] = B_55+ B_92;
        jcb[408] = 0.6*B_126;
        jcb[409] = B_80;
        jcb[410] = B_128;
        jcb[411] = 0;
        jcb[412] = 0;
        jcb[413] = 0;
        jcb[414] = 0.6*B_127;
        jcb[415] = 0;
        jcb[416] = B_42+ B_56+ B_58;
        jcb[417] = B_81+ B_93;
        jcb[418] = - B_47;
        jcb[419] = - B_48- B_82- B_84;
        jcb[420] = 0;
        jcb[421] = - B_83- B_85;
        jcb[422] = B_4;
        jcb[423] = - B_131;
        jcb[424] = 0;
        jcb[425] = B_118;
        jcb[426] = 0.4*B_126;
        jcb[427] = B_17+ B_28;
        jcb[428] = 0;
        jcb[429] = 0;
        jcb[430] = B_43+ B_78;
        jcb[431] = 0;
        jcb[432] = B_122;
        jcb[433] = B_6- B_12+ B_18+ B_119+ B_123+ 0.4*B_127;
        jcb[434] = B_7- B_10;
        jcb[435] = B_29+ B_44;
        jcb[436] = - B_63+ B_79;
        jcb[437] = - B_30+ B_47;
        jcb[438] = B_48+ B_84;
        jcb[439] = - B_11- B_13- 2*B_14- B_31- B_64- B_65- B_132;
        jcb[440] = - B_66+ B_85;
        jcb[441] = - B_73;
        jcb[442] = - B_75;
        jcb[443] = B_77;
        jcb[444] = 0;
        jcb[445] = - B_137;
        jcb[446] = 0;
        jcb[447] = B_59;
        jcb[448] = 0;
        jcb[449] = B_60;
        jcb[450] = - B_104- B_106- B_108;
        jcb[451] = - B_82- B_84;
        jcb[452] = - B_65;
        jcb[453] = - 2*B_61- 2*B_62- B_66- B_74- B_76- B_83- B_85- B_105- B_107- B_109- B_138;
    }

__device__ void Fun(REAL *var, const REAL * __restrict__ fix, const REAL * __restrict__ rconst, REAL *varDot, int &Nfun, const int VL_GLO){
    int index = blockIdx.x*blockDim.x+threadIdx.x;

    Nfun++;

 REAL dummy, A_0, A_1, A_2, A_3, A_4, A_5, A_6, A_7, A_8, A_9, A_10, A_11, A_12, A_13, A_14, A_15, A_16, A_17, A_18, A_19, A_20, A_21, A_22, A_23, A_24, A_25, A_26, A_27, A_28, A_29, A_30, A_31, A_32, A_33, A_34, A_35, A_36, A_37, A_38, A_39, A_40, A_41, A_42, A_43, A_44, A_45, A_46, A_47, A_48, A_49, A_50, A_51, A_52, A_53, A_54, A_55, A_56, A_57, A_58, A_59, A_60, A_61, A_62, A_63, A_64, A_65, A_66, A_67, A_68, A_69, A_70, A_71, A_72, A_73, A_74, A_75, A_76;

    {
        A_0 = rconst(index,0)*var[41]*fix(index,0);
        A_1 = rconst(index,1)*var[28]*fix(index,0);
        A_2 = rconst(index,2)*var[31]*fix(index,0);
        A_3 = rconst(index,3)*var[65]*var[66];
        A_4 = rconst(index,4)*var[40]*var[65];
        A_5 = rconst(index,5)*var[66]*var[71];
        A_6 = rconst(index,6)*var[65]*var[71];
        A_7 = rconst(index,7)*var[71]*var[71];
        A_8 = rconst(index,8)*var[41]*var[47];
        A_9 = 1.8e-12*var[52]*var[65];
        A_10 = rconst(index,10)*var[47]*var[47];
        A_11 = 1e+06*var[29];
        A_12 = rconst(index,12)*var[66]*var[67];
        A_13 = rconst(index,13)*var[69]*var[69];
        A_14 = rconst(index,14)*var[69]*var[69];
        A_15 = rconst(index,15)*var[69]*var[69];
        A_16 = rconst(index,16)*var[69]*var[69];
        A_17 = rconst(index,17)*var[30];
        A_18 = rconst(index,18)*var[52]*var[67];
        A_19 = rconst(index,19)*var[69]*var[71];
        A_20 = rconst(index,20)*var[58]*var[65];
        A_21 = rconst(index,21)*var[48]*var[69];
        A_22 = rconst(index,22)*var[55]*var[69];
        A_23 = rconst(index,23)*var[63];
        A_24 = rconst(index,24)*var[63]*var[67];
        A_25 = rconst(index,25)*var[42]*var[67];
        A_26 = rconst(index,26)*var[62]*var[67];
        A_27 = 5.9e-11*var[51]*var[67];
        A_28 = rconst(index,28)*var[69]*var[70];
        A_29 = rconst(index,29)*var[39]*var[65];
        A_30 = rconst(index,30)*var[45]*var[67];
        A_31 = 8e-11*var[46]*var[67];
        A_32 = rconst(index,32)*var[49]*var[67];
        A_33 = rconst(index,33)*var[43]*var[67];
        A_34 = rconst(index,34)*var[66]*var[68];
        A_35 = 2.7e-12*var[72]*var[72];
        A_36 = rconst(index,36)*var[72]*var[72];
        A_37 = rconst(index,37)*var[68]*var[71];
        A_38 = rconst(index,38)*var[71]*var[72];
        A_39 = rconst(index,39)*var[57]*var[65];
        A_40 = rconst(index,40)*var[59]*var[65];
        A_41 = 4.9e-11*var[61]*var[68];
        A_42 = rconst(index,42)*var[48]*var[72];
        A_43 = rconst(index,43)*var[55]*var[72];
        A_44 = rconst(index,44)*var[61];
        A_45 = rconst(index,45)*var[62]*var[68];
        A_46 = rconst(index,46)*var[51]*var[68];
        A_47 = rconst(index,47)*var[70]*var[72];
        A_48 = rconst(index,48)*var[70]*var[72];
        A_49 = rconst(index,49)*var[32]*var[65];
        A_50 = rconst(index,50)*var[45]*var[68];
        A_51 = rconst(index,51)*var[46]*var[68];
        A_52 = rconst(index,52)*var[49]*var[68];
        A_53 = rconst(index,53)*var[37]*var[65];
        A_54 = rconst(index,54)*var[36]*var[65];
        A_55 = 3.32e-15*var[56]*var[68];
        A_56 = 1.1e-15*var[54]*var[68];
        A_57 = rconst(index,57)*var[59]*var[67];
        A_58 = rconst(index,58)*var[69]*var[72];
        A_59 = rconst(index,59)*var[69]*var[72];
        A_60 = rconst(index,60)*var[69]*var[72];
        A_61 = 1.45e-11*var[56]*var[67];
        A_62 = rconst(index,62)*var[33]*var[65];
        A_63 = rconst(index,63)*var[34]*var[65];
        A_64 = rconst(index,64)*var[35]*var[65];
        A_65 = rconst(index,65)*var[44]*var[65];
        A_66 = rconst(index,66)*var[64]*var[65];
        A_67 = rconst(index,67)*var[64]*var[65];
        A_68 = rconst(index,68)*var[53]*var[64];
        A_69 = 1e-10*var[50]*var[65];
        A_70 = rconst(index,70)*var[60];
        A_71 = 3e-13*var[60]*var[66];
        A_72 = 5e-11*var[38]*var[71];
        A_73 = 3.3e-10*var[64]*var[67];
        A_74 = rconst(index,74)*var[64]*var[68];
        A_75 = 4.4e-13*var[64]*var[72];
        A_76 = rconst(index,76)*var[0];
        varDot[0] = - A_76;
        varDot[1] = A_68;
        varDot[2] = A_26+ A_45;
        varDot[3] = A_31+ A_51;
        varDot[4] = A_19;
        varDot[5] = A_15+ A_58;
        varDot[6] = A_38+ A_40+ A_47;
        varDot[7] = A_65;
        varDot[8] = 0.4*A_69+ A_72;
        varDot[9] = 2*A_76;
        varDot[10] = 2*A_76;
        varDot[11] = 0.666667*A_30+ 0.666667*A_50;
        varDot[12] = A_3;
        varDot[13] = A_5;
        varDot[14] = A_7;
        varDot[15] = A_5;
        varDot[16] = A_8;
        varDot[17] = A_3+ A_5+ A_8+ 2*A_13+ 2*A_14+ A_15+ A_28+ 2*A_35+ 2*A_36+ A_48+ A_58+ 2*A_59+ 2*A_60+ A_65+ A_68+ A_71 + A_75+ 3*A_76;
        varDot[18] = 2*A_35+ 2*A_36+ A_42+ A_48+ A_59+ A_60;
        varDot[19] = A_21+ A_28+ A_59+ A_60;
        varDot[20] = A_3+ A_5;
        varDot[21] = A_21+ A_42;
        varDot[22] = A_8;
        varDot[23] = A_3;
        varDot[24] = A_49;
        varDot[25] = 3*A_29;
        varDot[26] = A_40+ 0.6*A_69+ A_70;
        varDot[27] = 3*A_53+ 2*A_54+ A_62+ 2*A_63+ A_64;
        varDot[28] = A_0- A_1;
        varDot[29] = A_10- A_11;
        varDot[30] = A_16- A_17;
        varDot[31] = - A_2+ A_4;
        varDot[32] = - A_49;
        varDot[33] = - A_62;
        varDot[34] = - A_63;
        varDot[35] = - A_64;
        varDot[36] = - A_54;
        varDot[37] = - A_53;
        varDot[38] = A_71- A_72;
        varDot[39] = - A_29;
        varDot[40] = - A_4;
        varDot[41] = - A_0- A_8;
        varDot[42] = - A_25;
        varDot[43] = - A_33;
        varDot[44] = - A_65+ 0.6*A_69+ A_70;
        varDot[45] = - A_30- A_50;
        varDot[46] = - A_31- A_51;
        varDot[47] = A_4+ A_6- A_8+ A_9- 2*A_10+ 2*A_11+ A_20+ A_29+ A_39+ A_49+ A_53+ A_54+ A_62+ A_63+ A_64;
        varDot[48] = - A_21- A_42;
        varDot[49] = - A_32- A_52;
        varDot[50] = A_67- A_69+ A_75;
        varDot[51] = - A_27- A_46;
        varDot[52] = A_7- A_9- A_18;
        varDot[53] = A_24+ A_41- A_68;
        varDot[54] = A_13+ A_24- A_56+ A_61;
        varDot[55] = A_21- A_22+ A_23+ A_42- A_43+ A_44;
        varDot[56] = - A_55+ A_56+ A_57+ A_60- A_61;
        varDot[57] = A_37- A_39+ A_45+ A_46+ A_50+ A_51+ A_52+ A_74;
        varDot[58] = A_18- A_20+ A_25+ A_26+ A_27+ A_30+ A_31+ A_32+ A_33+ A_73;
        varDot[59] = A_36- A_40+ A_41+ A_55- A_57;
        varDot[60] = A_66+ A_68- A_70- A_71+ A_73+ A_74;
        varDot[61] = - A_41+ A_43- A_44;
        varDot[62] = - A_26+ A_27+ A_28- A_45+ A_47+ A_48+ A_66+ A_68+ A_69+ A_73+ A_74;
        varDot[63] = A_22- A_23- A_24;
        varDot[64] = - A_66- A_67- A_68- A_73- A_74- A_75;
        varDot[65] = - A_3- A_4+ A_5- A_6+ 2*A_8- A_9- A_20+ A_27- A_29- A_39- A_40- A_49- A_53- A_54- A_62- A_63- A_64- A_65 - A_66- A_67- A_69;
        varDot[66] = A_1- A_3- A_5- A_12- A_34- A_71;
        varDot[67] = - A_12+ 2*A_14+ A_15- A_18+ A_20+ A_21- A_24- A_25- A_26- A_27+ A_28+ 3*A_29- A_30- A_31- A_32- A_33+ A_55 + A_56- A_57+ A_59- A_61- A_73;
        varDot[68] = - A_34+ 2*A_35- A_37+ A_39+ A_40- A_41+ A_42- A_45- A_46+ A_48+ A_49- A_50- A_51- A_52+ 3*A_53+ 2*A_54 - A_55- A_56+ A_57+ A_58+ A_59+ A_61+ A_62+ 2*A_63+ A_64- A_74+ A_75;
        varDot[69] = A_12- 2*A_13- 2*A_14- 2*A_15- 2*A_16+ 2*A_17- A_19- A_21- A_22+ A_23- A_28- A_58- A_59- A_60;
        varDot[70] = A_25- A_28+ A_32+ A_33+ A_46- A_47- A_48+ A_52+ 0.6*A_69+ A_70;
        varDot[71] = A_2+ A_3- A_5- A_6- 2*A_7+ A_9+ A_18- A_19+ A_26+ A_28- A_37- A_38+ A_45+ A_48+ A_65+ A_67+ 0.4*A_69 - A_72;
        varDot[72] = A_34- 2*A_35- 2*A_36- A_38- A_42- A_43+ A_44- A_47- A_48- A_58- A_59- A_60- A_75;
    }
}

__device__ void ros_FunTimeDerivative(const REAL T, REAL roundoff, REAL * __restrict__ var, const REAL * __restrict__ fix,
                                      const REAL * __restrict__ rconst, REAL *dFdT, REAL *Fcn0, int &Nfun,
                                      const REAL * __restrict__ khet_st, const REAL * __restrict__ khet_tr,
                                      const REAL * __restrict__ jx,
                                      const int VL_GLO)
{
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    const REAL DELTAMIN = 1.0E-6;
    REAL delta,one_over_delta;

    delta = sqrt(roundoff)*fmax(DELTAMIN,fabs(T));
    one_over_delta = 1.0/delta;

    Fun(var, fix, rconst, dFdT, Nfun, VL_GLO);

    for (int i=0; i < NVAR; i++){
        dFdT(index,i) = (dFdT(index,i) - Fcn0(index,i)) * one_over_delta;
    }
}


__device__  static  int ros_Integrator(REAL * __restrict__ var, const REAL * __restrict__ fix, const REAL Tstart, const REAL Tend, REAL &T,
        //  Rosenbrock method coefficients
        const int ros_S, const REAL * __restrict__ ros_M, const REAL * __restrict__ ros_E, const REAL * __restrict__ ros_A, const REAL * __restrict__  ros_C,
        const REAL * __restrict__ ros_Alpha, const REAL * __restrict__ ros_Gamma, const REAL ros_ELO, const int * ros_NewF,
        //  Integration parameters
        const int autonomous, const int vectorTol, const int Max_no_steps,
        const REAL roundoff, const REAL Hmin, const REAL Hmax, const REAL Hstart, REAL &Hexit,
        const REAL FacMin, const REAL FacMax, const REAL FacRej, const REAL FacSafe,
        //  Status parameters
        int &Nfun, int &Njac, int &Nstp, int &Nacc, int &Nrej, int &Ndec, int &Nsol, int &Nsng,
        //  cuda global mem buffers              
        const REAL * __restrict__ rconst,  const REAL * __restrict__ absTol, const REAL * __restrict__ relTol, REAL * __restrict__ varNew, REAL * __restrict__ Fcn0,
        REAL * __restrict__ K, REAL * __restrict__ dFdT, REAL * __restrict__ jac0, REAL * __restrict__ Ghimj, REAL * __restrict__ varErr,
        // for update_rconst
        const REAL * __restrict__ khet_st, const REAL * __restrict__ khet_tr,
        const REAL * __restrict__ jx,
        // VL_GLO
        const int VL_GLO)
{
    int index = blockIdx.x*blockDim.x+threadIdx.x;

    REAL H, Hnew, HC, HG, Fac; // Tau - not used
    REAL Err; //*varErr;
    int direction;
    int rejectLastH, rejectMoreH;
    const REAL DELTAMIN = 1.0E-5;

    //   ~~~>  Initial preparations
    T = Tstart;
    Hexit = 0.0;
    H = fmin(Hstart,Hmax);
    if (fabs(H) <= 10.0*roundoff)
        H = DELTAMIN;

    if (Tend  >=  Tstart)
    {
        direction = + 1;
    }
    else
    {
        direction = - 1;
    }

    rejectLastH=0;
    rejectMoreH=0;
    //   ~~~> Time loop begins below

    // TimeLoop: 
    while((direction > 0) && ((T- Tend)+ roundoff <= ZERO) || (direction < 0) && ((Tend-T)+ roundoff <= ZERO))
    {
        if (Nstp > Max_no_steps) //  Too many steps
            return -6;
        //  Step size too small
        if (H <= roundoff){  //  Step size too small
            //if (((T+ 0.1*H) == T) || (H <= roundoff)) {
            return -7;
        }

        //   ~~~>  Limit H if necessary to avoid going beyond Tend
        Hexit = H;
        H = fmin(H,fabs(Tend-T));

        //   ~~~>   Compute the function at current time
        Fun(var, fix, rconst, Fcn0, Nfun, VL_GLO);	/// VAR READ - Fcn0 Write

        //   ~~~>  Compute the function derivative with respect to T
        if (!autonomous)
            ros_FunTimeDerivative(T, roundoff, var, fix, rconst, dFdT, Fcn0, Nfun, khet_st, khet_tr, jx,  VL_GLO); /// VAR READ - fcn0 read

        //   ~~~>   Compute the Jacobian at current time
        Jac_sp(var, fix, rconst, jac0, Njac, VL_GLO);   /// VAR READ 

        //   ~~~>  Repeat step calculation until current step accepted
        // UntilAccepted: 
        while(1)
        {
            ros_PrepareMatrix(H, direction, ros_Gamma[0], jac0, Ghimj, Nsng, Ndec, VL_GLO);
            //   ~~~>   Compute the stages
            // Stage: 
            for (int istage=0; istage < ros_S; istage++)
            {
                //   For the 1st istage the function has been computed previously
                if (istage == 0)
                {
                    for (int i=0; i<NVAR; i++){
                        varNew(index,i) = Fcn0(index,i);				// FCN0 Read
                    }
                }
                else if(ros_NewF[istage])
                {
                        for (int i=0; i<NVAR; i++){		
                            varNew(index,i) = var(index,i);
                        }

                    for (int j=0; j < (istage); j++){
                        for (int i=0; i<NVAR; i++){		
                            varNew(index,i) = K(index,j,i)*ros_A[(istage)*(istage-1)/2 + j]  + varNew(index,i);
                        }
                    }
                    Fun(varNew, fix, rconst, varNew, Nfun,VL_GLO); // FCN <- varNew / not overlap 
		} 

		for (int i=0; i<NVAR; i++)		
			K(index,istage,i)  = varNew(index,i);

		for (int j=0; j<(istage); j++)
		{
			HC = ros_C[(istage)*(istage-1)/2 + j]/(direction*H);
			for (int i=0; i<NVAR; i++){
				REAL tmp = K(index,j,i);
				K(index,istage,i) += tmp*HC;
			}
		}

                if ((!autonomous) && (ros_Gamma[istage] ))
                {
                    HG = direction*H*ros_Gamma[istage];
                    for (int i=0; i<NVAR; i++){
                        K(index,istage,i) += dFdT(index,i)*HG;
		     }
                }
		//	   R   ,RW, RW,  R,        R 
                ros_Solve(Ghimj, K, Nsol, istage, ros_S);


            } // Stage

            //  ~~~>  Compute the new solution
	    for (int i=0; i<NVAR; i++){
		    REAL tmpNew  = var(index,i); 					/// VAR READ
		    REAL tmpErr  = ZERO;

		    for (int j=0; j<ros_S; j++){
		    	    REAL tmp = K(index,j,i);

#ifdef DEBUG
			    if (isnan(tmp)){
			    	printf("Solver detected NAN!");
			    	tmp = 0;
			    }
#endif
			    tmpNew += tmp*ros_M[j];
			    tmpErr += tmp*ros_E[j];
		    }
		    varNew(index,i) = tmpNew;			// varNew is killed
		    varErr(index,i) = tmpErr;
	    }

            Err = ros_ErrorNorm(var, varNew, varErr, absTol, relTol, vectorTol);   /// VAR-varNew READ


//  ~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
            Fac  = fmin(FacMax,fmax(FacMin,FacSafe/pow(Err,ONE/ros_ELO)));
            Hnew = H*Fac;

//  ~~~>  Check the error magnitude and adjust step size
            Nstp = Nstp+ 1;
            if((Err <= ONE) || (H <= Hmin)) // ~~~> Accept step
            {
                Nacc = Nacc + 1;
                for (int j=0; j<NVAR ; j++)
                    var(index,j) =  fmax(varNew(index,j),ZERO);  /////////// VAR WRITE - last VarNew read

                T = T +  direction*H;
                Hnew = fmax(Hmin,fmin(Hnew,Hmax));
                if (rejectLastH)   // No step size increase after a rejected step
                    Hnew = fmin(Hnew,H);
                rejectLastH = 0;
                rejectMoreH = 0;
                H = Hnew;

            	break;  //  EXIT THE LOOP: WHILE STEP NOT ACCEPTED
            }
            else      // ~~~> Reject step
            {
                if (rejectMoreH)
                    Hnew = H*FacRej;
                rejectMoreH = rejectLastH;
                rejectLastH = 1;
                H = Hnew;
                if (Nacc >= 1)
                    Nrej += 1;
            } //  Err <= 1
        } // UntilAccepted
    } // TimeLoop
//  ~~~> Succesful exit
    return 0; //  ~~~> The integration was successful
}

typedef struct {
 REAL ros_A[15];
 REAL ros_C[15];
 int   ros_NewF[8];
 REAL ros_M[6];
 REAL ros_E[6];
 REAL ros_Alpha[6];
 REAL ros_Gamma[6];
 REAL ros_ELO;
 int    ros_S;
} ros_t;

/*
 * Lookup tables for different ROS for branch elimination. It is much faster in GPU.
 */
__device__ __constant__  ros_t ros[5] = {
    {       
        {.58578643762690495119831127579030,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* ros_A */
        {-1.17157287525380990239662255158060,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* ros_C */
        {1,1,0,0,0,0,0,0}, /* ros_NewF */
        {.87867965644035742679746691368545,.29289321881345247559915563789515,0,0,0,0}, /* ros_M */
        {.29289321881345247559915563789515,.29289321881345247559915563789515,0,0,0,0}, /* ros_E */
        {0,1.0,0,0,0,0}, /* ros_Alpha */
        {1.70710678118654752440084436210485,-1.70710678118654752440084436210485,0,0,0,0},  /* ros_Gamma */
        2.0, /* ros_ELO */
        2, /* ros_S*/
    }, /* Ros2 */
    {       
        {1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* ros_A */
        {-0.10156171083877702091975600115545E+01, 0.40759956452537699824805835358067E+01,0.92076794298330791242156818474003E+01,0,0,0,0,0,0,0,0,0,0,0,0}, /* ros_C */
        {1,1,0,0,0,0,0,0}, /* ros_NewF */
        {0.1E+01,0.61697947043828245592553615689730E+01,-0.42772256543218573326238373806514E+00,0,0,0}, /* ros_M */
        {0.5E+00,- 0.29079558716805469821718236208017E+01,0.22354069897811569627360909276199E+00,0,0,0}, /* ros_E */
        {0.0E+00,0.43586652150845899941601945119356E+00,0.43586652150845899941601945119356E+00,0,0,0}, /* ros_Alpha */
        {0.43586652150845899941601945119356E+00,0.24291996454816804366592249683314E+00,0.21851380027664058511513169485832E+01,0,0,0},  /* ros_Gamma */
        3.0, /* ros_ELO */
        3
    }, /* Ros3 */
    {       
        {0.2000000000000000E+01, 0.1867943637803922E+01, 0.2344449711399156E+00, 0.1867943637803922E+01, 0.2344449711399156E+00,0,0,0,0,0,0,0,0,0,0}, /* ros_A */
        {-0.7137615036412310E+01,0.2580708087951457E+01,0.6515950076447975E+00, - 0.2137148994382534E+01, - 0.3214669691237626E+00, - 0.6949742501781779E+00 ,0,0,0,0,0,0,0,0,0}, /* ros_C */
        {1,1,1,0,0,0,0,0}, /* ros_NewF */
        {0.2255570073418735E+01, 0.2870493262186792E+00, 0.4353179431840180E+00, 0.1093502252409163E+01,0,0}, /* ros_M */
        { -0.2815431932141155E+00, -0.7276199124938920E-01, -0.1082196201495311E+00, -0.1093502252409163E+01, 0, 0}, /* ros_E */
        {0.0, 0.1145640000000000E+01, 0.6552168638155900E+00, 0.6552168638155900E+00,0,0}, /* ros_Alpha */
        { 0.5728200000000000E+00, -0.1769193891319233E+01, 0.7592633437920482E+00, -0.1049021087100450E+00,0,0},  /* ros_Gamma */
        4.0, /* ros_ELO */
        4
    }, /* Ros4 */
    {       
        { 0.0E+00, 2.0E+00, 0.0E+00, 2.0E+00, 0.0E+00, 1.0E+00, 0,0,0,0,0,0,0,0,0}, /* ros_A */
        { 4.0E+00, 1.0E+00, - 1.0E+00,  1.0E+00, - 1.0E+00, - 2.66666666666666666666666666666666, 0,0,0,0,0,0,0,0,0}, /* ros_C */
        {1,0,1,1,0,0,0,0}, /* ros_NewF */
        {2.0,0,1.0,1.0,0,0}, /* ros_M */
        {0,0,0,1.0,0,0}, /* ros_E */
        {0,0,1.0,1.0,0,0}, /* ros_Alpha */
        {0.5,1.5,0,0,0,0},  /* ros_Gamma */
        3.0, /* ros_ELO */
        4
    }, /* Rodas3 */

    { 
        {
            0.1544000000000000E+01,  0.9466785280815826E+00, 0.2557011698983284E+00, 0.3314825187068521E+01,
            0.2896124015972201E+01,  0.9986419139977817E+00, 0.1221224509226641E+01, 0.6019134481288629E+01,
            0.1253708332932087E+02, -0.6878860361058950E+00, 0.1221224509226641E+01, 0.6019134481288629E+01,
            0.1253708332932087E+02, -0.6878860361058950E+00, 1.0E+00},  /* ros_A */ 

        {
            -0.5668800000000000E+01, -0.2430093356833875E+01, -0.2063599157091915E+00, -0.1073529058151375E+00,  
            -0.9594562251023355E+01, -0.2047028614809616E+02,  0.7496443313967647E+01, -0.1024680431464352E+02,  
            -0.3399990352819905E+02,  0.1170890893206160E+02,  0.8083246795921522E+01, -0.7981132988064893E+01,  
            -0.3152159432874371E+02,  0.1631930543123136E+02, -0.6058818238834054E+01}, /* ros_C */
        {1,1,1,1,1,1,0,0}, /* ros_NewF */
        {0.1221224509226641E+01,0.6019134481288629E+01,0.1253708332932087E+02,- 0.6878860361058950E+00,1,1}, /* ros_M */
        {0,0,0,0,0,1.0}, /* ros_E */
        {0.000,  0.386,  0.210,  0.630,  1.000, 1.000}, /* ros_Alpha */
        {0.2500000000000000E+00,  -0.1043000000000000E+00,  0.1035000000000000E+00,  0.3620000000000023E-01, 0, 0},  /* ros_Gamma */
        4.0, /* ros_ELO */
        6
    } /* Rodas4 */



};



//__device__ double rconst_local[MAX_VL_GLO*NREACT];

/* Initialize rconst local  */
//__device__ double * rconst_local;


__device__ REAL k_3rd(REAL temp, REAL cair, REAL k0_300K, REAL n, REAL kinf_300K, REAL m, REAL fc)
    /*
 *    
 * temp        temperature [K]
 * cair        air concentration [molecules/cm3]
 * k0_300K     low pressure limit at 300 K
 * n           exponent for low pressure limit
 * kinf_300K   high pressure limit at 300 K
 * m           exponent for high pressure limit
 * fc          broadening factor (usually fc=0.6)
 * 
 */
{

    REAL zt_help, k0_T, kinf_T, k_ratio, k_3rd_r;

    zt_help = 300.0/temp;
    k0_T    = k0_300K   *pow(zt_help,n) *cair;
    kinf_T  = kinf_300K *pow(zt_help,m);
    k_ratio = k0_T/kinf_T;
    k_3rd_r   = k0_T/(1.0+ k_ratio)*pow(fc,1.0/(1.0+ pow(log10(k_ratio),2)));
    return k_3rd_r;
}


__device__ REAL k_3rd_iupac(REAL temp, REAL cair, REAL k0_300K, REAL n, REAL kinf_300K, REAL m, REAL fc)
/*
 *    
 * temp        temperature [K]
 * cair        air concentration [molecules/cm3]
 * k0_300K     low pressure limit at 300 K
 * n           exponent for low pressure limit
 * kinf_300K   high pressure limit at 300 K
 * m           exponent for high pressure limit
 * fc          broadening factor (e.g. 0.45 or 0.6...)
 * nu          N
 * 
 */
{

    REAL zt_help, k0_T, kinf_T, k_ratio, nu, k_3rd_iupac_r;
    zt_help = 300.0/temp;
    k0_T    = k0_300K   *pow(zt_help,n) *cair;
    kinf_T  = kinf_300K *pow(zt_help,m);
    k_ratio = k0_T/kinf_T;
    nu      = 0.75- 1.27*log10(fc);
    k_3rd_iupac_r = k0_T/(1.0+ k_ratio)*pow(fc,1.0/(1.0+ pow(log10(k_ratio)/nu,2)));
    return k_3rd_iupac_r;
}



double * temp_gpu;
double * press_gpu;
double * cair_gpu;

#if defined(__SINGLEPREC)
float * temp_gpu_s;
float * press_gpu_s;
float * cair_gpu_s;

__global__ void doubleToFloat(float *out, double* in, int n)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    for (; i < n; i += gridDim.x * blockDim.x)
        out[i] = in[i];
}

__global__ void floatToDouble(double *out, float* in, int n)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    for (; i < n; i += gridDim.x * blockDim.x)
        out[i] = in[i];
}
#endif

__device__ void  update_rconst(const REAL * __restrict__ var, 
 			       const REAL * __restrict__ khet_st, const REAL * __restrict__ khet_tr,
 			       const REAL * __restrict__ jx, REAL * __restrict__ rconst, 
			       const REAL * __restrict__ temp_gpu, 
			       const REAL * __restrict__ press_gpu, 
			       const REAL * __restrict__ cair_gpu, 
			       const int VL_GLO)
{
    int index = blockIdx.x*blockDim.x+threadIdx.x;

    /* Set local buffer */

    {
        const REAL temp_loc  = temp_gpu[index];
        const REAL press_loc = press_gpu[index];
        const REAL cair_loc  = cair_gpu[index];

        REAL alpha_NO_HO2, beta_NO_HO2, k0_NO_HO2, k2d_NO_HO2, k1d_NO_HO2, k2w_NO_HO2, k1w_NO_HO2, k_PrO2_HO2, k_PrO2_NO, k_PrO2_CH3O2, k_HO2_HO2, k_NO3_NO2, k_NO2_HO2, k_HNO3_OH, k_CH3O2, k_CH3OOH_OH, k_CH3CO3_NO2, k_PAN_M, k_ClO_ClO, k_BrO_NO2, k_I_NO2, k_DMS_OH, G7402a_yield, k_O3s, beta_null_CH3NO3, beta_inf_CH3NO3, beta_CH3NO3, k_NO2_CH3O2, k_G4138, k_G9408, KRO2NO, KRO2HO2, KAPHO2, KAPNO, KRO2NO3, KNO3AL, J_IC3H7NO3, J_ACETOL, RO2;

        alpha_NO_HO2 = var[ind_H2O] *6.6E-27 *temp_loc *exp(3700. / temp_loc);
        beta_NO_HO2 = max(((530. / temp_loc)+ (press_loc *4.8004E-6)- 1.73) *0.01 , 0.);
        k0_NO_HO2 = 3.5E-12 *exp(250. / temp_loc);
        k2d_NO_HO2 = (beta_NO_HO2 *k0_NO_HO2) / (1.+ beta_NO_HO2);
        k1d_NO_HO2 = k0_NO_HO2 - k2d_NO_HO2;
        k2w_NO_HO2 = (beta_NO_HO2 *k0_NO_HO2 *(1.+ 42. *alpha_NO_HO2))/ ((1.+ alpha_NO_HO2) *(1.+ beta_NO_HO2));
        k1w_NO_HO2 = k0_NO_HO2 - k2w_NO_HO2;
        k_PrO2_HO2 = 1.9E-13 *exp(1300. / temp_loc);
        k_PrO2_NO = 2.7E-12 *exp(360. / temp_loc);
        k_PrO2_CH3O2 = 9.46E-14 *exp(431. / temp_loc);
        k_HO2_HO2 = (1.5E-12 *exp(19. / temp_loc)+ 1.7E-33 *exp(1000. / temp_loc) *cair_loc) * (1.+ 1.4E-21 *exp(2200. / temp_loc) *var[ind_H2O]);
        k_NO3_NO2 = k_3rd(temp_loc , cair_loc , 2.0E-30 , 4.4 , 1.4E-12 , 0.7 , 0.6);
        k_NO2_HO2 = k_3rd(temp_loc , cair_loc , 2.0E-31 , 3.4 , 2.9E-12 , 1.1 , 0.6);
        k_HNO3_OH = 2.4E-14 *exp(460. / temp_loc) + 1. / (1. / (6.5E-34 *exp(1335. / temp_loc) *cair_loc) + 1. / (2.7E-17 *exp(2199. / temp_loc)));
        k_CH3O2 = 1.03E-13 *exp(365. / temp_loc);
        k_CH3OOH_OH = 5.3E-12 *exp(190. / temp_loc);
        k_CH3CO3_NO2 = k_3rd(temp_loc , cair_loc , 9.7E-29 , 5.6 , 9.3E-12 , 1.5 , 0.6);
        k_PAN_M = k_CH3CO3_NO2 / (9.0E-29 *exp(14000. / temp_loc));
        k_ClO_ClO = k_3rd_iupac(temp_loc , cair_loc , 2.0E-32 , 4.0 , 1.0E-11 , 0.0 , 0.45);
        k_BrO_NO2 = k_3rd_iupac(temp_loc , cair_loc , 4.7E-31 , 3.1 , 1.8E-11 , 0.0 , 0.4);
        k_I_NO2 = k_3rd_iupac(temp_loc , cair_loc , 3.0E-31 , 1.0 , 6.6E-11 , 0.0 , 0.63);
        k_DMS_OH = 1.E-9 *exp(5820. / temp_loc) *var[ind_O2] / (1.E30+ 5. *exp(6280. / temp_loc) *var[ind_O2]);
        G7402a_yield = 0.8 / 1.1;
        k_O3s = (1.7E-12 *exp(- 940. / temp_loc)) *var[ind_OH] + (1.E-14 *exp(- 490. / temp_loc)) *var[ind_HO2] + jx(index,ip_O1D) *2.2E-10 *var[ind_H2O] / (3.2E-11 *exp(70. / temp_loc) *var[ind_O2] + 1.8E-11 *exp(110. / temp_loc) *var[ind_N2] + 2.2E-10 *var[ind_H2O]);
        beta_null_CH3NO3 = 0.00295 + 5.15E-22 *cair_loc * pow(temp_loc / 298, 7.4);
        beta_inf_CH3NO3 = 0.022;
        beta_CH3NO3 = (beta_null_CH3NO3 *beta_inf_CH3NO3) / (beta_null_CH3NO3 + beta_inf_CH3NO3);
        k_NO2_CH3O2 = k_3rd(temp_loc , cair_loc , 1.0E-30 , 4.8 , 7.2E-12 , 2.1 , 0.6);
        k_G4138 = 4.25E-12;
        k_G9408 = 3.66E-11;
        KRO2NO = 2.54E-12 *exp(360. / temp_loc);
        KRO2HO2 = 2.91E-13 *exp(1300. / temp_loc);
        KAPHO2 = 4.30E-13 *exp(1040. / temp_loc);
        KAPNO = 8.10E-12 *exp(270. / temp_loc);
        KRO2NO3 = 2.50E-12;
        KNO3AL = 1.4E-12 *exp(- 1900. / temp_loc);
        J_IC3H7NO3 = 3.7 *jx(index,ip_PAN);
        J_ACETOL = 0.65 *0.11 *jx(index,ip_CHOH);
        RO2 = 0.;
        if (ind_LISOPACO2>0) RO2 = RO2 + var[ind_LISOPACO2];
        if (ind_ISOPBO2>0) RO2 = RO2 + var[ind_ISOPBO2];
        if (ind_ISOPDO2>0) RO2 = RO2 + var[ind_ISOPDO2];
        if (ind_NISOPO2>0) RO2 = RO2 + var[ind_NISOPO2];
        if (ind_LHC4ACCO3>0) RO2 = RO2 + var[ind_LHC4ACCO3];
        if (ind_LC578O2>0) RO2 = RO2 + var[ind_LC578O2];
        if (ind_C59O2>0) RO2 = RO2 + var[ind_C59O2];
        if (ind_LNISO3>0) RO2 = RO2 + var[ind_LNISO3];
        if (ind_CH3O2>0) RO2 = RO2 + var[ind_CH3O2];
        if (ind_HOCH2O2>0) RO2 = RO2 + var[ind_HOCH2O2];
        if (ind_CH3CO3>0) RO2 = RO2 + var[ind_CH3CO3];
        if (ind_C2H5O2>0) RO2 = RO2 + var[ind_C2H5O2];
        if (ind_HOCH2CO3>0) RO2 = RO2 + var[ind_HOCH2CO3];
        if (ind_HYPROPO2>0) RO2 = RO2 + var[ind_HYPROPO2];
        if (ind_HCOCO3>0) RO2 = RO2 + var[ind_HCOCO3];
        if (ind_CO2H3CO3>0) RO2 = RO2 + var[ind_CO2H3CO3];
        if (ind_LHMVKABO2>0) RO2 = RO2 + var[ind_LHMVKABO2];
        if (ind_MACO3>0) RO2 = RO2 + var[ind_MACO3];
        if (ind_MACRO2>0) RO2 = RO2 + var[ind_MACRO2];
        if (ind_LMVKOHABO2>0) RO2 = RO2 + var[ind_LMVKOHABO2];
        if (ind_PRONO3BO2>0) RO2 = RO2 + var[ind_PRONO3BO2];
        if (ind_HOCH2CH2O2>0) RO2 = RO2 + var[ind_HOCH2CH2O2];
        if (ind_CH3COCH2O2>0) RO2 = RO2 + var[ind_CH3COCH2O2];
        if (ind_IC3H7O2>0) RO2 = RO2 + var[ind_IC3H7O2];
        if (ind_LC4H9O2>0) RO2 = RO2 + var[ind_LC4H9O2];
        if (ind_LMEKO2>0) RO2 = RO2 + var[ind_LMEKO2];
        if (ind_LAPINABO2> 0) RO2 = RO2 + var[ind_LAPINABO2];
        if (ind_C96O2> 0) RO2 = RO2 + var[ind_C96O2];
        if (ind_C97O2> 0) RO2 = RO2 + var[ind_C97O2];
        if (ind_C98O2> 0) RO2 = RO2 + var[ind_C98O2];
        if (ind_C85O2> 0) RO2 = RO2 + var[ind_C85O2];
        if (ind_C86O2> 0) RO2 = RO2 + var[ind_C86O2];
        if (ind_PINALO2> 0) RO2 = RO2 + var[ind_PINALO2];
        if (ind_C96CO3> 0) RO2 = RO2 + var[ind_C96CO3];
        if (ind_C89CO3> 0) RO2 = RO2 + var[ind_C89CO3];
        if (ind_C85CO3> 0) RO2 = RO2 + var[ind_C85CO3];
        if (ind_ROO6R1O2> 0) RO2 = RO2 + var[ind_ROO6R1O2];
        if (ind_RO6R1O2> 0) RO2 = RO2 + var[ind_RO6R1O2];
        if (ind_OHMENTHEN6ONEO2> 0) RO2 = RO2 + var[ind_OHMENTHEN6ONEO2];
        if (ind_C511O2> 0) RO2 = RO2 + var[ind_C511O2];
        if (ind_C106O2> 0) RO2 = RO2 + var[ind_C106O2];
        if (ind_CO235C6CO3> 0) RO2 = RO2 + var[ind_CO235C6CO3];
        if (ind_CHOC3COCO3> 0) RO2 = RO2 + var[ind_CHOC3COCO3];
        if (ind_CO235C6O2> 0) RO2 = RO2 + var[ind_CO235C6O2];
        if (ind_C716O2> 0) RO2 = RO2 + var[ind_C716O2];
        if (ind_C614O2> 0) RO2 = RO2 + var[ind_C614O2];
        if (ind_HCOCH2CO3> 0) RO2 = RO2 + var[ind_HCOCH2CO3];
        if (ind_BIACETO2> 0) RO2 = RO2 + var[ind_BIACETO2];
        if (ind_CO23C4CO3> 0) RO2 = RO2 + var[ind_CO23C4CO3];
        if (ind_C109O2> 0) RO2 = RO2 + var[ind_C109O2];
        if (ind_C811CO3> 0) RO2 = RO2 + var[ind_C811CO3];
        if (ind_C89O2> 0) RO2 = RO2 + var[ind_C89O2];
        if (ind_C812O2> 0) RO2 = RO2 + var[ind_C812O2];
        if (ind_C813O2> 0) RO2 = RO2 + var[ind_C813O2];
        if (ind_C721CO3> 0) RO2 = RO2 + var[ind_C721CO3];
        if (ind_C721O2> 0) RO2 = RO2 + var[ind_C721O2];
        if (ind_C722O2> 0) RO2 = RO2 + var[ind_C722O2];
        if (ind_C44O2> 0) RO2 = RO2 + var[ind_C44O2];
        if (ind_C512O2> 0) RO2 = RO2 + var[ind_C512O2];
        if (ind_C513O2> 0) RO2 = RO2 + var[ind_C513O2];
        if (ind_CHOC3COO2> 0) RO2 = RO2 + var[ind_CHOC3COO2];
        if (ind_C312COCO3> 0) RO2 = RO2 + var[ind_C312COCO3];
        if (ind_HOC2H4CO3> 0) RO2 = RO2 + var[ind_HOC2H4CO3];
        if (ind_LNAPINABO2> 0) RO2 = RO2 + var[ind_LNAPINABO2];
        if (ind_C810O2> 0) RO2 = RO2 + var[ind_C810O2];
        if (ind_C514O2> 0) RO2 = RO2 + var[ind_C514O2];
        if (ind_CHOCOCH2O2> 0) RO2 = RO2 + var[ind_CHOCOCH2O2];
        if (ind_ROO6R1O2> 0) RO2 = RO2 + var[ind_ROO6R1O2];
        if (ind_ROO6R3O2> 0) RO2 = RO2 + var[ind_ROO6R3O2];
        if (ind_RO6R1O2> 0) RO2 = RO2 + var[ind_RO6R1O2];
        if (ind_RO6R3O2> 0) RO2 = RO2 + var[ind_RO6R3O2];
        if (ind_BPINAO2> 0) RO2 = RO2 + var[ind_BPINAO2];
        if (ind_C8BCO2> 0) RO2 = RO2 + var[ind_C8BCO2];
        if (ind_NOPINDO2> 0) RO2 = RO2 + var[ind_NOPINDO2];
        if (ind_LNBPINABO2> 0) RO2 = RO2 + var[ind_LNBPINABO2];

        rconst(index,0) = (3.3E-11 *exp(55. / temp_loc));
        rconst(index,1) = (6.E-34 *( pow(temp_loc / 300., - 2.4) )*cair_loc);
        rconst(index,2) = (k_3rd(temp_loc , cair_loc , 4.4E-32 , 1.3 , 7.5E-11 , - 0.2 , 0.6));
        rconst(index,3) = (1.7E-12 *exp(- 940. / temp_loc));
        rconst(index,4) = (2.8E-12 *exp(- 1800. / temp_loc));
        rconst(index,5) = (1.E-14 *exp(- 490. / temp_loc));
        rconst(index,6) = (4.8E-11 *exp(250. / temp_loc));
        rconst(index,7) = (k_HO2_HO2);
        rconst(index,8) = (1.63E-10 *exp(60. / temp_loc));
        rconst(index,10) = (6.521E-26 *temp_loc *exp(1851.09 / temp_loc) *exp(- 5.10485E-3 *temp_loc) *1.E6);
        rconst(index,12) = (2.8E-11 *exp(- 250. / temp_loc));
        rconst(index,13) = (1.0E-12 *exp(- 1590. / temp_loc));
        rconst(index,14) = (3.0E-11 *exp(- 2450. / temp_loc));
        rconst(index,15) = (3.5E-13 *exp(- 1370. / temp_loc));
        rconst(index,16) = (k_ClO_ClO);
        rconst(index,17) = (k_ClO_ClO / (1.72E-27 *exp(8649. / temp_loc)));
        rconst(index,18) = (1.1E-11 *exp(- 980. / temp_loc));
        rconst(index,19) = (2.2E-12 *exp(340. / temp_loc));
        rconst(index,20) = (1.7E-12 *exp(- 230. / temp_loc));
        rconst(index,21) = (6.2E-12 *exp(295. / temp_loc));
        rconst(index,22) = (k_3rd_iupac(temp_loc , cair_loc , 1.6E-31 , 3.4 , 7.E-11 , 0. , 0.4));
        rconst(index,23) = (6.918E-7 *exp(- 10909. / temp_loc) *cair_loc);
        rconst(index,24) = (6.2E-12 *exp(145. / temp_loc));
        rconst(index,25) = (6.6E-12 *exp(- 1240. / temp_loc));
        rconst(index,26) = (8.1E-11 *exp(- 34. / temp_loc));
        rconst(index,28) = (3.3E-12 *exp(- 115. / temp_loc));
        rconst(index,29) = (1.64E-12 *exp(- 1520. / temp_loc));
        rconst(index,30) = (k_3rd_iupac(temp_loc , cair_loc , 1.85E-29 , 3.3 , 6.0E-10 , 0.0 , 0.4));
        rconst(index,32) = (k_3rd_iupac(temp_loc , cair_loc , 6.1e-30 , 3.0 , 2.0e-10 , 0. , 0.6));
        rconst(index,33) = (8.3E-11 *exp(- 100. / temp_loc));
        rconst(index,34) = (1.7E-11 *exp(- 800. / temp_loc));
        rconst(index,36) = (2.9E-14 *exp(840. / temp_loc));
        rconst(index,37) = (7.7E-12 *exp(- 450. / temp_loc));
        rconst(index,38) = (4.5E-12 *exp(500. / temp_loc));
        rconst(index,39) = (6.7E-12 *exp(155. / temp_loc));
        rconst(index,40) = (2.0E-11 *exp(240. / temp_loc));
        rconst(index,42) = (8.7E-12 *exp(260. / temp_loc));
        rconst(index,43) = (k_BrO_NO2);
        rconst(index,44) = (k_BrO_NO2 / (5.44E-9 *exp(14192. / temp_loc) *1.E6 *R_gas *temp_loc / (atm2Pa *N_A)));
        rconst(index,45) = (7.7E-12 *exp(- 580. / temp_loc));
        rconst(index,46) = (2.6E-12 *exp(- 1600. / temp_loc));
        rconst(index,47) = (G7402a_yield *5.7E-12);
        rconst(index,48) = ((1.- G7402a_yield) *5.7E-12);
        rconst(index,49) = (2.35E-12 *exp(- 1300. / temp_loc));
        rconst(index,50) = (2.8E-13 *exp(224. / temp_loc) / (1.+ 1.13E24 *exp(- 3200. / temp_loc) / var[ind_O2]));
        rconst(index,51) = (1.8e-11 *exp(- 460. / temp_loc));
        rconst(index,52) = (6.35e-15 *exp(440. / temp_loc));
        rconst(index,53) = (1.35E-12 *exp(- 600. / temp_loc));
        rconst(index,54) = (2.0E-12 *exp(- 840. / temp_loc));
        rconst(index,57) = (2.3E-10 *exp(135. / temp_loc));
        rconst(index,58) = (1.6E-12 *exp(430. / temp_loc));
        rconst(index,59) = (2.9E-12 *exp(220. / temp_loc));
        rconst(index,60) = (5.8E-13 *exp(170. / temp_loc));
        rconst(index,62) = (2.0E-12 *exp(- 840. / temp_loc));
        rconst(index,63) = (2.0E-12 *exp(- 840. / temp_loc));
        rconst(index,64) = (2.4E-12 *exp(- 920. / temp_loc));
        rconst(index,65) = (k_3rd(temp_loc , cair_loc , 3.3E-31 , 4.3 , 1.6E-12 , 0. , 0.6));
        rconst(index,66) = (1.13E-11 *exp(- 253. / temp_loc));
        rconst(index,67) = (k_DMS_OH);
        rconst(index,68) = (1.9E-13 *exp(520. / temp_loc));
        rconst(index,70) = (1.8E13 *exp(- 8661. / temp_loc));
        rconst(index,74) = (9.E-11 *exp(- 2386. / temp_loc));
        rconst(index,76) = (khet_tr(index,iht_N2O5));
        rconst(index,(10)-1) = 1.8e-12;
        rconst(index,(12)-1) = 1e+06;
        rconst(index,(28)-1) = 5.9e-11;
        rconst(index,(32)-1) = 8e-11;
        rconst(index,(36)-1) = 2.7e-12;
        rconst(index,(42)-1) = 4.9e-11;
        rconst(index,(56)-1) = 3.32e-15;
        rconst(index,(57)-1) = 1.1e-15;
        rconst(index,(62)-1) = 1.45e-11;
        rconst(index,(70)-1) = 1e-10;
        rconst(index,(72)-1) = 3e-13;
        rconst(index,(73)-1) = 5e-11;
        rconst(index,(74)-1) = 3.3e-10;
        rconst(index,(76)-1) = 4.4e-13;
    }
}

__global__
void Rosenbrock(REAL * __restrict__ conc, const REAL Tstart, const REAL Tend, REAL * __restrict__ rstatus, int * __restrict__ istatus,
                // values calculated from icntrl and rcntrl at host
                const int autonomous, const int vectorTol, const int UplimTol, const int method, const int Max_no_steps,
                REAL * __restrict__ d_jac0, REAL * __restrict__ d_Ghimj, REAL * __restrict__ d_varNew, REAL * __restrict__ d_K, REAL * __restrict__ d_varErr,REAL * __restrict__ d_dFdT ,REAL * __restrict__ d_Fcn0, REAL * __restrict__ d_var, REAL * __restrict__ d_fix, REAL * __restrict__ d_rconst,
                const REAL Hmin, const REAL Hmax, const REAL Hstart, const REAL FacMin, const REAL FacMax, const REAL FacRej, const REAL FacSafe, const REAL roundoff,
                // cuda global mem buffers              
                const REAL * __restrict__ absTol, const REAL * __restrict__ relTol,
                // for update_rconst
                const REAL * __restrict__ khet_st, const REAL * __restrict__ khet_tr,
                const REAL * __restrict__ jx,
                // global input
                const REAL * __restrict__ temp_gpu,
                const REAL * __restrict__ press_gpu,
                const REAL * __restrict__ cair_gpu,
                // extra
                const int VL_GLO)
{
    int index = blockIdx.x*blockDim.x+threadIdx.x;

    /* 
     *  In theory someone can aggregate accesses together,
     *  however due to algorithm, threads access 
     *  different parts of memory, making it harder to
     *  optimize accesses. 
     *
     */
    REAL *Ghimj  = &d_Ghimj[index*LU_NONZERO];
    REAL *K      = &d_K[index*NVAR*6];
    REAL *varNew = &d_varNew[index*NVAR];
    REAL *Fcn0   = &d_Fcn0[index*NVAR];
    REAL *dFdT   = &d_dFdT[index*NVAR];
    REAL *jac0   = &d_jac0[index*LU_NONZERO];
    REAL *varErr = &d_varErr[index*NVAR];
    REAL *var    = &d_var[index*NSPEC];
    REAL *fix    = &d_fix[index*NFIX];
    REAL *rconst = &d_rconst[index*NREACT];



    if (index < VL_GLO)
    {

        int Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng;
        REAL Texit, Hexit;

        Nfun = 0;
        Njac = 0;
        Nstp = 0;
        Nacc = 0;
        Nrej = 0;
        Ndec = 0;
        Nsol = 0;
        Nsng = 0;

        /* FIXME: add check for method */
        const REAL *ros_A     = &ros[method-1].ros_A[0];
        const REAL *ros_C     = &ros[method-1].ros_C[0];
        const REAL *ros_M     = &ros[method-1].ros_M[0];
        const REAL *ros_E     = &ros[method-1].ros_E[0];
        const REAL *ros_Alpha = &ros[method-1].ros_Alpha[0];
        const REAL *ros_Gamma = &ros[method-1].ros_Gamma[0];
        const int    *ros_NewF  = &ros[method-1].ros_NewF[0];
        const int     ros_S     =  ros[method-1].ros_S;
        const REAL  ros_ELO   =  ros[method-1].ros_ELO;



        /* Copy data from global memory to temporary array */
        /*
         * Optimization note: if we ever have enough constant
         * memory, we could use it for storing the data.
         * In current architectures if we use constant memory
         * only a few threads will be able to run on the fly.
         *
         */
        for (int i=0; i<NSPEC; i++)
            var(index,i) = conc(index,i);

        for (int i=0; i<NFIX; i++)
            fix(index,i) = conc(index,NVAR+i);


        update_rconst(var, khet_st, khet_tr, jx, rconst, temp_gpu, press_gpu, cair_gpu, VL_GLO); 

        ros_Integrator(var, fix, Tstart, Tend, Texit,
                //  Rosenbrock method coefficients
                ros_S, ros_M, ros_E, ros_A, ros_C, 
                ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, 
                //  Integration parameters
                autonomous, vectorTol, Max_no_steps, 
                roundoff, Hmin, Hmax, Hstart, Hexit, 
                FacMin, FacMax, FacRej, FacSafe,
                //  Status parameters
                Nfun, Njac, Nstp, Nacc, Nrej, Ndec, Nsol, Nsng,
                //  cuda global mem buffers              
                rconst, absTol, relTol, varNew, Fcn0,  
                K, dFdT, jac0, Ghimj,  varErr, 
                // For update rconst
                khet_st, khet_tr, jx,
                VL_GLO
                );

        for (int i=0; i<NVAR; i++)
            conc(index,i) = var(index,i); 


        /* Statistics */
        istatus(index,ifun) = Nfun;
        istatus(index,ijac) = Njac;
        istatus(index,istp) = Nstp;
        istatus(index,iacc) = Nacc;
        istatus(index,irej) = Nrej;
        istatus(index,idec) = Ndec;
        istatus(index,isol) = Nsol;
        istatus(index,isng) = Nsng;
        // Last T and H
        rstatus(index,itexit) = Texit;
        rstatus(index,ihexit) = Hexit; 
    }
}




                                                        // no int8 in CUDA :(
__global__ void reduce_istatus_1(int *istatus, int4 *tmp_out_1, int4 *tmp_out_2, int VL_GLO, int *xNacc, int *xNrej)
{
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    int idx_1 = threadIdx.x;
    int global_size = blockDim.x*gridDim.x;
    
    int foo;
    //no int8 in CUDA :(
    int4 accumulator_1 = make_int4(0,0,0,0);
    int4 accumulator_2 = make_int4(0,0,0,0);
    while (index < VL_GLO)
    {
        accumulator_1.x += istatus(index,0);
        accumulator_1.y += istatus(index,1);
        accumulator_1.z += istatus(index,2);
        //some dirty work on the side...
        foo = istatus(index,3);
        xNacc[index] = foo;
        accumulator_1.w += foo;
        foo = istatus(index,4);
        xNrej[index] = foo;
        accumulator_2.x += foo;
        accumulator_2.y += istatus(index,5);
        accumulator_2.z += istatus(index,6);
        accumulator_2.w += istatus(index,7);
        index += global_size;
    }
    //no int8 in CUDA :(
    __shared__ int4 buffer_1[REDUCTION_SIZE_1];
    __shared__ int4 buffer_2[REDUCTION_SIZE_1];
    
    buffer_1[idx_1] = accumulator_1;
    buffer_2[idx_1] = accumulator_2;
    __syncthreads();
    
    int idx_2, active_threads = blockDim.x;
    int4 tmp_1, tmp_2;
    while (active_threads != 1)
    {
        active_threads /= 2;
        if (idx_1 < active_threads)
        {
            idx_2 = idx_1+active_threads;
            
            tmp_1 = buffer_1[idx_1];
            tmp_2 = buffer_1[idx_2];
            
            tmp_1.x += tmp_2.x;
            tmp_1.y += tmp_2.y;
            tmp_1.z += tmp_2.z;
            tmp_1.w += tmp_2.w;
            
            buffer_1[idx_1] = tmp_1;
            
            
            tmp_1 = buffer_2[idx_1];
            tmp_2 = buffer_2[idx_2];
            
            tmp_1.x += tmp_2.x;
            tmp_1.y += tmp_2.y;
            tmp_1.z += tmp_2.z;
            tmp_1.w += tmp_2.w;
            
            buffer_2[idx_1] = tmp_1;
            
        }
        __syncthreads();
    }
    if (idx_1 == 0)
    {
        tmp_out_1[blockIdx.x] = buffer_1[0];
        tmp_out_2[blockIdx.x] = buffer_2[0];
    }
}            

__global__ void reduce_istatus_2(int4 *tmp_out_1, int4 *tmp_out_2, int *out)
{
    int idx_1 = threadIdx.x;
    //no int8 in CUDA :(
    __shared__ int4 buffer_1[REDUCTION_SIZE_2];
    __shared__ int4 buffer_2[REDUCTION_SIZE_2];
    
    buffer_1[idx_1] = tmp_out_1[idx_1];
    buffer_2[idx_1] = tmp_out_2[idx_1]; 
    __syncthreads();
    
    int idx_2, active_threads = blockDim.x;
    int4 tmp_1, tmp_2;
    while (active_threads != 1)
    {
        active_threads /= 2;
        if (idx_1 < active_threads)
        {
            idx_2 = idx_1+active_threads;
            
            tmp_1 = buffer_1[idx_1];
            tmp_2 = buffer_1[idx_2];
            
            tmp_1.x += tmp_2.x;
            tmp_1.y += tmp_2.y;
            tmp_1.z += tmp_2.z;
            tmp_1.w += tmp_2.w;
            
            buffer_1[idx_1] = tmp_1;
            
            
            tmp_1 = buffer_2[idx_1];
            tmp_2 = buffer_2[idx_2];
            
            tmp_1.x += tmp_2.x;
            tmp_1.y += tmp_2.y;
            tmp_1.z += tmp_2.z;
            tmp_1.w += tmp_2.w;
            
            buffer_2[idx_1] = tmp_1;
            
        }
        __syncthreads();
    }
    if (idx_1 == 0)
    {
        tmp_1 = buffer_1[0];
        tmp_2 = buffer_2[0];
        out[0] = tmp_1.x;
        out[1] = tmp_1.y;
        out[2] = tmp_1.z;
        out[3] = tmp_1.w;
        out[4] = tmp_2.x;
        out[5] = tmp_2.y;
        out[6] = tmp_2.z;
        out[7] = tmp_2.w;        
    }
}            

/* Assuming different processes */
enum { TRUE=1, FALSE=0 } ;
double *d_conc, *d_temp, *d_press, *d_cair, *d_khet_st, *d_khet_tr, *d_jx, *d_jac0, *d_Ghimj, *d_varNew, *d_K, *d_varErr, *d_dFdT, *d_Fcn0, *d_var, *d_fix, *d_rconst;
#if defined(__SINGLEPREC)
float *d_conc_s, *d_temp_s, *d_press_s, *d_cair_s, *d_khet_st_s, *d_khet_tr_s, *d_jx_s, *d_jac0_s, *d_Ghimj_s, *d_varNew_s, *d_K_s, *d_varErr_s, *d_dFdT_s, *d_Fcn0_s, *d_var_s, *d_fix_s, *d_rconst_s;
float *d_rstatus_s, *d_absTol_s, *d_relTol_s;
#endif
int initialized = FALSE;

/* Device pointers pointing to GPU */
double *d_rstatus, *d_absTol, *d_relTol;
int *d_istatus, *d_istatus_rd, *d_xNacc, *d_xNrej;
int4 *d_tmp_out_1, *d_tmp_out_2;

/* Allocate arrays on device for Rosenbrock */
__host__ void init_first_time(int pe, int VL_GLO, int size_khet_st, int size_khet_tr, int size_jx ){

    /* Select the proper GPU CARD */
    int deviceCount, device;
    gpuErrchk( cudaGetDeviceCount(&deviceCount) );
    device = pe % deviceCount;
    gpuErrchk( cudaSetDevice(device) );

    printf("PE[%d]: selected %d of total %d\n",pe,device,deviceCount);
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

    gpuErrchk( cudaMalloc ((void **) &d_conc   , sizeof(double)*VL_GLO*(NSPEC))        );
    gpuErrchk( cudaMalloc ((void **) &d_khet_st, sizeof(double)*VL_GLO*size_khet_st) );
    gpuErrchk( cudaMalloc ((void **) &d_khet_tr, sizeof(double)*VL_GLO*size_khet_tr) );
    gpuErrchk( cudaMalloc ((void **) &d_jx     , sizeof(double)*VL_GLO*size_jx)      );

    gpuErrchk( cudaMalloc ((void **) &d_rstatus    , sizeof(double)*VL_GLO*2)          );
    gpuErrchk( cudaMalloc ((void **) &d_istatus    , sizeof(int)*VL_GLO*8)             );
    gpuErrchk( cudaMalloc ((void **) &d_absTol     , sizeof(double)*NVAR)              );
    gpuErrchk( cudaMalloc ((void **) &d_relTol     , sizeof(double)*NVAR)              );

    /* Allocate input arrays */
    gpuErrchk( cudaMalloc ((void **) &temp_gpu     , sizeof(double)*VL_GLO)              );
    gpuErrchk( cudaMalloc ((void **) &press_gpu     , sizeof(double)*VL_GLO)              );
    gpuErrchk( cudaMalloc ((void **) &cair_gpu     , sizeof(double)*VL_GLO)              );

    /* Allocate arrays on device for reducing metrics */
    gpuErrchk( cudaMalloc ((void **) &d_istatus_rd  , sizeof(int)*8));
    gpuErrchk( cudaMalloc ((void **) &d_tmp_out_1   , sizeof(int4)*64));
    gpuErrchk( cudaMalloc ((void **) &d_tmp_out_2   , sizeof(int4)*64));
    gpuErrchk( cudaMalloc ((void **) &d_xNacc   , sizeof(int)*VL_GLO));
    gpuErrchk( cudaMalloc ((void **) &d_xNrej   , sizeof(int)*VL_GLO));

    /* Allocate arrays for solvers on device global memory to reduce the stack usage */
    gpuErrchk( cudaMalloc ((void **) &d_jac0, sizeof(double)*VL_GLO*LU_NONZERO)   );
    gpuErrchk( cudaMalloc ((void **) &d_Ghimj, sizeof(double)*VL_GLO*LU_NONZERO)   );
    gpuErrchk( cudaMalloc ((void **) &d_varNew, sizeof(double)*VL_GLO*NVAR)       );
    gpuErrchk( cudaMalloc ((void **) &d_Fcn0, sizeof(double)*VL_GLO*NVAR)       );
    gpuErrchk( cudaMalloc ((void **) &d_dFdT, sizeof(double)*VL_GLO*NVAR)       );

    gpuErrchk( cudaMalloc ((void **) &d_K, sizeof(double)*VL_GLO*NVAR*6)       );  // TODO: Change size according to solver steps
    gpuErrchk( cudaMalloc ((void **) &d_varErr, sizeof(double)*VL_GLO*NVAR)       );
    gpuErrchk( cudaMalloc ((void **) &d_var, sizeof(double)*VL_GLO*NSPEC)       );
    gpuErrchk( cudaMalloc ((void **) &d_fix, sizeof(double)*VL_GLO*NFIX)       );
    gpuErrchk( cudaMalloc ((void **) &d_rconst, sizeof(double)*VL_GLO*NREACT)       );

    #if defined(__SINGLEPREC)
    gpuErrchk( cudaMalloc ((void **) &d_conc_s   , sizeof(float)*VL_GLO*(NSPEC))        );
    gpuErrchk( cudaMalloc ((void **) &d_khet_st_s, sizeof(float)*VL_GLO*size_khet_st) );
    gpuErrchk( cudaMalloc ((void **) &d_khet_tr_s, sizeof(float)*VL_GLO*size_khet_tr) );
    gpuErrchk( cudaMalloc ((void **) &d_jx_s     , sizeof(float)*VL_GLO*size_jx)      );

    gpuErrchk( cudaMalloc ((void **) &d_rstatus_s    , sizeof(float)*VL_GLO*2)          );
    gpuErrchk( cudaMalloc ((void **) &d_absTol_s     , sizeof(float)*NVAR)              );
    gpuErrchk( cudaMalloc ((void **) &d_relTol_s     , sizeof(float)*NVAR)              );

    /* Allocate input arrays */
    gpuErrchk( cudaMalloc ((void **) &temp_gpu_s     , sizeof(float)*VL_GLO)              );
    gpuErrchk( cudaMalloc ((void **) &press_gpu_s     , sizeof(float)*VL_GLO)              );
    gpuErrchk( cudaMalloc ((void **) &cair_gpu_s     , sizeof(float)*VL_GLO)              );

    /* Allocate arrays for solvers on device global memory to reduce the stack usage */
    gpuErrchk( cudaMalloc ((void **) &d_jac0_s, sizeof(float)*VL_GLO*LU_NONZERO)   );
    gpuErrchk( cudaMalloc ((void **) &d_Ghimj_s, sizeof(float)*VL_GLO*LU_NONZERO)   );
    gpuErrchk( cudaMalloc ((void **) &d_varNew_s, sizeof(float)*VL_GLO*NVAR)       );
    gpuErrchk( cudaMalloc ((void **) &d_Fcn0_s, sizeof(float)*VL_GLO*NVAR)       );
    gpuErrchk( cudaMalloc ((void **) &d_dFdT_s, sizeof(float)*VL_GLO*NVAR)       );

    gpuErrchk( cudaMalloc ((void **) &d_K_s, sizeof(float)*VL_GLO*NVAR*6)       );  // TODO: Change size according to solver steps
    gpuErrchk( cudaMalloc ((void **) &d_varErr_s, sizeof(float)*VL_GLO*NVAR)       );
    gpuErrchk( cudaMalloc ((void **) &d_var_s, sizeof(float)*VL_GLO*NSPEC)       );
    gpuErrchk( cudaMalloc ((void **) &d_fix_s, sizeof(float)*VL_GLO*NFIX)       );
    gpuErrchk( cudaMalloc ((void **) &d_rconst_s, sizeof(float)*VL_GLO*NREACT)       );
    #endif

    initialized = TRUE;
}

/*
 * TODO: We should call it in some point..
 */
extern "C" void finalize_cuda(){
    /* Free memory on the device */
    gpuErrchk( cudaFree(d_conc        ) );
    gpuErrchk( cudaFree(d_temp        ) );
    gpuErrchk( cudaFree(d_press       ) );
    gpuErrchk( cudaFree(d_cair        ) );
    gpuErrchk( cudaFree(d_khet_st     ) );
    gpuErrchk( cudaFree(d_khet_tr     ) );
    gpuErrchk( cudaFree(d_jx          ) );
    gpuErrchk( cudaFree(d_rstatus     ) );
    gpuErrchk( cudaFree(d_istatus     ) );
    gpuErrchk( cudaFree(d_absTol      ) );
    gpuErrchk( cudaFree(d_relTol      ) );
    gpuErrchk( cudaFree(d_istatus_rd  ) );
    gpuErrchk( cudaFree(d_tmp_out_1   ) );
    gpuErrchk( cudaFree(d_tmp_out_2   ) );
    gpuErrchk( cudaFree(d_xNacc       ) );
    gpuErrchk( cudaFree(d_xNrej       ) );
    gpuErrchk( cudaFree(temp_gpu      ) );
    gpuErrchk( cudaFree(press_gpu     ) );
    gpuErrchk( cudaFree(cair_gpu      ) );

    gpuErrchk( cudaFree(d_jac0        ) );
    gpuErrchk( cudaFree(d_Ghimj       ) );
    gpuErrchk( cudaFree(d_varNew      ) );
    gpuErrchk( cudaFree(d_Fcn0        ) );
    gpuErrchk( cudaFree(d_dFdT        ) );
    gpuErrchk( cudaFree(d_K           ) );
    gpuErrchk( cudaFree(d_varErr      ) );
    gpuErrchk( cudaFree(d_var         ) );
    gpuErrchk( cudaFree(d_fix         ) );
    gpuErrchk( cudaFree(d_rconst      ) );

    #if defined(__SINGLEPREC)
    gpuErrchk( cudaFree(d_conc_s        ) );
    gpuErrchk( cudaFree(d_temp_s        ) );
    gpuErrchk( cudaFree(d_press_s       ) );
    gpuErrchk( cudaFree(d_cair_s        ) );
    gpuErrchk( cudaFree(d_khet_st_s     ) );
    gpuErrchk( cudaFree(d_khet_tr_s     ) );
    gpuErrchk( cudaFree(d_jx_s          ) );
    gpuErrchk( cudaFree(d_rstatus_s     ) );
    gpuErrchk( cudaFree(d_absTol_s      ) );
    gpuErrchk( cudaFree(d_relTol_s      ) );
    gpuErrchk( cudaFree(temp_gpu_s      ) );
    gpuErrchk( cudaFree(press_gpu_s     ) );
    gpuErrchk( cudaFree(cair_gpu_s      ) );

    gpuErrchk( cudaFree(d_jac0_s        ) );
    gpuErrchk( cudaFree(d_Ghimj_s       ) );
    gpuErrchk( cudaFree(d_varNew_s      ) );
    gpuErrchk( cudaFree(d_Fcn0_s        ) );
    gpuErrchk( cudaFree(d_dFdT_s        ) );
    gpuErrchk( cudaFree(d_K_s           ) );
    gpuErrchk( cudaFree(d_varErr_s      ) );
    gpuErrchk( cudaFree(d_var_s         ) );
    gpuErrchk( cudaFree(d_fix_s         ) );
    gpuErrchk( cudaFree(d_rconst_s      ) );

    #endif
}



extern "C" void kpp_integrate_cuda_( int *pe_p, int *sizes, double *time_step_len_p, double *conc, double *temp, double *press, double *cair, 
                                    double *khet_st, double *khet_tr, double *jx, double *absTol, double *relTol, int *ierr, int *istatus, 
                                    int *xNacc, int *xNrej, double *rndoff, int *icntrl=NULL, double *rcntrl=NULL
				    ) 
/*  // TODO
 *  Parameters:
 *          pe_p: scalar int - processor element
 *        VL_GLO: scalar int - size of the system
 *         NSPEC: scalar int - number of species
 *        NREACT: scalar int - number of reactions
 *          NVAR: scalar int - 
 *
 *  Input data:
 *          conc: 2D array of doubles - size: vl_glo x number of species
 *          temp: 1D array of doubles - size: vl_glo
 *         press: 1D array of doubles - size: vl_glo
 *          cair: 1D array of doubles - size: vl_glo
 *       khet_st: 2D array of doubles - size: vl_glo x number of species
 *       khet_tr: 2D array of doubles - size: vl_glo x number of species 
 *            jx: 2D array of doubles - size: vl_glo x number of species
 *        absTol: 1D array of doubles - size: number of species
 *        relTol: 1D array of doubles - size: number of species
 *     Control:
 *        icntrl: 1D array of ints   - size: 4
 *         sizes: 1D array of ints   - size: 4
 *        rcntrl: 1D array of doubles - size: 7
 * 
 * 
 */
{

    const REAL DELTAMIN = 1.0E-5;


    
    int VL_GLO       = sizes[0];
    int size_khet_st = sizes[1];
    int size_khet_tr = sizes[2];
    int size_jx      = sizes[3];
    REAL roundoff  = *rndoff; 
    
    REAL Tstart,Tend;
    Tstart = ZERO;
    Tend   =  *time_step_len_p;
    int pe = *pe_p;
    
    // variables from rcntrl and icntrl
    int autonomous, vectorTol, UplimTol, method, Max_no_steps;
    REAL Hmin, Hmax, Hstart, FacMin, FacMax, FacRej, FacSafe;
    
    //int rcntrl_bool = 0, icntrl_bool=0;
    if (rcntrl == NULL)
    {
        rcntrl = new double[7];
        for (int i=0; i < 7; i++)
            rcntrl[i] = 0.0;
    }
    if (icntrl == NULL)
    {
        icntrl = new int[4];
        for (int i=0; i < 4; i++)
            icntrl[i] = 0;
    }

    /* Allocate arrays on device for update_rconst kernel*/        
    if (initialized == FALSE)   init_first_time(pe, VL_GLO, size_khet_st, size_khet_tr, size_jx);

    /* Copy data from host memory to device memory */
    gpuErrchk( cudaMemcpy(d_conc   , conc   	, sizeof(double)*VL_GLO*NSPEC        , cudaMemcpyHostToDevice) );

    gpuErrchk( cudaMemcpy(temp_gpu   , temp   	, sizeof(double)*VL_GLO  , cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(press_gpu  , press  	, sizeof(double)*VL_GLO  , cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(cair_gpu   , cair   	, sizeof(double)*VL_GLO  , cudaMemcpyHostToDevice) );

    gpuErrchk( cudaMemcpy(d_khet_st, khet_st	, sizeof(double)*VL_GLO*size_khet_st , cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(d_khet_tr, khet_tr	, sizeof(double)*VL_GLO*size_khet_tr , cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(d_jx     , jx     	, sizeof(double)*VL_GLO*size_jx      , cudaMemcpyHostToDevice) );

    /* Copy arrays from host memory to device memory for Rosenbrock */    
    gpuErrchk( cudaMemcpy(d_absTol, absTol, sizeof(double)*NVAR, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(d_relTol, relTol, sizeof(double)*NVAR, cudaMemcpyHostToDevice) );


    /* Compute execution configuration for update_rconst */
    int block_size, grid_size;
    
    block_size = BLOCKSIZE;
    grid_size = (VL_GLO + block_size - 1)/block_size;  
    dim3 dimBlock(block_size);
    dim3 dimGrid(grid_size);


    #if defined(__SINGLEPREC)
    doubleToFloat<<<dimGrid,dimBlock>>>(d_conc_s,d_conc,VL_GLO*NSPEC);
    doubleToFloat<<<dimGrid,dimBlock>>>(temp_gpu_s,temp_gpu,VL_GLO);
    doubleToFloat<<<dimGrid,dimBlock>>>(press_gpu_s,press_gpu,VL_GLO);
    doubleToFloat<<<dimGrid,dimBlock>>>(cair_gpu_s,cair_gpu,VL_GLO);
    doubleToFloat<<<dimGrid,dimBlock>>>(d_khet_st_s,d_khet_st,VL_GLO*size_khet_st);
    doubleToFloat<<<dimGrid,dimBlock>>>(d_khet_tr_s,d_khet_tr,VL_GLO*size_khet_tr);
    doubleToFloat<<<dimGrid,dimBlock>>>(d_jx_s,d_jx,VL_GLO*size_jx);
    doubleToFloat<<<dimGrid,dimBlock>>>(d_absTol_s,d_absTol,NVAR);
    doubleToFloat<<<dimGrid,dimBlock>>>(d_relTol_s,d_relTol,NVAR);
    

 //   gpuErrchk( cudaFree(d_temp        ) );
 //   gpuErrchk( cudaFree(d_press       ) );
//    gpuErrchk( cudaFree(d_cair        ) );
    // gpuErrchk( cudaFree(d_khet_st     ) );
    // gpuErrchk( cudaFree(d_khet_tr     ) );
    // gpuErrchk( cudaFree(d_jx          ) );
  //  gpuErrchk( cudaFree(d_rstatus     ) );
  //  gpuErrchk( cudaFree(d_absTol      ) );
  //  gpuErrchk( cudaFree(d_relTol      ) );
 //   gpuErrchk( cudaFree(d_istatus_rd  ) );
 //   gpuErrchk( cudaFree(d_tmp_out_1   ) );
 //   gpuErrchk( cudaFree(d_tmp_out_2   ) );
 //   gpuErrchk( cudaFree(d_xNacc       ) );
 //  gpuErrchk( cudaFree(d_xNrej       ) );
 //   gpuErrchk( cudaFree(temp_gpu      ) );
 //   gpuErrchk( cudaFree(press_gpu     ) );
 //   gpuErrchk( cudaFree(cair_gpu      ) );

  //  gpuErrchk( cudaFree(d_jac0        ) );
  //  gpuErrchk( cudaFree(d_Ghimj       ) );
 //   gpuErrchk( cudaFree(d_varNew      ) );
 //   gpuErrchk( cudaFree(d_Fcn0        ) );
 //   gpuErrchk( cudaFree(d_dFdT        ) );
  //  gpuErrchk( cudaFree(d_K           ) );
 //   gpuErrchk( cudaFree(d_varErr      ) );
 //   gpuErrchk( cudaFree(d_var         ) );
//    gpuErrchk( cudaFree(d_fix         ) );
//    gpuErrchk( cudaFree(d_rconst      ) );

    #endif

    /* Execute the kernel */
    //update_rconst<<<dimGrid,dimBlock>>>(d_conc, d_khet_st, d_khet_tr, d_jx, VL_GLO); 

    GPU_DEBUG();
 
//  *------------------------------------------------------*
//  |    Default values vs input settings (icntrl, rcntrl) |
//  *------------------------------------------------------*
    int ierr_tmp=0;
    {
    //  autonomous or time dependent ODE. Default is time dependent.
        autonomous = !(icntrl[0] == 0);

    //  For Scalar tolerances (icntrl[1].NE.0)  the code uses absTol(0) and relTol(0)
    //  For Vector tolerances (icntrl[1] == 0) the code uses absTol(0:NVAR) and relTol(0:NVAR)
        if (icntrl[1] == 0)
        {
            vectorTol = 1; //bool
            UplimTol  = NVAR;
        }
        else
        {
            vectorTol = 0;
            UplimTol  = 1;
        }

    //  The particular Rosenbrock method chosen
        if (icntrl[2] == 0) 
        {
            method = 4;
        }
        else if ((icntrl[2] >= 1) && (icntrl[2] <= 5))
        {
            method = icntrl[2];
        }
        else
        {
            printf("User-selected Rosenbrock method: icntrl[2]=%d\n",method);
            ierr_tmp = -2;
        }
    //  The maximum number of steps admitted
        if (icntrl[3] == 0)
        {
            Max_no_steps = 100000;
        }
        else if (icntrl[3] > 0) 
        {
            Max_no_steps=icntrl[3];
        }
        else
        {
            printf("User-selected max no. of steps: icntrl[3]=%d\n",icntrl[3]);
            ierr_tmp = -1;
        }
    //  Unit roundoff (1+ roundoff>1)
        roundoff = machine_eps_flt(); 

    //  Lower bound on the step size: (positive value)
        if (rcntrl[0] == ZERO)
        {
            Hmin = ZERO;
        }
        else if (rcntrl[0] > ZERO) 
        {
            Hmin = rcntrl[0];
        }
        else
        {
            printf("User-selected Hmin: rcntrl[0]=%f\n",rcntrl[0]);
            ierr_tmp = -3;
        }
    //  Upper bound on the step size: (positive value)
        if (rcntrl[1] == ZERO) 
        {
            Hmax = fabs(Tend-Tstart);
        }
        else if (rcntrl[1] > ZERO) 
        {
            Hmax = fmin(fabs(rcntrl[1]),fabs(Tend-Tstart));
        }
        else
        {
            printf("User-selected Hmax: rcntrl[1]=%f\n",rcntrl[1]);
            ierr_tmp = -3;
        }
    //  Starting step size: (positive value)
        if (rcntrl[2] == ZERO) 
        {
            Hstart = fmax(Hmin,DELTAMIN);
        }
        else if (rcntrl[2] > ZERO) 
        {
            Hstart = fmin(fabs(rcntrl[2]),fabs(Tend-Tstart));
        }
        else
        {
            printf("User-selected Hstart: rcntrl[2]=%f\n",rcntrl[2]);
            ierr_tmp = -3;
        }
    //  Step size can be changed s.t.  FacMin < Hnew/Hexit < FacMax
        if (rcntrl[3] == ZERO)
        {
            FacMin = 0.2;
        }
        else if (rcntrl[3] > ZERO) 
        {
            FacMin = rcntrl[3];
        }
        else
        {
            printf("User-selected FacMin: rcntrl[3]=%f\n",rcntrl[3]);
            ierr_tmp = -4;
        }
        if (rcntrl[4] == ZERO) 
        {
            FacMax = 6.0;
        }
        else if (rcntrl[4] > ZERO) 
        {
            FacMax = rcntrl[4];
        }
        else
        {
            printf("User-selected FacMax: rcntrl[4]=%f\n",rcntrl[4]);
            ierr_tmp = -4; 
        }
    //  FacRej: Factor to decrease step after 2 succesive rejections
        if (rcntrl[5] == ZERO) 
        {
            FacRej = 0.1;
        }
        else if (rcntrl[5] > ZERO) 
        {
            FacRej = rcntrl[5];
        }
        else
        {
            printf("User-selected FacRej: rcntrl[5]=%f\n",rcntrl[5]);
            ierr_tmp = -4;
        }
    //  FacSafe: Safety Factor in the computation of new step size
        if (rcntrl[6] == ZERO) 
        {
            FacSafe = 0.9;
        }
        else if (rcntrl[6] > ZERO)
        {
            FacSafe = rcntrl[6];
        }
        else
        {
            printf("User-selected FacSafe: rcntrl[6]=%f\n",rcntrl[6]);
            ierr_tmp = -4;
        }
    //  Check if tolerances are reasonable
        for (int i=0; i < UplimTol; i++)
        {
            if ((absTol[i] <= ZERO) || (relTol[i] <= 10.0*roundoff) || (relTol[i] >= 1.0))
            {
                printf("CCC absTol(%d) = %f \n",i,absTol[i]);
                printf("CCC relTol(%d) = %f \n",i,relTol[i]);
                ierr_tmp = -5;
            }
        }
    }


#if defined(__SINGLEPREC) 
     Rosenbrock<<<dimGrid,dimBlock>>>(d_conc_s, Tstart, Tend, d_rstatus_s, d_istatus, 
                    // values calculated from icntrl and rcntrl at host 
                    autonomous, vectorTol, UplimTol, method, Max_no_steps, 
                    d_jac0_s, d_Ghimj_s,d_varNew_s, d_K_s, d_varErr_s, d_dFdT_s, d_Fcn0_s, d_var_s, d_fix_s, d_rconst_s,
                    Hmin, Hmax, Hstart, FacMin, FacMax, FacRej, FacSafe, roundoff,
                    //  cuda global mem buffers               
                    d_absTol_s, d_relTol_s, 
                    d_khet_st_s, d_khet_tr_s, d_jx_s, 
                    // Global input arrays 
                    temp_gpu_s, press_gpu_s, cair_gpu_s, 
                    // extra - vector lenght and processor 
                    VL_GLO); 
     #else 
     Rosenbrock<<<dimGrid,dimBlock>>>(d_conc, Tstart, Tend, d_rstatus, d_istatus, 
                    // values calculated from icntrl and rcntrl at host 
                    autonomous, vectorTol, UplimTol, method, Max_no_steps, 
                    d_jac0, d_Ghimj,d_varNew, d_K, d_varErr, d_dFdT, d_Fcn0, d_var, d_fix, d_rconst, 
                    Hmin, Hmax, Hstart, FacMin, FacMax, FacRej, FacSafe, roundoff, 
                    //  cuda global mem buffers             
                    d_absTol, d_relTol, 
                    d_khet_st, d_khet_tr, d_jx, 
                    // Global input arrays 
                    temp_gpu, press_gpu, cair_gpu, 
                    // extra - vector lenght and processor 
                    VL_GLO); 
     #endif
    GPU_DEBUG();

    
    reduce_istatus_1<<<REDUCTION_SIZE_2,REDUCTION_SIZE_1>>>(d_istatus, d_tmp_out_1, d_tmp_out_2, VL_GLO, d_xNacc, d_xNrej);


    GPU_DEBUG();

    reduce_istatus_2<<<1,REDUCTION_SIZE_2>>>(d_tmp_out_1, d_tmp_out_2, d_istatus_rd);

    GPU_DEBUG();
    
    /* Copy the result back */
    #if defined(__SINGLEPREC)
    floatToDouble<<<dimGrid,dimBlock>>>(d_conc,d_conc_s,VL_GLO*NSPEC);
    #endif
    
    gpuErrchk( cudaMemcpy( conc      , d_conc       , sizeof(double)*VL_GLO*NSPEC, cudaMemcpyDeviceToHost) );  
    gpuErrchk( cudaMemcpy( xNacc      , d_xNacc      , sizeof(int)*VL_GLO         , cudaMemcpyDeviceToHost) );  
    gpuErrchk( cudaMemcpy( xNrej      , d_xNrej      , sizeof(int)*VL_GLO         , cudaMemcpyDeviceToHost) ); 

    return;

}





/*
 *   Example input data for testing the verification of the transformation.
 *
 *   Copyright 2016 The Cyprus Institute 
 *
 *   Developers: Michail Alvanos - m.alvanos@cyi.ac.cy
 *               Theodoros Christoudias - christoudias@cyi.ac.cy
 *
 * */

#include <time.h>
#include <sys/time.h>

#define VL_GLO 5760
#define IPMAX 115

double conc[VL_GLO*NSPEC];
double temp[VL_GLO];
double press[VL_GLO];
double cair[VL_GLO];
double jx[VL_GLO*NSPEC];


int xNacc[VL_GLO];
int xNrej[VL_GLO];

double conc_cell[NSPEC] = {
0.000000000000000E+000,
0.000000000000000E+000,
1.130030837133365E-006,
2161.17681825926,
1.469481417859824E-004,
2.894067546497780E-004,
0.000000000000000E+000,
0.000000000000000E+000,
6.377486492629032E-031,
2.774602114035594E-004,
9.159068418074058E-022,
1.681545841334171E-030,
6.587848965925121E-036,
4.057130203198298E-031,
7.556675262619906E-006,
5.625822089563362E-006,
7.248546508346980E-010,
7.771754415762507E-039,
1.672965892516881E-032,
5.778276640099593E-029,
2.169623196599310E-031,
4.449685524913890E-029,
9.236991853178721E-028,
1.731254847935413E-009,
6.419363370200839E-028,
4.035724058634079E-029,
6234.08726448302,
25802.7788132849,
1.33974252411334,
11.1514176946459,
8.023966161170008E-032,
1.405402576145367E-030,
2.416365419045456E-029,
3.763980220765519E-033,
3.687747273615521E-004,
4.400695805857555E-030,
8.096351349854847E-009,
1.605777396541510E-008,
8.424266813161654E-005,
1.275728897910597E-029,
36780.6069067007,
44.2802185584881,
5.485594561042764E-010,
3.418234885986840E-032,
1.808885697309332E-008,
2.295321288609202E-030,
7.186736555958003E-032,
667193926.549068,
9.443976722997098E-030,
2.065479750965850E-030,
658798139.717353,
5013220.82927210,
6.594652607797343E-013,
4.779051920325237E-033,
0.241330392051758,
2.657031589287186E-030,
1.166890334972386E-014,
337.069782231658,
126494.977205691,
891.196915201611,
222.557367243832,
1.22451624669813,
4845.02754823106,
535329.616196368,
3.077774956209536E-002,
989833722.937206,
38527.6291432442,
1.857293910861109E-007,
5035616002.44018,
26824247.3107905,
211466.239175163,
60638129767802.7,
225227339137553.,
87651408241.1165
};

double abstol[NSPEC] = {
    0.0
};

double reltol[NSPEC] = {

    0.0
};

double khet_st[VL_GLO*NSPEC] = {
    0.0
};

double khet_tr[VL_GLO*NSPEC] = {
    0.0
};


#define COPY_DATA()\
    for (i=0;i<VL_GLO;i++){\
        for (j=0;j<NSPEC;j++){\
              conc[j*VL_GLO + i] = conc_cell[j];\
        }\
    }


#define PRINT_DATA()\
        printf("Results:");    \
        for (j=0;j<NSPEC;j++){\
            printf("   %.12e  ",conc[j*VL_GLO]);\
        }\
        printf("\n");    



int main(int argc, char **argv){
    
    int n = 1; // proccess element
   


    int istatus;
    int ierr;
    int i,j;

    int sizes[4] = {VL_GLO,NSPEC,NSPEC,IPMAX}; 
    int icntrl[20] = {0,0,2,0};

    double roundoff = 2.220446049250313E-016;
    double timestep = 720.0;

    for (i=0;i<VL_GLO;i++){
        for (j=0;j<NSPEC;j++){
              conc[i*NSPEC + j] = conc_cell[j];
        }
        temp[i] = 240.995971972245;
        press[i] = 0.994591236114502; 
        cair[i] = 298914285136738.0;

        khet_tr[i*4 + 0] = 7.408371201503456E-008;
        khet_tr[i*4 + 1] = 4.849455570110968E-007;
        khet_tr[i*4 + 2] =  0.000000000000000E+000;  
        khet_tr[i*4 + 3] = 2.718003287797325E-007;
    }
       

    for (i=0;i<NSPEC;i++){
        abstol[i] = 10.0; 
        reltol[i] = 0.5; 
    }


    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1); 

    COPY_DATA();

    kpp_integrate_cuda_( &n, sizes, &timestep, conc, temp, press, cair, khet_st, khet_tr, jx, abstol, reltol, &ierr, &istatus, xNacc, xNrej, &roundoff, icntrl);



        




    struct timeval start, end;

    if (argc==2){
        icntrl[2] = atoi(argv[1]);
        COPY_DATA()
        gettimeofday(&start, NULL);
        kpp_integrate_cuda_( &n, sizes, &timestep, conc, temp, press, cair, khet_st, khet_tr, jx, abstol, reltol, &ierr, &istatus, xNacc, xNrej, &roundoff, icntrl);
        gettimeofday(&end, NULL);
        printf("%d: %ld (ms)\n", icntrl[2],((end.tv_sec * 1000 + end.tv_usec/1000) - (start.tv_sec * 1000 + start.tv_usec/1000)));
        PRINT_DATA();

        return 0;
    }



    icntrl[2] = 1;

restart:



    COPY_DATA();
    gettimeofday(&start, NULL);

    kpp_integrate_cuda_( &n, sizes, &timestep, conc, temp, press, cair, khet_st, khet_tr, jx, abstol, reltol, &ierr, &istatus, xNacc, xNrej, &roundoff, icntrl);

    gettimeofday(&end, NULL);

    PRINT_DATA();

    printf("%d: %ld (ms)\n", icntrl[2],((end.tv_sec * 1000 + end.tv_usec/1000)
                - (start.tv_sec * 1000 + start.tv_usec/1000)));
    icntrl[2]++;
    if ( icntrl[2] >5) return;
    goto restart;





}










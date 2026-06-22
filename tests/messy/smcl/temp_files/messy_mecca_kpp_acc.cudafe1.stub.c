#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wcast-qual"
#define __NV_CUBIN_HANDLE_STORAGE__ static
#if !defined(__CUDA_INCLUDE_COMPILER_INTERNAL_HEADERS__)
#define __CUDA_INCLUDE_COMPILER_INTERNAL_HEADERS__
#endif
#include "crt/host_runtime.h"
#include "messy_mecca_kpp_acc.fatbin.c"
extern void __device_stub__Z13doubleToFloatPfPdi(float *, double *, int);
extern void __device_stub__Z13floatToDoublePdPfi(double *, float *, int);
extern void __device_stub__Z10RosenbrockPfffS_PiiiiiiS_S_S_S_S_S_S_S_S_S_ffffffffPKfS2_S2_S2_S2_S2_S2_S2_i(REAL *__restrict__, const REAL, const REAL, REAL *__restrict__, int *__restrict__, const int, const int, const int, const int, const int, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, const REAL, const REAL, const REAL, const REAL, const REAL, const REAL, const REAL, const REAL, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const int);
extern void __device_stub__Z16reduce_istatus_1PiP4int4S1_iS_S_(int *, struct int4 *, struct int4 *, int, int *, int *);
extern void __device_stub__Z16reduce_istatus_2P4int4S0_Pi(struct int4 *, struct int4 *, int *);
static void __nv_cudaEntityRegisterCallback(void **);
static void __sti____cudaRegisterAll(void) __attribute__((__constructor__));
void __device_stub__Z13doubleToFloatPfPdi(float *__par0, double *__par1, int __par2){__cudaLaunchPrologue(3);__cudaSetupArgSimple(__par0, 0UL);__cudaSetupArgSimple(__par1, 8UL);__cudaSetupArgSimple(__par2, 16UL);__cudaLaunch(((char *)((void ( *)(float *, double *, int))doubleToFloat)));}
# 3853 "messy_mecca_kpp_acc.cu"
void doubleToFloat( float *__cuda_0,double *__cuda_1,int __cuda_2)
# 3854 "messy_mecca_kpp_acc.cu"
{__device_stub__Z13doubleToFloatPfPdi( __cuda_0,__cuda_1,__cuda_2);




}
# 1 "./temp_files/messy_mecca_kpp_acc.cudafe1.stub.c"
void __device_stub__Z13floatToDoublePdPfi( double *__par0,  float *__par1,  int __par2) {  __cudaLaunchPrologue(3); __cudaSetupArgSimple(__par0, 0UL); __cudaSetupArgSimple(__par1, 8UL); __cudaSetupArgSimple(__par2, 16UL); __cudaLaunch(((char *)((void ( *)(double *, float *, int))floatToDouble))); }
# 3861 "messy_mecca_kpp_acc.cu"
void floatToDouble( double *__cuda_0,float *__cuda_1,int __cuda_2)
# 3862 "messy_mecca_kpp_acc.cu"
{__device_stub__Z13floatToDoublePdPfi( __cuda_0,__cuda_1,__cuda_2);




}
# 1 "./temp_files/messy_mecca_kpp_acc.cudafe1.stub.c"
void __device_stub__Z10RosenbrockPfffS_PiiiiiiS_S_S_S_S_S_S_S_S_S_ffffffffPKfS2_S2_S2_S2_S2_S2_S2_i( REAL *__restrict__ __par0,  const REAL __par1,  const REAL __par2,  REAL *__restrict__ __par3,  int *__restrict__ __par4,  const int __par5,  const int __par6,  const int __par7,  const int __par8,  const int __par9,  REAL *__restrict__ __par10,  REAL *__restrict__ __par11,  REAL *__restrict__ __par12,  REAL *__restrict__ __par13,  REAL *__restrict__ __par14,  REAL *__restrict__ __par15,  REAL *__restrict__ __par16,  REAL *__restrict__ __par17,  REAL *__restrict__ __par18,  REAL *__restrict__ __par19,  const REAL __par20,  const REAL __par21,  const REAL __par22,  const REAL __par23,  const REAL __par24,  const REAL __par25,  const REAL __par26,  const REAL __par27,  const REAL *__restrict__ __par28,  const REAL *__restrict__ __par29,  const REAL *__restrict__ __par30,  const REAL *__restrict__ __par31,  const REAL *__restrict__ __par32,  const REAL *__restrict__ __par33,  const REAL *__restrict__ __par34,  const REAL *__restrict__ __par35,  const int __par36) {  REAL *__T125;
 REAL *__T126;
 int *__T127;
 REAL *__T128;
 REAL *__T129;
 REAL *__T130;
 REAL *__T131;
 REAL *__T132;
 REAL *__T133;
 REAL *__T134;
 REAL *__T135;
 REAL *__T136;
 REAL *__T137;
 const REAL *__T138;
 const REAL *__T139;
 const REAL *__T140;
 const REAL *__T141;
 const REAL *__T142;
 const REAL *__T143;
 const REAL *__T144;
 const REAL *__T145;
__cudaLaunchPrologue(37); __T125 = __par0; __cudaSetupArgSimple(__T125, 0UL); __cudaSetupArgSimple(__par1, 8UL); __cudaSetupArgSimple(__par2, 12UL); __T126 = __par3; __cudaSetupArgSimple(__T126, 16UL); __T127 = __par4; __cudaSetupArgSimple(__T127, 24UL); __cudaSetupArgSimple(__par5, 32UL); __cudaSetupArgSimple(__par6, 36UL); __cudaSetupArgSimple(__par7, 40UL); __cudaSetupArgSimple(__par8, 44UL); __cudaSetupArgSimple(__par9, 48UL); __T128 = __par10; __cudaSetupArgSimple(__T128, 56UL); __T129 = __par11; __cudaSetupArgSimple(__T129, 64UL); __T130 = __par12; __cudaSetupArgSimple(__T130, 72UL); __T131 = __par13; __cudaSetupArgSimple(__T131, 80UL); __T132 = __par14; __cudaSetupArgSimple(__T132, 88UL); __T133 = __par15; __cudaSetupArgSimple(__T133, 96UL); __T134 = __par16; __cudaSetupArgSimple(__T134, 104UL); __T135 = __par17; __cudaSetupArgSimple(__T135, 112UL); __T136 = __par18; __cudaSetupArgSimple(__T136, 120UL); __T137 = __par19; __cudaSetupArgSimple(__T137, 128UL); __cudaSetupArgSimple(__par20, 136UL); __cudaSetupArgSimple(__par21, 140UL); __cudaSetupArgSimple(__par22, 144UL); __cudaSetupArgSimple(__par23, 148UL); __cudaSetupArgSimple(__par24, 152UL); __cudaSetupArgSimple(__par25, 156UL); __cudaSetupArgSimple(__par26, 160UL); __cudaSetupArgSimple(__par27, 164UL); __T138 = __par28; __cudaSetupArgSimple(__T138, 168UL); __T139 = __par29; __cudaSetupArgSimple(__T139, 176UL); __T140 = __par30; __cudaSetupArgSimple(__T140, 184UL); __T141 = __par31; __cudaSetupArgSimple(__T141, 192UL); __T142 = __par32; __cudaSetupArgSimple(__T142, 200UL); __T143 = __par33; __cudaSetupArgSimple(__T143, 208UL); __T144 = __par34; __cudaSetupArgSimple(__T144, 216UL); __T145 = __par35; __cudaSetupArgSimple(__T145, 224UL); __cudaSetupArgSimple(__par36, 232UL); __cudaLaunch(((char *)((void ( *)(REAL *__restrict__, const REAL, const REAL, REAL *__restrict__, int *__restrict__, const int, const int, const int, const int, const int, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, const REAL, const REAL, const REAL, const REAL, const REAL, const REAL, const REAL, const REAL, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const int))Rosenbrock))); }
# 4085 "messy_mecca_kpp_acc.cu"
void Rosenbrock( REAL *__restrict__ __cuda_0,const REAL __cuda_1,const REAL __cuda_2,REAL *__restrict__ __cuda_3,int *__restrict__ __cuda_4,const int __cuda_5,const int __cuda_6,const int __cuda_7,const int __cuda_8,const int __cuda_9,REAL *__restrict__ __cuda_10,REAL *__restrict__ __cuda_11,REAL *__restrict__ __cuda_12,REAL *__restrict__ __cuda_13,REAL *__restrict__ __cuda_14,REAL *__restrict__ __cuda_15,REAL *__restrict__ __cuda_16,REAL *__restrict__ __cuda_17,REAL *__restrict__ __cuda_18,REAL *__restrict__ __cuda_19,const REAL __cuda_20,const REAL __cuda_21,const REAL __cuda_22,const REAL __cuda_23,const REAL __cuda_24,const REAL __cuda_25,const REAL __cuda_26,const REAL __cuda_27,const REAL *__restrict__ __cuda_28,const REAL *__restrict__ __cuda_29,const REAL *__restrict__ __cuda_30,const REAL *__restrict__ __cuda_31,const REAL *__restrict__ __cuda_32,const REAL *__restrict__ __cuda_33,const REAL *__restrict__ __cuda_34,const REAL *__restrict__ __cuda_35,const int __cuda_36)
# 4101 "messy_mecca_kpp_acc.cu"
{__device_stub__Z10RosenbrockPfffS_PiiiiiiS_S_S_S_S_S_S_S_S_S_ffffffffPKfS2_S2_S2_S2_S2_S2_S2_i( __cuda_0,__cuda_1,__cuda_2,__cuda_3,__cuda_4,__cuda_5,__cuda_6,__cuda_7,__cuda_8,__cuda_9,__cuda_10,__cuda_11,__cuda_12,__cuda_13,__cuda_14,__cuda_15,__cuda_16,__cuda_17,__cuda_18,__cuda_19,__cuda_20,__cuda_21,__cuda_22,__cuda_23,__cuda_24,__cuda_25,__cuda_26,__cuda_27,__cuda_28,__cuda_29,__cuda_30,__cuda_31,__cuda_32,__cuda_33,__cuda_34,__cuda_35,__cuda_36);
# 4204 "messy_mecca_kpp_acc.cu"
}
# 1 "./temp_files/messy_mecca_kpp_acc.cudafe1.stub.c"
void __device_stub__Z16reduce_istatus_1PiP4int4S1_iS_S_( int *__par0,  struct int4 *__par1,  struct int4 *__par2,  int __par3,  int *__par4,  int *__par5) {  __cudaLaunchPrologue(6); __cudaSetupArgSimple(__par0, 0UL); __cudaSetupArgSimple(__par1, 8UL); __cudaSetupArgSimple(__par2, 16UL); __cudaSetupArgSimple(__par3, 24UL); __cudaSetupArgSimple(__par4, 32UL); __cudaSetupArgSimple(__par5, 40UL); __cudaLaunch(((char *)((void ( *)(int *, struct int4 *, struct int4 *, int, int *, int *))reduce_istatus_1))); }
# 4210 "messy_mecca_kpp_acc.cu"
void reduce_istatus_1( int *__cuda_0,struct int4 *__cuda_1,struct int4 *__cuda_2,int __cuda_3,int *__cuda_4,int *__cuda_5)
# 4211 "messy_mecca_kpp_acc.cu"
{__device_stub__Z16reduce_istatus_1PiP4int4S1_iS_S_( __cuda_0,__cuda_1,__cuda_2,__cuda_3,__cuda_4,__cuda_5);
# 4283 "messy_mecca_kpp_acc.cu"
}
# 1 "./temp_files/messy_mecca_kpp_acc.cudafe1.stub.c"
void __device_stub__Z16reduce_istatus_2P4int4S0_Pi( struct int4 *__par0,  struct int4 *__par1,  int *__par2) {  __cudaLaunchPrologue(3); __cudaSetupArgSimple(__par0, 0UL); __cudaSetupArgSimple(__par1, 8UL); __cudaSetupArgSimple(__par2, 16UL); __cudaLaunch(((char *)((void ( *)(struct int4 *, struct int4 *, int *))reduce_istatus_2))); }
# 4285 "messy_mecca_kpp_acc.cu"
void reduce_istatus_2( struct int4 *__cuda_0,struct int4 *__cuda_1,int *__cuda_2)
# 4286 "messy_mecca_kpp_acc.cu"
{__device_stub__Z16reduce_istatus_2P4int4S0_Pi( __cuda_0,__cuda_1,__cuda_2);
# 4342 "messy_mecca_kpp_acc.cu"
}
# 1 "./temp_files/messy_mecca_kpp_acc.cudafe1.stub.c"
static void __nv_cudaEntityRegisterCallback( void **__T212) {  __nv_dummy_param_ref(__T212); __nv_save_fatbinhandle_for_managed_rt(__T212); __cudaRegisterEntry(__T212, ((void ( *)(struct int4 *, struct int4 *, int *))reduce_istatus_2), _Z16reduce_istatus_2P4int4S0_Pi, (-1)); __cudaRegisterEntry(__T212, ((void ( *)(int *, struct int4 *, struct int4 *, int, int *, int *))reduce_istatus_1), _Z16reduce_istatus_1PiP4int4S1_iS_S_, (-1)); __cudaRegisterEntry(__T212, ((void ( *)(REAL *__restrict__, const REAL, const REAL, REAL *__restrict__, int *__restrict__, const int, const int, const int, const int, const int, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, REAL *__restrict__, const REAL, const REAL, const REAL, const REAL, const REAL, const REAL, const REAL, const REAL, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const REAL *__restrict__, const int))Rosenbrock), _Z10RosenbrockPfffS_PiiiiiiS_S_S_S_S_S_S_S_S_S_ffffffffPKfS2_S2_S2_S2_S2_S2_S2_i, (-1)); __cudaRegisterEntry(__T212, ((void ( *)(double *, float *, int))floatToDouble), _Z13floatToDoublePdPfi, (-1)); __cudaRegisterEntry(__T212, ((void ( *)(float *, double *, int))doubleToFloat), _Z13doubleToFloatPfPdi, (-1)); __cudaRegisterVariable(__T212, __shadow_var(ros,::ros), 0, 1280UL, 1, 0); }
static void __sti____cudaRegisterAll(void) {  __cudaRegisterBinary(__nv_cudaEntityRegisterCallback);  }

#pragma GCC diagnostic pop

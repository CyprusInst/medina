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

=#=#=#=#=#=#=#=#=#=#=defines_vars_2=#=#=#=#=#=#=#=#=#=#=

#define BLOCKSIZE 64


#define REDUCTION_SIZE_1 64
#define REDUCTION_SIZE_2 32

=#=#=#=#=#=#=#=#=#=#=defines_vars_1=#=#=#=#=#=#=#=#=#=#=

=#=#=#=#=#=#=#=#=#=#=defines_ind_1=#=#=#=#=#=#=#=#=#=#=

=#=#=#=#=#=#=#=#=#=#=defines_ind_2=#=#=#=#=#=#=#=#=#=#=

=#=#=#=#=#=#=#=#=#=#=defines_ind_3=#=#=#=#=#=#=#=#=#=#=

=#=#=#=#=#=#=#=#=#=#=defines_ind_4=#=#=#=#=#=#=#=#=#=#=

=#=#=#=#=#=#=#=#=#=#=defines_ind_5=#=#=#=#=#=#=#=#=#=#=

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

#define ASSOC 0
#define DISSOC 1
#define rconst(i,j)  rconst[(j)]

#ifdef REDUCE

/* reduced mem variant efficient memory access */
#define Ghimj(i,j) Ghimj[(j)*(VL_GLO)+(i)]
#define jac0(i,j) jac0[(j)*(VL_GLO)+(i)]
#define jcb(i,j) jcb[(j)*(VL_GLO)+(i)]
#define K(i,j,k) K[(j)*(VL_GLO)*(NVAR)+(k)*(VL_GLO)+(i)]

#else

#define jcb(i,j)     jcb[(j)]
#define K(i,j,k) K[(j)*(NVAR)+(k)]
#define jac0(i,j)    jac0[(j)]    
#define Ghimj(i,j)   Ghimj[(j)]

#endif


/* Temporary arrays allocated in stack */
#define var(i,j)     var[(j)]
#define fix(i,j)     fix[(j)]
#define varDot(i,j)  varDot[j]
#define varNew(i,j) varNew[(j)]
#define Fcn0(i,j)   Fcn0[(j)]
#define Fcn(i,j)    Fcn[(j)]
#define dFdT(i,j)   dFdT[(j)]
#define varErr(i,j) varErr[(j)]


/* Enable debug flags for GPU */
#define DEBUG

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

__device__ void  update_rconst(const REAL * __restrict__ var,
                               const REAL * __restrict__ khet_st, const REAL * __restrict__ khet_tr,
                               const REAL * __restrict__ jx, REAL * __restrict__ rconst,
                               const REAL * __restrict__ temp_gpu,
                               const REAL * __restrict__ press_gpu,
                               const REAL * __restrict__ cair_gpu,
			       const int VL_GLO);

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

__device__ static REAL k_arr(const REAL k_298,const REAL tdep, const REAL temp){
    // Arrhenius function
    return k_298 * exp(tdep*(1./temp-3.3540E-3));
}

__device__ static REAL k_limited(const REAL k3rd, const REAL cHp){
    // diffusion limitation caps 3rd order rate coefficients
    REAL const DiffLimit=1E10;
    return 1./(1./k3rd + cHp/DiffLimit);
}


__device__ static REAL k_siv_h2o2(const REAL k_298, const REAL tdep, const REAL chp, const REAL temp){
    return k_298 * exp(tdep * 1./temp-3.3540E-3) * chp / ( chp + 0.1);
}

__device__ static REAL * k_3rd_jpl_activation(const REAL temp, const REAL cair, const REAL k0_298K, const REAL n, const REAL kinf_298K, const REAL m, const REAL A, const REAL B){
    REAL zt_help, k0_TM, kinf_T, k_ratio, k_f, k_int, k_fCA;
    REAL static k_3rd_jpl[2];
    zt_help = 298./temp;
    k0_TM = k0_298K * pow(zt_help,n) * cair;
    kinf_T = kinf_298K * pow(zt_help, m);
    k_ratio = k0_TM/kinf_T;
    k_f = k0_TM/(1. + k_ratio)*pow(0.6,1./(1.+pow(log10(k_ratio),2)));
    k_int = A * exp(-B/temp);
    k_3rd_jpl[0] = k_f;
    k_3rd_jpl[1] = k_int * (1.-k_f/kinf_T);
    return k_3rd_jpl;
}

__device__ static REAL k_n2_o(const REAL temp, const REAL temp_ion){
    return 1.4E-10*pow(300./((temp_ion+temp)/2),0.44);
}

__device__ static REAL k_op_o2(const REAL temp, const REAL temp_ion){
    REAL temp_mean = 0.667*temp_ion + 0.333*temp;
    return 2.82E-11 - 7.74E-12*(temp_mean/300.) + 1.073E-12*pow(temp_mean/300.,2) - 5.17E-14*pow(temp_mean/300.,3) + 9.65E-16*pow(temp_mean/300.,4); 
}

__device__ static REAL k_op_n2(const REAL temp, const REAL temp_ion){
    REAL temp_mean = 0.6363*temp_ion + 0.3637*temp;
    return 2.91E-13 - 5.92E-13*(temp_mean/300.) + 8.6E-14*pow(temp_mean/300.,2);
}

__device__ static REAL uef(const REAL temp){
    return 8.4096792e+000 - 6.8484593e-002*temp + 2.3044184e-004*pow(temp,2) - 2.7257885e-007*pow(temp,3);
}

__device__ static REAL k_md (const REAL A, const REAL B, const REAL C1, const REAL C2, const REAL C3, const REAL pH, const REAL temp, const int SEL){
    REAL fac1,fac2,fac3,fac_back;
    fac1 = C1 * pow(10,-pH);
    fac2 = C2 - C3/temp * pow(10,pH);
    fac3 = A-(B/temp);
    fac_back = exp(1.449E-2+(5.609E2/temp));
    if (SEL == 1){ 
        return (1. + fac1 + fac2) * fac3;
    }
    else if (SEL == 2){
        return ((1. + fac1 + fac2) * fac3) / fac_back;
    }
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
            varMax = max(fabs(var[i]),fabs(varNew[i]));
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
            varMax = max(fabs(var[i]),fabs(varNew[i]));

            scale = absTol[0]+ relTol[0]*varMax;
            err += pow((REAL)varErr[i]/scale,2.0);
        }
        err  = sqrt((REAL) err/NVAR);
    }

    return err;


}

=#=#=#=#=#=#=#=#=#=#=kppSolve=#=#=#=#=#=#=#=#=#=#=

__device__ void ros_Solve(REAL * __restrict__ Ghimj, REAL * __restrict__ K, int &Nsol, const int istage, const int ros_S, int VL_GLO)
{

    int index = blockIdx.x*blockDim.x+threadIdx.x;
    #pragma unroll 4 
    for (int i=0;i<LU_NONZERO-16;i+=16){
        prefetch_ll1(&Ghimj(index,i));
    }

    kppSolve(Ghimj, K, istage, ros_S, VL_GLO);
    Nsol++;
}

=#=#=#=#=#=#=#=#=#=#=kppDecomp=#=#=#=#=#=#=#=#=#=#=

__device__ void ros_Decomp(REAL * __restrict__ Ghimj, int &Ndec, int VL_GLO)
{
    kppDecomp(Ghimj, VL_GLO);
    Ndec++;
}


=#=#=#=#=#=#=#=#=#=#=ros_PrepareMatrix=#=#=#=#=#=#=#=#=#=#=

=#=#=#=#=#=#=#=#=#=#=Jac_sp=#=#=#=#=#=#=#=#=#=#=

=#=#=#=#=#=#=#=#=#=#=Fun=#=#=#=#=#=#=#=#=#=#=

__device__ void ros_FunTimeDerivative(const REAL T, REAL roundoff, REAL * __restrict__ var, const REAL * __restrict__ fix,
                                      REAL * __restrict__ rconst, REAL *dFdT, REAL *Fcn0, int &Nfun,
                                      const REAL * __restrict__ khet_st, const REAL * __restrict__ khet_tr,
                                      const REAL * __restrict__ jx, const REAL * __restrict__ temp_gpu, const REAL * __restrict__ press_gpu, const REAL * __restrict__ cair_gpu,
                                      const int VL_GLO)
{
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    const REAL DELTAMIN = 1.0E-6;
    REAL delta,one_over_delta;

    delta = sqrt(roundoff)*max(DELTAMIN,fabs(T));
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
        const REAL FacMin, const REAL FacMax, const REAL FacRej, const REAL FacSafe, const REAL kk, const REAL bb,
        //  Status parameters
        int &Nfun, int &Njac, int &Nstp, int &Nacc, int &Nrej, int &Ndec, int &Nsol, int &Nsng,
        //  cuda global mem buffers
        REAL * __restrict__ rconst,  const REAL * __restrict__ absTol, const REAL * __restrict__ relTol, REAL * __restrict__ varNew, REAL * __restrict__ Fcn0,
        REAL * __restrict__ K, REAL * __restrict__ dFdT, REAL * __restrict__ jac0, REAL * __restrict__ Ghimj, REAL * __restrict__ varErr,
        // for update_rconst
        const REAL * __restrict__ khet_st, const REAL * __restrict__ khet_tr,
        const REAL * __restrict__ jx, const REAL * __restrict__ temp_gpu, const REAL * __restrict__ press_gpu, const REAL * __restrict__ cair_gpu,
        // VL_GLO
        const int VL_GLO)
{
    int index = blockIdx.x*blockDim.x+threadIdx.x;

    REAL H, Hnew, HC, HG, Fac; // Tau - not used
    REAL Err, ErrOld = 1., FacOld = 1.; //*varErr;
    int direction;
    int rejectLastH, rejectMoreH;
    const REAL DELTAMIN = 1.0E-5;

    //   ~~~>  Initial preparations
    T = Tstart;
    Hexit = 0.0;
    H = min(Hstart,Hmax);
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
        H = min(H,fabs(Tend-T));

        //   ~~~>   Compute the function at current time
        Fun(var, fix, rconst, Fcn0, Nfun, VL_GLO);	/// VAR READ - Fcn0 Write

        //   ~~~>  Compute the function derivative with respect to T
        if (!autonomous)
            ros_FunTimeDerivative(T, roundoff, var, fix, rconst, dFdT, Fcn0, Nfun, khet_st, khet_tr, jx, temp_gpu, press_gpu,cair_gpu,  VL_GLO); /// VAR READ - fcn0 read

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
                ros_Solve(Ghimj, K, Nsol, istage, ros_S, VL_GLO);


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
            Fac  = min(FacMax,max(FacMin,FacSafe/pow(Err,ONE/ros_ELO)));
            Hnew = H*Fac;

//  ~~~>  Check the error magnitude and adjust step size
            Nstp = Nstp+ 1;
            if((Err <= ONE) || (H <= Hmin)) // ~~~> Accept step
            {
                Nacc = Nacc + 1;
                for (int j=0; j<NVAR ; j++)
                    var(index,j) =  max(varNew(index,j),ZERO);  /////////// VAR WRITE - last VarNew read

                T = T +  direction*H;
                Hnew = max(Hmin,min(Hnew,Hmax));
                if (rejectLastH)   // No step size increase after a rejected step
                    Hnew = min(Hnew,H);
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
                //if (Nacc >= 1)
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

=#=#=#=#=#=#=#=#=#=#=update_rconst=#=#=#=#=#=#=#=#=#=#=


__global__ 
void Rosenbrock(REAL * __restrict__ conc, const REAL Tstart, const REAL Tend, REAL * __restrict__ rstatus, int * __restrict__ istatus,
                // values calculated from icntrl and rcntrl at host
                const int autonomous, const int vectorTol, const int UplimTol, const int method, const int Max_no_steps,
#ifdef REDUCE
                REAL * __restrict__ d_jac0, REAL * __restrict__ d_Ghimj, REAL * __restrict__ d_varNew, REAL * __restrict__ d_K, REAL * __restrict__ d_varErr,REAL * __restrict__ d_dFdT ,REAL * __restrict__ d_Fcn0,
#endif 
                const REAL Hmin, const REAL Hmax, const REAL Hstart, const REAL FacMin, const REAL FacMax, const REAL FacRej, const REAL FacSafe, const REAL roundoff, const REAL kk, const REAL bb,
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

#ifdef REDUCE
    REAL *Ghimj  = d_Ghimj;
    REAL *K      = d_K;
    REAL *varNew = &d_varNew[index*NVAR];
    REAL *Fcn0   = &d_Fcn0[index*NVAR];
    REAL *dFdT   = &d_dFdT[index*NVAR];
    REAL *jac0   = d_jac0;
    REAL *varErr = &d_varErr[index*NVAR];
#else 
    REAL varNew_stack[NVAR];
    REAL varErr_stack[NVAR];
    REAL Fcn0_stack[NVAR];
    REAL jac0_stack[LU_NONZERO];
    REAL dFdT_stack[NVAR];
    REAL Ghimj_stack[LU_NONZERO];
    REAL K_stack[6*NVAR];

    /* Allocated in stack */
    REAL *Ghimj  = Ghimj_stack;
    REAL *K      = K_stack;
    REAL *varNew = varNew_stack;
    REAL *Fcn0   = Fcn0_stack;
    REAL *dFdT   = dFdT_stack;
    REAL *jac0   = jac0_stack;
    REAL *varErr = varErr_stack;
#endif

    /* Temporary arrays allocated in stack */
    REAL var_stack[NSPEC];
    REAL fix_stack[NFIX];
    REAL rconst_stack[NREACT];
    REAL *var    = var_stack;
    REAL *fix    = fix_stack;
    REAL *rconst = rconst_stack;

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
                FacMin, FacMax, FacRej, FacSafe, kk, bb,
                //  Status parameters
                Nfun, Njac, Nstp, Nacc, Nrej, Ndec, Nsol, Nsng,
                //  cuda global mem buffers              
                rconst, absTol, relTol, varNew, Fcn0,  
                K, dFdT, jac0, Ghimj,  varErr, 
                // For update rconst
                khet_st, khet_tr, jx, temp_gpu, press_gpu, cair_gpu,
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


=#=#=#=#=#=#=#=#=#=#=special_ros=#=#=#=#=#=#=#=#=#=#=


/* Assuming different processes */
enum { TRUE=1, FALSE=0 } ;
double *d_conc, *d_temp, *d_press, *d_cair, *d_khet_st, *d_khet_tr, *d_jx; 
#ifdef REDUCE
REAL *d_jac0, *d_Ghimj, *d_varNew, *d_K, *d_varErr, *d_dFdT, *d_Fcn0;
#endif
#if defined(__SINGLEPREC)
float *d_conc_s, *d_temp_s, *d_press_s, *d_cair_s, *d_khet_st_s, *d_khet_tr_s, *d_jx_s;
float *d_rstatus_s, *d_absTol_s, *d_relTol_s;
#endif
int initialized = FALSE;

/* Device pointers pointing to GPU */
double *d_rstatus, *d_absTol, *d_relTol;
int *d_istatus;

/* Allocate arrays on device for Rosenbrock */
__host__ void init_first_time(int pe, int VL_GLO, int size_khet_st, int size_khet_tr, int size_jx ){

    /* Select the proper GPU CARD */
    int deviceCount, device;
    gpuErrchk( cudaGetDeviceCount(&deviceCount) );
    device = pe % deviceCount;
    gpuErrchk( cudaSetDevice(device) );

    printf("PE[%d]: selected %d of total %d\n",pe,device,deviceCount);
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

#ifndef _OPEN_ACC
    gpuErrchk( cudaMalloc ((void **) &d_conc   , sizeof(double)*VL_GLO*(NSPEC))        );
    gpuErrchk( cudaMalloc ((void **) &d_khet_st, sizeof(double)*VL_GLO*size_khet_st) );
    gpuErrchk( cudaMalloc ((void **) &d_khet_tr, sizeof(double)*VL_GLO*size_khet_tr) );
    gpuErrchk( cudaMalloc ((void **) &d_jx     , sizeof(double)*VL_GLO*size_jx)      );
    gpuErrchk( cudaMalloc ((void **) &d_istatus    , sizeof(int)*VL_GLO*8)             );
#endif

    gpuErrchk( cudaMalloc ((void **) &d_rstatus    , sizeof(double)*VL_GLO*2)          );

#ifndef _OPEN_ACC
    gpuErrchk( cudaMalloc ((void **) &d_absTol     , sizeof(double)*NVAR)              );
    gpuErrchk( cudaMalloc ((void **) &d_relTol     , sizeof(double)*NVAR)              );

    /* Allocate input arrays */
    gpuErrchk( cudaMalloc ((void **) &temp_gpu     , sizeof(double)*VL_GLO)              );
    gpuErrchk( cudaMalloc ((void **) &press_gpu     , sizeof(double)*VL_GLO)              );
    gpuErrchk( cudaMalloc ((void **) &cair_gpu     , sizeof(double)*VL_GLO)              );
#endif


#ifdef REDUCE
    /* Allocate arrays for solvers on device global memory to reduce the stack usage */
    gpuErrchk( cudaMalloc ((void **) &d_jac0, sizeof(REAL)*VL_GLO*LU_NONZERO)   );
    gpuErrchk( cudaMalloc ((void **) &d_Ghimj, sizeof(REAL)*VL_GLO*LU_NONZERO)   );
    gpuErrchk( cudaMalloc ((void **) &d_varNew, sizeof(REAL)*VL_GLO*NVAR)       );
    gpuErrchk( cudaMalloc ((void **) &d_Fcn0, sizeof(REAL)*VL_GLO*NVAR)       );
    gpuErrchk( cudaMalloc ((void **) &d_dFdT, sizeof(REAL)*VL_GLO*NVAR)       );

    gpuErrchk( cudaMalloc ((void **) &d_K, sizeof(REAL)*VL_GLO*NVAR*6)       );  // TODO: Change size according to solver steps
    gpuErrchk( cudaMalloc ((void **) &d_varErr, sizeof(REAL)*VL_GLO*NVAR)       );
#endif
#if defined(__SINGLEPREC)
    gpuErrchk( cudaMalloc ((void **) &d_conc_s   , sizeof(float)*VL_GLO*(NSPEC))        );
    gpuErrchk( cudaMalloc ((void **) &d_khet_st_s, sizeof(float)*VL_GLO*size_khet_st) );
    gpuErrchk( cudaMalloc ((void **) &d_khet_tr_s, sizeof(float)*VL_GLO*size_khet_tr) );
    gpuErrchk( cudaMalloc ((void **) &d_jx_s     , sizeof(float)*VL_GLO*size_jx)      );

    gpuErrchk( cudaMalloc ((void **) &d_rstatus_s    , sizeof(float)*VL_GLO*2)          );
    gpuErrchk( cudaMalloc ((void **) &d_absTol_s     , sizeof(float)*NVAR)              );
    gpuErrchk( cudaMalloc ((void **) &d_relTol_s     , sizeof(float)*NVAR)              );
    gpuErrchk( cudaMalloc ((void **) &temp_gpu_s     , sizeof(float)*VL_GLO)              );
    gpuErrchk( cudaMalloc ((void **) &press_gpu_s     , sizeof(float)*VL_GLO)              );
    gpuErrchk( cudaMalloc ((void **) &cair_gpu_s     , sizeof(float)*VL_GLO)              );
#endif


    initialized = TRUE;
}

/*
 * TODO: We should call it in some point..
 */
extern "C" void finalize_cuda(){
    /* Free memory on the device */
#ifndef _OPEN_ACC
    gpuErrchk( cudaFree(d_conc        ) );
    gpuErrchk( cudaFree(d_khet_st     ) );
    gpuErrchk( cudaFree(d_khet_tr     ) );
    gpuErrchk( cudaFree(d_jx          ) );
    gpuErrchk( cudaFree(d_absTol      ) );
    gpuErrchk( cudaFree(d_relTol      ) );
    gpuErrchk( cudaFree(temp_gpu      ) );
    gpuErrchk( cudaFree(press_gpu     ) );
    gpuErrchk( cudaFree(cair_gpu      ) );
#endif
    gpuErrchk( cudaFree(d_temp        ) );
    gpuErrchk( cudaFree(d_press       ) );
    gpuErrchk( cudaFree(d_cair        ) );
    gpuErrchk( cudaFree(d_rstatus     ) );
    gpuErrchk( cudaFree(d_istatus     ) );

#ifdef REDUCE
    gpuErrchk( cudaFree(d_jac0        ) );
    gpuErrchk( cudaFree(d_Ghimj       ) );
    gpuErrchk( cudaFree(d_varNew      ) );
    gpuErrchk( cudaFree(d_Fcn0        ) );
    gpuErrchk( cudaFree(d_dFdT        ) );
    gpuErrchk( cudaFree(d_K           ) );
    gpuErrchk( cudaFree(d_varErr      ) );
#endif
#if defined(__SINGLEPREC)
    gpuErrchk( cudaFree(d_conc_s        ) );
    gpuErrchk( cudaFree(d_khet_st_s     ) );
    gpuErrchk( cudaFree(d_khet_tr_s     ) );
    gpuErrchk( cudaFree(d_jx_s          ) );
    gpuErrchk( cudaFree(d_rstatus_s     ) );
    gpuErrchk( cudaFree(d_absTol_s      ) );
    gpuErrchk( cudaFree(d_relTol_s      ) );
#endif
}



extern "C" void kpp_integrate_cuda_( int *pe_p, int *sizes, double *time_step_len_p, double *conc, double *temp, double *press, double *cair, 
                                    double *khet_st, double *khet_tr, double *jx, double *absTol, double *relTol, int *ierr, int *istatus, 
                                    double *rndoff, int *icntrl=NULL, double *rcntrl=NULL
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
 *        icntrl: 1D array of ints   - size: 20
 *         sizes: 1D array of ints   - size: 4
 *        rcntrl: 1D array of doubles - size: 20
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
    REAL Hmin, Hmax, Hstart, FacMin, FacMax, FacRej, FacSafe, kk, bb;
    
    //int rcntrl_bool = 0, icntrl_bool=0;
    if (rcntrl == NULL)
    {
        rcntrl = new double[10];
        for (int i=0; i < 10; i++)
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
#ifdef _OPEN_ACC

    /* Pointers to device memory */
    d_conc = conc;
    temp_gpu = temp;
    press_gpu = press;
    cair_gpu = cair;
    d_khet_st = khet_st;
    d_khet_tr = khet_tr;
    d_jx = jx;

    /* Pointers to device memory for Rosenbrock */
    d_absTol = absTol;
    d_relTol = relTol;
    d_istatus = istatus
#else
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
#endif

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
            Hmin = (REAL) rcntrl[0];
        }
        else
        {
            printf("User-selected Hmin: rcntrl[0]=%f\n",rcntrl[0]);
            ierr_tmp = -3;
        }
    //  Upper bound on the step size: (positive value)
        if (rcntrl[1] == ZERO) 
        {
            Hmax = (REAL) fabs(Tend-Tstart);
        }
        else if (rcntrl[1] > ZERO) 
        {
            Hmax = (REAL) fmin(fabs(rcntrl[1]),fabs(Tend-Tstart));
        }
        else
        {
            printf("User-selected Hmax: rcntrl[1]=%f\n",rcntrl[1]);
            ierr_tmp = -3;
        }
    //  Starting step size: (positive value)
        if (rcntrl[2] == ZERO) 
        {
            Hstart = (REAL) fmax(Hmin,DELTAMIN);
        }
        else if (rcntrl[2] > ZERO) 
        {
            Hstart = (REAL) fmin(fabs(rcntrl[2]),fabs(Tend-Tstart));
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
            FacMin = (REAL) rcntrl[3];
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
            FacMax = (REAL) rcntrl[4];
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
            FacRej = (REAL) rcntrl[5];
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
            FacSafe = (REAL) rcntrl[6];
        }
        else
        {
            printf("User-selected FacSafe: rcntrl[6]=%f\n",rcntrl[6]);
            ierr_tmp = -4;
        }
    // Values for aggrasive time stepping
        if (rcntrl[7] == ZERO)
        {
            bb = 1.0;
        }
        else if (rcntrl[7] > ZERO)
        {
            bb = (REAL) rcntrl[7];
        }
        else
        {
            printf("User-selected bb: rcntrl[7]=%f\n",rcntrl[7]);
            ierr_tmp = -4;
        }
        if (rcntrl[8] == ZERO)
        {
            kk = 2.5;
        }
        else if (rcntrl[8] > ZERO)
        {
            kk = (REAL) rcntrl[8];
        }
        else
        {
            printf("User-selected kk: rcntrl[8]=%f\n",rcntrl[8]);
            ierr_tmp = -4;
        }

#ifndef _OPEN_ACC
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
#endif
    }
    // set return value of error code
    *ierr = ierr_tmp;

    =#=#=#=#=#=#=#=#=#=#=call_kernel=#=#=#=#=#=#=#=#=#=#=

    GPU_DEBUG();

/* Copy the result back */
#if defined(__SINGLEPREC)
    floatToDouble<<<dimGrid,dimBlock>>>(d_conc,d_conc_s,VL_GLO*NSPEC);
#endif
#ifndef _OPEN_ACC
    gpuErrchk( cudaMemcpy( conc      , d_conc       , sizeof(double)*VL_GLO*NVAR, cudaMemcpyDeviceToHost) );  
    gpuErrchk( cudaMemcpy( istatus   , d_istatus    , sizeof(int)*VL_GLO*8      , cudaMemcpyDeviceToHost) ); 
#else
    gpuErrchk( cudaDeviceSynchronize() );
#endif
    
    return;

}






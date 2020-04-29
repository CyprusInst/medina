
__device__  static  int ros_Integrator_rodas3(double * __restrict__ var, const double * __restrict__ fix, const double Tstart, const double Tend, double &T,

        //  Integration parameters
        const int autonomous, const int vectorTol, const int Max_no_steps, 
        const double roundoff, const double Hmin, const double Hmax, const double Hstart, double &Hexit, 
        const double FacMin, const double FacMax, const double FacRej, const double FacSafe, 
        //  Status parameters
        int &Nfun, int &Njac, int &Nstp, int &Nacc, int &Nrej, int &Ndec, int &Nsol, int &Nsng,
        //  cuda global mem buffers              
        const double * __restrict__ rconst,  const double * __restrict__ absTol, const double * __restrict__ relTol, double * __restrict__ varNew, double * __restrict__ Fcn0, 
        double * __restrict__ K, double * __restrict__ dFdT, double * __restrict__ jac0, double * __restrict__ Ghimj, double * __restrict__ varErr,
        // for update_rconst
        const double * __restrict__ khet_st, const double * __restrict__ khet_tr,
        const double * __restrict__ jx,
        // VL_GLO
        const int VL_GLO)
{
    int index = blockIdx.x*blockDim.x+threadIdx.x;

    double H, Hnew, HC, HG, Fac; // Tau - not used
    double Err; //*varErr;
    int direction;
    int rejectLastH, rejectMoreH;
    const double DELTAMIN = 1.0E-5;
    const int ros_S = 4; 


    /* Not sure if it worth it for shared */
    double ros_A[15];
    double ros_C[15];
    int   ros_NewF[8];
    double ros_M[6];
    double ros_E[6];
    double ros_Alpha[6];
    double ros_Gamma[6];


    ros_A[0] = 0.0E+00;
    ros_A[1] = 2.0E+00;
    ros_A[2] = 0.0E+00;
    ros_A[3] = 2.0E+00;
    ros_A[4] = 0.0E+00;
    ros_A[5] = 1.0E+00;

    ros_C[0] = 4.0E+00;
    ros_C[1] = 1.0E+00;
    ros_C[2] =- 1.0E+00;
    ros_C[3] = 1.0E+00;
    ros_C[4] =- 1.0E+00;
    ros_C[5] =- (8.0E+00/3.0E+00);

    ros_NewF[0]  = 1;
    ros_NewF[1]  = 0;
    ros_NewF[2]  = 1;
    ros_NewF[3]  = 1;

    ros_M[0] = 2.0E+00;
    ros_M[1] = 0.0E+00;
    ros_M[2] = 1.0E+00;
    ros_M[3] = 1.0E+00;

    ros_E[0] = 0.0E+00;
    ros_E[1] = 0.0E+00;
    ros_E[2] = 0.0E+00;
    ros_E[3] = 1.0E+00;

    ros_Alpha[0] = 0.0E+00;
    ros_Alpha[1] = 0.0E+00;
    ros_Alpha[2] = 1.0E+00;
    ros_Alpha[3] = 1.0E+00;

    ros_Gamma[0] = 0.5E+00;
    ros_Gamma[1] = 1.5E+00;
    ros_Gamma[2] = 0.0E+00;
    ros_Gamma[3] = 0.0E+00;


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
            ros_PrepareMatrix(H, direction, 0.5E+00, jac0, Ghimj, Nsng, Ndec, VL_GLO);
            //   ~~~>   Compute the stages
            // Stage:
            {
                int istage = 0;
		for (int i=0; i<NVAR; i++)		
			K(index,istage,i)  = varNew(index,i) = Fcn0(index,i);

                if ((!autonomous))
                {
                    HG = direction*H*ros_Gamma[istage];
                    for (int i=0; i<NVAR; i++){
                        K(index,istage,i) += dFdT(index,i)*HG;
		     }
                }
		//	   R   ,RW, RW,  R,        R 
                ros_Solve(Ghimj, K, Nsol, istage, ros_S);
            } // Stage

           #pragma unroll
           for (int istage=1; istage < ros_S; istage++)
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
				double tmp = K(index,j,i);
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
		    double tmpNew  = var(index,i); 					/// VAR READ
		    double tmpErr  = ZERO;

                    #pragma unroll
		    for (int j=0; j<ros_S; j++){
		    	    double tmp = K(index,j,i);
			    tmpNew += tmp*ros_M[j];
			    tmpErr += tmp*ros_E[j];
		    }
		    varNew(index,i) = tmpNew;			// varNew is killed
		    varErr(index,i) = tmpErr;
	    }

            Err = ros_ErrorNorm(var, varNew, varErr, absTol, relTol, vectorTol);   /// VAR-varNew READ


//  ~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
            Fac  = fmin(FacMax,fmax(FacMin,FacSafe/pow(Err,ONE/3.0)));
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

__global__ 
void Rosenbrock_rodas3(double * __restrict__ conc, const double Tstart, const double Tend, double * __restrict__ rstatus, int * __restrict__ istatus,
                const int autonomous, const int vectorTol, const int UplimTol, const int Max_no_steps,
                double * __restrict__ d_jac0, double * __restrict__ d_Ghimj, double * __restrict__ d_varNew, double * __restrict__ d_K, double * __restrict__ d_varErr,double * __restrict__ d_dFdT ,double * __restrict__ d_Fcn0, double * __restrict__ d_var, double * __restrict__ d_fix, double * __restrict__ d_rconst,
                const double Hmin, const double Hmax, const double Hstart, const double FacMin, const double FacMax, const double FacRej, const double FacSafe, const double roundoff,
                const double * __restrict__ absTol, const double * __restrict__ relTol,
    	        const double * __restrict__ khet_st, const double * __restrict__ khet_tr,
		const double * __restrict__ jx,
                const double * __restrict__ temp_gpu,
                const double * __restrict__ press_gpu,
                const double * __restrict__ cair_gpu,
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
    double *Ghimj  = &d_Ghimj[index*LU_NONZERO];    
    double *K      = &d_K[index*NVAR*3];
    double *varNew = &d_varNew[index*NVAR];
    double *Fcn0   = &d_Fcn0[index*NVAR];
    double *dFdT   = &d_dFdT[index*NVAR];
    double *jac0   = &d_jac0[index*LU_NONZERO];
    double *varErr = &d_varErr[index*NVAR];
    double *var    = &d_var[index*NSPEC];
    double *fix    = &d_fix[index*NFIX];
    double *rconst = &d_rconst[index*NREACT];

    if (index < VL_GLO)
    {

        int Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng;
        double Texit, Hexit;

        Nfun = 0;
        Njac = 0;
        Nstp = 0;
        Nacc = 0;
        Nrej = 0;
        Ndec = 0;
        Nsol = 0;
        Nsng = 0;
    
        const int method = 4; 

        /* FIXME: add check for method */
        const double *ros_A     = &ros[method-1].ros_A[0]; 
        const double *ros_C     = &ros[method-1].ros_C[0];
        const double *ros_M     = &ros[method-1].ros_M[0]; 
        const double *ros_E     = &ros[method-1].ros_E[0];
        const double *ros_Alpha = &ros[method-1].ros_Alpha[0]; 
        const double *ros_Gamma = &ros[method-1].ros_Gamma[0]; 
        const int    *ros_NewF  = &ros[method-1].ros_NewF[0];
        const int     ros_S     =  ros[method-1].ros_S; 
        const double  ros_ELO   =  ros[method-1].ros_ELO; 


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

        ros_Integrator_rodas3(var, fix, Tstart, Tend, Texit,
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





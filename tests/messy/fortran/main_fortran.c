/*
 *   Example input data for testing the verification of the transformation.
 *
 *   Copyright 2016 The Cyprus Institute 
 *
 *   Developers: Michail Alvanos - m.alvanos@cyi.ac.cy
 *               Theodoros Christoudias - christoudias@cyi.ac.cy
 *
 * */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#define NSPEC 74
#define NREACT 77
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


//          ,-- 2 Leading underscores to start
//          | ,-- then the module name
//          | |               ,-- then _MOD_
//          | |               |    ,-- then the subroutine name
//          V V               V    V
extern void __messy_mecca_kpp_MOD_kpp_integrate();
extern void __messy_mecca_kpp_MOD_kpp_integrate_init();

#define INIT_TABLES(){\
    for (i=0;i<VL_GLO;i++){\
        temp[i] = 240.995971972245;\
        press[i] = 0.994591236114502; \
        cair[i] = 298914285136738.0;\
        khet_tr[i*4 + 0] = 7.408371201503456E-008;\
        khet_tr[i*4 + 1] = 4.849455570110968E-007;\
        khet_tr[i*4 + 2] =  0.000000000000000E+000;  \
        khet_tr[i*4 + 3] = 2.718003287797325E-007;\
    }\
    for (i=0;i<NSPEC;i++){\
        abstol[i] = 10.0; \
        reltol[i] = 0.5; \
    }\
    for (i=0;i<VL_GLO;i++){\
        for (j=0;j<NSPEC;j++){\
              conc[j*VL_GLO + i] = conc_cell[j];\
        }\
    }\
    __messy_mecca_kpp_MOD_kpp_integrate_init( &vl_glo, temp, press, cair, jx, khet_tr, abstol, reltol, icntrl);\
}

#define PRINT_DATA()\
        printf("Results:");    \
        for (j=0;j<NSPEC;j++){\
            printf("   %.12e  ",conc[j*VL_GLO]);\
        }\
        printf("\n");    


int main(int argc, char **argv){
    
    int n = 1; // proccess element
    int istatus[20];
    int ierr;
    int i,j;
    int PE = 0;
    int l_debug = 0;
    int vl_glo = VL_GLO;
    int sizes[4] = {VL_GLO,NSPEC,NSPEC,IPMAX}; 
    int icntrl[20] = {0,0,2,0};
    double roundoff = 2.220446049250313E-016;
    double timestep = 720.0;

    struct timeval start, end;




    if (argc==2){
        icntrl[2] = atoi(argv[1]);
        INIT_TABLES()
        gettimeofday(&start, NULL);
        __messy_mecca_kpp_MOD_kpp_integrate( &timestep, conc, &ierr, xNacc, xNrej, istatus, &l_debug, &PE);
        gettimeofday(&end, NULL);
        printf("%d: %ld (ms)\n", icntrl[2],((end.tv_sec * 1000 + end.tv_usec/1000) - (start.tv_sec * 1000 + start.tv_usec/1000)));
        PRINT_DATA();
        return 0;
    }



    icntrl[2] = 1;

restart:

    INIT_TABLES()

        gettimeofday(&start, NULL);

    __messy_mecca_kpp_MOD_kpp_integrate( &timestep, conc, &ierr, xNacc, xNrej, istatus, &l_debug, &PE);

    gettimeofday(&end, NULL);

    printf("%d: %ld (ms)\n", icntrl[2],((end.tv_sec * 1000 + end.tv_usec/1000)
                - (start.tv_sec * 1000 + start.tv_usec/1000)));

    PRINT_DATA();

    icntrl[2]++;
    if ( icntrl[2] >5) return 0;
    goto restart;



}








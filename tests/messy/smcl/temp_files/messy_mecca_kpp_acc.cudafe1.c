# 1 "messy_mecca_kpp_acc.cu"
# 3844 "messy_mecca_kpp_acc.cu"
double *temp_gpu = 0;
double *press_gpu = 0;
double *cair_gpu = 0;


float *temp_gpu_s = 0;
float *press_gpu_s = 0;
float *cair_gpu_s = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_conc = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_temp = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_press = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_cair = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_khet_st = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_khet_tr = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_jx = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_jac0 = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_Ghimj = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_varNew = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_K = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_varErr = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_dFdT = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_Fcn0 = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_var = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_fix = 0;
# 4346 "messy_mecca_kpp_acc.cu"
double *d_rconst = 0;

float *d_conc_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_temp_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_press_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_cair_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_khet_st_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_khet_tr_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_jx_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_jac0_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_Ghimj_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_varNew_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_K_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_varErr_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_dFdT_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_Fcn0_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_var_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_fix_s = 0;
# 4348 "messy_mecca_kpp_acc.cu"
float *d_rconst_s = 0;
float *d_rstatus_s = 0;
# 4349 "messy_mecca_kpp_acc.cu"
float *d_absTol_s = 0;
# 4349 "messy_mecca_kpp_acc.cu"
float *d_relTol_s = 0;

extern int initialized;


double *d_rstatus = 0;
# 4354 "messy_mecca_kpp_acc.cu"
double *d_absTol = 0;
# 4354 "messy_mecca_kpp_acc.cu"
double *d_relTol = 0;
int *d_istatus = 0;
# 4355 "messy_mecca_kpp_acc.cu"
int *d_istatus_rd = 0;
# 4355 "messy_mecca_kpp_acc.cu"
int *d_xNacc = 0;
# 4355 "messy_mecca_kpp_acc.cu"
int *d_xNrej = 0;
struct int4 *d_tmp_out_1 = 0;
# 4356 "messy_mecca_kpp_acc.cu"
struct int4 *d_tmp_out_2 = 0;
# 4879 "messy_mecca_kpp_acc.cu"
double conc[426240];
double temp[5760];
double press[5760];
double cair[5760];
double jx[426240];


int xNacc[5760];
int xNrej[5760];

extern double conc_cell[74];
# 4966 "messy_mecca_kpp_acc.cu"
extern double abstol[74];



extern double reltol[74];




extern double khet_st[426240];



extern double khet_tr[426240];
# 4351 "messy_mecca_kpp_acc.cu"
int initialized = 0;
# 4889 "messy_mecca_kpp_acc.cu"
double conc_cell[74] = {(0.0),(0.0),(1.130030837133365E-6),(2161.17681825926),(1.469481417859824E-4),(2.89406754649778E-4),(0.0),(0.0),(6.377486492629032E-31),(2.774602114035594E-4),(9.159068418074058E-22),(1.681545841334171E-30),(6.587848965925121E-36),(4.057130203198298E-31),(7.556675262619906E-6),(5.625822089563362E-6),(7.24854650834698E-10),(7.771754415762507E-39),(1.672965892516881E-32),(5.778276640099593E-29),(2.16962319659931E-31),(4.44968552491389E-29),(9.23699185317872E-28),(1.731254847935413E-9),(6.419363370200839E-28),(4.035724058634079E-29),(6234.08726448302),(25802.7788132849),(1.33974252411334),(11.1514176946459),(8.023966161170008E-32),(1.405402576145367E-30),(2.416365419045456E-29),(3.763980220765519E-33),(3.687747273615521E-4),(4.400695805857555E-30),(8.096351349854847E-9),(1.60577739654151E-8),(8.424266813161654E-5),(1.275728897910597E-29),(36780.6069067007),(44.2802185584881),(5.485594561042764E-10),(3.41823488598684E-32),(1.808885697309332E-8),(2.295321288609202E-30),(7.186736555958003E-32),(6.67193926549068E8),(9.443976722997098E-30),(2.06547975096585E-30),(6.58798139717353E8),(5013220.8292721),(6.594652607797343E-13),(4.779051920325237E-33),(0.241330392051758),(2.657031589287186E-30),(1.166890334972386E-14),(337.069782231658),(126494.977205691),(891.196915201611),(222.557367243832),(1.22451624669813),(4845.02754823106),(535329.616196368),(0.03077774956209536),(9.89833722937206E8),(38527.6291432442),(1.857293910861109E-7),(5.03561600244018E9),(2.68242473107905E7),(211466.239175163),(6.06381297678027E13),(2.25227339137553E14),(8.76514082411165E10)};
# 4966 "messy_mecca_kpp_acc.cu"
double abstol[74] = {(0.0)};



double reltol[74] = {(0.0)};




double khet_st[426240] = {(0.0)};



double khet_tr[426240] = {(0.0)};

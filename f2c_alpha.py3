#!/usr/bin/env python3
#########################################################################################################
#
# MECCA - KPP Fortran to CUDA parser
#
# Copyright 2016-2020 The Cyprus Institute
#
# Developers: Michail Alvanos - m.alvanos@cyi.ac.cy
#             Theodoros Christoudias - christoudias@cyi.ac.cy
#             Giannis Ashiotis
#
#########################################################################################################


import os
import shutil
import re
import subprocess, string
import argparse

smcl = "../../smcl/"

def remove_comments(source):
    print("Removing comments...")
    lines = source[:]
    lines.reverse()
    out = []
    while True:
        if lines == []:
            break
        line = lines.pop()
        if "!" in line:
            if line[0] == "!":
                continue
            line = line[:line.find("!")-1]+"\n"
            line = line.strip()
        if (line != ''):
            out.append(line)
    return out

def strip_and_unroll_lines(source):
    lines = source[:]
    lines.reverse()
    out = []
    while True:
        if lines == []:
            break
        line = lines.pop()
        line = line.strip()
        if line != "":
            if line[-1] == "&":
                line = line[0:-1].strip() + " "
                while True:
                    next_line = lines.pop()
                    next_line = next_line.strip()
                    if next_line!= "" and next_line[0] == "&":
                        next_line = " " + next_line[1:].strip()
                    if "&" in next_line:
                        line = line + next_line[:-1]
                    else:
                        line = line + next_line
                        break
        line = line + "\n"
        out.append(line)
    return out

def find_subroutines(file_in, subroutine_names):
    subroutines = {}
    for subroutine in subroutine_names:
        subroutine = subroutine.lower()
        file_in.seek(0)
        lines = file_in.readlines()
        lines.reverse()
        source = []
        while True:
            if lines == []:
                break
            line = lines.pop().lower()
            if ( (("subroutine "+subroutine + " ") in line.lower()) or (("subroutine "+subroutine + "(") in line.lower())  ):
                while True:
                    if lines == []:
                        break
                    line = lines.pop()
                    if ( ("subroutine "+subroutine) in line.lower()):
                        break
                    source.append(line)
                break
        subroutines[subroutine.strip()] = source
    return subroutines

def decapitalize_vars(source,keys):
    fixed=[]
    for i in range(len(source)):
        line = source[i]
        for key in keys:
            key = key.lower()
            key_len = len(key)
            if key in line.lower():
                index = 0
                while True:
                    index = line.lower().find(key, index)
                    if index == -1:
                        break
                    line = line[:index]+key+line[index+key_len:]
                    index = index + key_len
        fixed.append(line)
    return fixed



def fix_power_op(source):
    operators = "*/+-<>=.,"
    # can use this in the future:
    #(\(?[0-9,a-z,A-Z,_,+,-,/,., ,)]+\)?)\*\*(\(?[-,\d,+,.,a-zA-Z_,/, ,\t]+\)?)
    #
    var_name = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_'
    precision_qualifiers = ["_sp","_dp","_qp"]
    fixed = []
    for line in source:
        if "**" in line:
            index = len(line)
            while True:
                index = line.rfind("**",0, index-1)
                if index == -1:
                    break
                left = line[:index].strip()
                right = line[index+2:].strip()
                if right[0]=="(":
                    left_bra=0
                    right_bra = 0
                    pos = 0
                    for i in right:
                        if i == "(":
                            left_bra = left_bra + 1
                        elif i == ")":
                            right_bra = right_bra + 1
                        if left_bra == right_bra:
                            break
                        pos = pos + 1
                    exponent = right[1:pos]
                    right = right[pos+1:].strip()
                else:
                    exponent = right
                    right = ""

                if left[-1]==")": #if it's (...) or a3_(...)
                    left_bra=0
                    right_bra = 0
                    pos = 0
                    for i in reversed(left):
                        if i == "(":
                            left_bra = left_bra + 1
                        elif i == ")":
                            right_bra = right_bra + 1
                        if left_bra == right_bra:
                            break
                        pos = pos - 1
                    base = left[pos:-1].strip()
                    left = left[:pos-1].strip()
                    if left[-1] in var_name: #if it's a3_(...)
                        pos = 0
                        for i in reversed(left):
                            if i not in var_name:
                                break
                            pos = pos - 1
                        base = left[pos:].strip()+"("+base+")"
                        left = left[:pos].strip()
                elif left[-3:].lower() in precision_qualifiers: # if it's 33._dp, 33.33E-33_dp or a_333_dp
                    if left[-4] == ".": # if it's 33._dp
                        pos=-4
                        for i in reversed(left[:-4]):
                            if i not in string.digits:
                                if pos == -4:
                                    print("ERROR 1 in \"parser\" \n")
                                break
                            pos = pos -1
                        base = left[pos:-3]
                        left = left[:pos].strip()
                    elif left[-4] in var_name:  # if it's 33.33E-33_dp or a_333_dp
                        pos=-3
                        scientif = False
                        isnumber = True
                        for i in reversed(left[:-3]):
                            if i not in string.digits+".":
                                if (i in "+-") and abs(pos-2)<=len(left) and (left[pos-2].lower() == "e"):
                                    scientif = True
                                elif scientif and i.lower() == "e":
                                    pass
                                elif i in string.letters+"_": # if it's a_333_dp
                                    pos=-3
                                    for i in reversed(left[:-3]):
                                        if i not in var_name:
                                            break
                                        pos = pos - 1
                                    isnumber = False
                                    break
                                else:
                                    break
                            pos = pos - 1
                        base = left[pos:-3]+((not isnumber) and (not scientif))*left[-3:]
                        left = left[:pos].strip()
                else: # if it's 3. , 3.3E-33 or a_3
                    if left[-1] == ".": # if it's 33.
                        pos=-1
                        for i in reversed(left[:-1]):
                            if i not in string.digits:
                                if pos == -1:
                                    print("ERROR 2 in \"parser\" \n")
                                break
                            pos = pos - 1
                        base = left[pos:]
                        left = left[:pos].strip()
                    elif left[-1] in var_name:   # if it's 33.33E-33 or a_3
                        pos=0
                        scientif = False
                        isnumber = True
                        for i in reversed(left):
                            if i not in string.digits+".":
                                if (i in "+-") and abs(pos-2)<=len(left) and (left[pos-2].lower() == "e"):
                                    scientif = True
                                elif scientif and i.lower() == "e":
                                    pass
                                elif i in string.letters+"_": # if it's a_3
                                    pos=0
                                    for i in reversed(left):
                                        if i not in var_name:
                                            break
                                        pos = pos - 1
                                    isnumber = False
                                    break
                                else:
                                    break
                            pos = pos - 1
                        #base = left[pos:]+((not isnumber) and (not scientif))*left[:]
                        base = left[pos:]
                        left = left[:pos].strip()
                    else:
                        print("OOPS! Missed something...")
                line = left+" pow("+base+", "+ exponent+") "+right+"\n"
        fixed.append(line)
    return fixed

def strip_lines(source):
    stripped = []
    for line in source:
        stripped.append(line.strip())
    return stripped

def split_beta(source,key):
    the_line = 0
    for line_num in range(len(source)):
        if source[line_num].strip()[0:len(key)] == key:
            the_line = line_num
            break
    main_body = []
    for line_num in range(the_line,len(source)):
        if "return" in source[line_num].lower():
            source[line_num] = ""
        main_body.append(source[line_num])
    return main_body

def fix_indices(source,keys):
    lines = source[:]
    lines.reverse()
    processed = []
    while True:
        if lines == []:
            break
        line = lines.pop()
        for var in keys:
            if var[0]+"(" in line:
                index = 0
                while True:
                    index = line.find(var[0]+"(", index)
                    if index == -1:
                        break
                    if len(var) == 2:
                        extra = ""
                    elif len(var) == 3:
                        extra = var[2]+","
                    else:
                        raise ValueError
                    if  ((line[index+len(var[0])+1:line.find(")",index)]).strip().isdigit()):
                        line = line[:index]+var[1]+"(index,"+extra+str(int(line[index+len(var[0])+1:line.find(")",index)])-1)+line[line.find(")",index):].replace(")=",") =")
                    else:

                        print("Value error : "+ str(line[index+len(var[0])+1:line.find(")",index)]))
                        raise ValueError

                    index = index + len(var[1])+1
        processed.append(line.strip())
    return processed


pass
#########################################################################################################
#########################################################################################################


def split_rconst(source):
    lines = source[:]
    lines.reverse()
    rconst_decls = []
    rconst_ops = []
    while True:
        if lines == []:
            break
        line = lines.pop()
        if line[0:6].lower() == "rconst":
            rconst_ops.append(line)
        elif ("=" in line) and (line.find("=") > 2):
            rconst_decls.append(line)
    return rconst_ops,rconst_decls

def rconst_preprocessor_1(source):
    lines = source[:]
    lines.reverse()
    file_temp = open("file_temp.c","w")
    file_temp.write("#define rconst(i) rconst(index,i)\n")
    file_temp.write("#define jx(i,j) jx(index,j)\n")
    file_temp.write("#define khet_st(i,j) khet_st(index,j)\n")
    file_temp.write("#define khet_tr(i,j) khet_tr(index,j)\n")
    file_temp.write("#define exp(i) exp(i)\n")
    file_temp.write("#define C(i) var[i]\n")
    file_temp.write("#define c(i) var[i]\n")
    file_temp.write("#define REAL( i, SP) (i)\n")
    file_temp.write("#define temp(i) temp_loc\n")
    file_temp.write("#define cair(i) cair_loc\n")
    file_temp.write("#define press(i) press_loc\n")
    file_temp.write("#define log(i) log(i)\n")
    while True:
        if lines == []:
            break
        line = lines.pop()
        if "rconst" in line:
            index = 0
            while True:
                index = line.find("rconst(", index)
                if index == -1:
                    break
                line = line[:index+7]+str(int(line[index+7:line.find(")",index)])-1)+line[line.find(")",index):]
                index = index + 7
        file_temp.write(line)
    file_temp.close()
    file_temp = open("file_temp2.c","w")
    p1 = subprocess.Popen(["gcc","-E","file_temp.c"], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["grep","-v","\#"], stdin=p1.stdout, stdout=file_temp)
    p2.wait()
    file_temp.close()
    file_temp = open("file_temp2.c","r")
    file_temp.seek(0)
    preprocessed = file_temp.readlines()
    file_temp.close()
    subprocess.call(["rm","-f","file_temp.c","file_temp2.c"])
    return preprocessed

def get_rconst_locals(source):
    rconst_locals=[]
    for line in source:
        if "=" in line:
            if "IF" not in line:
                rconst_locals.append(line.split("=",1)[0].strip())
    return rconst_locals


def create_rconst_init(source):
    rconst_init=[]
    for line in source:
        if "rconst" in line:
            eline = line.replace("\n","\n")
            eline = re.sub(r"rconst(\([0-9]+\))",r"rconst(index,\1-1)",eline)
            rconst_init.append( eline )
    return rconst_init

def generate_update_rconst(rconst_ops,rconst_decls,locals,rcint):
    update_rconst = []
    rename_tmp = False

    # these are for renaming the rconst ops
    for line in rconst_ops:
        if ( "temp)" in line):
            rename_tmp = True
            break

    if (rename_tmp == True):
        rconst_ops = [w.replace('temp', 'temp_loc') for w in rconst_ops]
        rconst_decls = [w.replace('temp', 'temp_loc') for w in rconst_decls]
        rconst_ops = [w.replace('press', 'press_loc') for w in rconst_ops]
        rconst_decls = [w.replace('press', 'press_loc') for w in rconst_decls]
        rconst_ops = [w.replace('cair', 'cair_loc') for w in rconst_ops]
        rconst_decls = [w.replace('cair', 'cair_loc') for w in rconst_decls]



    update_rconst.append( \
    "__device__ void  update_rconst(const double * __restrict__ var, \n \
			       const double * __restrict__ khet_st, const double * __restrict__ khet_tr,\n \
			       const double * __restrict__ jx, double * __restrict__ rconst, \n\
			       const double * __restrict__ temp_gpu, \n\
			       const double * __restrict__ press_gpu, \n\
			       const double * __restrict__ cair_gpu, \n\
			       const int VL_GLO)\n")
    update_rconst.append("{\n")
    update_rconst.append("    int index = blockIdx.x*blockDim.x+threadIdx.x;\n\n")
    update_rconst.append("    /* Set local buffer */\n")
    update_rconst.append("\n")
    update_rconst.append("    {\n")
    update_rconst.append("        const double temp_loc  = temp_gpu[index];\n")
    update_rconst.append("        const double press_loc = press_gpu[index];\n")
    update_rconst.append("        const double cair_loc  = cair_gpu[index];\n")
    update_rconst.append("\n")
    line = "        double"
    for i in locals:
        line = line+" "+i+","
    line = line[:-1]+";\n"
    update_rconst.append(line)
    update_rconst.append("\n")
    for line in rconst_decls:
        line = re.sub(r"IF \(",r"if (",line)
        update_rconst.append("        "+line.strip()+";\n")
    update_rconst.append("\n")
    for line in rconst_ops:
        update_rconst.append("        "+line.strip()+";\n")
    for line in rcint:
        update_rconst.append("        "+line.strip()+";\n")
    update_rconst.append("    }\n")
    update_rconst.append("}\n")
    return update_rconst

pass
#########################################################################################################
#########################################################################################################

def generate_kppsolve(source):
    kppsolve=[]
    kppsolve.append("__device__ void kppSolve(const double * __restrict__ Ghimj, double * __restrict__ K, \n\
                         const int istage, const int ros_S )")
    kppsolve.append("{\n")
    kppsolve.append("    int index = blockIdx.x*blockDim.x+threadIdx.x;\n")
    kppsolve.append("\n")
    kppsolve.append("       K = &K[istage*NVAR];\n")
    kppsolve.append("\n")
    for line in source:
        line = line.strip()
        if line != "":
            line = re.sub(r"Ghimj\(index,([0-9]+)\)",r"Ghimj[\1]",line)
            line = re.sub(r"K\(index,istage,([0-9]+)\)",r"K[\1]",line)
            kppsolve.append("        "+line+";\n")
    kppsolve.append("}\n")
    return kppsolve

pass
#########################################################################################################
#########################################################################################################

def generate_kppDecomp(source,NSPEC,lu_diag,lu_crow,lu_icol):
    kppDecomp = []
    kppDecomp.append("__device__ void kppDecomp(double *Ghimj, int VL_GLO)\n")
    kppDecomp.append("{\n")
    kppDecomp.append("    double a=0.0;\n")

    kppDecomp.append("\n")
    kppDecomp.append(" double dummy")
    for var in range(NSPEC):
        kppDecomp.append(", W_" + str(var))
    kppDecomp.append(";\n\n")

    for line in source:
        line = line.strip()
        if line != "":
            line = re.sub(r"Ghimj\(index,([0-9]+)\)",r"Ghimj[\1]",line)
            line = re.sub(r"W\(index,([0-9]+)\)",r"W_\1",line)
            kppDecomp.append("        "+line+";\n")
    kppDecomp.append("}\n")
    return kppDecomp

pass


#########################################################################################################
#########################################################################################################

def generate_kppDecompIndirect(source,NSPEC,lu_diag,lu_crow,lu_icol):
    kppDecomp = []

    kppDecomp.append("\n")

    s_lu =  "__device__ const int LU_DIAG[" + str(len(lu_diag)+1) + "] =  { "
    for i in range(len(lu_diag)):
        s_lu = s_lu + str(lu_diag[i]) + ","
    s_lu = s_lu + "0  };\n"
    kppDecomp.append(s_lu)

    s_lu =  "__device__ const int LU_CROW[" + str(len(lu_crow)+1) + "] =  { "
    for i in range(len(lu_crow)):
        s_lu = s_lu + str(lu_crow[i]) + ","
    s_lu = s_lu + "0  };\n"
    kppDecomp.append(s_lu)

    s_lu =  "__device__ const int LU_ICOL[" + str(len(lu_icol)+1) + "] =  { "
    for i in range(len(lu_icol)):
        s_lu = s_lu + str(lu_icol[i]) + ","
    s_lu = s_lu + "0  };\n"
    kppDecomp.append(s_lu)


    kppDecomp.append("\n")
    kppDecomp.append("__device__ void kppDecomp(double *Ghimj, const int VL_GLO)\n")
    kppDecomp.append("{\n")
    kppDecomp.append("    double a=0.0;\n")
    kppDecomp.append("    int k, kk, j, jj;\n")
    kppDecomp.append("    double W[" + str(NSPEC) +"];\n")

    kppDecomp.append("\n")

    loop = "\n\
    for (k=0;k<NVAR;k++){ \n\
        for ( kk = LU_CROW[k]; kk< (LU_CROW[k+1]-1); kk++){ \n\
            W[ LU_ICOL[kk] ]= Ghimj[kk];\n\
        }\n\
        for ( kk = LU_CROW[k]; kk < (LU_DIAG[k]- 1); k++){\n\
            j = LU_ICOL[kk];\n\
            a = - W[j] / Ghimj[ LU_DIAG[j]];\n\
            W[j] = - a;\n\
            for ( jj = LU_DIAG[j]+1; jj< (LU_CROW[j+ 1]- 1); jj++) {\n\
                W[ LU_ICOL[jj] ] = W[ LU_ICOL[jj]]+  a*Ghimj[jj];\n\
            }\n\
        }\n\
        for (kk = LU_CROW[k]; kk< (LU_CROW[k+ 1]- 1); kk++ ) {\n\
            Ghimj[kk] = W[ LU_ICOL[kk]];\n\
        }\n\
    }\n"


    kppDecomp.append(loop)
    kppDecomp.append("}\n")
    return kppDecomp

pass

#########################################################################################################
#########################################################################################################

def generate_jac_sp(source,NBSIZE):
    jac_sp = []
    jac_sp.append("__device__ void Jac_sp(const double * __restrict__ var, const double * __restrict__ fix,\n\
                 const double * __restrict__ rconst, double * __restrict__ jcb, int &Njac, const int VL_GLO)\n")
    jac_sp.append("{\n")
    jac_sp.append("    int index = blockIdx.x*blockDim.x+threadIdx.x;\n")

    jac_sp.append("\n")
    jac_sp.append(" double dummy")
    for var in range(NBSIZE):
        jac_sp.append(", B_" + str(var))
    jac_sp.append(";\n\n")

    jac_sp.append("\n")
    jac_sp.append("    Njac++;\n")
    jac_sp.append("\n")



    for line in source:
        line = line.strip()
        if line != "":
            line = re.sub(r"B\(index,([0-9]+)\)",r"B_\1",line)
            line = re.sub(r"jcb\(index,([0-9]+)\)",r"jcb[\1]",line)
            line = re.sub(r"var\(index,([0-9]+)\)",r"var[\1]",line)
            line = re.sub(r"fix\(index,([0-9]+)\)",r"fix[\1]",line)
            jac_sp.append("        "+line+";\n")
    jac_sp.append("    }\n")
    return jac_sp

pass
#########################################################################################################
#########################################################################################################
def generate_fun(source,NREACT):
    fun = []
    fun.append("__device__ void Fun(double *var, const double * __restrict__ fix, const double * __restrict__ rconst, double *varDot, int &Nfun, const int VL_GLO)")
    fun.append("{\n")
    fun.append("    int index = blockIdx.x*blockDim.x+threadIdx.x;\n")
    fun.append("\n")
    fun.append("    Nfun++;\n")
    fun.append("\n")

    fun.append(" double dummy")
    for var in range(NREACT):
        fun.append(", A_" + str(var))
    fun.append(";\n\n")


    fun.append("    {\n")
    for line in source:
        line = line.strip()
        if line != "":
            line = re.sub(r"A\(index,([0-9]+)\)",r"A_\1",line)
            line = re.sub(r"var\(index,([0-9]+)\)",r"var[\1]",line)
            line = re.sub(r"varDot\(index,([0-9]+)\)",r"varDot[\1]",line)
            fun.append("        "+line+";\n")
    fun.append("    }\n")
    fun.append("}\n")
    return fun

pass
#########################################################################################################
#########################################################################################################



def find_LU_DIAG(file_in, NVAR):
    file_in.seek(0)
    source = file_in.readlines()
    the_line = 0
    glu_diag = []
    long_tables = False

    for line_num in range(len(source)):
        if "lu_diag_0" in source[line_num].lower():
            print("Detected long tables!")
            long_tables = True
            the_line = line_num
            break


    for line_num in range(len(source)):
        if (long_tables == True):
            if "lu_diag " in source[line_num].lower():
                end_line = line_num
        else:
            if "lu_diag" in source[line_num].lower():
                the_line = line_num
                break

    lu_diag = []
    if (long_tables == True):
        for line_num in range(the_line,end_line):
            lu_diag.append(source[line_num])
    else:
        for line_num in range(the_line,len(source)):
            lu_diag.append(source[line_num])
            if "/)" in source[line_num]:
                break;


    lu_diag = remove_comments(lu_diag)
    lu_diag = strip_and_unroll_lines(lu_diag)
    lu_diag = "".join(lu_diag)
    lu_diag = lu_diag.lower()
    lu_diag = lu_diag[lu_diag.find("(/")+2:lu_diag.rfind("/)")]

    lu_diag = re.sub(r"\/\)\ninteger, parameter, dimension\([0-9]+\)?\s::?\slu_diag_[0-9]?\s=?\s\(/",r",",lu_diag)


    # if failes break it in smaller
    lu_diag = re.sub(r"dimension\([0-9]+\)::lu_diag_[0-9]\s?=\s?\(\/",r",",lu_diag)
    lu_diag = re.sub(r"dimension\([0-9]+\)", r"",lu_diag)
    lu_diag = re.sub(r"::", r"",lu_diag)
    lu_diag = re.sub(r"lu_diag_[0-9]+\s?=\s?",r"",lu_diag)
    lu_diag = re.sub(r"\(/",r"",lu_diag)
    lu_diag = lu_diag.replace("/)\ninteger","")
    lu_diag = lu_diag.replace("parameter,","")


    lu_diag = lu_diag.replace(" ","")
    lu_diag = lu_diag.split(",")
    for line_num in range(len(lu_diag)):
        lu_diag[line_num] = str(int(lu_diag[line_num])-1)


    return lu_diag




def find_LU_CROW(file_in, NVAR):
    file_in.seek(0)
    source = file_in.readlines()
    the_line = 0
    glu_diag = []
    long_tables = False

    for line_num in range(len(source)):
        if "lu_crow_0" in source[line_num].lower():
            print("Detected long tables!")
            long_tables = True
            the_line = line_num
            break


    for line_num in range(len(source)):
        if (long_tables == True):
            if "lu_crow " in source[line_num].lower():
                end_line = line_num
        else:
            if "lu_crow" in source[line_num].lower():
                the_line = line_num
                break

    lu_diag = []
    if (long_tables == True):
        for line_num in range(the_line,end_line):
            lu_diag.append(source[line_num])
    else:
        for line_num in range(the_line,len(source)):
            lu_diag.append(source[line_num])
            if "/)" in source[line_num]:
                break;


    lu_diag = remove_comments(lu_diag)
    lu_diag = strip_and_unroll_lines(lu_diag)
    lu_diag = "".join(lu_diag)
    lu_diag = lu_diag.lower()
    lu_diag = lu_diag[lu_diag.find("(/")+2:lu_diag.rfind("/)")]

    lu_diag = re.sub(r"\/\)\ninteger, parameter, dimension\([0-9]+\)?\s::?\slu_crow_[0-9]?\s=?\s\(/",r",",lu_diag)

    # if failes break it in smaller
    lu_diag = re.sub(r"dimension\([0-9]+\)::lu_crow_[0-9]\s?=\s?\(\/",r",",lu_diag)
    lu_diag = re.sub(r"dimension\([0-9]+\)", r"",lu_diag)
    lu_diag = re.sub(r"::", r"",lu_diag)
    lu_diag = re.sub(r"lu_crow_[0-9]\s?=\s?",r"",lu_diag)
    lu_diag = re.sub(r"\(/",r"",lu_diag)
    lu_diag = lu_diag.replace("/)\ninteger","")
    lu_diag = lu_diag.replace("parameter,","")



    lu_diag = lu_diag.replace(" ","")
    lu_diag = lu_diag.split(",")
    for line_num in range(len(lu_diag)):
        lu_diag[line_num] = str(int(lu_diag[line_num])-1)


    return lu_diag

def find_LU_ICOL(file_in, NVAR):
    file_in.seek(0)
    source = file_in.readlines()
    the_line = 0
    glu_diag = []
    long_tables = False

    for line_num in range(len(source)):
        if "lu_icol_0" in source[line_num].lower():
            print("Detected long tables!")
            long_tables = True
            the_line = line_num
            break


    for line_num in range(len(source)):
        if (long_tables == True):
            if "lu_icol " in source[line_num].lower():
                end_line = line_num
        else:
            if "lu_icol" in source[line_num].lower():
                the_line = line_num
                break

    lu_diag = []
    if (long_tables == True):
        for line_num in range(the_line,end_line):
            lu_diag.append(source[line_num])
    else:
        for line_num in range(the_line,len(source)):
            lu_diag.append(source[line_num])
            if "/)" in source[line_num]:
                break;


    lu_diag = remove_comments(lu_diag)
    lu_diag = strip_and_unroll_lines(lu_diag)
    lu_diag = "".join(lu_diag)
    lu_diag = lu_diag.lower()
    lu_diag = lu_diag[lu_diag.find("(/")+2:lu_diag.rfind("/)")]

    lu_diag = re.sub(r"\/\)\ninteger, parameter, dimension\([0-9]+\)?\s::?\slu_icol_[0-9]?\s=?\s\(/",r",",lu_diag)


    # if failes break it in smaller
    lu_diag = re.sub(r"dimension\([0-9]+\)::lu_icol_[0-9]\s?=\s?\(\/",r",",lu_diag)
    lu_diag = re.sub(r"dimension\([0-9]+\)", r"",lu_diag)
    lu_diag = re.sub(r"::", r"",lu_diag)
    lu_diag = re.sub(r"lu_icol_[0-9]+\s?=\s?",r"",lu_diag)
    lu_diag = re.sub(r"\(/",r"",lu_diag)
    lu_diag = lu_diag.replace("/)\ninteger","")
    lu_diag = lu_diag.replace("parameter,","")

    lu_diag = lu_diag.replace(" ","")
    lu_diag = lu_diag.split(",")
    for line_num in range(len(lu_diag)):
        lu_diag[line_num] = str(int(lu_diag[line_num])-1)


    return lu_diag


#########################################################################################################

def generate_prepareMatrix(lu_diag):
    prepareMatrix = []
    prepareMatrix.append("__device__ void ros_PrepareMatrix(double &H, int direction, double gam, double *jac0, double *Ghimj,  int &Nsng, int &Ndec, int VL_GLO)\n")
    prepareMatrix.append("{\n")
    prepareMatrix.append("    int index = blockIdx.x*blockDim.x+threadIdx.x;\n")
    prepareMatrix.append("    int ising, nConsecutive;\n")
    prepareMatrix.append("    double ghinv;\n")
    prepareMatrix.append("    \n")
    prepareMatrix.append("        ghinv = ONE/(direction*H*gam);\n")
    prepareMatrix.append("        for (int i=0; i<LU_NONZERO; i++)\n")
    prepareMatrix.append("            Ghimj[i] = -jac0[i];\n\n")
    for i in lu_diag:
        prepareMatrix.append("        Ghimj["+i+"] += ghinv;\n")
    prepareMatrix.append("        ros_Decomp(Ghimj, Ndec, VL_GLO);\n")
    prepareMatrix.append("}\n")
    return prepareMatrix


pass

#########################################################################################################

def generate_special_ros(ros,inject_rconst):



    if ( ros == '2'):
        file_ros = open("./source/ros2.cu","r")
    elif (ros == '3'):
        file_ros = open("./source/ros3.cu","r")
    elif (ros == '4'):
        file_ros = open("./source/ros4.cu","r")
    elif (ros == '5'):
        file_ros = open("./source/rodas3.cu","r")
    elif (ros == '6'):
        file_ros = open("./source/rodas4.cu","r")
    else:
        return ''


    rosfunc = []
    source = file_ros.readlines()
    for line in source:
        if ( inject_rconst is True ):
            line = line.replace("Jac_sp(var, fix, rconst, jac0, Njac, VL_GLO)","update_rconst(var, khet_st, khet_tr, jx, VL_GLO);   \n        Jac_sp(var, fix, rconst, jac0, Njac, VL_GLO)")
            line = line.replace("Fun(varNew, fix, rconst, varNew, Nfun,VL_GLO);","update_rconst(var, khet_st, khet_tr, jx, VL_GLO);   \n            Fun(varNew, fix, rconst, varNew, Nfun,VL_GLO);")
            line = line.replace("Fun(var, fix, rconst, Fcn0, Nfun, VL_GLO);","update_rconst(var, khet_st, khet_tr, jx, VL_GLO);   \n           Fun(var, fix, rconst, Fcn0, Nfun, VL_GLO);")
        rosfunc.append(line)




    return rosfunc


pass
#########################################################################################################


def generate_special_ros_caller(ros):

    roscall = []

    default_call = '      Rosenbrock<<<dimGrid,dimBlock>>>(d_conc, Tstart, Tend, d_rstatus, d_istatus,\n\
                    // values calculated from icntrl and rcntrl at host\n\
                    autonomous, vectorTol, UplimTol, method, Max_no_steps,\n\
                    d_jac0, d_Ghimj,d_varNew, d_K, d_varErr, d_dFdT, d_Fcn0, d_var, d_fix, d_rconst,\n\
                    Hmin, Hmax, Hstart, FacMin, FacMax, FacRej, FacSafe, roundoff,\n\
                    //  cuda global mem buffers              \n\
                    d_absTol, d_relTol,   \n\
                    d_khet_st, d_khet_tr, d_jx, \n\
                    // Global input arrays\n\
                    temp_gpu, press_gpu, cair_gpu, \n\
                    // extra - vector lenght and processor\n\
                    VL_GLO); '

    if ( ros == '2'):
        rosscall = '     switch (method){\n\
        case 1:\n\
            Rosenbrock_ros2<<<dimGrid,dimBlock>>>(d_conc, Tstart, Tend, d_rstatus, d_istatus,\n\
                    autonomous, vectorTol, UplimTol, Max_no_steps,\n\
                    d_jac0, d_Ghimj,d_varNew, d_K, d_varErr, d_dFdT, d_Fcn0, d_var, d_fix, d_rconst,\n\
                    Hmin, Hmax, Hstart, FacMin, FacMax, FacRej, FacSafe, roundoff,\n\
                    d_absTol, d_relTol,\n\
                    d_khet_st, d_khet_tr, d_jx, \n\
                    temp_gpu, press_gpu, cair_gpu, \n\
                    VL_GLO);\n\
            break;\n\
        default: \n' + default_call + '\n\
        \n\
                    break;\n\
    }\n'

    elif (ros == '3'):
        rosscall = '      switch (method){\n\
        case 2:\n\
            Rosenbrock_ros3<<<dimGrid,dimBlock>>>(d_conc, Tstart, Tend, d_rstatus, d_istatus,\n\
                    autonomous, vectorTol, UplimTol, Max_no_steps,\n\
                    d_jac0, d_Ghimj,d_varNew, d_K, d_varErr, d_dFdT, d_Fcn0, d_var, d_fix, d_rconst,\n\
                    Hmin, Hmax, Hstart, FacMin, FacMax, FacRej, FacSafe, roundoff,\n\
                    d_absTol, d_relTol,\n\
                    d_khet_st, d_khet_tr, d_jx, \n\
                    temp_gpu, press_gpu, cair_gpu, \n\
                    VL_GLO);\n\
            break;\n\
        default: \n' + default_call + '\n\
        \n\
                    break;\n\
    }\n'


    elif (ros == '4'):
        rosscall = '      switch (method){\n\
        case 3:\n\
            Rosenbrock_ros4<<<dimGrid,dimBlock>>>(d_conc, Tstart, Tend, d_rstatus, d_istatus,\n\
                    autonomous, vectorTol, UplimTol, Max_no_steps,\n\
                    d_jac0, d_Ghimj,d_varNew, d_K, d_varErr, d_dFdT, d_Fcn0, d_var, d_fix, d_rconst,\n\
                    Hmin, Hmax, Hstart, FacMin, FacMax, FacRej, FacSafe, roundoff,\n\
                    d_absTol, d_relTol,\n\
                    d_khet_st, d_khet_tr, d_jx, \n\
                    temp_gpu, press_gpu, cair_gpu, \n\
                    VL_GLO);\n\
            break;\n\
        default: \n' + default_call + '\n\
        \n\
                    break;\n\
    }\n'


    elif (ros == '5'):
        rosscall = '      switch (method){\n\
        case 4:\n\
            Rosenbrock_rodas3<<<dimGrid,dimBlock>>>(d_conc, Tstart, Tend, d_rstatus, d_istatus,\n\
                    autonomous, vectorTol, UplimTol, Max_no_steps,\n\
                    d_jac0, d_Ghimj,d_varNew, d_K, d_varErr, d_dFdT, d_Fcn0, d_var, d_fix, d_rconst,\n\
                    Hmin, Hmax, Hstart, FacMin, FacMax, FacRej, FacSafe, roundoff,\n\
                    d_absTol, d_relTol,\n\
                    d_khet_st, d_khet_tr, d_jx, \n\
                    temp_gpu, press_gpu, cair_gpu, \n\
                    VL_GLO);\n\
            break;\n\
        default: \n' + default_call + '\n\
        \n\
                    break;\n\
    }\n'


    elif (ros == '6'):
        rosscall = '      switch (method){\n\
        case 5:\n\
            Rosenbrock_rodas4<<<dimGrid,dimBlock>>>(d_conc, Tstart, Tend, d_rstatus, d_istatus,\n\
                    autonomous, vectorTol, UplimTol, Max_no_steps,\n\
                    d_jac0, d_Ghimj,d_varNew, d_K, d_varErr, d_dFdT, d_Fcn0, d_var, d_fix, d_rconst,\n\
                    Hmin, Hmax, Hstart, FacMin, FacMax, FacRej, FacSafe, roundoff,\n\
                    d_absTol, d_relTol,\n\
                    d_khet_st, d_khet_tr, d_jx, \n\
                    temp_gpu, press_gpu, cair_gpu, \n\
                    VL_GLO);\n\
            break;\n\
        default: \n' + default_call + '\n\
        \n\
                    break;\n\
    }\n'
    else:
        return default_call



    return rosscall


pass



#########################################################################################################
#########################################################################################################

def generate_define_indices_one_line(file_in,prefix):
    file_in.seek(0)
    source = file_in.readlines()
    source = remove_comments(source)
    source = strip_and_unroll_lines(source)
    the_line = 0
    for line in source:
        if prefix+"_" in line:
            the_line = line
            break
    the_line = the_line[the_line.find("::")+2:].replace("=","").replace("  "," ").replace("  "," ").replace("  "," ").replace("  "," ").replace("  "," ").replace("\n"," ")
    the_line = the_line.split(",")
    for i in range(len(the_line)):
        if (len(the_line[i])<3):
            continue;
        #    return;
        the_line[i] = the_line[i].strip()
        name, value = the_line[i].split(" ")
        the_line[i] = "#define "+name+" "+str(int(value)-1)+"\n"
    return the_line



def generate_define_indices_many_lines(file_in,prefix):
    file_in.seek(0)
    source = file_in.readlines()
    source = remove_comments(source)
    source = strip_and_unroll_lines(source)
    the_line = 0
    for i in range(len(source)):
        if " "+prefix+"_" in source[i]:
            the_line = i
            break
    ind_s = []
    for i in range(the_line,len(source)):
        source[i] = source[i].strip()
        if prefix+"_" not in source[i] and source[i] != "":
            break
        if prefix+"_" in source[i]:
            ind_s.append(source[i])
    for i in range(len(ind_s)):
        ind_s[i] = ind_s[i][ind_s[i].find("::")+2:].strip().replace("=","").replace("  "," ").replace("  "," ").replace("  "," ").replace("  "," ").replace("  "," ")
        name,value = ind_s[i].split(" ")
        ind_s[i] = "#define "+name+" "+str(int(value)-1)+"\n"
    return ind_s

def generate_define_vars(file_in,var_names):
    file_in.seek(0)
    source = file_in.readlines()
    source = remove_comments(source)
    source = strip_and_unroll_lines(source)
    out = []
    for line in source:
        if var_names == []:
            break
        for var_name in var_names:
            if var_name in line and "=" in line:
                var = line[line.find("::")+2:].strip()
                name,value = var.split("=")
                name = name.strip()
                value = value.strip().replace("_dp","").replace("_DP","").replace("_Dp","")
                if "real" in line.lower():
                    value = float(value)
                elif "integer" in line.lower():
                    value = int(value)
                else:
                    value = int(value)
                out.append("#define "+name+" "+str(value)+"\n")
                var_names.remove(var_name)
    if var_names != []:
        print("Warning: variables "+str(var_names)+" were not found")
    return out

#
# Takes prefix of variables as input and the file
# Returns definitions using index
#
def generate_definitions_global(file_in,var_prefix):
    file_in.seek(0)
    source = file_in.readlines()
    source = remove_comments(source)
    source = strip_and_unroll_lines(source)
    out = []

    for var_name in var_prefix:
        for line in source:

            # ignore some definitions that are not double
            if "INTEGER" in line:
                continue

            # we reached after the definitions
            if "interface" in line:
                break

            allvars = re.findall(r'(' + var_name  + '(\w+)(\s+)?)=\s+(([0-9,E,\-,.])+(\s+)?)[,&,\n]',line)

            if ( len(allvars)  > 0):
                for definition in allvars:
                    out.append("#define "+definition[0]+"  ("+str(definition[3])+")\n")

    return out

pass
#########################################################################################################
#########################################################################################################

def gen_kpp_integrate_cuda(file_prototype, source_cuda, inject_rconst):
    file_prototype.seek(0)
    lines_prototype = file_prototype.readlines()
    file_out = open(smcl + "messy_mecca_kpp_acc.cu","w")
    for line in lines_prototype:
        if "=#=#=#=#=#=#=#=#=#=#=" in line:
            chunk_name = line.replace("=#=#=#=#=#=#=#=#=#=#=","").replace("=#=#=#=#=#=#=#=#=#=#=","").strip().lower()
            chunk = source_cuda[chunk_name]
            if chunk is not None:
                for chunk_line in chunk:
                    chunk_line = remove_precision_qualifiers(chunk_line)
                    file_out.write(chunk_line)
        else:

            if ( inject_rconst is True ):
                line = line.replace("Jac_sp(var, fix, rconst, jac0, Njac, VL_GLO)","update_rconst(var, khet_st, khet_tr, jx, VL_GLO);   \n        Jac_sp(var, fix, rconst, jac0, Njac, VL_GLO)")
                line = line.replace("Fun(varNew, fix, rconst, varNew, Nfun,VL_GLO);","update_rconst(var, khet_st, khet_tr, jx, VL_GLO);   \n            Fun(varNew, fix, rconst, varNew, Nfun,VL_GLO);")
                line = line.replace("Fun(var, fix, rconst, Fcn0, Nfun, VL_GLO);","update_rconst(var, khet_st, khet_tr, jx, VL_GLO);   \n           Fun(var, fix, rconst, Fcn0, Nfun, VL_GLO);")
                line = line.replace("Fun(var, fix, rconst, dFdT, Nfun, VL_GLO);","update_rconst(var, khet_st, khet_tr, jx, VL_GLO);   \n       Fun(var, fix, rconst, dFdT, Nfun, VL_GLO);")



            file_out.write(line)
    file_out.close()




pass
#########################################################################################################
#########################################################################################################

def generate_define_NBSIZE(source):
    lines = source[:]
    lines = remove_comments(lines)
    lines = strip_and_unroll_lines(lines)
    nbsize = ""
    for line in lines:
        if ":: b(" in line.lower():
            index = line.lower().find(":: b(")+5
            nbsize = "#define NBSIZE "+line[index:line.find(")",index)].strip()
            break
    return nbsize


pass
#########################################################################################################
#########################################################################################################

def remove_precision_qualifiers(line):
    operators = "*/+-<>=.,"
    var_name = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_'
    precision_qualifiers = ["_sp","_dp","_qp"]
    truth_mat = []
    for qual in precision_qualifiers:
        truth_mat.append((qual in line.lower())*1)
    if any(truth_mat):
        for truth_id in range(len(truth_mat)):
            qual = 0
            if truth_mat[truth_id]:
                qual = precision_qualifiers[truth_id]
                index = len(line)
                while True:
                    index = line.rfind(qual,0, index-1)
                    if index == -1:
                        break
                    left = line[:index]
                    right = line[index:]
                    #if left[-3:].lower() in precision_qualifiers: # if it's 33._dp, 33.33E-33_dp or a_333_dp
                    number = 0
                    if left[-1] == ".": # if it's 33._dp
                        number = 1
                    elif left[-1] in var_name:  # if it's 33.33E-33_dp or a_333_dp
                        pos=0
                        scientif = False
                        isnumber = True
                        for i in reversed(left[:]):
                            if i not in string.digits+".":
                                if (i in "+-") and abs(pos-2)<=len(left) and (left[pos-2].lower() == "e"):
                                    scientif = True
                                elif scientif and i.lower() == "e":
                                    pass
                                elif i in string.ascii_letters+"_": # if it's a_333_dp
                                    pos=0
                                    for i in reversed(left[:]):
                                        if i not in var_name:
                                            break
                                        pos = pos - 1
                                    isnumber = False
                                    break
                                else:
                                    break
                            pos = pos - 1
                        number = not ((not isnumber) and (not scientif))
                    line = left + (not number)*right[:3] + right[3:]
    return line

pass
#########################################################################################################
#########################################################################################################

def generate_c2f_interface(file_in):
    file_in.seek(0)
    source = file_in.readlines()
    start = -1
    stop = -1
    for i in range(len(source)):
        if 'subroutinekpp_integrate' in source[i].lower().replace(" ",""):
            if start == -1:
                start = i
            elif "end" in source[i].lower().replace(" ",""):
                stop = i
                break
            else:
                print("Something went wrong in generate c2f_interface")
                return
    file_out = open(smcl + "messy_mecca_kpp.f90","w")
    for i in range(start):
        file_out.write(source[i])
    out = "SUBROUTINE kpp_integrate (time_step_len,Conc,ierrf,xNacc,xNrej,istatus,l_debug,PE) \n\
                                                                    \n\
  IMPLICIT NONE                                                     \n\
                                                                    \n\
  REAL(dp),INTENT(IN)                   :: time_step_len           \n\
  REAL(dp),INTENT(INOUT),dimension(:,:) :: Conc                    \n\
  INTEGER, INTENT(OUT),OPTIONAL        :: ierrf(:)                \n\
  INTEGER, INTENT(OUT),OPTIONAL        :: xNacc(:)                \n\
  INTEGER, INTENT(OUT),OPTIONAL        :: xNrej(:)                \n\
  INTEGER, INTENT(INOUT),OPTIONAL      :: istatus(:)              \n\
  INTEGER, INTENT(IN),OPTIONAL         :: PE                      \n\
  LOGICAL, INTENT(IN),OPTIONAL         :: l_debug                 \n\
  \n\
  integer, save                        :: counter = 0  ! For debuging\n\
  integer, save                        :: newcounter = 0  ! For debuging\n\
\n\
                                                                    \n\
  INTEGER                                :: k   ! loop variable     \n\
  REAL(dp)                               :: dt                      \n\
  REAL(dp)                               :: roundoff\n\
  integer,dimension(:),allocatable      :: xNacc_u, xNrej_u, ierr_u2\n\
  integer                                :: ierr_u                  \n\
  integer                                :: i\n\
  integer                                :: j\n\
  integer                                :: istep\n\
  integer,dimension(20)                  :: istatus_u               \n\
  integer,dimension(5)                   :: sizes\n\
                                                                    \n\
  LOGICAL :: file_exists\n\
\n\
  character(len=10) :: filename\n\
  CHARACTER(LEN=3000) :: rowfmt\n\
                                                                    \n\
  if (present (istatus)) istatus = 0                              \n\
\n\
\n\
  allocate(xNacc_u(VL_GLO))\n\
  allocate(xNrej_u(VL_GLO))\n\
  allocate(ierr_u2(VL_GLO))\n\
\n\
\n\
  sizes(1) = VL_glo\n\
  sizes(2) = size(khet_st,2)\n\
  sizes(3) = size(khet_tr,2)\n\
  sizes(4) = size(jx,2)\n\
  roundoff = WLAMCH('E')\n\
\n\
\n\
#if 1\n\
            \n\
  CALL kpp_integrate_cuda(PE, sizes, time_step_len, Conc, temp, press, cair, &\n\
        khet_st, khet_tr, jx, aTol, rTol, ierr_u2, istatus_u, xNacc_u, xNrej_u, roundoff, icntrl, rcntrl)\n\
\n\
  DO k=1,VL_glo,VL_DIM                                              \n\
    is = k                                                          \n\
    ie = min(k+VL_DIM-1,VL_glo)                                     \n\
    vl = ie-is+1                                                    \n\
                                                                    \n\
    ! Return Diagnostic Information                                 \n\
                                                                    \n\
    if(Present(ierrf))    ierrf(is) = IERR_U                      \n\
    if(Present(xNacc))    xNacc(is) = istatus_u(4)                \n\
    if(Present(xNrej))    xNrej(is) = istatus_u(5)                \n\
                                                                    \n\
    if (present (istatus)) then                                   \n\
      istatus(1:8) = istatus(1:8) + istatus_u(1:8)                 \n\
    end if                                                          \n\
                                                                    \n\
  END DO      \n\
\n\
#endif\n\
\n\
#if 0\n\
\n\
  DO k=1,VL_glo,VL_DIM                                              \n\
    is = k                                                          \n\
    ie = min(k+VL_DIM-1,VL_glo)                                     \n\
    vl = ie-is+1                                                    \n\
                                                                    \n\
    C(:) = Conc(is,:)                                             \n\
                                                                    \n\
    CALL update_rconst                                              \n\
                                                                    \n\
    dt = time_step_len                                              \n\
                                                                    \n\
    ! integrate from t=0 to t=dt                                    \n\
    CALL integrate(0._dp,dt,icntrl,rcntrl,istatus_u = istatus_u,ierr_u=ierr_u)\n\
                       \n\
                                                                    \n\
    IF (PRESENT(l_debug) .AND. PRESENT(PE)) THEN                       \n\
      IF (l_debug) CALL error_output(Conc(is,:),ierr_u,PE)           \n\
    ENDIF           \n\
    Conc(is,:) = C(:)                                                 \n\
\n\
    if(Present(ierrf))    ierrf(is) = IERR_U                      \n\
    if(Present(xNacc))    xNacc(is) = istatus_u(4)                \n\
    if(Present(xNrej))    xNrej(is) = istatus_u(5)                \n\
                                                                    \n\
    if (present (istatus)) then                                   \n\
      istatus(1:8) = istatus(1:8) + istatus_u(1:8)                 \n\
    end if             \n\
\n\
  END DO  \n\
#endif\n\
\n\
! Deallocate input arrays                                           \n\
                                                                    \n\
  if (allocated(TEMP))   deallocate(TEMP)   \n\
  if (allocated(cair))   deallocate(cair)   \n\
  if (allocated(press))   deallocate(press)   \n\
  if (allocated(temp_ion))   deallocate(temp_ion)   \n\
  if (allocated(temp_elec))   deallocate(temp_elec)   \n\
  if (allocated(xaer))   deallocate(xaer)   \n\
  if (allocated(cvfac))   deallocate(cvfac)   \n\
  if (allocated(lwc))   deallocate(lwc)   \n\
  if (allocated(k_exf))   deallocate(k_exf)   \n\
  if (allocated(k_exb))   deallocate(k_exb)   \n\
  if (allocated(k_exf_N2O5))   deallocate(k_exf_N2O5)   \n\
  if (allocated(k_exf_ClNO3))   deallocate(k_exf_ClNO3)   \n\
  if (allocated(k_exf_BrNO3))   deallocate(k_exf_BrNO3)   \n\
  if (allocated(jx))   deallocate(jx)   \n\
  if (allocated(khet_Tr))   deallocate(khet_Tr)   \n\
  if (allocated(khet_St))   deallocate(khet_St)   \n\
  if (allocated(mcexp))   deallocate(mcexp)   \n\
\n\
  deallocate(xNacc_u)\n\
  deallocate(xNrej_u)\n\
  deallocate(ierr_u2)      \n\
  data_loaded = .false.                                             \n\
                                                                    \n\
  return                                                            \n\
END SUBROUTINE kpp_integrate\n"
    file_out.write(out)


    for i in range(stop+1,len(source)):
        file_out.write(source[i])


pass

#########################################################################################################
#########################################################################################################

def add_cuda_compilation(file_specific,file_makefile,arch):

    file_makefile.seek(0)
    out = "\nmessy_mecca_kpp_acc.o: messy_mecca_kpp_acc.cu specific.mk \n\
\tnvcc  -v  --ptxas-options=-v  " + arch +"  --ftz=false --prec-div=true --prec-sqrt=true --fmad=false    -O3  -g   -c  $<"
    file_specific.write(out)


    temp = open('__temp', 'wb')
    for line in file_makefile:
        if line.startswith('$(LIB): depend $(OBJS)'):
            line = line.strip() + ' $(OBJSCUDA)\n'
        if line.startswith('\t$(AR) $(ARFLAGS) $(LIB) $(OBJS)'):
            line = line.rstrip() + ' $(OBJSCUDA)\n'
        if line.startswith('depend $(MAKEFILE_INC): $(SRCS)'):
            line = line.rstrip() + ' $(ACC_SRCS)\n'
        if line.startswith('OBJS  := $(SRCS:.f90=.o)'):
            line ='OBJSCUDA  := $(SRCS_ACC:.cu=.o)\n' +  line
        if line.startswith('SRCS  := $(filter-out F%.f90, $(SRCS0)'):
            line ='SRCS_ACC  := $(wildcard *.cu) \n' +  line

        if line.startswith('.SUFFIXES: $(SUFFIXES) .f90 .md5'):
            line ='.SUFFIXES: $(SUFFIXES) .f90 .md5 .cu\n'

        temp.write(line.encode())

    temp.close()
    os.rename('__temp', smcl + "Makefile.m")



pass


#########################################################################################################
#########################################################################################################

# Based on the input files, select the proper flags
def get_transformation_flags():

    multifile = False
    vectorize = False
    indirect  = False
    inject_rconst = False

    # Check if kpp created indirect indexing
    if ('LU_CROW(k+1)' in open(smcl + "messy_mecca_kpp.f90").read()) or ('LU_CROW(k+ 1)' in open(smcl + "messy_mecca_kpp.f90").read()):
        print("Warning: Can't convert indirect indexing of file.")
        print("--> Change the decomp in the conf file or modify the output file.\n")
        indirect = True


    # Check if kpp created vector length chemistry
    if '= C(1:VL,:)' in open(smcl + "messy_mecca_kpp.f90").read():
        print("Can't convert vectorized version of file.")
        print("--> Change the rosenbrock_vec to reosenbrock_mz in the conf file.\n")
        print("Exiting... \n")
        vectorized = True
        exit(-1)


    # Check if kpp created the multiple files version.
    if ( os.path.isfile(smcl + "messy_mecca_kpp_global.f90") == True             and
         os.path.isfile(smcl + "messy_mecca_kpp_jacobian.f90") == True
        ):
        print("Multifile version detected!")
        multifile = True

    if (multifile == True):
        file_messy_mecca_kpp = open(smcl + "messy_mecca_kpp_linearalgebra.f90","r")
        subroutines = find_subroutines(file_messy_mecca_kpp, ["KppDecomp","KppDecompCmplx"])
        infile = " ".join(subroutines["kppdecomp"])
        if 'LU_ICOL(kk)' in infile:
            print("Multiple files with indirect indexing detected.\n")
            indirect = True

    if (multifile == True):
        file_messy_mecca_kpp = open(smcl + "messy_mecca_kpp_integrator.f90","r")
        lines = file_messy_mecca_kpp.readlines()
        infile = " ".join(lines)
        if '!     CALL Update_RCONST()' not in infile:
            inject_rconst = True;
    else:
        file_messy_mecca_kpp = open(smcl + "messy_mecca_kpp.f90","r")
        lines = file_messy_mecca_kpp.readlines()
        infile = " ".join(lines)
        if '!     CALL Update_RCONST()' not in infile:
            inject_rconst = True;

    return multifile, vectorize, indirect, inject_rconst

pass

#########################################################################################################
#########################################################################################################

def print_warning():
    print('\033[1m' + "\n####################################################################################################")
    print("## WARNING!! BETA VERSION ! PLEASE REPORT TO PACKAGE MAINTAINERS ANY BUGS OR UNEXPECTED BEHAVIOUR.")
    print("####################################################################################################\n")
    print('\033[0m')
pass

#########################################################################################################
#########################################################################################################


def select_architecture(ans):
    if ans=="1":
        arch = "--gpu-architecture=compute_20 -maxrregcount=128 "
    elif ans=="2":
        arch = "--gpu-architecture=compute_35"
    elif ans=="3":
        arch = "--gpu-architecture=compute_50"
    elif ans=="4":
        arch = "--gpu-architecture=compute_60"
    elif ans=="5":
        arch = "--gpu-architecture=compute_70"
    else:
        arch = "--gpu-architecture=compute_20"

    return arch

def print_menu_make_selection(ros,gpu):

    if gpu is None:
        print ("""

    Select CUDA architecture (1-5): 

                1. CUDA 2.0 ( FERMI GPU architecture   )
                2. CUDA 3.5 ( KEPLER GPU architecture  )
                3. CUDA 5.2 ( MAXWELL GPU architecture )
                4. CUDA 6.0 ( PASCAL GPU architecture  )
                5. CUDA 7.0 ( VOLTA GPU architecture   )

                """)

        gpu = input("Option (Default 1): ")

    arch = select_architecture(gpu)


    if ros is None:
        print ("""

Select Rosenbrock solver (1-6): 

            1. All    ( Selects based on the runtime option  )
            2. Ros2   ( 2-stage L-stable - FASTEST           )
            3. Ros3   ( 3-stage L-stable - RECOMMENDED       )
            4. Ros4   ( 4-stage L-stable                     )
            5. Rodas3 ( 4-stage stiffly accurate             )
            6. Rodas4 ( 6-stage stiffly accurate - SLOWEST   )

            """)

        ros = input("Option (Default 1): ")

    if ros not in ['1','2','3','4','5','6']:
        ros = "0"

    print("Selected options: " + arch + " with ros: "  + ros + "\n")
    return ros,arch

#########################################################################################################
#########################################################################################################
##
##  Main program
##
#########################################################################################################
#########################################################################################################

# Set variables for checking
multifile = False
vectorize = False
indirect  = False
inject_rconst = False


###############################################

# check if we have the arguments
parser = argparse.ArgumentParser(description='MEDINA: FORTRAN to CUDA KPP for EMAC Preprocessor.')
parser.add_argument('-r', '--ros', help='An integer value of the Rosenbrock solver produced [1: all (selected at runtime), 2: Ros2, 3: Ros3, 4: Rodas3, 5: Rodas4]')
parser.add_argument('-g', '--gpu', help='An integer value of the architecture [1: FERMI, 2: KEPLER, 3: MAXWELL, 4: PASCAL]')
parser.add_argument('-s', '--smcl', help='smcl folder location, default: "../../smcl/"')
args = parser.parse_args()

if args.smcl:
  smcl = args.smcl
ros = args.ros
gpu = args.gpu


# get the options for the architecture and the rosenbrock kernel
ros,arch = print_menu_make_selection(ros,gpu)


###############################################
# Print generic information - header
print("\n+===================================================================+ ")
print("| KPP Fortran to CUDA praser - Copyright 2016 The Cyprus Institute  |")
print("+===================================================================+ \n")

print_warning()

# First check if the files exist
# In the future, we will check also  the path of the binary

if ( os.path.isfile(smcl + "messy_mecca_kpp.f90") == False             or
     os.path.isfile(smcl + "messy_cmn_photol_mem.f90") == False        or
     os.path.isfile(smcl + "messy_main_constants_mem.f90") == False    or
     os.path.isfile("./source/kpp_integrate_cuda_prototype.cu") == False or
     os.path.isfile(smcl + "specific.mk") == False or
     os.path.isfile(smcl + "Makefile.m") == False
     ):
    print("Can't find one or more files. \n")
    print("--> Run the script at ./messy/util directory of messy. \n")
    print("Exiting... \n")
    exit(-1)


multifile, vectorize, indirect, inject_rconst = get_transformation_flags()


###############################################
### Backup files
print("==> Step 0: Backup files.\n")

shutil.copyfile(smcl + "specific.mk", smcl + "specific.mk.old")
shutil.copyfile(smcl + "Makefile.m", smcl + "Makefile.m.old")
shutil.copyfile(smcl + "messy_mecca_kpp.f90", smcl + "messy_mecca_kpp.f90.old")
os.remove(smcl + "messy_mecca_kpp.f90")


# Open the files
file_messy_mecca_kpp = open(smcl + "messy_mecca_kpp.f90.old","r")
file_messy_cmn_photol_mem = open(smcl + "messy_cmn_photol_mem.f90","r")
file_messy_main_constants_mem = open(smcl + "messy_main_constants_mem.f90","r")
file_prototype = open("./source/kpp_integrate_cuda_prototype.cu","r")
file_specific = open(smcl + "specific.mk","a")
file_makefile = open(smcl + "Makefile.m","r+")


###############################################
print("==> Step 1: Detect subroutines in the file.")

subroutine_names = ["ros_PrepareMatrix","kppSolve","kppDecomp","Jac_sp","Fun","update_rconst","Initialize"]

subroutines = {}
source_cuda = {}

# if multiple files then we have to extract the functions from multiple files
if (multifile == True):
    file_messy = open("messy_mecca_kpp_linearalgebra.f90","r")
    subroutines1 = find_subroutines(file_messy, ["KppSolve","kppDecomp"])

    file_messy = open("messy_mecca_kpp_integrator.f90","r")
    subroutines2 = find_subroutines(file_messy, ["ros_PrepareMatrix"])

    file_messy = open("messy_mecca_kpp_jacobian.f90","r")
    subroutines3 = find_subroutines(file_messy,  ["Jac_SP"])

    file_messy = open("messy_mecca_kpp_function.f90","r")
    subroutines4 = find_subroutines(file_messy, ["Fun"])

    file_messy = open("messy_mecca_kpp_rates.f90","r")
    subroutines5 = find_subroutines(file_messy, ["Update_RCONST"])

    file_messy = open("messy_mecca_kpp_initialize.f90","r")
    subroutines6 = find_subroutines(file_messy, ["Initialize"])


    subroutines = dict(  list(subroutines1.items()) + list(subroutines2.items()) + list(subroutines3.items()) + list(subroutines4.items()) + list(subroutines5.items())  + list(subroutines6.items()) )

else:
    subroutines = find_subroutines(file_messy_mecca_kpp, subroutine_names)




###############################################

print("\n==> Step 2: Replacing variables.")

source_cuda["defines_vars_1"] = generate_define_vars(file_messy_main_constants_mem,["R_gas","atm2Pa","N_A"])
source_cuda["defines_ind_1"] = generate_define_indices_one_line(file_messy_cmn_photol_mem,"ip")




if (multifile == True):
    file_messy_mecca_kpp_global = open("messy_mecca_kpp_global.f90","r")
    file_messy_mecca_kpp_parameters = open("messy_mecca_kpp_parameters.f90","r")
    source_cuda["defines_vars_2"] = generate_define_vars(file_messy_mecca_kpp_parameters,["NSPEC","NVAR","NFIX","NREACT","LU_NONZERO","NBSIZE"])
    source_cuda["defines_vars_2"].append(generate_define_NBSIZE(subroutines["jac_sp"]))
    source_cuda["defines_ind_2"] = generate_define_indices_many_lines(file_messy_mecca_kpp_parameters,"ind")
    source_cuda["defines_ind_3"] = generate_define_indices_one_line(file_messy_mecca_kpp_global,"ihs")
    source_cuda["defines_ind_4"] = generate_define_indices_one_line(file_messy_mecca_kpp_global,"iht")
    source_cuda["defines_ind_5"] = generate_definitions_global(file_messy_mecca_kpp_global ,["k_","f_","a_"])
else:
    source_cuda["defines_vars_2"] = generate_define_vars(file_messy_mecca_kpp,["NSPEC","NVAR","NFIX","NREACT","LU_NONZERO","NBSIZE"])
    source_cuda["defines_vars_2"].append(generate_define_NBSIZE(subroutines["jac_sp"]))
    source_cuda["defines_ind_2"] = generate_define_indices_many_lines(file_messy_mecca_kpp,"ind")
    source_cuda["defines_ind_3"] = generate_define_indices_one_line(file_messy_mecca_kpp,"ihs")
    source_cuda["defines_ind_4"] = generate_define_indices_one_line(file_messy_mecca_kpp,"iht")
    source_cuda["defines_ind_5"] = generate_definitions_global(file_messy_mecca_kpp,["k_","f_","a_"])



# read the values
NSPEC = int(source_cuda["defines_vars_2"][0].split(" ")[2].strip())
NVAR = int(source_cuda["defines_vars_2"][1].split(" ")[2].strip())
NFIX = int(source_cuda["defines_vars_2"][2].split(" ")[2].strip())
NREACT = int(source_cuda["defines_vars_2"][3].split(" ")[2].strip())
LU_NONZERO = int(source_cuda["defines_vars_2"][4].split(" ")[2].strip())
NBSIZE = int(source_cuda["defines_vars_2"][5].split(" ")[2].strip())


# read the tables
if (multifile == True):
    file_messy_jacobian = open("messy_mecca_kpp_jacobiansp.f90","r")
    lu_diag = find_LU_DIAG(file_messy_jacobian, NVAR)
    lu_crow = find_LU_CROW(file_messy_jacobian, NVAR)
    lu_icol = find_LU_ICOL(file_messy_jacobian, NVAR)
else:
    lu_diag = find_LU_DIAG(file_messy_mecca_kpp, NVAR)
    lu_crow = find_LU_CROW(file_messy_mecca_kpp, NVAR)
    lu_icol = find_LU_ICOL(file_messy_mecca_kpp, NVAR)


###############################################
print("\n==> Step 3: Parsing function update_rconst.")

source = subroutines['update_rconst']
source = remove_comments(source)
source = strip_and_unroll_lines(source)
source = fix_power_op(source)
source = decapitalize_vars(source,["rconst","jx","khet_st","khet_tr","cair","press","temp","exp","log","max","min"])

# These are fixes with multifile: jx and khet are 2d
if (multifile == True):
    for i in range(len(source)):
        source[i] = source[i].replace("USE messy_main_constants_mem","")
        source[i] = source[i].replace("USE messy_cmn_photol_mem","k = is")
        source[i] = source[i].replace("jx(","jx(k,")
        source[i] = source[i].replace("khet_st(","khet_st(k,")
        source[i] = source[i].replace("khet_tr(","khet_tr(k,")


source = rconst_preprocessor_1(source)
rconst_ops,rconst_decls = split_rconst(source)
flocals = get_rconst_locals(rconst_decls)

source = subroutines['initialize']
source = remove_comments(source)
source = decapitalize_vars(source,["rconst"])
rinit  = create_rconst_init(source)

source_cuda["update_rconst"] = generate_update_rconst(rconst_ops,rconst_decls,flocals,rinit)

###############################################
print("\n==> Step 4: Parsing function kppsolve.")

source = subroutines['kppsolve']
source = remove_comments(source)
source = strip_and_unroll_lines(source)
source = fix_power_op(source)
source = split_beta(source,"X(")

source = fix_indices(source,[("X","K","istage"),("JVS","Ghimj")])
source = strip_lines(source)
source_cuda["kppsolve"] = generate_kppsolve(source)

###############################################
print("\n==> Step 5: Parsing function kppdecomp.")

source = subroutines['kppdecomp']
source = remove_comments(source)
source = strip_and_unroll_lines(source)
source = fix_power_op(source)


if ( indirect == True):
    source = split_beta(source,"DO k=1,NVAR")
    print("Indirect transformation.")
    source_cuda["kppdecomp"] = generate_kppDecompIndirect(source,NSPEC,lu_diag,lu_crow,lu_icol)
else:
    source = split_beta(source,"W(")
    source = fix_indices(source,[("W","W"),("JVS","Ghimj")])
    source_cuda["kppdecomp"] = generate_kppDecomp(source,NSPEC,lu_diag,lu_crow,lu_icol)




###############################################
print("\n==> Step 6: Parsing function jac_sp.")

source = subroutines["jac_sp"]
source = remove_comments(source)
source = strip_and_unroll_lines(source)
source = fix_power_op(source)
source = split_beta(source, "B(")
source = fix_indices(source,[("B","B"),("RCT","rconst"),("F","fix"),("V","var"),("JVS","jcb")])
source_cuda["jac_sp"] = generate_jac_sp(source, NBSIZE)


###############################################
print("\n==> Step 7: Parsing function fun.")

source = subroutines["fun"]
source = remove_comments(source)
source = strip_and_unroll_lines(source)
source = fix_power_op(source)
source = split_beta(source, "A(")
source = fix_indices(source,[("A","A"),("RCT","rconst"),("F","fix"),("V","var"),("Vdot","varDot")])
source_cuda["fun"] = generate_fun(source,NREACT)

###############################################
print("\n==> Step 8: Parsing and preparing diagonal.")


source_cuda["ros_preparematrix"] = generate_prepareMatrix(lu_diag)

###############################################
print("\n==> Step 9: Generating customized solver.")

source_cuda["special_ros"] = generate_special_ros(ros,inject_rconst)



###############################################
print("\n==> Step 10: Generating calls to customized solver.")

source_cuda["call_kernel"] = generate_special_ros_caller(ros)



###############################################

print("\n==> Step 11: Generating kpp_integrate_cuda.")

gen_kpp_integrate_cuda(file_prototype, source_cuda, inject_rconst)



###############################################

print("\n==> Step 12: Generating messy_mecca_kpp replacement.")
generate_c2f_interface(file_messy_mecca_kpp)

###############################################


print("\n==> Step 13: Modifying specific.mk and Makefile")

add_cuda_compilation(file_specific,file_makefile,arch)


###############################################

print("\n##################################################################\n")
print("Don't forget to add the '-lcudart' in the linking options during configuration")
print("For example, it can be added to the SPEC_NETCDF_LIB variable:")
print("SPEC_NETCDF_LIB = -L$EBROOTNETCDFMINFORTRAN/lib -lnetcdff   -lcudart  -lstdc++")


print_warning()


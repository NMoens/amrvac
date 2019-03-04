import numpy as np
import matplotlib.pyplot as plt
import os

#READ IN STRUCTURE FILE
################################################################################

struc = os.getcwd()+ '/init_struc_total'
struc = np.loadtxt(struc,dtype = 'string')

I,J = np.shape(struc)
dum = np.ones((I,J))

struc_rho = []
struc_tr = []
struc_pg = []
struc_tau = []
struc_mc = []
struc_gammar = []
struc_er = []
struc_y = []

for i in range(I):
    struc_rho.append(float(struc[i][0].replace('D', 'e'))) #CGS
    struc_tr.append(float(struc[i][1].replace('D', 'e'))) #CGS
    struc_pg.append(float(struc[i][2].replace('D', 'e'))) #CGS
    struc_tau.append(float(struc[i][3].replace('D', 'e'))) #CGS
    struc_mc.append(float(struc[i][4].replace('D', 'e'))) #CGS
    struc_gammar.append(float(struc[i][5].replace('D', 'e'))) #CGS
    struc_er.append(float(struc[i][6].replace('D', 'e'))) #CGS
    struc_y.append(float(struc[i][7].replace('D', 'e'))) #CGS

# READ IN PARAMETERS
################################################################################

params = os.getcwd()+ '/init_params'
dum = open(params, 'r')
dum1 = dum.read()
dum.close()
params = dum1

print params[1]

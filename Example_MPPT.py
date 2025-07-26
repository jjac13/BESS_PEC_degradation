# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a 
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 1.0

###############################################################################
###########################   Imported modules   ##############################

from MPPT_PO_controller import MPPT

###############################################################################
#################################   Main code  ################################
    
f_MPPT = 4e3
dt = 1/60/60#/10


V_MPPT = [0, 44.9]
P_MPPT = [0]
I_MPPT = [0]
G_list = [1, 1, 0.4, 1, 1] # 0.2, 0.6, 1, 0.5, 0.8, 0.9, 0.4 # 1, 0.2, 0.4, 0.6, 0.8, 1, 0.8, 0.6, 0.4, 0.2 # 1, 1, 1, 1, 1, 1, 1, 1
T_list = [20+273 for i in range(len(G_list))]
P_curtailment_list = [False, 300, 300, 300, False] #  [False, 200, 300, 400, False, 400, 300, 200] # [False for  i in range(len(G_list))] #
P_ref_plot = []
G_ref_plot = []


for G, T, P_curtailment in zip(G_list, T_list, P_curtailment_list):
    [V_i, I_i, P_i] = MPPT(V_MPPT[-1], G, T, P_curtailment, V_MPPT[-2], I_MPPT[-1], P_MPPT[-1], dt = dt)

    V_MPPT+=V_i
    P_MPPT+=P_i
    I_MPPT+=I_i
    
    for i in V_i:
        G_ref_plot.append(G)
        if P_curtailment:
            P_ref_plot.append(P_curtailment)
        else:
            P_ref_plot.append(416)

             
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family": "Times New Roman",
    'font.size': 16
})    
            

t = [i*1/f_MPPT for i in range(len(P_MPPT)-1)]
    
    
fig, (axP, axV, axI) = plt.subplots(3, sharex=True)
axV.plot(t, V_MPPT[2:]) # label = ('$P_{grid}^{min} = $'+str(round(ts, 2)))
axV.set_xlim([0, max(t)])
axV.set_ylabel('Voltage, $V_{PV}$, [V]') 
axV.grid()

axI.plot(t, I_MPPT[1:]) # label = ('$P_{grid}^{min} = $'+str(round(ts, 2)))
axI.set_xlim([0, max(t)])         
axI.set_ylabel('Current, $I_{PV}$, [I]') 
axI.set_xlabel('Time, t, [s]') 
axI.grid()        


axP.plot(t, P_ref_plot, 'r', label = "$P_{ref}$") # label = ('$P_{grid}^{min} = $'+str(round(ts, 2)))
axP.plot(t, P_MPPT[1:], 'b', label = "$P_{PV}$") # label = ('$P_{grid}^{min} = $'+str(round(ts, 2)))
axP.set_xlim([0, max(t)])
axP.set_ylabel('Power, $P$, [W]') 

axP2 = axP.twinx()
axP2.plot(t, G_ref_plot, 'c', label = "$G_{ref}$")
axP2.set_ylabel('Irradiance, $G$, [W/$m^{2}$]') 

axP.legend(loc='lower left', fancybox=True, shadow=False, ncol=1)
axP.grid()  
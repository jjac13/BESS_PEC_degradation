# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a 
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 1.0

###############################################################################
###########################   Imported modules   ##############################

import csvreader
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family": "Times New Roman",
    'font.size': 16
})  

###############################################################################
###########################   Functions definition  ###########################

def PV_degradation(year, p0 = 0.975, pf = 0.8, dt_years = 25):
    
    return year*(pf - p0)/dt_years + p0

###############################################################################

def PV_Current_Power(V, G, T):
    from math import exp

    # (5) (PDF) Evaluating MPPT converter topologies using a Matlab PV Model. Available from: https://www.researchgate.net/publication/37630747_Evaluating_MPPT_converter_topologies_using_a_Matlab_PV_Model [accessed Mar 05 2024].    
    
    # https://www.canadiansolar.com/wp-content/uploads/2019/12/Canadian_Solar-Datasheet-HiDM_CS1U-MS_EN.pdf
    # current given voltage, illumination and temperature
    # Ia,Va = array current,voltage
    # G = num of Suns (1 Sun = 1000 W/mˆ2)
    # T = Temp in Deg C
    
    k = 1.38e-23 # Boltzman’s const
    q = 1.60e-19 # charge on an electron
    
    # enter the following constants here, and the model will be
    # calculated based on these. for 1000W/mˆ2
    A = 1.2 # "diode quality" factor, =2 for crystaline, <2 for amorphous
    Vg = 1.12 # band gap voltage, 1.12eV for xtal Si,  ̃1.75 for amorphous Si.
    Ns = 40 # number of series connected cells (diodes)
    
    
    T1 = 273 + 25
    Voc_T1 = 53.4 /Ns # open cct voltage per cell at temperature T1
    Isc_T1 = 9.60 # short cct current per cell at temp T1
    
    T2 = 273 + 65
    Voc_T2 = 47.1 /Ns # open cct voltage per cell at temperature T2
    Isc_T2 = 9.9 # short cct current per cell at temp T2
    
#    TaK = 273 + T # array working temp
    TrK = 273 + 25 # reference temp    
    
    
    # when Va = 0, light generated current Iph_T1 = array short cct current
    # constant "a" can be determined from Isc vs T
    
    Iph_T1 = Isc_T1 * G#/1000
    a = (Isc_T2 - Isc_T1)/Isc_T1 * 1/(T2 - T1)
    Iph = Iph_T1 * (1 + a*(T - T1))
    Vt_T1 = k * T1 / q # = A * kT/q
    Ir_T1 = Isc_T1 / (exp(Voc_T1/(A*Vt_T1))-1)
    Ir_T2 = Isc_T2 / (exp(Voc_T2/(A*Vt_T1))-1)
    
    b = Vg * q/(A*k)
    Ir = Ir_T1 * (T/T1)**(3/A)*exp(-b*(1/T - 1/T1))
    
    X2v = Ir_T1/(A*Vt_T1) * exp(Voc_T1/(A*Vt_T1))
    dVdI_Voc = - 1.15/Ns / 2   # dV/dI at Voc per cell --
                                # from manufacturers graph
    Rs = - dVdI_Voc - 1/X2v # series resistance per cell
    
    # Ia = 0:0.01:Iph
    Vt_Ta = A * 1.38e-23 * T / 1.60e-19
    
    # = A * kT/q# Ia1 = Iph - Ir.*( exp((Vc+Ia.*Rs)./Vt_Ta) -1)
    # solve for Ia: f(Ia) = Iph - Ia - Ir.*( exp((Vc+Ia.*Rs)./Vt_Ta) -1) = 0
    # Newton’s method: Ia2 = Ia1 - f(Ia1)/f’(Ia1)
    
    Vc = V/Ns
    Ia = 0 #  zeros(len(Vc)).tolist()
    # Iav = Ia
    for j in range(5):
        Ia = Ia - (Iph - Ia - Ir*( exp((Vc+Ia*Rs)/Vt_Ta) -1))/ (-1 - (Ir*( exp((Vc+Ia*Rs)/Vt_Ta) -1))*Rs/Vt_Ta)
        # Iav = [IavIa] 
        # to observe convergence for debugging.

    
    return [Ia, V*Ia]
    
###############################################################################

def PV_curves(G, T, Vmin = 0, Vmax = 55, dV = 0.1, reuturn_lists = False, plots = False):

    V_range = [i*dV for i in range(int(Vmin/dV), int(Vmax/dV), 1)]

    IV_curve = []
    PV_curve = []    
    
    for V in V_range:
        [I, P] = PV_Current_Power(V, G, T)
        
        if I < 0 or P < 0:
            break
        
        IV_curve.append(I)
        PV_curve.append(P)

    if plots:
        adjust_factor = 1
    
        plt.figure(constrained_layout=True)
        plt.plot(V_range[0:len(IV_curve)], IV_curve)
        plt.grid()
        plt.xlim([0, Vmax])
        plt.xlabel('Voltage, $V_{PV}$, [V]')    
        plt.ylabel('Current, $I_{PV}$, [A]')
        plt.show()
            
        plt.figure(constrained_layout=True)
        plt.plot(V_range[0:len(PV_curve)], [i*adjust_factor for i in PV_curve])
        plt.grid()
        plt.xlim([0, Vmax])
        plt.xlabel('Voltage, $V_{PV}$, [V]')            
        plt.ylabel('Power, $P_{PV}$, [W]')    
        plt.show()

        fig, (axI, axP) = plt.subplots(2, sharex=True)
        axI.plot(V_range[0:len(IV_curve)], IV_curve) # label = ('$P_{grid}^{min} = $'+str(round(ts, 2)))
        axI.set_xlim([0, Vmax])
        axI.set_ylabel('Current, $I_{PV}$, [A]') 
        axI.grid()


        axP.plot(V_range[0:len(PV_curve)], [i*adjust_factor for i in PV_curve]) # label = ('$P_{grid}^{min} = $'+str(round(ts, 2)))
        axP.set_xlim([0, Vmax])          
        axP.set_ylabel('Power, $P_{PV}$, [W]') 
        axP.set_xlabel('Voltage, $V_{PV}$, [V]') 
        axP.grid()     
        
    if not reuturn_lists:
        return [V_range[PV_curve.index(max(PV_curve))], IV_curve[PV_curve.index(max(PV_curve))], PV_curve[PV_curve.index(max(PV_curve))]]
        
    else:
            
        return [IV_curve, PV_curve]
    
###############################################################################    
    
def Curtailment_voltage(P_curtailed, G, T, Vmin = 0, Vmax = 60, dV = 0.1, reuturn_lists = True, plots = False):
    
    V_range = [i*dV for i in range(int(Vmin/dV), int(Vmax/dV), 1)]
    
    [IV_curve, PV_curve] = PV_curves(G, T, Vmin, Vmax, dV, reuturn_lists, plots)

    for P_PV in PV_curve:
        if P_curtailed <= P_PV:
            return V_range[PV_curve.index(P_PV)]
        
    if P_curtailed > max(PV_curve):        
        return V_range[PV_curve.index(max(PV_curve))]          

###############################################################################

def PO_algorithm(V_module, G, T, V_ref, I_ref, P_ref, Vmin = 0, Vmax = 55, dV = 0.1):
    
    [I_PV, P_PV] = PV_Current_Power(V_module, G, T)
         
    if round(P_PV, 2) == round(P_ref, 2):
        return [V_module, I_PV, P_PV]
    
    elif P_PV > P_ref:
        if V_module > V_ref:
            return [V_module + dV, I_PV, P_PV]
        else:
            return [V_module - dV, I_PV, P_PV]
    
    else:
        if V_module > V_ref:
            return [V_module - dV, I_PV, P_PV]
        else:
            return [V_module + dV, I_PV, P_PV]
    
###############################################################################

def PO_algorithm_curtailment(V_module, G, T, V_ref, I_ref, P_ref, P_curtailment, Vmin = 0, Vmax = 55, dV = 0.1):
    
    [I_PV, P_PV] = PV_Current_Power(V_module, G, T)
     
    if round(P_PV, 2) == round(P_curtailment, 2):
        return [V_module, I_PV, P_PV]    
    
    elif round(abs(P_PV - P_curtailment), 1) <= round(abs(P_ref - P_curtailment), 1):
        if V_module > V_ref:
            return [V_module + dV, I_PV, P_PV]
        else:
            return [V_module - dV, I_PV, P_PV]
    else:
        if V_module > V_ref:
            return [V_module - dV, I_PV, P_PV]
        else:
            return [V_module + dV, I_PV, P_PV]        

###############################################################################
    
def MPPT(V_module, G, T, P_curtailment = 0, V_ref = 0, I_ref = 0, P_ref = 0, dV = 0.01, f = 4e3, dt = 1/60, Vmin = 0, Vmax = 55, plots = False): # , dt = 0.25
    dt *= 3600
    
    V = []
    I = []
    P = []
    t = []

    if not P_curtailment:
    
        for ts in range(int(f*dt)):
            [V_MPPT, I_MPTT, P_MPPT] = PO_algorithm(V_module, G, T, V_ref, I_ref, P_ref, Vmin, Vmax, dV)
            V.append(V_MPPT)
            I.append(I_MPTT)
            P.append(P_MPPT)
            t.append(ts/f)
            V_ref = V_module
            
            V_module = V[-1]
            P_ref = P[-1]
    else:
    
        for ts in range(int(f*dt)):
            [V_MPPT, I_MPTT, P_MPPT] = PO_algorithm_curtailment(V_module, G, T, V_ref, I_ref, P_ref, P_curtailment, Vmin, Vmax, dV)
            V.append(V_MPPT)
            I.append(I_MPTT)
            P.append(P_MPPT)
            t.append(ts/f)
            V_ref = V_module
            
            V_module = V[-1]
            P_ref = P[-1]            
        
        
    if plots:

        fig, (axV, axI, axP) = plt.subplots(3, sharex=True)
        axV.plot(t, V) # label = ('$P_{grid}^{min} = $'+str(round(ts, 2)))
        axV.set_xlim([0, dt])         
        axV.set_ylabel('Voltage, $V_{PV}$, [V]') 
        axV.grid()
        
        axI.plot(t, I) # label = ('$P_{grid}^{min} = $'+str(round(ts, 2)))
        axI.set_xlim([0, dt])         
        axI.set_ylabel('Current, $I_{PV}$, [I]') 
        axI.grid()        


        axP.plot(t, P) # label = ('$P_{grid}^{min} = $'+str(round(ts, 2)))
        axP.set_xlim([0, dt])       
        axP.set_ylabel('Power, $P_{PV}$, [W]') 
        axP.set_xlabel('Time, t, [s]') 
        axP.grid()  
    
    return [V, I, P]
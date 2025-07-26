# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a 
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 1.0

###############################################################################
###########################   Imported modules   ##############################

from math import floor, ceil
from MPPT_PO_controller import *

###############################################################################
###########################   Functions definition  ###########################

def Yearly_Curtailed_Power(P_PV_registry, Sampling = True, sampling_frequency = 4, kW = True):

    if Sampling:
        sampled_PV_registry = [] 
        if kW:
            for year in P_PV_registry:      
                sampled_PV_registry.append([year[0][i*sampling_frequency]*1000 for i in range(int(len(year[0])/sampling_frequency))])
                
            return sampled_PV_registry        
        else:
                
            for year in P_PV_registry:      
                sampled_PV_registry.append([year[0][i*sampling_frequency] for i in range(int(len(year[0])/sampling_frequency))])
                
            return sampled_PV_registry

###############################################################################

def Yearly_Voltage(P_PV_registry, G, T_amb, Sampling = True, sampling_frequency = 4, n_modules = 5, n_series = 1):    
    Yearly_voltage_registry = []
    
    if Sampling:
        sampled_PV_registry = []
        l = 0
        for year in P_PV_registry:      
            sampled_PV_registry.append([year[0][i*sampling_frequency]/n_modules for i in range(int(len(year[0])/sampling_frequency))])
            
        for threshold_registry in sampled_PV_registry: 
            print('Threshold ', l)                  
            Yearly_voltage_registry.append([n_modules*Curtailment_voltage(threshold_registry[i]*1000, G[i*sampling_frequency], T_amb[i*sampling_frequency])/n_series for i in range(len(threshold_registry))])

            l+=1
        return Yearly_voltage_registry
    
    else:
        for year in P_PV_registry:
            Yearly_voltage_registry.append([n_modules*Curtailment_voltage(year[0][i]*1000/n_modules, G[i], T_amb[i])/n_series for i in range(len(year[0]))])
        
        return Yearly_voltage_registry

###############################################################################

def duty_cycle(V_PV, V_bus = 400):
    return 1 - V_PV/V_bus

###############################################################################

def i_av_switch(V_PV, P_PV, V_bus = 400):
    if V_PV:
        return P_PV*duty_cycle(V_PV)/V_PV
    else:
        return 0

###############################################################################

def i_av_diode(V_PV, P_PV, V_bus = 400):
    if V_PV:    
        return P_PV*(1-duty_cycle(V_PV))/V_PV
    else:
        return 0
  
###############################################################################    
  
def i_RMS_switch(V_PV, P_PV, V_bus = 400, f = 20*1e3, L = 1.45*1e-3):
    if V_PV:    
        from math import sqrt
        return sqrt(duty_cycle(V_PV)*(P_PV/V_PV)**2 + (1/3)*(V_PV*duty_cycle(V_PV)/(2*f*L))**2)
    else:
        return 0

###############################################################################

def i_RMS_diode(V_PV, P_PV, V_bus = 400, f = 20*1e3, L = 1.45*1e-3):
    if V_PV:    
        from math import sqrt
        return sqrt((1-duty_cycle(V_PV))*(P_PV/V_PV)**2 + (1/3)*(V_PV*duty_cycle(V_PV)/(2*f*L))**2)
    else:
        return 0

###############################################################################

def get_currents(V_PV_registry, P_PV_registry, n_series = 1, V_bus = 400, f = 20*1e3, L = 1.45*1e-3):
    i_av_switch_registry = []
    i_RMS_switch_registry = []
    i_av_diode_registry = []
    i_RMS_diode_registry = []
    
    for j in range(len(V_PV_registry)):
        i_av_switch_registry.append([i_av_switch(V_PV_registry[j][i], P_PV_registry[j][i]) for i in range(len(V_PV_registry[j]))])
        i_RMS_switch_registry.append([i_RMS_switch(V_PV_registry[j][i], P_PV_registry[j][i]) for i in range(len(V_PV_registry[j]))])
        i_av_diode_registry.append([i_av_diode(V_PV_registry[j][i], P_PV_registry[j][i]) for i in range(len(V_PV_registry[j]))])
        i_RMS_diode_registry.append([i_RMS_diode(V_PV_registry[j][i], P_PV_registry[j][i]) for i in range(len(V_PV_registry[j]))])

    return [i_av_switch_registry, i_RMS_switch_registry, i_av_diode_registry, i_RMS_diode_registry]

###############################################################################

def conduction_loss_IGBT(i_av, i_RMS, R_CE = 0.0856, VT = 1.198):
    
    return VT*i_av + R_CE*i_RMS**2

###############################################################################

def energy_losses_IGBT(i_av, i_RMS, dynamics = [0.0195, 0.011, 0.0005]):
    
    return (dynamics[0] + dynamics[1]*i_av + dynamics[2]*i_RMS**2)/1000

###############################################################################

def switching_losses_IGBT(i_av, i_RMS, f = 20*1e3, V_bus = 400, V_nom = 600):
    return f*energy_losses_IGBT(i_av, i_RMS)*V_bus/V_nom

###############################################################################

def Calc_HeatSink_R(IGBT_losses_Registry, Diode_losses_Registry, T_amb, R_jc_IGBT = 1.7, R_jc_diode = 4, T_j_max = 175):
    
    T_a = max(T_amb) - 273
    
    Q_IGBT = max(max(IGBT_losses_Registry))
    Q_diode = max(max(Diode_losses_Registry))
    
    R_ch_IGBT = (T_j_max - Q_IGBT*R_jc_IGBT - T_a)/(Q_IGBT + Q_diode)
    
    R_ch_diode = (T_j_max - Q_diode*R_jc_diode - T_a)/(Q_IGBT + Q_diode)
    
    
    return min([R_ch_IGBT, R_ch_diode])

###############################################################################
    
def Foster_impedance(R_i = [0.29566, 0.25779, 0.19382, 0.05279], tau_i = [6.478*1e-2, 6.12*1e-3, 4.679*1e-4, 6.45*1e-5], t = 100):
    from math import exp
    
    Z_jc = 0
    
    for R, tau in zip(R_i, tau_i):
        Z_jc += R*(1-exp(-t/tau))
        
    return Z_jc

###############################################################################

def switch_total_losses_IGBT(i_av_switch_registry, i_RMS_switch_registry, R_CE = 0.0856, VT = 1.198, f = 20*1e3, V_bus = 400, V_nom = 600):
    
    return [[conduction_loss_IGBT(av, RMS, R_CE, VT) + switching_losses_IGBT(av, RMS, f, V_bus, V_nom) for av, RMS in zip(i_av, i_RMS)] for i_av, i_RMS in zip(i_av_switch_registry, i_RMS_switch_registry)]

###############################################################################

def switch_temperature_IGBT(T_amb, i_av, i_RMS, R_CE = 0.0856, VT = 1.198, f = 20*1e3, V_bus = 400, V_nom = 600, R_th = [1.7, 8], R_tch = 60, R_i = [0.29566, 0.25779, 0.19382, 0.05279], tau_i = [6.478*1e-2, 6.12*1e-3, 4.679*1e-4, 6.45*1e-5], t = 100):

    if not R_th:    
        return (conduction_loss_IGBT(i_av, i_RMS, R_CE, VT) + switching_losses_IGBT(i_av, i_RMS, f, V_bus, V_nom))*(Foster_impedance(R_i, tau_i, t) + R_tch) + T_amb
    else:
        return (conduction_loss_IGBT(i_av, i_RMS, R_CE, VT) + switching_losses_IGBT(i_av, i_RMS, f, V_bus, V_nom))*(R_th[0] + R_th[1]) + T_amb

###############################################################################

def get_switch_temperature_IGBT(T_amb, i_av_switch_registry, i_RMS_switch_registry, R_CE = 0.0856, VT = 1.198, f = 20*1e3, V_bus = 400, V_nom = 600, R_th = [1.7, 8], R_tch = 39, R_i = [0.29566, 0.25779, 0.19382, 0.05279], tau_i = [6.478*1e-2, 6.12*1e-3, 4.679*1e-4, 6.45*1e-5], t = 100):

    return [[switch_temperature_IGBT(T, i_av, i_RMS, R_CE, VT, f, V_bus, V_nom, R_th, t = t) for T, i_av, i_RMS in  zip(T_amb, i_av_switch, i_RMS_switch)] for i_av_switch,i_RMS_switch in zip(i_av_switch_registry, i_RMS_switch_registry)]

###############################################################################

def conduction_loss_MOSFET(i_RMS, R_DS = 0.7):
    
    return R_DS*i_RMS**2

###############################################################################

def switching_losses_MOSFET(i_RMS, V_DS = 400,  f = 20*1e3, t_r = 40*1e-9, t_f = 35*1e-9):

    return V_DS*i_RMS*f*(t_r + t_f)/2    

###############################################################################

def switch_temperature_MOSFET(T_amb, i_RMS, V_DS = 400,  f = 20*1e3, R_DS = 0.7, t_r = 40*1e-9, t_f = 35*1e-9, R_chc = 2.5, R_cha = 62.5):
    
    return (conduction_loss_MOSFET(i_RMS, R_DS) + switching_losses_MOSFET(i_RMS, V_DS, f, t_r, t_f))*(R_chc + R_cha) + T_amb

###############################################################################

def get_switch_temperature_MOSFET(T_amb, i_RMS_switch_registry, V_DS = 400,  f = 20*1e3, R_DS = 0.7, t_r = 40*1e-9, t_f = 35*1e-9, R_chc = 2.5, R_cha = 62.5):
    
    return [[switch_temperature_MOSFET(T, i_RMS, V_DS,  f, R_DS, t_r, t_f, R_chc, R_cha) for T, i_RMS in zip(T_amb, i_RMS_switch)] for i_RMS_switch in i_RMS_switch_registry]

###############################################################################

def get_t_on(V_PV, V_bus = 400, f = 20*1e3):
    
    return (V_bus - V_PV)/(f*V_bus)

###############################################################################

def Bayerer_method(DT, T_min, t_on):
    from math import exp
    # Bayerer Lifetime model
    A = 9.34*1e14
    b1 = -4.416
    b2 = 1285
    b3 = -0.463
    b4 = -0.716
    b5 = -0.761
    b6 = -0.5
    I_b = 10
    V_b = 600/100
    D = 0.45*1e-3
      
    return A*pow(DT, b1)*exp(b2/(T_min))*pow(t_on, b3)*pow(I_b, b4)*pow(V_b, b5)*pow(D, b6)
       
###############################################################################

def damage(T_Switch, t_on_list = 0, plot_histogram = True, dt = 3600):
    import rainflow
    
    delta_T = []
    T_min = []
    
    if not t_on_list:
        t_on_list = []
        
    n = 0
    
    for rng, mean, count, i_start, i_end in rainflow.extract_cycles(T_Switch):
        delta_T.append(abs(T_Switch[i_end] - T_Switch[i_start]))
        T_min.append(min([T_Switch[i_end], T_Switch[i_start]]))
        t_on_list.append(dt*(i_end - i_start))
        
        n += 1
        
    D = 0        
    for DT, T, t_on in zip(delta_T, T_min, t_on_list):
        D +=  1/Bayerer_method(DT, T, t_on)
        
    if plot_histogram:
        import matplotlib.pyplot as plt
    
        plt.rcParams.update({
            "font.family": "Times New Roman",
            'font.size': 16            
        })       
      
    
        plt.figure(constrained_layout=True)
        plt.hist(delta_T, bins = 50)
        plt.xlabel('Thermal_cycle, $\Delta$T, [K]')    
        plt.ylabel('Frequency')           
        plt.show()

        plt.figure(constrained_layout=True)
        plt.hist(T_min, bins = 50)
        plt.xlabel('Minimum temperature, T$_{min}$, [K]')    
        plt.ylabel('Frequency')           
        plt.show() 

        plt.figure(constrained_layout=True)
        plt.hist(t_on_list, bins = 50)
        plt.xlabel('Time on, t$_{on}$, [s]')    
        plt.ylabel('Frequency')           
        plt.show()           
        
    return D

###############################################################################

def get_lifetime(T_Switch_registry, P_Grid_min_list = 0, t_on_registry = 0, plotting = True, return_Lifetime = True):
    
    if not t_on_registry:
        D = [damage(T_Switch, plot_histogram = False) for T_Switch in T_Switch_registry]        
    else:
        
        D = [damage(T_Switch, T_on_list, False) for T_Switch, T_on_list in zip(T_Switch_registry, t_on_registry)]
    
    L = [1/i for i in D]
    
    if plotting:
        import matplotlib.pyplot as plt
        plt.rcParams.update({
            "font.family": "Times New Roman",
            'font.size': 16
        })     
        
        plt.figure(constrained_layout=True)
        plt.plot(P_Grid_min_list, D)
        plt.xlim([int(min(P_Grid_min_list)), int(max(P_Grid_min_list))])    
        plt.xlabel('Curtailment threshold, $P_{grid}^{min}$, [kW]')    
        plt.ylabel('Damage, $D$ [1/year]')     
        plt.grid()
        plt.show()     
    
        plt.figure(constrained_layout=True)
        plt.plot(P_Grid_min_list, [i/max(L) for i in L])
        plt.xlim([int(min(P_Grid_min_list)), int(max(P_Grid_min_list))])   
        plt.xlabel('Curtailment threshold, $P_{grid}^{min}$, [kW]')    
        plt.ylabel('Relative lifetime, $L$, [-]')     
        plt.grid()
        plt.show() 
        
    if return_Lifetime:
        return L

###############################################################################

def plot_I_T_day(i_switch_registry, T_IGBT, P_Grid_min_list, i_plot = [0, 5, 10, 15, 20], day = 172, dt = 1): # 15, 172
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "Times New Roman",
        'font.size': 16
    })        

    
    fig, (axi, axT) = plt.subplots(2, sharex=True)
    for i in i_plot:
        ts = P_Grid_min_list[i]
        axi.plot(i_switch_registry[i], label = (str(round(ts, 2))+' kW')) # label = ('$P_{grid}^{min} = $'+str(round(ts, 2)))
        axi.set_xlim([day*24, (day+1)*24])
        axi.set_ylim([0, 10])                   
        axi.set_ylabel(r'$\bar{i}_{s}$, [A]') 
        axi.grid()

    for i in i_plot:
        ts = P_Grid_min_list[i]
        axT.plot([j-273 for j in T_IGBT[i]], label = (str(round(ts, 2))+' kW')) # label = ('$P_{grid}^{min} = $'+str(round(ts, 2))) 
        axT.set_xlim([day*24, (day+1)*24])   
        axT.set_ylim([0, 200])     
        axT.set_xlabel('Time, [hour]') 
        axT.set_ylabel('$T_{j}$, [°C]') 
        axT.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=5)
        axT.grid()
        
###############################################################################    

def BESS_degradation_NMC(SoC_0, SoC_f, Capacity_BESS, Capacity_BESS_BOL = 3.36, T = 20 + 273, dt = 0.25):
    from numpy import exp
    
    SoC_0*=100
    SoC_f*=100
    
    DOD = abs(SoC_f - SoC_0)
    mSOC = (SoC_f + SoC_0)/2
    N = Capacity_BESS*DOD/(2*Capacity_BESS_BOL)
    
    if SoC_0 > SoC_f:   # Discharge
        Cch = 0
        Cdc = Capacity_BESS*(DOD/100)/(Capacity_BESS_BOL*dt)

    elif SoC_0 < SoC_f:   # Charge
        Cch = Capacity_BESS*(DOD/100)/(Capacity_BESS_BOL*dt)
        Cdc = 0

    else:               # Rest
        Cch = 0
        Cdc = 0  
        

    beta = 0.001673
    k_t = 21.6745
    k_dod = 0.022
    k_ch = 0.2553
    k_dc = 0.1571
    k_soc = -0.0212
    alpha = 0.915
    T_ref = 293
    soc_ref = 42
    
    delta = beta*exp(k_t*(T - T_ref)/T + k_dod*DOD + k_ch*Cch + k_dc*Cdc)*(1+k_soc*mSOC*(1-(mSOC/(2*soc_ref))))
    
    soh = 100 - delta*pow(N,alpha)
    return soh*Capacity_BESS_BOL/100
 
###############################################################################    
 
def BESS_degradation_LFP(SoC_0, SoC_f, Capacity_BESS, Capacity_BESS_BOL = 3.36, T = 20 + 273, V = 62.7, model = 'Wang', dt = 0.25):
    #Modelling the cycling degradation of Li-ion batteries: Chemistry influenced stress factors - 10.1016/j.est.2021.102765    
    from numpy import exp

    if model == 'Olmos':    
        SoC_0*=100
        SoC_f*=100
    
        DOD = abs(SoC_f - SoC_0)
        mSOC = (SoC_f + SoC_0)/2
        N = Capacity_BESS*DOD/(2*Capacity_BESS_BOL)
        
        if SoC_0 > SoC_f:   # Discharge
            Cch = 0
            Cdc = Capacity_BESS*(DOD/100)/(Capacity_BESS_BOL*dt)
    
        elif SoC_0 < SoC_f:   # Charge
            Cch = Capacity_BESS*(DOD/100)/(Capacity_BESS_BOL*dt)
            Cdc = 0
    
        else:               # Rest
            Cch = 0
            Cdc = 0       
            
        
        beta = 0.003414
        k_t = 5.8755
        k_dod = -0.0046
        k_ch = 0.1038
        k_dc = 0.296
        k_soc = 0.0513
        alpha = 0.869
        T_ref = 293
        soc_ref = 42
        
        delta = beta*exp(k_t*(T - T_ref)/T + k_dod*DOD + k_ch*Cch + k_dc*Cdc)*(1+k_soc*mSOC*(1-(mSOC/(2*soc_ref))))
        
        soh = 100 - delta*pow(N,alpha)
        return soh*Capacity_BESS/100


#########
    elif model == 'Vermeer':
        from math import sqrt
        
        isa = 1000*Capacity_BESS*abs(SoC_f - SoC_0)/(V*dt)
    
        c=[0.0008/3600, 0.39, 1.035, 50, 14.876/sqrt(24*3600)] # aging parameters
        R = 8.314 # Regnault constant 

        dt *= 3600
        
        ilossCycle = c[0]*c[2]/c[3]*exp(c[1]*abs(isa))*(1-SoC_0)*abs(isa)*dt # cyclic aging   
        ilossCal = c[4]*sqrt(dt)*exp(-24e3/R/T) # calendar aging
        iloss = ilossCycle + ilossCal # total aging
        
        dt /= 3600        
        
        return Capacity_BESS - iloss*V*dt
    
    elif model == 'Wang':
        from math import sqrt
        
        a = 8.61*1e-6
        b = -5.13*1e-3
        c = 7.63*1e-1
        d = -6.7*1e-3
        e = 2.35
        f = 14.876/sqrt(24)
        E_a = -24.5*1000
        R = 8.314 # Regnault constant         
        
        DOD = abs(SoC_f - SoC_0)        
        i_rate = DOD/dt # Capacity_BESS*(DOD)/(Capacity_BESS_BOL*dt)
        Ah_throughput = 1000*Capacity_BESS*DOD/(V) #1000*Capacity_BESS*DOD/(V*dt)

    
        
        ilossCycle = (a*T**2 + b*T + c) * exp((d*T+e)*i_rate) * Ah_throughput # cyclic aging   
        ilossCal = f * sqrt(dt) * exp(E_a/(R*T)) # calendar aging
        iloss = ilossCycle - ilossCal # total aging

        return Capacity_BESS*(100+iloss)/100

###############################################################################
        
def BESS_perm_min(SoC, Capacity_BESS = 3.36, SoCmax = 0.9, P_BESS_max = 1.28, dt = 0.25):
    from numpy import clip
    
    return clip(Capacity_BESS*(SoC - SoCmax)/dt, -P_BESS_max, P_BESS_max)
    
###############################################################################
    
def BESS_perm_max(SoC, Capacity_BESS = 3.36, SoCmin = 0.2, P_BESS_max = 1.28, dt = 0.25):
    from numpy import clip
    
    return clip(Capacity_BESS*(SoC - SoCmin)/dt, -P_BESS_max, P_BESS_max)

###############################################################################

# Enphase IQ3 https://enphase.com/download/iq-battery-3-data-sheet
def update_BESS_PeakShaving(SoC_0, P_Load, P_PV, SoCmax = 0.9, SoCmin = 0.2, P_BESS_max = 1.28, P_Grid_max = 1.5, P_Grid_min = -1.5, Capacity_BESS = 3.36, charge_efficiency = 0.943, discharge_efficiency = 0.943, P_SD = 0, dt = 0.25):
    
    
    E_BESS_0 = Capacity_BESS*SoC_0
    
    if SoC_0 >= SoCmax:                  # Battery can only discharge

        if (P_Load-P_PV) <= P_Grid_max:        # No peakshaving needed
            P_BESS = 0
            E_BESS = E_BESS_0*(1-P_SD)
            SoC_BESS = E_BESS/Capacity_BESS
            P_Grid = P_Load - P_PV 
        
        else:                                                       # Peakshaving needed
            if P_Load - P_PV - P_Grid_max  <= BESS_perm_max(SoC_0, P_BESS_max = P_BESS_max, Capacity_BESS = Capacity_BESS):    # Below the BESS max power
                P_BESS = P_Load - P_PV - P_Grid_max
                E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS/Capacity_BESS
                P_Grid = P_Grid_max
                
            else:                                                       # Above the BESS max power
                P_BESS = BESS_perm_max(SoC_0, P_BESS_max = P_BESS_max, Capacity_BESS = Capacity_BESS)
                E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS/Capacity_BESS
                P_Grid = -BESS_perm_max(SoC_0, P_BESS_max = P_BESS_max, Capacity_BESS = Capacity_BESS) - P_PV + P_Load     

    
    elif SoC_0 > SoCmin:                  # Battery can charge and discharge
        
        if (P_Load-P_PV) <= P_Grid_max:        # No peakshaving needed
            
            if P_Load >= P_PV:                    # PV below demanded load
                P_BESS = 0
                E_BESS = E_BESS_0*(1-P_SD)
                SoC_BESS = E_BESS/Capacity_BESS
                P_Grid = P_Load - P_PV     

            else:                                                   # Surplus of PV power
                if P_PV - P_Load <= -BESS_perm_min(SoC_0, P_BESS_max = P_BESS_max, Capacity_BESS = Capacity_BESS):    # Below the BESS max power
                    P_BESS = P_Load - P_PV
                    E_BESS = E_BESS_0 - P_BESS*dt*charge_efficiency
                    E_BESS = E_BESS*(1-P_SD)
                    SoC_BESS = E_BESS/Capacity_BESS
                    P_Grid = 0
                    
                else:                                                       # Above the BESS max power
                    P_BESS = BESS_perm_min(SoC_0, P_BESS_max = P_BESS_max, Capacity_BESS = Capacity_BESS)
                    E_BESS = E_BESS_0 - P_BESS*dt*charge_efficiency
                    E_BESS = E_BESS*(1-P_SD)
                    SoC_BESS = E_BESS /  Capacity_BESS
                    P_Grid = -BESS_perm_min(SoC_0, P_BESS_max = P_BESS_max, Capacity_BESS = Capacity_BESS) - P_PV + P_Load

                    
        else:                                                       # Peakshaving needed
            if P_Load - P_PV - P_Grid_max <= BESS_perm_max(SoC_0, P_BESS_max = P_BESS_max, Capacity_BESS = Capacity_BESS):    # Below the BESS max power
                P_BESS = P_Load - P_PV - P_Grid_max
                E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS /  Capacity_BESS
                P_Grid = P_Grid_max

                
            else:                                                       # Above the BESS max power
                P_BESS = BESS_perm_max(SoC_0, P_BESS_max = P_BESS_max, Capacity_BESS = Capacity_BESS)
                E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS /  Capacity_BESS
                P_Grid = -BESS_perm_max(SoC_0, P_BESS_max = P_BESS_max, Capacity_BESS = Capacity_BESS) - P_PV + P_Load  


    
    else: # self.BESS.SoC[1,i] <= SoCmin:                 # Battery can only charge     
        
        if P_Load >= P_PV:                        # PV below demanded load
            P_BESS = 0
            E_BESS = E_BESS_0*(1-P_SD)
            SoC_BESS = E_BESS/Capacity_BESS
            P_Grid = P_Load - P_PV 
             
        else:                                                       # Surplus of PV power
            
            if P_PV - P_Load <= -BESS_perm_min(SoC_0, P_BESS_max = P_BESS_max, Capacity_BESS = Capacity_BESS):    # Below the BESS max power
                P_BESS = P_Load - P_PV
                E_BESS = E_BESS_0 - P_BESS*dt*charge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS /  Capacity_BESS
                P_Grid = 0

                
            else:                                                       # Above the BESS max power
                P_BESS = BESS_perm_min(SoC_0, P_BESS_max = P_BESS_max, Capacity_BESS = Capacity_BESS)
                E_BESS = E_BESS_0 - P_BESS*dt*charge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS /  Capacity_BESS
                P_Grid = -BESS_perm_min(SoC_0, P_BESS_max = P_BESS_max, Capacity_BESS = Capacity_BESS) - P_PV + P_Load      
                
    return [SoC_BESS, P_Grid, P_BESS, E_BESS]   

###############################################################################
    

def yearly_PeakShaving(PV_degradation, P_Load, P_PV_av, P_BESS_0, SoC_BESS_0, P_Grid_0, E_BESS_0, t, SoCmax = 0.9, SoCmin = 0.2, P_BESS_max = 1.28, P_Grid_max = 1.5, P_Grid_min = -1.5, Capacity_BESS = 3.36, Capacity_BESS_BOL = 3.36, charge_efficiency = 0.943, discharge_efficiency = 0.943, P_SD = 0, dt = 60*15, t_final = 365*24*4, replace_BESS = False):

    P_BESS = [P_BESS_0]    
    SoC_BESS = [SoC_BESS_0]
    P_Grid = [P_Grid_0]
    E_BESS = [E_BESS_0]
    Capacity_BESS = [Capacity_BESS]
    
    T_TESS = [0]
    Qdot_TESS = [0]
    Qdot_HP = [0]
    Qdot_HP_TESS = [0]
    
    
    P_PV = [i*PV_degradation for i in P_PV_av]

    for step in range(t_final-1):

        # BESS
        [SoC_BESS_state, P_Grid_state, P_BESS_state, E_BESS_state] = update_BESS_PeakShaving(SoC_BESS[-1], P_Load[step], P_PV[step], P_Grid_max = P_Grid_max, P_BESS_max = P_BESS_max, Capacity_BESS = Capacity_BESS[-1])
        
        P_BESS.append(P_BESS_state)    
        SoC_BESS.append(SoC_BESS_state)
        P_Grid.append(P_Grid_state)
        E_BESS.append(E_BESS_state)

        T_TESS.append(0)
        Qdot_TESS.append(0)
        Qdot_HP.append(0)
        Qdot_HP_TESS.append(0)        
        
        if not replace_BESS:
            Capacity_BESS.append(BESS_degradation_LFP(SoC_BESS[-2], SoC_BESS[-1], Capacity_BESS[-1], Capacity_BESS_BOL))
        else:
            if Capacity_BESS[-1] > 0.8*Capacity_BESS_BOL:
                Capacity_BESS.append(BESS_degradation_LFP(SoC_BESS[-2], SoC_BESS[-1], Capacity_BESS[-1], Capacity_BESS_BOL))
            else:
                Capacity_BESS.append(Capacity_BESS_BOL)
    
        t.append(dt*step/3600)
        
    return [P_PV, P_BESS, SoC_BESS, P_Grid, E_BESS, Capacity_BESS, T_TESS, Qdot_TESS, Qdot_HP, Qdot_HP_TESS, t]
    
###############################################################################

def PeakShaving_Projection(P_Load, P_PV_av, P_BESS_max, P_Grid_max, Capacity_BESS_0, Capacity_BESS_BOL, P_BESS_0 = 0, SoC_BESS_0 = 0.5, P_Grid_0 = 0, t = [0], replace_BESS = False, simulation_years = 26):

    E_BESS_0 = SoC_BESS_0*Capacity_BESS_0
    
    P_PV = []
    P_BESS = []    
    SoC_BESS = []
    P_Grid = []
    E_BESS = []
    Capacity_BESS = []
    T_TESS = []
    Qdot_TESS = []
    Qdot_HP = []
    Qdot_HP_TESS = []
   
    
    for year in range(simulation_years):
    
        [P_PV_year, P_BESS_year, SoC_BESS_year, P_Grid_year, E_BESS_year, Capacity_BESS_year, T_TESS_year, Qdot_TESS_year, Qdot_HP_year, Qdot_HP_TESS_year, t_year] = yearly_PeakShaving(PV_degradation(year), P_Load, P_PV_av, P_BESS_0, SoC_BESS_0, P_Grid_0, E_BESS_0, t, P_BESS_max = P_BESS_max, P_Grid_max = P_Grid_max, Capacity_BESS = Capacity_BESS_0, Capacity_BESS_BOL = Capacity_BESS_BOL, replace_BESS = replace_BESS)
        
        P_BESS_0 = P_BESS_year[-1]
        SoC_BESS_0 = SoC_BESS_year[-1]
        P_Grid_0 = P_Grid_year[-1]
        E_BESS_0 = E_BESS_year[-1]
        Capacity_BESS_0 = Capacity_BESS_year[-1]
    
        P_PV.append(P_PV_year)
        P_BESS.append(P_BESS_year)    
        SoC_BESS.append(SoC_BESS_year)
        P_Grid.append(P_Grid_year)
        E_BESS.append(E_BESS_year)  
        Capacity_BESS.append(Capacity_BESS_year)
        T_TESS.append(T_TESS_year)
        Qdot_TESS.append(Qdot_TESS_year)
        Qdot_HP.append(Qdot_HP_year)
        Qdot_HP_TESS.append(Qdot_HP_TESS_year)
        
    return [t_year, P_PV, SoC_BESS, E_BESS, P_BESS, P_Grid, Capacity_BESS, T_TESS, Qdot_TESS, Qdot_HP, Qdot_HP_TESS]

###############################################################################

def SoH_estimation(Capacity_BESS):
    
    return [[100*ts/Capacity_BESS[0][0] for ts in year] for year in Capacity_BESS]

###############################################################################

def EoL_BESS(SoH, dt = 0.25):

    l = []
    
    for year in SoH:
        l += year
        
    return sum([i>80 for i in l])*dt/(24*365)   
    

###############################################################################

def Energy_EoL(P_BESS, EoL, dt = 0.25, accumulated = True):
    years, days = divmod(EoL, 1)
    P_BESS_EoL = []
    
    for year in P_BESS[0:int(years)]:
        P_BESS_EoL += year
    
    if years < len(P_BESS)-1:
    
        P_BESS_EoL += P_BESS[int(years)+1][0:int(days*365*24/dt)]
        
    if accumulated:
        return -sum([i for i in P_BESS_EoL if i<0])*dt
    else:
        return P_BESS_EoL

###############################################################################

def LCoS_calc(Discharged_energy, CAPEX = 6500, OPEX = 0):
    
    return (CAPEX + OPEX)/(Discharged_energy)

###############################################################################    
    
def plot_BESS_PeakShaving(t, P_Load, P_PV_av, SoC_BESS, E_BESS, P_BESS, P_Grid, Capacity_BESS, SoH, start_day = 0, end_day = 365, year = 0, dt = 60*15):
    import matplotlib.pyplot as plt
    from numpy import arange

    plt.rcParams.update({
        "font.family": "Times New Roman",
        'font.size': 16
    })
    
    plt.figure(constrained_layout=True)
    plt.plot(t[int(start_day*24*3600/dt): int(end_day*24*3600/dt)], P_PV_av[year][int(start_day*24*3600/dt): int(end_day*24*3600/dt)])
    plt.grid()
    plt.xlim([start_day*24, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    plt.ylabel('PV power, $P_{PV}$, [kW]')
    plt.show()  
    
    plt.figure(constrained_layout=True)
    plt.plot(t[int(start_day*24*3600/dt): int(end_day*24*3600/dt)], [i*100 for i in SoC_BESS[year]])
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.ylabel('State-of-charge of the BESS, $SoC_{BESS}$ [%]')
    plt.show()
    
    plt.figure(constrained_layout=True)
    plt.plot(t[int(start_day*24*3600/dt): int(end_day*24*3600/dt)], P_BESS[year], 'b', label='Power delivered by the BESS')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    plt.ylabel('Power of the BESS, $P_{BESS}$, [kW]')
    plt.show()


    plt.figure(constrained_layout=True)
    plt.plot(t[int(start_day*24*3600/dt): int(end_day*24*3600/dt)], [(Power/(dt/3600))/Cap if Power!=0 else 0 for Cap,Power in zip(Capacity_BESS[year],P_BESS[year])], 'b', label='BESS C-rate')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    plt.ylabel('BESS C-rate, $c_{BESS}$, [-]')
    plt.show()
    
    plt.figure(constrained_layout=True)
    plt.plot(t[int(start_day*24*3600/dt): int(end_day*24*3600/dt)], P_Grid[year], 'b', label='Power consumed from the grid')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    plt.ylabel('Power from the grid, $P_{G}$, [kW]')
    plt.show()
    
    imbalance = [Load - PV for Load, PV in zip(P_Load, P_PV_av[year])]
    plt.figure(constrained_layout=True)
    plt.plot(t[int(start_day*24*3600/dt): int(end_day*24*3600/dt)], imbalance[int(start_day*24*3600/dt): int(end_day*24*3600/dt)], 'b', label='Power imbalance between the load and the PV')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    plt.ylabel('Power imbalance, $P_{i}$, [kW]')
    plt.show()

    plt.figure(constrained_layout=True)
    plt.plot(t[int(start_day*24*3600/dt): int(end_day*24*3600/dt)], Capacity_BESS[year][int(start_day*24*3600/dt): int(end_day*24*3600/dt)], 'b', label='Capacity of the BESS')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    plt.ylabel('Capacity of the BESS, $C_{BESS}$, [kWh]')
    plt.show()    

    E_Grid_net = [sum(year)*dt/3600 for year in P_Grid]
    E_Grid_ret = [-sum([i for i in year if i<0])*dt/3600 for year in P_Grid]
    plt.figure(constrained_layout=True)
    plt.plot(arange(len(P_Grid)), E_Grid_net, 'b', label='Net')
    plt.plot(arange(len(P_Grid)), E_Grid_ret, 'r', label='Returned')
    plt.grid()
    plt.xlim([0, 25])
    plt.ylabel('Anualized energy exchange with the grid, $E_{G}$, [kWh]')    
    plt.xlabel('Time [years]')
    plt.legend(loc='lower right')       
    plt.show()
    
    SoH = [100*year[0]/Capacity_BESS[0][0] for year in Capacity_BESS]
    plt.figure(constrained_layout=True)
    plt.plot(SoH)
    plt.grid()
    plt.xlim([0, 25])
    plt.ylim([0, 100])    
    plt.xlabel('Time [years]')    
    plt.ylabel('State-of-Health of the BESS [%]')    
    plt.show()    

###############################################################################

def PeakShaving_limit_violations(single_profile_P_Grid_EoL, P_Grid_max_list, single_profile_P_BESS_EoL, plotting = True, return_violations = False):
    from statistics import mean
    
    violations = []
    histogram_P_Grid_max_list = []

    for i in range(len(single_profile_P_BESS_EoL)):
        for j in single_profile_P_Grid_EoL[i]:
            if (j - P_Grid_max_list[i]) > 0:
                histogram_P_Grid_max_list.append(P_Grid_max_list[i])                
                violations.append(j - P_Grid_max_list[i])
    
    average_violations = [mean([(P_Grid - threshold) for P_Grid in P_Grid_threshold[0]  if P_Grid > threshold]) if len([(P_Grid - threshold) for P_Grid in P_Grid_threshold[0]  if P_Grid > threshold]) > 0 else 0 for threshold, P_Grid_threshold in zip(P_Grid_max_list, P_Grid_registry)]
    max_violations = [max([(P_Grid - threshold) for P_Grid in P_Grid_threshold[0]  if P_Grid > threshold]) if len([(P_Grid - threshold) for P_Grid in P_Grid_threshold[0]  if P_Grid > threshold]) > 0 else 0 for threshold, P_Grid_threshold in zip(P_Grid_max_list, P_Grid_registry)]

    if plotting:
        import matplotlib.pyplot as plt
    
        plt.rcParams.update({
            "font.family": "Times New Roman",
            'font.size': 16
        })        
        plt.figure(constrained_layout=True)
        plt.plot(P_Grid_max_list, average_violations, 'b', label='Average power demanded above the peak-shaving threshold')
        plt.plot(P_Grid_max_list, max_violations, 'r', label='Maximum power demanded above the peak-shaving threshold')      
        plt.grid()
        plt.xlim([int(min(P_Grid_max_list)), ceil(max(P_Grid_max_list))])
        plt.ylim(bottom=0)        
        plt.xlabel('Peak-shaving threshold, $P_{grid}^{max}$, [kW]')  
        plt.ylabel('Peak-shaving threshold non-compliance, $P_{grid}^{nc}$, [kW]')   
        plt.legend(loc='upper right')          
        plt.show()   
    

        plt.figure(constrained_layout=True)
        plt.hist2d(histogram_P_Grid_max_list, violations, bins=100, cmap='gist_heat_r', vmax = 10)
        plt.xlim([min(P_Grid_max_list), max(P_Grid_max_list)])
        plt.ylim([0, max(violations)])
        plt.xlabel('Peak-shaving threshold, $P_{grid}^{max}$, [kW]')   
        plt.ylabel('Peak-shaving non-compliance, $P_{grid}^{nc}$, [kW]')        
        plt.show()      
        
    if return_violations:
        
        return [average_violations, max_violations]
    
###############################################################################
    
def PowerCurtailment(P_Load, P_PV_av, P_Grid_min = -1.5, P_Grid_max = 1.5, dt = 0.25):
    
    P_PV_max = P_Load - P_Grid_min
    
    if P_PV_av <= P_PV_max:     # No curtailment needed
        
        P_PV = P_PV_av
        P_Grid = P_Load - P_PV_av
        P_PV_curtailed = 0
        
    else:                       # Curtailment is needed
        P_PV = P_PV_max
        P_Grid = P_Load - P_PV_max
        P_PV_curtailed = P_PV_av - P_PV_max       
        
    return [P_Grid, P_PV, P_PV_curtailed]

###############################################################################

def yearly_Curtailment(PV_degradation, P_Load, P_PV_av, P_Grid_0, t, P_Grid_min = -1.5, dt = 60*15, t_final = 365*24*4):
    
    P_PV_av = [i*PV_degradation for i in P_PV_av]
    P_Grid = [P_Grid_0]
    P_PV = [0]
    P_PV_curtailed = [0]
    
    for step in range(t_final-1):
    
        # BESS
        [P_Grid_state, P_PV_state, P_curtailed_state] = PowerCurtailment(P_Load[step], P_PV_av[step], P_Grid_min)
        
        P_PV.append(P_PV_state)
        P_PV_curtailed.append(P_curtailed_state)
        P_Grid.append(P_Grid_state)
    
        t.append(dt*step/3600)
        
    return [P_PV, P_PV_curtailed, P_Grid, t]

###############################################################################

def Curtailment_Projection(P_Load, P_PV_av, P_Grid_0, P_Grid_min, t = [0], simulation_years = 26):    
    P_PV = []
    P_PV_curtailed = []    
    P_Grid = [] 
    
    for year in range(simulation_years):
        
        [P_PV_year, P_PV_curtailed_year, P_Grid_year, t_year] = yearly_Curtailment(PV_degradation(year), P_Load, P_PV_av, P_Grid_0, t, P_Grid_min)
        
        P_Grid_0 = P_Grid_year[-1]   
        
        P_PV.append(P_PV_year)
        P_PV_curtailed.append(P_PV_curtailed_year)
        P_Grid.append(P_Grid_year)
    
    P_PV_av = [[i + j for i, j in zip(P_PV[year_plot], P_PV_curtailed[year_plot])] for year_plot in range(len(P_PV))]
        
    return(t_year, P_PV, P_PV_curtailed, P_PV_av, P_Grid)
    
###############################################################################
    
def LCoE_calc(P_PV, Investment, Operational = 0, degradation = 1, Maintenance = 0, Fuel = 0, Period = 25, discount_rate = 0.08, inflation_rate = 0.02):
    # Using Fedd-in tariff
    
    Total_Cost = Investment + replacement_OPEX(degradation)
    Inflation = 1
    
    for year in range(Period):
        Total_Cost += (Operational[year] + Maintenance + Fuel)*Inflation/((1 + discount_rate)**year)
        
        Inflation *= (1 + inflation_rate)
        
    return Total_Cost/sum([sum(i)*0.25 for i in P_PV])

###############################################################################

def replacement_OPEX(degradation = 1, replacement_cost = 1250, discount_rate = 0.08, Expected_life = 15, Period = 25):
    from math import floor
    
    OPEX = 0
    replacement_period = floor(degradation*Expected_life)    
    replacement_year = replacement_period


    while replacement_year <= Period:

        OPEX += replacement_cost/((1+discount_rate)**replacement_year)
        replacement_year += replacement_period
        
        
    return OPEX

###############################################################################

def replacement_LCoE(degradation, P_PV, Investment, Operational = 0, Maintenance = 0, Fuel = 0, Period = 25, discount_rate = 0.08, inflation_rate = 0.02, replacement_cost = 1250):
    Total_Cost = Investment + replacement_OPEX(degradation, replacement_cost, discount_rate, Period)
    Inflation = 1
    
    for year in range(Period):
        Total_Cost += (Operational[year] + Maintenance[year] + Fuel[year])*Inflation/((1 + discount_rate)**year)
        
        Inflation *= (1 + inflation_rate)
        
    return Total_Cost/sum([sum(i)*0.25 for i in P_PV])

###############################################################################

def plot_LCoE_comparison(LCoE_list, Lifetime_norm, P_PV_registry, P_Grid_min_list, Investment = 4000, Operational = 0, Maintenance = 0, Fuel = 0, Period = 25):
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "Times New Roman",
        'font.size': 16
    })
    
    if not Operational:
        Operational = [0 for i in range(Period)]

    if not Maintenance:
        Maintenance = [0 for i in range(Period)]
        
    if not Fuel:
        Fuel = [0 for i in range(Period)]    

    LCoE_degraded = [LCoE_calc(P_PV, Investment, Operational, round(degradation, 8)) for degradation, P_PV in zip(Lifetime_norm, P_PV_registry)]
    
    plt.figure(constrained_layout=True)
    plt.plot(P_Grid_min_list, LCoE_list, 'b', label='LCoE - base')    
    plt.plot(P_Grid_min_list, LCoE_degraded, 'r', label='LCoE - replacements')        
    plt.grid()
    plt.xlim([int(min(P_Grid_min_list)), int(max(P_Grid_min_list))])
    plt.ylim(bottom=0.05)          
    plt.xlabel('Curtailment threshold, $P_{grid}^{min}$, [kW]')       
    plt.ylabel('Levelized cost of energy, $LCoE_{PV}$, [€/kWh]')  
    plt.legend(loc='upper left')
    plt.show()

    plt.figure(constrained_layout=True)
    plt.plot(P_Grid_min_list, [100*i/j for i,j in zip(LCoE_degraded, LCoE_list)])
    plt.grid()   
    plt.xlim([int(min(P_Grid_min_list)), int(max(P_Grid_min_list))]) 
    plt.xlabel('Curtailment threshold, $P_{grid}^{min}$, [kW]')       
    plt.ylabel('Relative difference between LCoEs, [%]')  
    plt.show()    
    
    return LCoE_degraded
    

###############################################################################

def plot_PowerCurtailment(t, P_Load, P_PV, P_PV_curtailed, P_Grid, start_day = 0, end_day = 365, year = 0, dt = 60*15):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter
    from numpy import arange, ones

    plt.rcParams.update({
        "font.family": "Times New Roman",
        'font.size': 16
    })
    
    P_PV_av = [[i + j for i, j in zip(P_PV[year_plot], P_PV_curtailed[year_plot])] for year_plot in range(len(P_PV))]
    
    plt.figure(constrained_layout=True)
    plt.grid()    
    plt.hist(P_Load, weights=ones(len(P_Load)) / len(P_Load), range=[0, int(max(P_Load))], bins=int(max(P_Load)/0.1), align='mid')    
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlim([0, int(max(P_Load))])
    plt.ylabel('Frequency')
    plt.xlabel('Demanded power, $P_{L}$, [kW]')      
    plt.show()
    
    plt.figure(constrained_layout=True)
    plt.plot(t[int(start_day*24*3600/dt): int(end_day*24*3600/dt)], P_PV_av[year][int(start_day*24*3600/dt): int(end_day*24*3600/dt)], 'b', label='PV power available, $P_{PV}^{av}$')
    plt.plot(t[int(start_day*24*3600/dt): int(end_day*24*3600/dt)], P_PV_curtailed[year][int(start_day*24*3600/dt): int(end_day*24*3600/dt)], 'r', label='PV power curtailed, $P_{PV}$')
    plt.grid()
    plt.xlim([start_day*24, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    plt.ylabel('PV power, $P_{PV}$, [kW]')
    plt.legend(loc='upper right')    
    plt.show()  
    
    
    plt.figure(constrained_layout=True)
    plt.plot(t[int(start_day*24*3600/dt): int(end_day*24*3600/dt)], P_Grid[year], 'b', label='Power consumed from the grid')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    plt.ylabel('Power from the grid, $P_{G}$, [kW]')
    plt.show()
    
    imbalance = [Load - PV for Load, PV in zip(P_Load, P_PV_av[year])]
    plt.figure(constrained_layout=True)
    plt.plot(t[int(start_day*24*3600/dt): int(end_day*24*3600/dt)], imbalance[int(start_day*24*3600/dt): int(end_day*24*3600/dt)], 'b', label='Power imbalance between the load and the PV')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    plt.ylabel('Power imbalance, $P_{i}$, [kW]')
    plt.show()
    
    plt.figure(constrained_layout=True)
    plt.grid()    
    plt.hist(imbalance, weights=ones(len(imbalance)) / len(imbalance), range=[int(min(imbalance))-1, int(max(imbalance))], bins=int(max(imbalance)/0.1), align='mid')    
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlim([int(min(imbalance))-1, int(max(imbalance))])
    plt.ylabel('Frequency')
    plt.xlabel('Power imbalance, $P_{i}$, [kW]')      
    plt.show()    

    E_Grid_net = [sum(year)*dt/3600 for year in P_Grid]
    E_Grid_ret = [-sum([i for i in year if i<0])*dt/3600 for year in P_Grid]
    plt.figure(constrained_layout=True)
    plt.plot(arange(len(P_Grid)), E_Grid_net, 'b', label='Net')
    plt.plot(arange(len(P_Grid)), E_Grid_ret, 'r', label='Returned')
    plt.grid()
    plt.xlim([0, 25])
    plt.xlabel('Time [years]')
    plt.ylabel('Anualized energy exchange with the grid, $E_{G}$, [kWh]')
    plt.legend(loc='lower right')       
    plt.show()
    
    E_curtailed = [sum(year)*dt/3600 for year in P_PV_curtailed]
    plt.figure(constrained_layout=True)
    plt.plot(E_curtailed)    
    plt.grid()
    plt.xlim([0, 25]) 
    plt.xlabel('Time [years]')    
    plt.ylabel('Annualized PV energy curtailed, $E_{PV}^{cur}$, [kWh]')    
    plt.show()

###############################################################################

def DoD_cycle_counting(SoC, return_DoD = False, plot_histogram = True):
    import rainflow
    
    DoD = []
    
    for rng, mean, count, i_start, i_end in rainflow.extract_cycles(SoC):
        if abs(SoC[i_end] - SoC[i_start]) > 0.01:
            DoD.append(int(100*(SoC[i_end] - SoC[i_start])))
        
    if plot_histogram:
        import matplotlib.pyplot as plt
    
        plt.rcParams.update({
            "font.family": "Times New Roman",
            'font.size': 16
        })       
    
        plt.figure(constrained_layout=True)
        plt.hist(DoD, bins = 50)
        plt.xlabel('Depth of discharge, DoD, [%]')    
        plt.ylabel('Frequency')           
        plt.show()
    
    if return_DoD:
        return DoD
        
###############################################################################

def parameters_degradation(SoC, plot_graphs = False):
    import rainflow
        
    range_registry = []
    mean_registry = []
    count_registry = []
    i_start_registry = []
    i_end_registry = []
    DoD_registry = []

    for rng, mean, count, i_start, i_end in rainflow.extract_cycles(SoC):
        if abs(SoC[i_end] - SoC[i_start]) > 0.01:
            range_registry.append(rng)
            mean_registry.append(mean)
            count_registry.append(count)
            i_start_registry.append(i_start)
            i_end_registry.append(i_end)
            DoD_registry.append(int(100*(SoC[i_end] - SoC[i_start])))
            
    if plot_graphs:
        import matplotlib.pyplot as plt
    
        plt.rcParams.update({
            "font.family": "Times New Roman",
            'font.size': 16
        })       
    
        plt.figure(constrained_layout=True)
        plt.plot(range(len(SoC)), SoC, 'c', label = "SoC")
        plt.scatter(i_end_registry, [SoC[i] for i in i_end_registry], color = 'blue', label="end")
        plt.scatter(i_start_registry, [SoC[i] for i in i_start_registry], color = 'red', label="start")
        
###############################################################################        
        
def P_SoC_2d_histogram(single_profile_P_BESS_EoL, single_profile_SoC_EoL, P_Grid_max_list, threshold, P_range = [-2.56, 2.56], SoC_range = [20, 90], vmax_lim = 10):
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "Times New Roman",
        'font.size': 16
    })     
    
    threshold_label = '$P_{grid}^{max}$ = ' + str(P_Grid_max_list[threshold]) + ' kW'
    

    plt.figure(constrained_layout=True)
    plt.hist2d(single_profile_P_BESS_EoL[threshold], [i*100 for i in single_profile_SoC_EoL[threshold]], bins=100, cmap='gist_heat_r', vmax = vmax_lim)
    plt.xlim(P_range)
    plt.ylim(SoC_range)
    plt.ylabel('State-of-charge, SoC, [%]')
    plt.xlabel('BESS power, $P_{BESS}$, [kW]')     
    plt.title(threshold_label)
    plt.show()   

###############################################################################

def Load_histogram(Load, xlabel = 'Demanded power, $P_{L}$, [kW]'):
    from numpy import ones
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter    

    plt.rcParams.update({
        "font.family": "Times New Roman",
        'font.size': 16
    })  
    
    plt.figure(constrained_layout=True)
    
    plt.grid()    
    plt.hist(Load, weights=ones(len(Load)) / len(Load), range=[floor(min(Load)), int(max(Load))], bins=int(max(Load)/0.1), align='mid')    
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlim([floor(min(Load)), int(max(Load))])
    plt.ylabel('Frequency')
    plt.xlabel(xlabel)      
    plt.show()    
    
###############################################################################    
    
def Save_CSV(variable, file_name = ''):
    import csv

    with open(file_name, 'w') as file:
        writer = csv.writer(file)
        writer.writerows(variable)      
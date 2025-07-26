# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a 
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 1.0

###############################################################################
###########################   Imported modules   ##############################

import csvreader
import time
from numpy import arange
from math import floor, ceil
from PeakShaving_Curtailment import *

import warnings
warnings.filterwarnings("ignore")

###############################################################################
#################################   Main code  ################################

start_day = 0
end_day = 365

scale = 1
scale_PV = 2
scale_BESS = 2
    
Replacement_years = []
Replacement_cost = 2500*scale_PV

Operational = [0 for i in range(25)] # replacement_OPEX()

CSVDataTamb = csvreader.read_data(csv='Tamb_15min.csv', address='')
CSVDataP_Load = csvreader.read_data(csv='P_Load_HP_TESS.csv', address='', delim=',') # P_Load_HP, Load_Profile_15min, P_Load_HP_TESS
CSVDataRad = csvreader.read_data(csv='Radiation_1min.csv', address='')
CSVDataTamb.data2array()
CSVDataP_Load.data2array()
CSVDataRad.data2array()
T_amb = [i[0]+273 for i in CSVDataTamb.ar[start_day*24*4:end_day*24*4]]
P_Load = [i*scale for i in CSVDataP_Load.ar[0]]
a = arange(0,len(CSVDataRad.ar),15)
G = [CSVDataRad.ar[i][0]/3600/1000 for i in a]

# PV
# https://www.canadiansolar.com/wp-content/uploads/2019/12/Canadian_Solar-Datasheet-HiDM_CS1U-MS_EN.pdf
n_modules = 5*scale_PV
module_power_ref = 0.315
module_power = 0.400  
n_STC = 0.194
A_module = 2.078*0.992
P_PV = [0]
P_PV_curtailed = [0]
CSVDataPV = csvreader.read_data(csv='PV_15min.csv', address='')
CSVDataPV.data2array()
P_PV_av_0 = [i*n_STC*n_modules*A_module for i in G]
Energy_accumulated = []
Energy_curtailed_accumulated = []
Energy_av_accumulated = []
LCoE_list = []

Investment = 4000*scale_PV


# Grid
P_Grid_0 = 0            # In kW
P_Grid_max_list = arange(0, 10, 0.1).tolist()
P_Grid_min_list = arange(floor(min([i - j for i,j in zip(P_Load, P_PV_av_0)])), 0.1, 0.1).tolist()
P_grid_accumulated = []
E_grid_accumulated = []


# Initial conditions
dt = 60*15                          # In s
t = [start_day*24*3600/dt]          # In s
t_final = int((end_day - start_day)*24*3600/dt)         # In s


mode = 'Curtailment' # modes = ['PeakShaving', 'Curtailment']
simulation_years = 26
year_plot = 0 
threshold_plot = 25

start = time.time()    # The timer is initializad.

if mode == 'PeakShaving':
    
    # Registries
    P_PV_registry = []
    SoC_BESS_registry = []
    E_BESS_registry = []
    P_BESS_registry = []
    P_Grid_registry = []
    Capacity_BESS_registry = []
    T_TESS_registry = []
    Qdot_TESS_registry = []
    Qdot_HP_registry = []
    Qdot_HP_TESS_registry = []    

    # BESS
    replace_BESS = False
    P_BESS_max = 1.28*2*scale_BESS            # In kW
    Capacity_BESS_BOL = 3.36*3*scale_BESS              # in kWh
    Capacity_BESS_0 = Capacity_BESS_BOL     # in kWh
    CAPEX_BESS = 6500*scale_BESS
    EoL = []
    Stored_Energy = []              # Energy required to be stored
    Stored_Energy_EoL = []          # Energy the BESS can store before its EoL
    
    
    for P_Grid_max in P_Grid_max_list:
        [t_year, P_PV, SoC_BESS, E_BESS, P_BESS, P_Grid, Capacity_BESS, T_TESS, Qdot_TESS, Qdot_HP, Qdot_HP_TESS] = PeakShaving_Projection(P_Load, P_PV_av_0, P_BESS_max, P_Grid_max, Capacity_BESS_0, Capacity_BESS_BOL, replace_BESS = replace_BESS, simulation_years = simulation_years)
        
        EoL.append(EoL_BESS(SoH_estimation(Capacity_BESS)))
        P_grid_accumulated.append(P_Grid)
        E_grid_accumulated.append(sum([sum(year) for year in P_Grid])*0.25)
        
        Stored_Energy.append(-sum(sum([i for i in year if i<0]) for year in P_BESS))
        Stored_Energy_EoL.append(Energy_EoL(P_BESS, EoL[-1]))

        P_PV_registry.append(P_PV)
        SoC_BESS_registry.append(SoC_BESS)
        E_BESS_registry.append(E_BESS)
        P_BESS_registry.append(P_BESS)
        P_Grid_registry.append(P_Grid)
        Capacity_BESS_registry.append(Capacity_BESS)
        T_TESS_registry.append(T_TESS)
        Qdot_TESS_registry.append(Qdot_TESS)
        Qdot_HP_registry.append(Qdot_HP)
        Qdot_HP_TESS_registry.append(Qdot_HP_TESS)
    
    replacements = [ceil(25/years) for years in EoL]
    LCoS_list = [LCoS_calc(i, CAPEX_BESS) for i in Stored_Energy_EoL]
    
    

    single_profile_SoC_EoL = [Energy_EoL(S,E, accumulated = False) for S,E in zip(SoC_BESS_registry, EoL)]
    single_profile_P_BESS_EoL = [Energy_EoL(P,E, accumulated = False) for P,E in zip(P_BESS_registry, EoL)]
    single_profile_P_Grid_EoL = [Energy_EoL(P,E, accumulated = False) for P,E in zip(P_Grid_registry, EoL)]
    
    single_profile_DoD_EoL = [DoD_cycle_counting(SoC, return_DoD = True, plot_histogram = False) for SoC in single_profile_SoC_EoL]
    
    
    plot_BESS_PeakShaving(t_year, P_Load, P_PV_registry[threshold_plot], SoC_BESS_registry[threshold_plot], E_BESS_registry[threshold_plot], P_BESS_registry[threshold_plot], P_Grid_registry[threshold_plot], Capacity_BESS_registry[threshold_plot], SoH_estimation(Capacity_BESS_registry[threshold_plot]), start_day, end_day, year_plot, dt)
    PeakShaving_limit_violations(single_profile_P_Grid_EoL, P_Grid_max_list, single_profile_P_BESS_EoL)
    
    
elif mode == 'Curtailment':

    # Registries
    P_PV_registry = []
    P_PV_curtailed_registry = []
    P_PV_av_registry = []
    P_Grid_registry = []
 
    
    for P_Grid_min in P_Grid_min_list:

        [t_year, P_PV, P_PV_curtailed, P_PV_av, P_Grid] = Curtailment_Projection(P_Load, P_PV_av_0, P_Grid_0, P_Grid_min, simulation_years = simulation_years)

        P_grid_accumulated.append(P_Grid)
        E_grid_accumulated.append(sum([sum(year) for year in P_Grid])*0.25)        
        Energy_av_accumulated.append(sum([sum(i)*0.25 for i in P_PV_av]))
        Energy_accumulated.append(sum([sum(i)*0.25 for i in P_PV]))
        Energy_curtailed_accumulated.append(sum([sum(i)*0.25 for i in P_PV_curtailed]))
        LCoE_list.append(LCoE_calc(P_PV, Investment, Operational))
        
        P_PV_registry.append(P_PV)
        P_PV_curtailed_registry.append(P_PV_curtailed)
        P_PV_av_registry.append(P_PV_av)
        P_Grid_registry.append(P_Grid)      
    
    V_PV_registry = Yearly_Voltage(P_PV_registry, G, T_amb, n_modules = n_modules, n_series = scale_PV)
    P_PV_registry_sampled = Yearly_Curtailed_Power(P_PV_registry)
    [i_av_switch_registry, i_RMS_switch_registry, i_av_diode_registry, i_RMS_diode_registry] = get_currents(V_PV_registry, P_PV_registry_sampled, scale_PV)
    Conduction_losses_IGBT_registry = [[conduction_loss_IGBT(i_av, i_RMS) for i_av, i_RMS in zip(threshold_i_av, threshold_i_RMS)] for threshold_i_av, threshold_i_RMS in zip(i_av_switch_registry, i_RMS_switch_registry)]
    Switching_losses_IGBT_registry = [[switching_losses_IGBT(i_av, i_RMS) for i_av, i_RMS in zip(threshold_i_av, threshold_i_RMS)] for threshold_i_av, threshold_i_RMS in zip(i_av_switch_registry, i_RMS_switch_registry)]
    IGBT_losses_Registry = [[Conduction_losses + Switching_losses for Conduction_losses, Switching_losses in zip(Conduction_threshold, Switching_threshold)] for Conduction_threshold, Switching_threshold in zip(Conduction_losses_IGBT_registry, Switching_losses_IGBT_registry)]

    T_IGBT_registry = get_switch_temperature_IGBT(T_amb, i_av_switch_registry, i_RMS_switch_registry, R_th = [1.7, 8+1.33])
    
    L = get_lifetime(T_IGBT_registry, P_Grid_min_list)
    
    LCoE_degraded = plot_LCoE_comparison(LCoE_list, [i/max(L) for i in L], P_PV_registry, P_Grid_min_list, Investment)
    
    plot_I_T_day(i_av_switch_registry, T_IGBT_registry, P_Grid_min_list, day = 15)
    plot_I_T_day(i_av_switch_registry, T_IGBT_registry, P_Grid_min_list, day = 172)

end = time.time()    # The timer is initializad.
totalelapsed = end - start  # The total time is calculated.


###############################################################################
#################################   Figures  ##################################

import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "Times New Roman",
    'font.size': 14
})


if mode == 'PeakShaving':
    plt.figure(constrained_layout=True)
    plt.plot(P_Grid_max_list, EoL, 'b', label='End of life')
    plt.grid()
    plt.xlim([int(min(P_Grid_max_list)), ceil(max(P_Grid_max_list))])
    plt.ylim(bottom=0)        
    plt.xlabel('Peak-shaving threshold, $P_{grid}^{max}$, [kW]')  
    plt.ylabel('End of life, EoL, [years]')   
    plt.show()    

    plt.figure(constrained_layout=True)
    plt.plot(P_Grid_max_list, Stored_Energy_EoL, 'b', label='Accumulated energy stored before EoL')
    plt.grid()
    plt.xlim([int(min(P_Grid_max_list)), ceil(max(P_Grid_max_list))])
    plt.ylim(bottom=0)        
    plt.xlabel('Peak-shaving threshold, $P_{grid}^{max}$, [kW]')  
    plt.ylabel('Accumulated energy stored before EoL, $E_{BESS}^{EoL}$, [kWh]')   
    plt.show() 

# It increases because as the energy cannot be returned to the grid, so the net possitive energy (PV to grid) is reduced.   
    plt.figure(constrained_layout=True)
    plt.plot(P_Grid_max_list, [i/1000  for i in E_grid_accumulated], 'b', label='Accumulated Energy')    
    plt.grid()
    plt.xlim([int(min(P_Grid_max_list)), ceil(max(P_Grid_max_list))])  
    plt.xlabel('Peak-shaving threshold, $P_{grid}^{max}$, [kW]')  
    plt.ylabel('Accumulated energy from the grid, $E_{grid}$, [MWh]')    
    plt.show()      

    plt.figure(constrained_layout=True)
    plt.plot(P_Grid_max_list, replacements, 'b', label='Replacements')    
    plt.grid()
    plt.xlim([int(min(P_Grid_max_list)), ceil(max(P_Grid_max_list))])
    plt.ylim(bottom=0)     
    plt.xlabel('Peak-shaving threshold, $P_{grid}^{max}$, [kW]')  
    plt.ylabel('Required replacements, $n_{BESS}$, [-]')    
    plt.show()   
    
    plt.figure(constrained_layout=True)
    plt.plot(P_Grid_max_list, LCoS_list, 'b', label='LCoS')    
    plt.grid()
    plt.xlim([int(min(P_Grid_max_list)), ceil(max(P_Grid_max_list))])
    plt.ylim(bottom=0.05)        
    plt.xlabel('Peak-shaving threshold, $P_{grid}^{max}$, [kW]')  
    plt.ylabel('Levelized cost of storage, $LCoS_{BESS}$, [€/kWh]')    
    plt.show()  



    plt.figure(constrained_layout=True)
    plt.plot(P_Grid_max_list, [i/max(EoL) for i in EoL], 'b', label='Normalized end of life')    
    plt.grid()
    plt.xlim([int(min(P_Grid_max_list)), ceil(max(P_Grid_max_list))])  
    plt.xlabel('Peak-shaving threshold, $P_{grid}^{max}$, [kW]')  
    plt.ylabel('Normalized EoL, $EoL_{norm}$')    
    plt.show()

    for i in range(10):
        P_SoC_2d_histogram(single_profile_P_BESS_EoL, single_profile_SoC_EoL, P_Grid_max_list, 10*i, P_range = [-6, 6])
      
    
elif mode == 'Curtailment':
    
    plt.figure(constrained_layout=True)
    plt.plot(P_Grid_min_list, [i/1000  for i in Energy_curtailed_accumulated], 'b', label='Total generation curtailed')    
    plt.plot(P_Grid_min_list, [i/1000  for i in Energy_accumulated], 'r', label='Total generation used')    
    plt.grid()
    plt.xlim([int(min(P_Grid_min_list)), int(max(P_Grid_min_list))])
    plt.ylim(bottom=0)           
    plt.xlabel('Curtailment threshold, $P_{grid}^{min}$, [kW]')    
    plt.ylabel('Accumulated PV energy, $E_{PV}$ [MWh]')    
    plt.legend(loc='center left')    
    plt.show()
    
    plt.figure(constrained_layout=True)
    plt.plot(P_Grid_min_list, LCoE_list, 'b', label='LCoE')    
    plt.grid()
    plt.xlim([int(min(P_Grid_min_list)), int(max(P_Grid_min_list))])
    plt.ylim(bottom=0.05)           
    plt.xlabel('Curtailment threshold, $P_{grid}^{min}$, [kW]')       
    plt.ylabel('Levelized cost of energy, $LCoE_{PV}$, [€/kWh]')    
    plt.show()


# It increases because as the energy cannot be returned to the grid, so the net possitive energy (PV to grid) is reduced.
    plt.figure(constrained_layout=True)
    plt.plot(P_Grid_min_list, [i/1000 for i in E_grid_accumulated], 'b', label='Energy from the Grid, $E_{grid}$')    
    plt.plot(P_Grid_min_list, [i/1000  for i in Energy_accumulated], 'r', label='PV generation used, $E_{PV}$')       
    plt.grid()
    plt.xlim([int(min(P_Grid_min_list)), int(max(P_Grid_min_list))])
    plt.ylim(bottom=0)          
    plt.xlabel('Curtailment threshold, $P_{grid}^{min}$, [MWh]')       
    plt.ylabel('Accumulated energy, [MWh]')    
    plt.legend(loc='lower left')     
    plt.show()  
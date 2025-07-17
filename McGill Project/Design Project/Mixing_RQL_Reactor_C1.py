import cantera as ct
import csv
from numpy import *
import pandas as pd
from matplotlib import pyplot as plt
from Blend_Mass_Fractions import BlendingMassFractions

def mixing(gas_a,gas_b,phi):
    res_a = ct.Reservoir(gas_a)
    res_b = ct.Reservoir(gas_b)
    downstream = ct.Reservoir(gas_b)

    mixer = ct.IdealGasReactor(gas_b)
    M_air = gas_a.density
    M_fuel = gas_b.density
    print('Air density is',M_air)
    print('Fuel density is',M_fuel)
    mfc1 = ct.MassFlowController(res_a, mixer, mdot=69/phi)
    mfc2 = ct.MassFlowController(res_b, mixer, mdot=118)

    outlet = ct.Valve(mixer, downstream, K=10.0)
    sim = ct.ReactorNet([mixer])
    sim.advance_to_steady_state()

    # view the state of the gas in the mixer
    return mixer

# Solve the first flame
def RQL(phi_first,phi_last,blending_ratio):
    Tin = 800
    P = 30*ct.one_atm
    grid_initial = 0.03*linspace(0,1,20)	# Flame length
    gas_init = ct.Solution('NOxC1.yaml')
    bl_ratio = blending_ratio
    phi = phi_first
    A,Y,Z,H = BlendingMassFractions(bl_ratio)
    reactants = {'POSF11498':A,'H2':H,'O2':Y/phi, 'N2':Z/phi}
    gas_init.TPY = Tin, P, reactants
    f1 = ct.FreeFlame(gas_init, grid_initial)
    f1.energy_enabled = True
    tol_ss = [1.0e-5, 1.0e-9] # [rtol atol] for steady-state problem
    tol_ts = [1.0e-5, 1.0e-9] # [rtol atol] for time stepping
    f1.flame.set_steady_tolerances(default=tol_ss)
    f1.flame.set_transient_tolerances(default=tol_ts)
    f1.set_max_jac_age(50, 50)
    f1.set_refine_criteria(ratio=5, slope=0.5,curve=0.5)
    f1.solve(loglevel=1, refine_grid=True, auto=False) 
    gas_c1_boundary = ''
    for i in range(0,gas_init.n_species):
        if i == 0:
            gas_c1_boundary += '%s:%.9f'%(gas_init.species_name(i),f1.Y[gas_init.species_index(gas_init.species_name(i))][-1])
        else:
            gas_c1_boundary += ',%s:%.9f'%(gas_init.species_name(i),f1.Y[gas_init.species_index(gas_init.species_name(i))][-1])
    #with open('Flame1_C1.in', 'w') as f:
    #    f.write(gas_c1_boundary)

    # Mix the gases
    gas_a = ct.Solution('NOxC1.yaml') # Gas 'a' is the air, right after the HPC
    gas_a.TPY = 800.0, 30*ct.one_atm, 'O2:0.21, N2:0.79'
    # Use conditions from Combustor 1 to initialize the gas to be mixed with HPC air
    gas_b = ct.Solution('NOxC1.yaml')
    gas_b.TP = f1.T[-1], 30*ct.one_atm
    gas_b.Y = gas_c1_boundary # Get the gas composition space from excel data of the first
    gas = mixing(gas_a,gas_b,phi_last)
    print('Mixed gas temperature is',gas.thermo.T)
    # Prepare the species from this reactor to the one with Lean burner calculation
    boundary_species_lb = ''
    for i in range(100000):
        try:
            name = gas.thermo.species_name(i)
            y_val = gas.thermo.Y[gas.thermo.species_index(name)]
            if i == 0:
                boundary_species_lb += '%s:%.5f'%(name,y_val)
            else:
                boundary_species_lb += ' ,%s:%.5f'%(name,y_val)
        except:
            break

    # Calculate combustor 2 in Well stirred reactor
    gas_f = ct.Solution('NOxC1.yaml')
    gas_f.TP = gas.thermo.T,gas.thermo.P
    gas_f.Y = boundary_species_lb
    inlet = ct.Reservoir(gas_f)
    gas_f.equilibrate('HP')
    combustor = ct.IdealGasReactor(gas_f)
    combustor.volume = 1.0
    exhaust = ct.Reservoir(gas_f)
    inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=1.0)
    outlet_mfc = ct.PressureController(combustor, exhaust, primary=inlet_mfc, K=0.01)
    sim = ct.ReactorNet([combustor])
    sim.advance_to_steady_state()
    Temp = combustor.T
    NO = 1000000*combustor.thermo.X[gas.thermo.species_index('NO')]
    CO = 1000000*combustor.thermo.X[gas.thermo.species_index('CO')]
    return Temp,NO,CO



phi_last = [0.2,0.3,0.4,0.5,0.6,0.7]
phi_first = [1.2,1.3,1.4,1.5,1.6,1.7,1.8]
for bl in [0,1,2,3]:
    Temp = []
    NOXs = []
    COs = []
    phi_v = []
    for phif in phi_first:
        T,NOx,CO = RQL(phif,0.4,bl)
        phi_v.append(phif)
        Temp.append(T)
        NOXs.append(NOx)
        COs.append(CO)
    data = {'phi_f':phi_v,'Temp':Temp,'NOx':NOXs,'COs':COs}
    df = pd.DataFrame(data)
    df.to_csv('RQL_%d_phif.csv'%(bl))

phi_last = [0.2,0.3,0.4,0.5,0.6,0.7,0.8]
for bl in [0,1,2,3]:
    Temp = []
    NOXs = []
    COs = []
    phi_v = []
    for phil in phi_last:
        T,NOx,CO = RQL(1.5,phil,bl)
        phi_v.append(phil)
        Temp.append(T)
        NOXs.append(NOx)
        COs.append(CO)
    data = {'phi_l':phi_v,'Temp':Temp,'NOx':NOXs,'COs':COs}
    df = pd.DataFrame(data)
    df.to_csv('RQL_%d_phil.csv'%(bl))

# %%
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct


# %%
gas = ct.Solution('Jerzembeck.yaml')

# %%
equiv_ratio = 0.39  # lean combustion
gas.TP = 800.0, ct.one_atm*30 
gas.set_equivalence_ratio(equiv_ratio, 'IXC8H18:1.0', 'O2:1.0, N2:3.76')
inlet = ct.Reservoir(gas)
gas.equilibrate('HP')
combustor = ct.IdealGasReactor(gas)
combustor.volume = 0.40
exhaust = ct.Reservoir(gas)

# %%
residence_time = 0.01
def mdot(t):
    return combustor.mass / residence_time


# %%
inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=mdot)
outlet_mfc = ct.PressureController(combustor, exhaust, primary=inlet_mfc, K=0.01)


# %%
sim = ct.ReactorNet([combustor])
states = ct.SolutionArray(gas, extra=['tres'])
gas.equilibrate("HP")
T_PSR = gas.T
print(T_PSR)

# %%
species = {S.name : S for S in ct.Species.list_from_file("Jerzembeck.yaml")}
species = ['CO','NO']
conc = []

for spec in species :
    conc.append(gas[spec].X[0])
print(conc)

# %%
while combustor.T > 800:
    sim.initial_time = 0.0  # reset the integrator
    sim.advance_to_steady_state()
    print('tres = {:.2e}; T = {:.1f}'.format(residence_time, combustor.T))
    states.append(combustor.thermo.state, tres=residence_time)
    residence_time *= 0.09  # decrease the residence time for the next iteration


# %%
#net = ct.ReactorNet([combustor])
#net.max_time_step = residence_time
#net.initialize()
#end_time = 5.0*residence_time
#time= []
#T=[]
#mdot=[]
#while net.time <= end_time:
    #time.append(net.time)
    #T.append(combustor.T)
    #mdot.append(inlet.mass_flow_rate)
    #net.step()
#print(time)


# %%
fig.ax = plt.subplots()
ax.plot(time, T, label='Temperature')
ax.plot(time,np.array(mdot)/2, label= 'Mass flow rate/2' )
ax.legend(loc='best')



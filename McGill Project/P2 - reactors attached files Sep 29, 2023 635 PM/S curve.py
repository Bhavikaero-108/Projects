# %%
import pandas as pd
import numpy as np
import math
import time
import cantera as ct
print('Runnning Cantera version: ' + ct.__version__)
import matplotlib.pyplot as plt

# %%
pressures = np.logspace (-3,4,num=20)
temperatures = np.linspace(650,900,num=31)
estimated_ignition_delay_time= 10 

# %%
for press in pressures :
    for temps in temperatures:

        time_bank = np.array([])
        temperature_history = np.array([])
        gas = ct.Solution('gri30.yaml')
        gas.TP = temps,press
        gas.set_equivalence_ratio(phi=1.0, fuel="H2", oxidizer={"o2": 1.0, "n2": 3.76})

        r= ct.IdealGasReactor(contents=gas, name="Batch Reactor")
        reactor_network = ct.ReactorNet([r])

        time_history = ct.SolutionArray(gas, extra="t")
        
        t= 0
        while t < estimated_ignition_delay_time:
            t= reactor_network.step()
            time_bank= np.append(time_bank,t)
            temperature_history = np.append(temperature_history, r.thermo.T)

        if np.max (temperature_history)-temps == 500 > 0 :
            ignition_temp = np.append(ignition_temp,temps)
            print(ignition_temp)
            break

        



# %%
plt.figure()
plt.semilogy(ignition_temp,pressures/ct.one_atm)
plt.xlabel("Ignition Temp")
plt.ylabel("Pressure")
plt.show ()



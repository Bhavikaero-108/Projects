import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option("display.precision", 4)

def BlendingMassFractions(i):
    A = np.array([0.06088,0.04914,0.0396,0.02925,0.02501,0.01932,0.0144,0.0101,0.0063,0.003,0]) # C1 fuel
    Y = np.array([0.2189,0.22,0.22156,0.2225,0.2234,0.224,0.2247,0.2252,0.2257,0.2261,0.2265]) # O2
    Z = np.array([0.72022,0.7254,0.72895,0.73218,0.73492,0.73736,0.7393,0.7410,0.7426,0.7440,0.7452]) # N2
    H = np.array([0,0.00546,0.00989,0.01357,0.01667,0.01932,0.0216,0.0236,0.0254,0.0269,0.0283]) # H2
    return A[i],Y[i],Z[i],H[i]

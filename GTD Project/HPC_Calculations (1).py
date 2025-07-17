import numpy as np
import matplotlib.pyplot as plt 
import random
from scipy.optimize import fsolve

#from sympy import symbols, Eq, solve
                                                                                           ##### INLET #####
# Given Conditions 
m_dot = 4.082    
Po1 = 201335
To1 = 416.4124
To2 = 597.056
To5 = To2
Po5 = 608243.5
r1_h_turbine = 0.097429
R_idealgas = 287
gamma = 1.4
cp_a = 1.005 
HPC_eta = 0.85
eta_i = 0.89
M5_abs = 0.1

# Inputs 
N_rpm = 29400 
r1_h = 0.722 * r1_h_turbine                          # Hub Radius at Inlet
Alpha_1 =   (10*np.pi)/180                           # 10 degrees in radians
r1_sh_values = np.arange(0.09301, 0.140010, 0.001)   # Range of r_sh values to try in m
i1_h = (7*np.pi)/180                                 ###
i1_m = (3.5*np.pi)/180                               # Incidence Angle at inlet h, m, sh converted to radians selected from the ranges given
i1_sh = (2*np.pi)/180                                ###
Alpha2r_star = -25                                   # Guessed metal angle (For corrected tip speed part)
sigma_ax = 0.000254                                  # tip clearance in m 
N_b = 32                                             # Blade count
B2_star = 0.19                                       # B2*. Range between 0.15 and 0.23

m_dot_3 = m_dot_2 = m_dot
   
## Shroud Inlet Calculations
# Initialize minimum relative Mach number
min_M1_rel_sh = []
rho_static_list = []
U1_sh_list = []
C1_w_list = []
V1_w_sh_list = []
V1_sh_list = []
Alpha1r_sh_list = []
C1_x_list = []
r1_m_list = []
U1_m_list = []

# Iterate over values of r1_sh
for r1_sh in r1_sh_values:

    # Calculate initial conditions
    Rho_tot_1 = (Po1/(R_idealgas*To1))
    Area_1 = np.pi*((r1_sh**2) - (r1_h**2))
    
    # First Iteration (with total density)
    C1_x = m_dot/(Rho_tot_1*Area_1)                                     
    C1 = C1_x/(np.cos(Alpha_1))                             
    M1_abs = np.sqrt((C1**2)/(To1*gamma*R_idealgas-(((gamma-1)/2)*(C1**2)))) 
    a = C1/M1_abs                                                               
    
    # Iterations for Static density at inlet
    for i in range(100):
        rho_static = Rho_tot_1/((1+((gamma-1)/2)*(M1_abs**2))**(1/(gamma-1)))
        rho_static_new = 0
        C1_x = m_dot/(rho_static*Area_1)
        C1 = C1_x/(np.cos(Alpha_1))
        M1_abs = np.sqrt((C1**2)/(To1*gamma*R_idealgas-(((gamma-1)/2)*(C1**2))))
        a = C1/M1_abs
        rho_static_new = rho_static
        if rho_static - rho_static_new == 0.0000001:
            break

    # Inlet shroud calculations
    U1_sh = 2*np.pi*r1_sh*(N_rpm/60)     
    C1_w = C1_x*np.tan(Alpha_1)  
    V1_w_sh = U1_sh - C1_w
    V1_sh = np.sqrt((C1_x**2) + (V1_w_sh**2))
    Alpha1r_sh = np.arctan(V1_w_sh/C1_x)                                     
    M1_rel_sh = V1_sh/a                                 
    min_M1_rel_sh.append(M1_rel_sh)
    rho_static_list.append(rho_static)
    U1_sh_list.append(U1_sh)
    C1_w_list.append(C1_w)
    V1_w_sh_list.append(V1_w_sh)
    V1_sh_list.append(V1_sh)
    Alpha1r_sh_list.append(Alpha1r_sh)
    C1_x_list.append(C1_x)
    r1_m = (1/2)*(r1_h + r1_sh)
    U1_m = 2*np.pi*r1_m*(N_rpm/60)  
    r1_m_list.append(r1_m)
    U1_m_list.append(U1_m)


# Printing corresponding values with minimum relative Mach Number at inlet Shroud
min_index = min_M1_rel_sh.index(min(min_M1_rel_sh))
print("----------------Impeller Inlet Shroud Conditions-----------------")
print("Static Density at Inlet:", rho_static_list[min_index], 'kg/m^3' )
print("U1_sh:", U1_sh_list[min_index], 'm/s')
print("V1_w_sh:", V1_w_sh_list[min_index], 'm/s')
print("V1_sh:", V1_sh_list[min_index], 'm/s')
print("Alpha1r_sh:", (180*Alpha1r_sh_list[min_index])/(np.pi), 'degrees')
print("r1_sh:", r1_sh_values[min_index], 'm')
print('r1_h:', r1_h, 'm')
print("C1_w:", C1_w_list[min_index], 'm/s')
print("Minimum M1_rel_sh:", min(min_M1_rel_sh))


# Assigning the corresponding values to the parameters for future calculations
M1_rel_sh = min(min_M1_rel_sh) 
rho_static = rho_static_list[min_index]
U1_sh = U1_sh_list[min_index]
V1_w_sh = V1_w_sh_list[min_index]
V1_sh = V1_sh_list[min_index]
Alpha1r_sh = Alpha1r_sh_list[min_index]
r1_sh = r1_sh_values[min_index] 
C1_w = C1_w_list[min_index]



Area_1 = np.pi*((r1_sh**2) - (r1_h**2))
C1_x = m_dot/(rho_static*Area_1)
C1 = C1_x/(np.cos(Alpha_1))
M1_abs = M1_abs = np.sqrt((C1**2)/(To1*gamma*R_idealgas-(((gamma-1)/2)*(C1**2))))
print('C1:', C1, 'm/s')
print('M1:', M1_abs)

print('Cx1:', C1_x)
print('inlet Area:', Area_1, 'm^2')

C1_x_U1_m = np.array(C1_x_list)/np.array(U1_m_list)
min_M1_rel_sh_resized = np.resize(min_M1_rel_sh, C1_x_U1_m.shape)

##### Plotting Mach vs Cx/Um
plt.plot(C1_x_U1_m, min_M1_rel_sh)
plt.xlabel('Cx1/U1m')
plt.ylabel('Relative Mach Number')
plt.title('Change in Inlet Relative Mach No. at Shroud vs Cx1/U1m')
plt.show()

## Mean Inlet Calculations
r1_m = r1_m_list[min_index]
U1_m = U1_m_list[min_index]   
V1_w_m = U1_m - C1_w
V1_m = np.sqrt((C1_x**2) + (V1_w_m**2))
Alpha1r_m = np.arctan(V1_w_m/C1_x)                                     
M1_rel_m = V1_m/a  

print("----------------Impeller Inlet Mean Conditions-----------------")
print("r1_m:", r1_m, 'm')
print("U1_m:", U1_m, 'm/s')
print("V1_w_m:", V1_w_m, 'm/s')
print("V1_m:", V1_m, 'm/s')
print("Alpha1r_m:", (180*Alpha1r_m)/(np.pi), 'degrees')
print("M1_rel_m:", M1_rel_m)

##Hub Inlet Calculations
U1_h = 2*np.pi*r1_h*(N_rpm/60)     
V1_w_h = U1_h - C1_w
V1_h = np.sqrt((C1_x**2) + (V1_w_h**2))
Alpha1r_h = np.arctan(V1_w_h/C1_x)                                     
M1_rel_h = V1_h/a 

print("----------------Impeller Inlet Hub Conditions-----------------")
print("r1_h:", r1_h, 'm')
print("U1_h:", U1_h, 'm/s')
print("V1_w_h:", V1_w_h, 'm/s')
print("V1_h:", V1_h, 'm/s')
print("Alpha1r_h:", (180*Alpha1r_h)/(np.pi), 'degrees')
print("M1_rel_h:", M1_rel_h)

T1 = (To1)/(1+(((gamma-1)/2)*M1_abs**2))
P1 = Po1 * ((T1/To1)**(gamma/(gamma-1)))
Q_rate = (m_dot/rho_static) * 35.3147               #(ft/s)
ho_is = (HPC_eta * cp_a * (To5 - To1)) * 0.429923   #(Btu/lb)
N_s = (N_rpm * np.sqrt(Q_rate))/(778.26 * ho_is)**(3/4)

Tr_1_sh = T1 + ((V1_sh**2)/(2*(cp_a*10**3)))
Pr_1_sh = P1 * (Tr_1_sh/T1)**(gamma/(gamma-1))
Tr_1_m = T1 + ((V1_m**2)/(2*(cp_a*10**3)))
Pr_1_m= P1 * (Tr_1_m/T1)**(gamma/(gamma-1))
Tr_1_h= T1 + ((V1_h**2)/(2*(cp_a*10**3)))
Pr_1_h = P1 * (Tr_1_h/T1)**(gamma/(gamma-1))

print("----------------Other Impeller Inlet Conditions-----------------")
print("Tr1sh:", Tr_1_sh, 'K')
print("Pr1sh:", Pr_1_sh, 'Pa')
print("Tr1m:", Tr_1_m, 'K')
print("Pr1m:", Pr_1_m, 'Pa')
print("Tr1h:", Tr_1_h, 'K')
print("Pr1h:", Pr_1_h, 'Pa')
print('T1 =', T1, 'K')
print('P1 =', P1, 'Pa')
print('Volumetric Flow Rate =', Q_rate, 'ft/s')
print('ho_is =', ho_is)
print ('Centrifugal Specific Speed =', N_s)


### Incidence Angle i in radians
Alpha1_star_h = Alpha1r_h - i1_h
Alpha1_star_m = Alpha1r_m - i1_m
Alpha1_star_sh = Alpha1r_sh - i1_sh

print('Alpha1r Star Hub:', (180*Alpha1_star_h)/(np.pi), 'deg')
print('Alpha1r Star Mean:', (180*Alpha1_star_m)/(np.pi), 'deg')
print('Alpha1r Star Shroud:', (180*Alpha1_star_sh)/(np.pi), 'deg')
################################################# Impeller Corrected tip speed  (Graph)
Theta_1 = (1.8*To1)/518.7
P_ratio = 3                  
if -30 <= Alpha2r_star <= 0:
    U2_corr1 = ((Alpha2r_star*(1348.8295-1226.8381))/(-30)) + 1226.8381
else:
    U2_corr1 = (((Alpha2r_star-(-45))*(1348.8295-1422.189))/(-30+45)) + 1422.189
                
U2_corr2 = (((180*Alpha_1)/(np.pi))*(1241.1689-1190.26975))/50 + (1190.26975)

U2_corr = U2_corr1*(U2_corr2/1190.26975)   #(ft/s)

print("---------------- Impeller Exit Conditions-----------------")
print("U2 corrected:", U2_corr, "ft/s")

################################################## Exit Geometry


U_2 = (U2_corr*np.sqrt(Theta_1))*0.3048   #(m/s)
C2_w = ((cp_a*1000)*(To2 - To1) + (U1_m*C1_w))/U_2
C2_w_inf = C2_w/(1-((0.63*np.pi)/N_b))
V2_w_inf = -(U_2 - C2_w_inf)
V2_inf = abs(V2_w_inf/(np.sin((Alpha2r_star*np.pi)/180)))

print("U_2:", U_2, 'm/s')
print("C2_w:", C2_w, 'm/s')
print("C2_w_inf:", C2_w_inf, 'm/s')
print("V2_w_inf:", V2_w_inf, 'm/s')
print("V2_inf:", V2_inf, 'm/s')
print('Slip Factor:', (1-((0.63*np.pi)/N_b)))


V_2 = 0.54 * V1_sh

# Calculate V2_w, Alpha2r, and C2_r
V2_w = -(U_2 - C2_w)
Alpha2r = np.arcsin(V2_w/V_2)
C2_r = V2_w/(np.tan(Alpha2r))
print('C2r:', C2_r, 'm/s')
print('V2w:', V2_w, 'm/s')
print('V2:', V_2, 'm/s')
#print('Cr2/Cx1: ', C2_r/C1_x)


max_iterations = 500
counter = 0
print('C2r/C1x:', C2_r/C1_x)
print('Alpha range: less or equal to:', (-40*np.pi)/180)
print('Alpha 2r: ', (180*Alpha2r)/(np.pi), 'degrees')


C2 = np.sqrt((C2_r**2)+(C2_w**2))
Alpha2 = np.arctan(C2_w/C2_r)
M2_abs = np.sqrt((C2**2)/(To2*gamma*R_idealgas-(((gamma-1)/2)*(C2**2))))
Po2 = Po1*(( 1 + (( eta_i * ( To2 - To1 )) / (To1)))**( gamma / (gamma-1)))
T2 = (C2**2)/(gamma*R_idealgas*(M2_abs**2))
P2 = Po2 * ((T2/To2)**(gamma/(gamma-1)))
Rho_static_2 = P2/(R_idealgas*T2)
M_rel_2 = V_2/(np.sqrt(gamma*R_idealgas*T2))
r_2 = (30*U_2)/(N_rpm*np.pi)
delta = Alpha2r - ((np.pi*Alpha2r_star)/180)
Tr_2 = T2 + ((V_2**2)/(2*(cp_a*10**3)))
Pr_2 = P2 * (Tr_2/T2)**(gamma/(gamma-1))


print("Tr2:", Tr_2, 'K')
print("Pr2:", Pr_2, 'Pa')
print('C2:', C2, 'm/s')
print('Alpha2:', (180*Alpha2)/(np.pi), 'degrees')
print('M2_abs:', M2_abs)
print('Po2:', Po2, 'Pa')
print('T2:', T2, 'K')
print('P2:', P2, 'Pa')
print('Static Density at Exit:', Rho_static_2, 'kg/m^3')
print('M2_rel: ',M_rel_2)
print('r2:', r_2, 'm')
print('Delta:',(180*delta)/(np.pi), 'degrees')


## Axial Length
Length_imp = (r_2 - r1_h)/1.05    # Range is specified

print('Length of Impeller:', Length_imp, 'm')

# Channel Height at Impeller Exit b_2
b2 = m_dot_2/(2 * np.pi * Rho_static_2 * r_2 * C2_r * (1-B2_star))
h_2 = b2 - sigma_ax                                                # Channel Height
Area_2 = 2 * np.pi * r_2 * b2
print('Area 2:', Area_2, 'm^2')
print('Channel Height:', b2, 'm')
print('Blade Height:', h_2, 'm')

### Channel Height at Vaneless Diffuser Check (Graph):

print('b2/r2:', b2/r_2)
print('Alpha2:', Alpha2)


### Radius Ratio at Vaneless Diffuser and Diffuser Loss: 
r_3 = 1.075 * r_2    # between 1.05 and 1.1
Po3 = Po2 * (1-0.01)
print("----------------Vaned Diffuser Leading Edge-----------------")
print('Po3:', Po3, 'Pa')

### Flow properties at Vaned Diffuser Leading Edge
b3 = b2

C3 = C2
C3_w = C2_w * (r_2/r_3)
C3_r_t = np.sqrt((C3**2)-(C3_w**2))
To3 = To2
M3_abs = np.sqrt((C3**2)/(To3*gamma*R_idealgas-(((gamma-1)/2)*(C3**2))))
T3 = (C3**2)/(gamma*R_idealgas*(M3_abs**2))
P3 = Po3 * ((T3/To3)**(gamma/(gamma-1)))
Rho_static_3 = P3/(R_idealgas*T3)
Area_3 = 2 * np.pi * r_3 * b3 
C3_r_t = (m_dot_3)/(Rho_static_3 * Area_3) 
Alpha3 = np.arctan(C3_r_t/C3_w)
C3_t = np.sqrt((C3_w**2)+(C3_r_t**2))

C3 = C3_t

converged = False
iteration_count = 0
tolerance = 1e-6
max_iterations = 400


for i in range (100):
    C3 = C3_t
    C3_w = C2_w * (r_2/r_3)
    C3_r_t = np.sqrt((C3**2)-(C3_w**2))
    To3 = To2
    M3_abs = np.sqrt((C3**2)/(To3*gamma*R_idealgas-(((gamma-1)/2)*(C3**2))))
    T3 = (C3**2)/(gamma*R_idealgas*(M3_abs**2))
    P3 = Po3 * ((T3/To3)**(gamma/(gamma-1)))
    Rho_static_3 = P3/(R_idealgas*T3)
    Area_3 = 2 * np.pi * r_3 * b3 
    C3_r_t = (m_dot_3)/(Rho_static_3 * Area_3) 
    Alpha3 = np.arctan(C3_w/C3_r_t)
    new_C3_t = np.sqrt((C3_w**2)+(C3_r_t**2))
    C3 = new_C3_t
    if C3 - new_C3_t == 0.0001:
        break
    

C3 = new_C3_t
print('C3:', C3, 'm')
print('Cw3:', C3_w, 'm/s')
print('C3rt:', C3_r_t, 'm/s')
print('To3:', To3, 'K')
print('P3:', P3, 'Pa')
print('Static Density at stage 3:', Rho_static_3, 'kg/m^3')
print('Area3:', Area_3, 'm^2')
print('Alpha3:', (180*Alpha3)/(np.pi), 'degrees')
print('T3:', T3, 'K')
print('M3:', M3_abs)

### At Throat
print("----------------At Throat-----------------")
To3_choke = To3
Po3_choke = Po3
M3_choke = 1
T3_choke = To3_choke/(1+(((gamma-1)/2)*(M3_choke**2)))
P3_choke = Po3_choke * ((T3_choke/To3_choke)**(gamma/(gamma-1)))
P3_rcy = (P3_choke - P3)/(Po3 - P3)

print('P3_rcy =', P3_rcy, 'Pa')
print('T3_choke:', T3_choke, 'K')
print('P3_choke:', P3_choke, 'Pa')
print('P3_choke:', P3_choke, 'Pa')
    ################################################################   Find B3_star using graph and P3_rcy

B3_star = 0.035
C3_choke = np.sqrt(gamma*R_idealgas*T3_choke)
m_dot_3_choke = (1+0.01)*(m_dot_3)
Area_3_choke = m_dot_3_choke/((1-B3_star)*(P3_choke/(R_idealgas*T3_choke))*(C3_choke))
N_v_values = np.arange(1, 32, 1)

print('m_dot_3_choke:', m_dot_3_choke, 'kg/s')
print('Area_3_choke:', Area_3_choke, 'm^2')
print('C3_choke:', C3_choke, 'm/s')

for N_v in N_v_values: 
    C3_choke = np.sqrt(gamma*R_idealgas*T3_choke)
    m_dot_3_choke = (1 + 0.01)*(m_dot_3)
    Area_3_choke = m_dot_3_choke/((1-B3_star)*(P3_choke/(R_idealgas*T3_choke))*(C3_choke))
    W3_star = Area_3_choke / (N_v * b3)
    if 0.8 <= b3/W3_star <= 1.2:
        break

print('Number of vanes that satisfy the condition =', N_v)
print('W3* = ', W3_star)
print('B*/W*:', b3/W3_star )


### Diffuser Vane
## From Diffuser Vane Graph:
i3 = (-4.7*np.pi)/180
Alpha3_star = Alpha3 - i3   ###radians
Phi = ((360 / N_v)*(np.pi/180))   ### radians              

print('Phi:', Phi, 'radians')
print ('Alpha3star:', (180*Alpha3_star)/(np.pi), 'degrees')
print('i3:', i3, 'radians')

#Dirty Method
print("----------------Diffuser Wedge Angle-----------------")

Omega = (34.925*np.pi)/180        #### radians  (GUESSED) 34.9
a2ss = -np.tan( (np.pi/2) - Alpha3_star + Phi - (Omega/2) )
b2ss = r_3 * np.cos(Phi) + r_3 * np.sin(Phi) * np.tan((np.pi/2) - Alpha3_star + Phi - (Omega/2))

x_values = np.arange (0, W3_star, 0.0000001)    

for x in x_values:
    y1 = a2ss * x + b2ss
    y2 = np.sqrt((W3_star**2) - x**2) + r_3
    y3 = -np.sqrt((W3_star**2) - x**2) + r_3
    if abs(y1 - y2) == 0.0001:
        break
    
print (x)
print('y1:', y1)
print('y2:', y2)
print('y3:', y3)

print("----------------Diffuser Exit-----------------")

### Throat to exit
T5 = To5/(1+(((gamma-1)/2)*(M5_abs**2)))
P5 = Po5 * ((T5/To5)**(gamma/(gamma-1)))
P5_rcy = (P5 - P3_choke)/(Po3 - P3_choke)


print('T5:', T5, 'K')
print('P5:', P5, 'Pa')
print('P5_rcy:', P5_rcy, 'Pa')


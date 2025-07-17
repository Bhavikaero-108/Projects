import numpy as np
import matplotlib.pyplot as plt

m_dot = 5.2163
T_combustor_exit = 1251.89
To3 = 662.75

LHV = 40007.2
Po4 = 1182.006

Eta_HPT = 0.84
Eta_m = 0.99
Vane = 0.03
Disk = 0.0165

Eta_LPT = 0.91
Eta_PT = 0.92

Me = 0.15

Cp = 1.148
Cpa = 1.005
Cpg = 1.148

gamma = 1.33

W_LPC = 873.47
W_HPC = 1046.925

Pa = 101.325

Temp_combustor = []
Mass_flow_bleeded = []
Mass_flow_bleeded_vane = []
Mass_flow_bleeded_disk = []
HPT_efficiency = []
SFC_calculation = []


def temp_raise(raise_of_T):

    nb_degrees = raise_of_T
    Temp_max = T_combustor_exit + nb_degrees
    Temp_combustor.append(Temp_max)

    #number_of_values = int((Temp_max - T_combustor_exit) // 310.928)
    #print(number_of_values)

    #for i in range(number_of_values + 1):

        #Temp_combustor.append(T_combustor_exit + i * 310.928)


def mass_flow_bleeded():

    for i in range(len(Temp_combustor)):

        Mass_flow_bleeded.append(m_dot - m_dot * (0.09 + i * 0.01))
        Mass_flow_bleeded_vane.append(Mass_flow_bleeded[i] * (0.03 + i * 0.01))
        Mass_flow_bleeded_disk.append(Mass_flow_bleeded[i] * (0.0165 + i * 0.01))

def efficiency_HPT():

    for i in range(len(Temp_combustor)):

        HPT_efficiency.append(Eta_HPT - i * 0.002)


def calculation(Temp_max):

    temp_raise(Temp_max)
    mass_flow_bleeded()
    efficiency_HPT()

    print('Temp =', Temp_combustor)
    print('ma_dot', Mass_flow_bleeded)
    print('mv_dot =', Mass_flow_bleeded_vane)
    print('md_dot =', Mass_flow_bleeded_disk)
    print('HPT efficiency', HPT_efficiency)

    for i in range(len(Temp_combustor)):

        print('-------------------------Temp = ', Temp_combustor[i],'Efficiency =', HPT_efficiency[i], '------------------------------')

        # Combustor

        mf_dot = Mass_flow_bleeded[i] * (Cpg * Temp_combustor[i] - Cpa * To3) / (LHV - Cpg * Temp_combustor[i])

        #mf_dot = Mass_flow_bleeded[i] * 0.02
        print('ma_dot =', Mass_flow_bleeded[i])
        print('mf_dot = ', mf_dot)

        # HPT

        mg1_dot = Mass_flow_bleeded[i] + mf_dot
        mv_dot = Mass_flow_bleeded_vane[i] #* Vane
        md_dot = Mass_flow_bleeded_disk[i] #* Disk

        To4a = (mg1_dot * Cp * Temp_combustor[i] + mv_dot * Cpa * To3) / ((mg1_dot + mv_dot) * Cpg)
        W_HPT = W_HPC / 0.99

        mg2_dot = mg1_dot + mv_dot
        Tg = To4a

        To5 = -W_HPT / ((mg2_dot + md_dot) * Cpg) + To4a

        Po5 = Po4 * (1 - (1 / HPT_efficiency[i]) * (1 - (To5 / To4a))) ** (gamma / (gamma - 1))

        mg3_dot = mg2_dot + mv_dot

        To5b = (mg2_dot * Cpg * To5 + md_dot * Cpa * To3) / ((mg2_dot + md_dot) * Cpg)

        print('----------------HPT calculation---------------------')
        print('mg1_dot = ', mg1_dot)
        print('mv_dot = ', mv_dot)
        print('md_dot = ', md_dot)
        print('To4a = ', To4a)
        print('W_HPT = ', W_HPT)
        print('mg2_dot = ', mg2_dot)
        print('Tg = ', Tg)
        print('To5 = ', To5)
        print('Po5 = ', Po5)
        print('mg3_dot = ', mg3_dot)
        print('To5b = ', To5b)



        # LPT

        mg3_dot = mg3_dot + md_dot

        md2_dot = 0.011 * (md_dot + mg2_dot)

        W_LPT = W_LPC / 0.99

        To6 = (-W_LPT / ((mg3_dot + md2_dot) * Cpg)) + To5b

        Po6 = Po5 * (1 - (1 / Eta_LPT) * (1 - (To6 / To5))) ** (gamma / (gamma - 1))

        To6a = (mg3_dot * Cpg * To6 + md2_dot * Cpa * To3) / ((mg3_dot + md2_dot) * Cpg)

        print('----------------LPT calculation---------------------')
        print('mg3_dot = ', mg3_dot)
        print('md2_dot = ', md2_dot)
        print('W_LPT = ', W_LPT)
        print('To6 = ', To6)
        print('Po6 = ', Po6)
        print('To6a = ', To6a)

        # PT

        Po7 = Po6 * (1 - 0.006)

        To7 = To6a

        Mg4 = mg3_dot - 0.0125 * mg3_dot

        print('-----------------PT calculation---------------------')
        print('Po7 = ', Po7)
        print('To7 = ', To7)
        print('Mg4 = ', Mg4)

        # Nozzle expansion

        Po9 = (Pa / (1 + ((gamma - 1) / 2) * Me ** 2) ** (-gamma / (gamma - 1)))

        Po8 = Po9 + 0.02 * Po9

        To8a = (Po8 / Po7) ** ((gamma - 1) / gamma) * To7

        To8 = -Eta_PT * (To7 - To8a) + To7

        W_PT = Mg4 * Cpg * (To7 - To8) * Eta_m

        print('-----------------Nozzle expansion---------------------')
        print('Po9 = ', Po9)
        print('Po8 = ', Po8)
        print('To8a = ', To8a)
        print('To8 = ', To8)
        print('W_PT = ', W_PT)

        # SFC

        #W_PT = 935
        SFC = mf_dot * 3600 / (Mass_flow_bleeded[i] * (W_PT / Mg4))
        SFC_impacted = SFC + SFC * (i * 0.02) * 1.011
        SFC_calculation.append(SFC_impacted)

        print('-----------------SFC---------------------')
        print('SFC =', SFC)


def graphics(Temp_max):

    calculation(Temp_max)

    plt.plot(Temp_combustor, SFC_calculation)
    plt.xlabel('Temperature at combustor exit (°K)')
    plt.ylabel('SFC (kg/kWh)')
    plt.show()

    for i in range(len(HPT_efficiency)):

        HPT_efficiency[i] = HPT_efficiency[i] * 100

    plt.plot(Temp_combustor, HPT_efficiency)
    plt.xlabel('Temperature at combustor exit (°K)')
    plt.ylabel('Efficiency of the turbine (%)')
    plt.show()


#graphics(310.928)
calculation(0)











try:
    from pointgroup import PointGroup
    pointgroup_install = 'YES'
except ImportError:
    print("WARRNING : Module 'pointgroup' not found. But it is not nessesary.\n")
    print("pointgourp module used for auto-estimation of RXN Symmetry Number.")
    print("Your can install module for auto via 'pip install pointgroup'")
    print("")
    pointgroup_install = 'NO'

import numpy as np
import math
from decimal import Decimal

# Define Basic conversion factor
wavenumber2frequency = 29979245800 #[cm/s] It's a speed of light
hartree2jmol = 2625.5002e3
kcalmol2jmol = 4.184e3
kjmol2jmol = 1e3

# more symmetry number to be add
rot_symmetry_number = {'C1' : 1,
                       'Cs' : 1,
                       'C2' : 2,
                       'C2v': 2,
                       'C3v': 3,
                       'D2h': 4,
                       'D3h': 6,
                       'D5h': 10,
                       'D3d': 6,
                       'Td' : 12,
                       'Oh' : 24}

def sec_to_time(seconds):
    if seconds < 1:
        out = str(round(seconds,4)) + ' seconds'
        return out 
    else:
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)

        # suppose 365 days = 1 year, 30days = 1 month
        years = days // 365
        months = (days % 365) // 30
        days = (days % 365) % 30

        result = []
        count = 0
        if years:
            result.append(f"{years} years")
            count += 1
        if months and count < 2:
            result.append(f"{int(months)} months")
            count += 1
        if days and count < 2:
            result.append(f"{int(days)} days")
            count += 1
        if hours and count < 2:
            result.append(f"{int(hours)} hours")
            count += 1
        if minutes and count < 2:
            result.append(f"{int(minutes)} minutes")
            count += 1
        if seconds and count < 2:
            result.append(f"{int(seconds)} seconds")
            count += 1

        return " ".join(result)

# wigner transmission correction factor
def kappa_wigner(v, T):    
    # constants
    h_plk = 6.62607015e-34 # [J/s]
    k_b = 1.380649e-23 # [J/K]
                        
    kappa = 1 + (1/24) * np.power((h_plk * np.abs(v) / (k_b * T)),2)
    return kappa

# Eyring equation
def rate_constant_tst(T, dG):
    # dG in J/mol
    k_b = 1.380649e-23 # [J/K]
    h_plk = 6.62607015e-34 # [J/s]
    R_gas = 8.31446261815324 # [J/K*mol]
    k_TST = (k_b * T / h_plk) * np.exp(-dG / (R_gas * T))
    return k_TST

# percent consumption time ( when percent =50, it is half life of enantiomeric excess) 
def consumption_time(percent, k_r):
    # unit
    # k_rac [1/s]
    # consumption_time [s]
    consumption_ratio = percent / 100.0
    consumption_time = - math.log(1.0-consumption_ratio) / k_r
    return consumption_time


##############################################

#####       to be implemented         ##### 
def pg_from_coordinates(xyz_coordi):
    ## to be implemented
    pg = 'none'
    return pg

##############################################

# Rxn info
# This code is witten for monomolecular rxn, A -> TS -> B

print("\n 1 : Hartree(default) \n 2 : kcal/mol \n 3 : kJ/mol")
Energy_unit = input("\n Unit type : ")
print("")

unit = {'1' : 'Hartree', '2' : 'kcal/mol', '3' : 'kJ/mol'}

conversion_factor = hartree2jmol # default
if Energy_unit == "1":
    conversion_factor = hartree2jmol
elif Energy_unit == "2":
    conversion_factor = kcalmol2jmol
elif Energy_unit == "3":
    conversion_factor = kjmol2jmol
else :
    Energy_unit = "1"
    print("Set Default(hartree)\n")

Gibbs_Rc_own = float(input(" Gibbs_Rc : "))   # [own unit]
Gibbs_TS_own = float(input(" Gibbs_TS : "))  # [own unit]
Barrier_own = Gibbs_TS_own - Gibbs_Rc_own  # [own unit]

Gibbs_Rc = Gibbs_Rc_own * conversion_factor  # [J/mol]
Gibbs_TS = Gibbs_TS_own * conversion_factor # [J/mol]
Barrier =  Barrier_own * conversion_factor # [J/mol]


# Rxn symmetry
if pointgroup_install == 'YES':
    print("\n # RXN symmetry number #\n")
    print(" 1 : auto-estimation mode")
    print(" 2 : manual mode")
    RSN_mode = input("\n Answer : ")
    if RSN_mode == '1':
        reactant_xyz_path = input(" Reactant xyz file path : ")
        TS_xyz = input(" TS xyz file path : ")
        with open(reactant_xyz_path, 'r') as file:
            reactant_xyz = file.read()
            reactant_pg = pg_from_coordinates(reactant_xyz)
            reactant_rsn = rot_symmetry_number[reactant_pg] # rotational symmetry number
        with open(TS_xyz, 'r') as file:
            TS_xyz = file.read()
            TS_pg = pg_from_coordinates(TS_xyz)
            TS_rsn = rot_symmetry_number[TS_pg] # rotational symmetry number
        Sym_num = float(reactant_rsn) / float(TS_rsn)
    elif RSN_mode == '2':
        Sym_num = float(input("\n RXN Symmetry Number : "))
    else:
        print(" Invalid input")
elif pointgroup_install == 'NO':
    Sym_num = float(input("\n RXN Symmetry Number : "))


TS_freq_cm = float(input(" Transision State Imaginary Frequency [ cm^-1 ] : "))
TS_freq_s = TS_freq_cm * wavenumber2frequency

print("\n---------------------------------------------------------------------------------")
print("\n                                  # Rxn info. #                                  \n")
print(f"   Reaction Symmetry :   {Sym_num}\n")
print(f"   Imaginary Freq.   :   {TS_freq_cm}  [cm^-1]         {TS_freq_s}  [s^-1] \n")
print(f"   Reactant Energy   :   {Gibbs_Rc}  [J/mol]         {Gibbs_Rc_own}  [{unit[Energy_unit]}] ")
print(f"   Product  Energy   :   {Gibbs_TS}  [J/mol]         {Gibbs_TS_own}  [{unit[Energy_unit]}] ")
print(f"   Reation Barrier   :   {Barrier}  [J/mol]         {Barrier_own}  [{unit[Energy_unit]}] ")
print("\n---------------------------------------------------------------------------------\n\n")

# Temperature input 

print("Type 'help' for help")
print("\n")
while True:
    Temp_input = input("Temperatures? [K] : ")
    if Temp_input == 'help':
        print("")
        print("# Temperature input Styles #")
        print('----------------------------------------------------------------------------------------------')
        print("Style #1  :  Temp_1   Temp_2   Temp_3   Temp_4 ...")
        print("Style #2  :  start_Temp   end_Temp   interval")
    else:
        try:
            if isinstance(float(Temp_input.split()[0]),float):
                Temp = [float(tmp) for tmp in Temp_input.split()]
                if len(Temp) == 1:
                    break
                elif Temp[-1] > Temp[-2]:
                    break
                elif Temp[-1] < Temp[-2] and len(Temp) == 3:
                    Temp = [tmp for tmp in np.arange(Temp[0], Temp[1]+Temp[2], Temp[2])]
                    break
                else:
                    print("Code may be fix...")
        except ValueError:
            print("Invalid Temp input")
print("")

while True:
    STemp_input = input("\n* No Special Temperature? Then type 'ENTER'\n* Recommended : 'auto'\n\nSpecial Temperature? [K] : ")
    if STemp_input == 'help':
        print("")
        print("# Special Temperature input Styles #")
        print('----------------------------------------------------------------------------------------------')
        print("Style #1  :  Stemp_1 Stemp_2 Stemp_3 ... \n e.g. 273.15 293.15 298.15")
        print("Style #2  :  'auto'  will  set  273.15(STP) 293.15(NTP) 298.15(SATP) ")
    else:
        try:
            if STemp_input.strip() == '':
                STemp = []
                break
            elif STemp_input == 'auto':
                STemp = [273.15, 293.15, 298.15]
                break
            else:
                for_ValueError = int(" ")
        except ValueError:
            STemp = [float(tmp) for tmp in STemp_input.split()]
            if isinstance(float(STemp[0]),float):
                break
            else:
                print("Invalid STemp input")

print("\n\n")
print('----------------------------------------------------------------------------------------------')
print("\n  #  Temperature and Special Temperature input  #  \n")
print("Temperature List : ", Temp)
print("Special Temperature List : ", STemp)
print("")
print('----------------------------------------------------------------------------------------------')
Temp.extend(STemp)
Temp_list = sorted(Temp)
print("\n  #  Temperature List  #  \n")
print("Total Temperature List", Temp_list)
print("")
print('----------------------------------------------------------------------------------------------')

# consumption time
# half life (default)
consumption_percent = 50 

# activate here!
"""
print("This input is for calculating the decay-time of your reactant.")
print("But meaning of decay-time is depends on your reacton system.")
print("Plz be aware that 'What is the decay-time for which concentration?'\n")
print("50% (default) for half-life")
consumption_percent = float(input(" Consumption_percent : "))
"""

# formatting
print('Temperature       kappa_wig             k_eyring               dacay time')
print('----------------------------------------------------------------------------------------------')
for temperature in Temp_list:
    raw = ''
    k_wigner = kappa_wigner(TS_freq_s, temperature)
    k_TST_crd = rate_constant_tst(temperature, Barrier) * Sym_num * k_wigner
    decay_time = sec_to_time(consumption_time(consumption_percent, k_TST_crd))
    raw = '            '.join(map(str, ['%.2f' %temperature, '%.4E' % Decimal(k_wigner), '%.4E' % Decimal(k_TST_crd), decay_time]))
    print(raw)

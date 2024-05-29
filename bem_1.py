import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from math import acos, exp, sqrt, atan, cos, sin, radians, degrees

def BEM(v0, omega, pitch):
    # Constants
    B = 3  # number of blades
    R = 63  # rotor radius
    hubrad = 1.5  # hub radius
    rou = 1.225  # density of air
    EPS = 0.00001  # iterative precision tolerance

    a = 0
    a_prime = 0

    # Read blade section data
    BS = pd.read_csv('Blade/Blade section/Blade section.csv').to_numpy()

    # Read aero data files
    import glob
    AD_files = glob.glob('Blade/Aero data/*.csv')
    AD = []
    for file in AD_files:
        AD.append(pd.read_csv(file).to_numpy())

    NBS = len(BS)  # Number of blade sections
    Rx = np.zeros(NBS)
    FN = np.zeros(NBS)
    FT = np.zeros(NBS)
    Mx = np.zeros(NBS)

    # Main loop
    for i in range(NBS):
        ADofBS = int(BS[i, 1]) - 1
        r = BS[i, 2]
        Rx[i] = r
        dr = BS[i, 3]
        Theta = BS[i, 4]
        chord = BS[i, 5]
        alpha = AD[ADofBS][:, 0]
        Cl = AD[ADofBS][:, 1]
        Cd = AD[ADofBS][:, 2]
        Sigma = chord * B / (2 * np.pi * r)
        ax = a
        ax_prime = a_prime
        a = ax - 10 * EPS
        a_prime = ax_prime - 10 * EPS

        numite = 0
        while abs(ax - a) >= EPS or abs(ax_prime - a_prime) >= EPS:
            numite += 1
            a = ax
            a_prime = ax_prime
            Phi = atan((1 - a) * v0 / ((1 + a_prime) * r * omega))
            Phi = degrees(Phi)
            Alpha = Phi - Theta - pitch
            interp_Cl = interp1d(alpha, Cl, fill_value="extrapolate")
            interp_Cd = interp1d(alpha, Cd, fill_value="extrapolate")
            Cla = interp_Cl(Alpha)
            Cda = interp_Cd(Alpha)
            Cn = Cla * cos(radians(Phi)) + Cda * sin(radians(Phi))
            Ct = Cla * sin(radians(Phi)) - Cda * cos(radians(Phi))
            f_tiploss = B / 2 * (R - r) / (r * sin(radians(Phi)))
            F_tiploss = (2 / np.pi) * acos(exp(-f_tiploss))
            f_hubloss = B / 2 * (r - hubrad) / (r * sin(radians(Phi)))
            F_hubloss = (2 / np.pi) * acos(exp(-f_hubloss))
            F = F_tiploss * F_hubloss
            ac = 0.2
            if ax > ac:
                K = 4 * F * sin(radians(Phi)) ** 2 / (Sigma * Cn)
                ax = 0.5 * (2 + K * (1 - 2 * ac) - sqrt((K * (1 - 2 * ac) + 2) ** 2 + 4 * (K * ac ** 2 - 1)))
            else:
                ax = 1 / (4 * F * (sin(radians(Phi))) ** 2 / (Sigma * Cn) + 1)
            ax_prime = 1 / (4 * F * sin(radians(Phi)) * cos(radians(Phi)) / (Sigma * Ct) - 1)
            if numite >= 100:
                ax = 0.3
                ax_prime = 0.1

        FN[i] = 0.5 * rou * ((r * omega * (1 + a_prime)) ** 2 + (v0 * (1 - a)) ** 2) * chord * Cn * dr
        FT[i] = 0.5 * rou * ((r * omega * (1 + a_prime)) ** 2 + (v0 * (1 - a)) ** 2) * chord * Ct * dr
        Mx[i] = FT[i] * r

    M = sum(Mx)  # rotor torque from one blade
    P = M * omega * 3 * 0.944  # Power calculation
    return Rx, FN, FT, P

# Example usage:
v0 = 10  # inflow wind speed in m/s
omega = 1  # rotational speed in rad/s
pitch = 2  # pitch angle in degrees
Rx, FN, FT, P = BEM(v0, omega, pitch)
print('Radius:', Rx)
print('Normal Force:', FN)
print('Tangential Force:', FT)
print('Power:', P)
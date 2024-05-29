import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from numpy import arccos, exp, sqrt, arctan, cos, sin, radians, degrees

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
    blade_sections = pd.read_csv('Blade/Blade section/Blade section.csv').to_numpy()

    # Read aero data files
    import glob
    airfoil_data_files = glob.glob('Blade/Aero data/*.csv')
    airfoil_data = []
    for file in airfoil_data_files:
        airfoil_data.append(pd.read_csv(file).to_numpy())

    n_sections = len(blade_sections)  # Number of blade sections
    Rx = np.zeros(n_sections)
    FN = np.zeros(n_sections)
    FT = np.zeros(n_sections)
    Mx = np.zeros(n_sections)

    # Main loop: from root to tip section
    for i, blade_section in enumerate(blade_sections):
        _, airfoil, r, dr, theta_deg, chord = blade_section
        airfoil_index = int(airfoil) - 1

        # Airfoil data
        alphas = airfoil_data[airfoil_index][:, 0]  # Degrees!
        lift_coefficients = airfoil_data[airfoil_index][:, 1]
        drag_coefficients = airfoil_data[airfoil_index][:, 2]
        interp_Cl = interp1d(alphas, lift_coefficients, fill_value="extrapolate")
        interp_Cd = interp1d(alphas, drag_coefficients, fill_value="extrapolate")

        Sigma = chord * B / (2 * np.pi * r)  # solidity
        ax = a  # change value
        ax_prime = a_prime  # change value
        a = ax - 10 * EPS  # generate error, active iteration
        a_prime = ax_prime - 10 * EPS  # generate error, active iteration

        numite = 0
        while abs(ax - a) >= EPS or abs(ax_prime - a_prime) >= EPS:
            numite += 1

            #  record results of last step
            a = ax
            a_prime = ax_prime
            phi = arctan((1 - a) * v0 / ((1 + a_prime) * r * omega))
            alpha_deg = np.degrees(phi) - theta_deg - pitch

            # find Cl and Cd
            lift_coefficient = interp_Cl(alpha_deg)
            drag_coefficient = interp_Cd(alpha_deg)

            # projection in and out of plane
            Cn = lift_coefficient * cos(phi) + drag_coefficient * sin(phi)
            Ct = lift_coefficient * sin(phi) - drag_coefficient * cos(phi)

            # Prandtl loss
            f_tiploss = B / 2 * (R - r) / (r * sin(phi))
            F_tiploss = (2 / np.pi) * arccos(exp(-f_tiploss))
            f_hubloss = B / 2 * (r - hubrad) / (r * sin(phi))
            F_hubloss = (2 / np.pi) * arccos(exp(-f_hubloss))
            F = F_tiploss * F_hubloss

            if np.isnan(F):
                F = 1

            # Glauert correction
            ac = 0.2
            if ax > ac:
                K = 4 * F * sin(phi) ** 2 / (Sigma * Cn)
                ax = 0.5 * (2 + K * (1 - 2 * ac) - sqrt((K * (1 - 2 * ac) + 2) ** 2 + 4 * (K * ac ** 2 - 1)))
            else:
                ax = 1 / (4 * F * (sin(phi)) ** 2 / (Sigma * Cn) + 1)
            ax_prime = 1 / (4 * F * sin(phi) * cos(phi) / (Sigma * Ct) - 1)

            # in case of iterative convergence failure
            if np.isnan(ax) or np.isnan(ax_prime):
                ax = 0.3
                ax_prime = 0.1

        # Update value
        a = ax
        a_prime = ax_prime

        # force in two directions and bending moment
        FN[i] = 0.5 * rou * ((r * omega * (1 + a_prime)) ** 2 + (v0 * (1 - a)) ** 2) * chord * Cn * dr
        FT[i] = 0.5 * rou * ((r * omega * (1 + a_prime)) ** 2 + (v0 * (1 - a)) ** 2) * chord * Ct * dr
        Mx[i] = FT[i] * r
        Rx[i] = r

    M = sum(Mx)  # rotor torque from one blade
    P = M * omega * 3 * 0.944  # Power calculation
    return Rx, FN, FT, P
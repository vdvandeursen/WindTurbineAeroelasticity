import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from numpy import arccos, exp, sqrt, arctan, cos, sin, radians, degrees
from Shape_function import ShapeFunction


def BEM(v0, omega, pitch, Vf="no input", Ve="no input", shape_functions="no input"):
    # Constants
    B = 3  # number of blades
    R = 63  # rotor radius
    hubrad = 1.5  # hub radius
    rou = 1.225  # density of air
    EPS = 0.001  # iterative precision tolerance
    MAX_ITER = 1000
    # Initial inductions
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

    # Set blade and flapwise arrays to 0 if there is no input
    if isinstance(Vf, str):
        Vf = 0
    if isinstance(Ve, str):
        Ve = 0
    if isinstance(shape_functions, str):
        shape_functions = []
        shape_functions.append(ShapeFunction([0], R))
        shape_functions.append(ShapeFunction([0], R))

    Rx = np.zeros(n_sections)
    FN = np.zeros(n_sections)
    FT = np.zeros(n_sections)
    Ff = np.zeros(n_sections)
    Fe = np.zeros(n_sections)
    Mx = np.zeros(n_sections)

    # Main loop: from root to tip section
    for i, blade_section in enumerate(blade_sections):
        _, airfoil, r, dr, theta_deg, chord = blade_section
        airfoil_index = int(airfoil) - 1

        # Determine flapwise and edgewise velocities at location r using shape function
        Vf_i = Ve * shape_functions[0].f(r)
        Ve_i = Vf * shape_functions[1].f(r)

        # Blade velocity (Translate to normal and tangential)
        Vtb = np.cos(np.radians(theta_deg + pitch)) * Ve_i + np.sin(np.radians(theta_deg + pitch)) * Vf_i
        Vnb = np.cos(np.radians(theta_deg + pitch)) * Vf_i - np.sin(np.radians(theta_deg + pitch)) * Ve_i

        # Airfoil data
        alphas = airfoil_data[airfoil_index][:, 0]  # Degrees!
        lift_coefficients = airfoil_data[airfoil_index][:, 1]
        drag_coefficients = airfoil_data[airfoil_index][:, 2]
        interp_Cl = interp1d(alphas, lift_coefficients, fill_value="extrapolate")
        interp_Cd = interp1d(alphas, drag_coefficients, fill_value="extrapolate")

        Sigma = chord * B / (2 * np.pi * r)  # Solidity
        ax = a  # Change value
        ax_prime = a_prime  # Change value
        a = ax - 10 * EPS  # Generate error, active iteration
        a_prime = ax_prime - 10 * EPS  # Generate error, active iteration
        numite = 0
        while (abs(ax - a) >= EPS or abs(ax_prime - a_prime) >= EPS) and numite < MAX_ITER:
            numite += 1

            # Record results of last step
            a = 0.8 * a + 0.2 * ax
            a_prime = 0.8 * a_prime + 0.2 * ax_prime
            phi = arctan(((1 - a) * v0 + Vnb) / ((1 + a_prime) * r * omega - Vtb))
            alpha_deg = np.degrees(phi) - theta_deg - pitch

            # Find Cl and Cd
            lift_coefficient = interp_Cl(alpha_deg)
            drag_coefficient = interp_Cd(alpha_deg)

            # Projection in and out of plane
            Cn = lift_coefficient * cos(phi) + drag_coefficient * sin(phi)
            Ct = lift_coefficient * sin(phi) - drag_coefficient * cos(phi)

            # Prandtl loss
            f_tiploss = B / 2 * (R - r) / (r * sin(phi))
            F_tiploss = (2 / np.pi) * arccos(exp(-f_tiploss))
            f_hubloss = B / 2 * (r - hubrad) / (r * sin(phi))
            F_hubloss = (2 / np.pi) * arccos(exp(-f_hubloss))
            F = F_tiploss * F_hubloss

            if np.isnan(F) or F < EPS:
                F = 1

            # Glauert correction
            ac = 0.2
            if ax > ac:
                K = 4 * F * sin(phi) ** 2 / (Sigma * Cn)
                ax = 0.5 * (2 + K * (1 - 2 * ac) - sqrt((K * (1 - 2 * ac) + 2) ** 2 + 4 * (K * ac ** 2 - 1)))
            else:
                ax = 1 / (4 * F * (sin(phi)) ** 2 / (Sigma * Cn) + 1)
            ax_prime = 1 / (4 * F * sin(phi) * cos(phi) / (Sigma * Ct) - 1)

            # Handle convergence failure
            if np.isnan(ax) or np.isnan(ax_prime):
                ax = 0.3
                ax_prime = 0.1
                break  # Exit the loop if NaN values persist

        # Update value
        a = ax
        a_prime = ax_prime
        if np.isnan(ax) or np.isnan(ax_prime):
            ax = 0.3
            ax_prime = 0.1

        # Force in two directions and bending moment
        FN[i] = (0.5 * rou * ((r * omega * (1 + a_prime)) ** 2 + (v0 * (1 - a)) ** 2) * chord * Cn * dr)
        FT[i] = (0.5 * rou * ((r * omega * (1 + a_prime)) ** 2 + (v0 * (1 - a)) ** 2) * chord * Ct * dr)
        Ff[i] = np.cos(np.radians(theta_deg + pitch)) * FN[i] + np.sin(np.radians(theta_deg + pitch)) * FT[i]
        Fe[i] = np.cos(np.radians(theta_deg + pitch)) * FT[i] - np.sin(np.radians(theta_deg + pitch)) * FN[i]
        Mx[i] = FT[i] * r
        Rx[i] = r

    M = sum(Mx)  # Rotor torque from one blade
    P = M * omega * 3 * 0.944  # Power calculation
    return Rx, Ff/B, Fe/B
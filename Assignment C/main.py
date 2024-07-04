import pandas as pd

from dynamic_stall import DynamicStallModel
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
from matplotlib.ticker import LinearLocator


if __name__ == "__main__":

    case1 = {
        'name': 'DYN1',
        'U1_U0': 1,
        'DeltaU': 0.5,
        'OmegaR_U0': 0
    }

    case2 = {
        'name': 'DYN2',
        'U1_U0': 1,
        'DeltaU': 0.5,
        'OmegaR_U0': 0.3
    }

    speed_wind_free_stream = 10
    span_blade = 50
    n_blades = 3
    ratio_tip_speed = 8
    ratio_blade_start = 0.2
    pitch_blade = np.radians(-2)

    def function_twist(r_R): return np.radians(14-14*r_R)
    def function_chord(r_R): return 3*(1-r_R) + 1

    def Uinf(t, U0, case_dict):
        omega_fun = case_dict['OmegaR_U0'] * U0 / span_blade
        omega_rotor = ratio_tip_speed * U0 / span_blade
        angle_azimuth = t * omega_rotor

        return case_dict['U1_U0'] * U0 + case_dict['DeltaU'] * np.cos(omega_fun * t) * np.cos(angle_azimuth)

    # Runs four cases: DYN2 with leading edge sep, DYN1 without, DYN1 with, and DYN1 with steady aero
    for ls_mod, case, steady_aero in zip([True, False, True, False], [case2, case1, case1, case1], [False, False, False, True]):
    # for ls_mod, case, steady_aero in zip([False], [case1], [True]):
        if case['name'] == 'DYN2':
            dur = 30
        else:
            dur = 10

        model = DynamicStallModel(
            chord_func=function_chord,
            twist_func=function_twist,
            span_blade=span_blade,
            ratio_blade_start=ratio_blade_start,
            n_blades=n_blades,
            airfoil='DU95W180'
        )

        res = model.run(
            time_start=0,
            time_end=dur,
            num_time=1000,
            pitch=pitch_blade,
            omega=ratio_tip_speed * speed_wind_free_stream / span_blade,
            u_inf_func=lambda t: Uinf(t, U0=speed_wind_free_stream, case_dict=case),
            leading_edge_separation=ls_mod,
            steady_aero=steady_aero
        )

        for quantity in ['alpha', 'F_t', 'F_n']:
            for orientation in ['left', 'right']:
                matplotlib.rc('font', size=12)
                fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8, 7))

                if orientation == 'left':
                    ax.view_init(elev=ax.elev, azim=ax.azim-45)

                # Make data.
                X = res['r']
                Y = res['t']
                X, Y = np.meshgrid(X, Y)
                # R = np.sqrt(X ** 2 + Y ** 2)
                Z = res[quantity].T

                # # Add color bar for U0
                # C = Uinf(Y, speed_wind_free_stream, case1)
                # norm = plt.Normalize(C.min(), C.max())
                # colors = matplotlib.cm.coolwarm(norm(C))
                # m = matplotlib.cm.ScalarMappable(cmap=matplotlib.cm.coolwarm, norm=norm)
                # m.set_array(C)

                surf = ax.plot_surface(X, Y, Z, cmap=matplotlib.cm.coolwarm, linewidth=0, antialiased=False)

                # Set axis labels
                ax.set_xlabel('$r$ [m]')
                ax.set_xlim(0, None)
                ax.set_ylabel(f't [s]')

                # Z label depends on the plotted quantity
                if 'F' in quantity:
                    z_label = f'${quantity}$ [N/m]'
                elif 'M' in quantity:
                    z_label = f'${quantity}$ [Nm]'
                else:
                    z_label = f'$\\{quantity}$ [deg]'

                ax.set_zlabel(z_label)
                fig.colorbar(surf, shrink=0.5, aspect=10, label=z_label)

                plt.tight_layout()

                if steady_aero:
                    plt.savefig(f'{case["name"]} {dur}s {quantity} {orientation} STEADY AERO.pdf')
                else:
                    plt.savefig(f'{case["name"]} {dur}s {quantity} {orientation} {"with" if ls_mod else "no"} LS.pdf')
                # plt.show()

    print('break')

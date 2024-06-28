from beddoes_leishman import DynamicStallModel
import numpy as np


if __name__ == "__main__":

    case1 = {
        'U1_U0': 1,
        'DeltaU': 0.5,
        'OmegaR_U0': 0
    }

    case2 = {
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
        time_end=10,
        pitch=pitch_blade,
        omega=ratio_tip_speed * speed_wind_free_stream / span_blade,
        u_inf_func=lambda t: Uinf(t, U0=speed_wind_free_stream, case_dict=case1),
    )

    print('break')

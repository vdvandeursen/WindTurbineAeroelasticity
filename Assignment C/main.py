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
    length_blade = 50
    n_blades = 3
    ratio_tip_speed = 8
    ratio_blade_start = 0.2
    pitch_blade = np.radians(-2)
    def function_twist(r): return np.radians(14-14*r/length_blade)
    def function_chord(r): return 3*(1-r/length_blade) + 1

    def Uinf(t, angle_azimuth, U0, case_dict):
        omega = case_dict['OmegaR_U0'] * U0 / length_blade

        return case_dict['U1_U0'] * U0 + case_dict['DeltaU'] * np.cos(omega * t) * np.cos(angle_azimuth)


    model = DynamicStallModel()

    print('break')

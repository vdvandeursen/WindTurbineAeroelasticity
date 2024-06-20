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



    model = DynamicStallModel()

    print('break')

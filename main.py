from structural_model import StructuralModel
from bem import BEM
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    v0 = np.array([3, 4, 5, 6, 7, 8, 9, 10, 11, 11.4, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25])
    rotational_frequancies = np.array(
        [5.7173, 6.9957, 7.470993355, 7.887136213, 8.396777409, 9.014936877, 10.14196013, 11.27405316, 11.85744186,
         12.1, 12.10207641, 12.10166113, 12.10111296, 12.10069767, 12.10004983, 12.09983389, 12.09961794, 12.09928571,
         12.09950166, 12.09960133, 12.09965116, 12.09975083, 12.09945183, 12.09956811])
    pitch_angles = np.array(
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.83, 6.6, 8.7, 10.45, 12.06, 13.54, 14.92, 16.23, 17.47, 18.68, 19.85, 21.02,
         22.12, 23.15])

    structural_model = StructuralModel(filepath='structural_data.csv')

    print(f'Nat freqs: {structural_model.natural_frequencies}')

    labels = ['flap displ', 'flap velocity', 'edge displ', 'edge velocity']
    plt.figure()

    for v, omega, pitch in zip(v0, rotational_frequancies, pitch_angles):
        if v < 10:
            continue

        r_R, FN, FT, P = BEM(v, omega, pitch)

        res = structural_model.calculate_time_response(
            time_span=(0, 0.1),
            f_flap=FN,
            f_edge=FT,
            r_R=r_R
        )

        for i in range(4):
            plt.plot(res.t, res.y[i, :], label=f'{labels[i]} for {v} m/s')

        print('Radius:', r_R)
        print('Normal Force:', FN)
        print('Tangential Force:', FT)
        print('Power:', P)

        plt.legend()
        plt.show()
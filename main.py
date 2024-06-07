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

    rotational_frequancies = rotational_frequancies / 60 * 2 * np.pi

    pitch_angles = np.array(
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.83, 6.6, 8.7, 10.45, 12.06, 13.54, 14.92, 16.23, 17.47, 18.68, 19.85, 21.02,
         22.12, 23.15])

    structural_model = StructuralModel(filepath='structural_data.csv')

    print(f'Nat freqs: {structural_model.natural_frequencies}')

    labels = ['flap displ', 'edge displ', 'flap velocity', 'edge velocity']
    plt.figure()
    flap_dis = []
    edge_dis = []
    timestamps = np.linspace(0, 500, 2000)
    initial_conditions = np.array([0,0,0,0])
    for v, omega, pitch in zip(v0, rotational_frequancies, pitch_angles):
        res = structural_model.calculate_time_response_static_load(timestamps,initial_conditions,v,omega,pitch)
        flap_dis.append(res.y[0,-1])
        edge_dis.append(res.y[1,-1])

    # for i in range(4):
    #     plt.plot(res.t, res.y[i, :], label=f'{labels[i]} for {v} m/s')
    # plt.legend()
    # plt.show()

    plt.plot(v0,flap_dis,label="Flapwise Displacement at the tip")
    plt.plot(v0,edge_dis,label="Edgewise Displacement at the tip")
    plt.legend()
    plt.grid()
    plt.xlabel("Wind Speed [m/s]")
    plt.ylabel("Displacement [m]")
    plt.show()

    res = structural_model.calculate_time_response_varying_velocity(timestamps=np.linspace(0, 200, 2000))
    for i in range(4):
        plt.plot(res.t, res.y[i, :], label=f'{labels[i]} for 15 m/s')
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement & Speed [m] / [m/s]")
    plt.show()
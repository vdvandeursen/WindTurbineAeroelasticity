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

    structural_model = StructuralModel(filepath='./Assignment A/structural_data.csv')

    print(f'Nat freqs: {structural_model.natural_frequencies}')

    labels = ['Flapwise displacement', 'Edgewise displacement', 'Flapwise velocity', 'Edgewise velocity']
    colors = ['b','r','g','y']
    plt.figure()
    flap_dis = []
    edge_dis = []
    timestamps = np.linspace(0, 400, 1000)
    N = len(timestamps)
    T = timestamps[1] - timestamps[0]
    frequencies = np.fft.fftfreq(N, T)
    initial_conditions = np.array([0,0,0,0])

    #Static displacement for varyind wind speed 3-25 m/s
    for v, omega, pitch in zip(v0, rotational_frequancies, pitch_angles):
        res, Mn, FF, FE = structural_model.calculate_time_response_static_load(timestamps,initial_conditions,v,omega,pitch)
        flap_dis.append(res.y[0,-1])
        edge_dis.append(res.y[1,-1])
    plt.plot(v0,flap_dis,label=f'{labels[0]}',color=colors[0])
    plt.plot(v0,edge_dis,label=f'{labels[2]}',color=colors[1])
    plt.title("Static tip displasment for varying wind speed")
    plt.legend()
    plt.grid()
    plt.xlabel("Wind Speed [m/s]")
    plt.ylabel("Displacement [m]")
    plt.show()

    #Constant velocity of 15 m/s without taking into account the blade velocity
    res, Mn, FF, FE = structural_model.calculate_time_response_static_load(timestamps, initial_conditions, v0[13],rotational_frequancies[13],pitch_angles[13])
    #Plotting Displacement and velocities
    for i in [2,0]:
        plt.plot(res.t, res.y[i, :], label=f'{labels[i]}',color=colors[i])
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement & Velocity [m] / [m/s]")
    plt.show()

    for i in [3,1]:
        plt.plot(res.t, res.y[i, :], label=f'{labels[i]}',color=colors[i])
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement & Velocity [m] / [m/s]")
    plt.show()

    #Constant velocity of 15 m/s with taking into account the blade velocity
    res, Mn_lst, FF_lst, FE_lst= structural_model.calculate_time_response_dynamic_load(timestamps,initial_conditions,v0[13],rotational_frequancies[13],pitch_angles[13],periodic="False",blade_velocities="True")
    #Time Plot velocities and displacement
    for i in [2,0]:
        plt.plot(timestamps, res[i, :], label=f'{labels[i]}',color=colors[i])
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement & Velocity [m] / [m/s]")
    plt.show()

    for i in [3,1]:
        plt.plot(timestamps, res[i, :], label=f'{labels[i]}',color=colors[i])
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement & Velocity [m] / [m/s]")
    plt.show()

    #Time plot Forces and Moments
    plt.plot(timestamps[1:],FF_lst, label = "Flapwise Force at the tip [N]")
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Force [N]")
    plt.show()

    plt.plot(timestamps[1:],FE_lst, label = "Edgewise Force at the tip [N]")
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Force [N]")
    plt.show()

    plt.plot(timestamps[1:],Mn_lst, label = "Root Bending Moment [N m]")
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Moment [Nm m]")
    plt.show()

    #Frequency Plot
    for i in range(2):
        fft_values = np.fft.fft(res[i, :])
        magnitude = np.abs(fft_values)
        plt.plot(frequencies[0:N//2-1], 2*magnitude[0:N//2-1], label=f'{labels[i]}')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude")
    plt.show()


    #Varying velocity around 15 m/s without taking into account the blade velocity
    res, Mn_lst, FF_lst, FE_lst = structural_model.calculate_time_response_dynamic_load(timestamps,initial_conditions,v0[13],rotational_frequancies[13],pitch_angles[13],periodic="True",blade_velocities="False")
    #Time Plot velocities and displacemen
    for i in [2,0]:
        plt.plot(timestamps, res[i, :], label=f'{labels[i]}',color=colors[i])
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement & Velocity [m] / [m/s]")
    plt.show()

    for i in [3,1]:
        plt.plot(timestamps, res[i, :], label=f'{labels[i]}',color=colors[i])
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement & Velocity [m] / [m/s]")
    plt.show()


    #Varying velocity around 15 m/s with taking into account the blade velocity
    res, Mn_lst, FF_lst, FE_lst = structural_model.calculate_time_response_dynamic_load(timestamps,initial_conditions,v0[13],rotational_frequancies[13],pitch_angles[13],periodic="True",blade_velocities="True")
    #Time Plot
    for i in [0,2]:
        plt.plot(timestamps, res[i, :], label=f'{labels[i]}',color=colors[i])
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement & Velocity [m] / [m/s]")
    plt.show()

    for i in [1,3]:
        plt.plot(timestamps, res[i, :], label=f'{labels[i]}',color=colors[i])
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement & Velocity [m] / [m/s]")
    plt.show()

    # Time plot Forces and Moments
    plt.plot(timestamps[1:], FF_lst, label="Flapwise Force at the tip [N]")
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Force [N]")
    plt.show()

    plt.plot(timestamps[1:], FE_lst, label="Edgewise Force at the tip [N]")
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Force [N]")
    plt.show()

    plt.plot(timestamps[1:], Mn_lst, label="Root Bending Moment [N m]")
    plt.legend()
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Moment [Nm m]")
    plt.show()

    #Frequency Plot
    for i in range(2):
        fft_values = np.fft.fft(res[i, :])
        magnitude = np.abs(fft_values)
        plt.plot(frequencies[0:N//2-1], 2*magnitude[0:N//2-1], label=f'{labels[i]}')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude")
    plt.show()
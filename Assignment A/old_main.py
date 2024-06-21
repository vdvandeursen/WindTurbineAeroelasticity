from structural_model import StructuralModel
from bem import BEM
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter,ScalarFormatter
import os
import itertools


save_plot = False
folder_path = r"C:\Users\olegr\Documents\0. AE Master Wind Energy\Q4 AE4W21-14 Wind Turbine Aeroelasticity\Assignements"
os.makedirs(folder_path, exist_ok=True)
def one_decimal(x, pos):
    return f'{x:.1f}'

def plot_saver(plot_name, save_plot):
    if save_plot:
        plt.tight_layout()
        plt.savefig(os.path.join(folder_path, f'{plot_name}.pdf'))
        plt.close()
    else:
        plt.title(plot_name)
        plt.tight_layout()
        plt.show()

def adjust_x_y_ticks(start_plot,end_time):
    plt.xlim(start_plot,end_time)
    plt.xticks(range(start_plot, end_time+ 1, 20),  # Tick positions
               range(0, end_time - start_plot + 1, 20))  # Tick labels
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.yaxis.get_offset_text().set_fontsize(10)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

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

    labels = ['Flapwise displacement', 'Edgewise displacement', 'Flapwise force', 'Edgewise force','Root bending moment' ]
    colors = ["tab:blue","tab:orange","tab:green","tab:red","tab:purple"]
    plt.figure()
    N_T = 1000/400
    start_time = 0
    start_plot = 400 # time from where you want to plot the output
    end_time = 500
    timestamps = np.linspace(0, end_time, int(end_time * N_T))

    N = len(timestamps[int(start_plot * N_T):int(end_time * N_T)])
    T = timestamps[1] - timestamps[0]

    frequencies = np.fft.fftfreq(N, T)
    initial_conditions = np.array([0,0,0,0])

    flap_dis_no_stiff = []
    edge_dis_no_stiff = []
    flap_vel_no_stiff = []
    edge_vel_no_stiff = []

    #Static displacement for varyind wind speed 3-25 m/s
    for v, omega, pitch in zip(v0, rotational_frequancies, pitch_angles):
        res, Mn, FF, FE = structural_model.calculate_time_response_static_load(timestamps,initial_conditions,v,omega,pitch,geometric_stiffness=False)
        flap_dis_no_stiff.append(res.y[0,-1])
        edge_dis_no_stiff.append(res.y[1,-1])
        flap_vel_no_stiff.append(res.y[2, -1])
        edge_vel_no_stiff.append(res.y[3, -1])

    plt.plot(v0,flap_dis_no_stiff,label=f'{labels[0]}',color=colors[0])
    plt.plot(v0,edge_dis_no_stiff,label=f'{labels[1]}',color=colors[1])
    plt.legend(loc='upper right')
    plt.grid()
    plt.gca().yaxis.set_major_formatter(FuncFormatter(one_decimal))
    plt.xlabel("Wind Speed [m/s]")
    plt.xlim(v0[0], v0[-1])
    plt.ylabel("Displacement [m]")
    plot_name = "Static_displacement_no_geo_stiffening"
    plot_saver(plot_name,save_plot)

    flap_dis_with_stiff = []
    edge_dis_with_stiff = []
    flap_vel_with_stiff = []
    edge_vel_with_stiff = []
    #Static displacement for varyind wind speed 3-25 m/s with geometric stiffness
    for v, omega, pitch in zip(v0, rotational_frequancies, pitch_angles):
        res, Mn, FF, FE = structural_model.calculate_time_response_static_load(timestamps,initial_conditions,v,omega,pitch,geometric_stiffness=True)
        flap_dis_with_stiff.append(res.y[0,-1])
        edge_dis_with_stiff.append(res.y[1,-1])
        flap_vel_with_stiff.append(res.y[2, -1])
        edge_vel_with_stiff.append(res.y[3, -1])

    plt.plot(v0,flap_dis_with_stiff,label=f'{labels[0]}',color=colors[0])
    plt.plot(v0,edge_dis_with_stiff,label=f'{labels[1]}',color=colors[1])
    plt.legend(loc='upper right')
    plt.grid()
    plt.gca().yaxis.set_major_formatter(FuncFormatter(one_decimal))
    plt.xlabel("Wind Speed [m/s]")
    plt.xlim(v0[0],v0[-1])
    plt.ylabel("Displacement [m]")
    plot_name = "Static_displacement_with_geo_stiffening"
    plot_saver(plot_name,save_plot)

    #extra plots static displacement
    max_flap_no_stiff = np.max(flap_dis_no_stiff)
    max_flap_with_stiff = np.max(flap_dis_with_stiff)
    plt.plot(v0, flap_dis_no_stiff/max_flap_no_stiff , label='Without geometrical stiffening', color=colors[0])
    plt.plot(v0,flap_dis_with_stiff/max_flap_with_stiff, label='With geometrical stiffening', color=colors[0],linestyle="--")
    plt.legend(loc='upper right')
    plt.grid()
    plt.gca().yaxis.set_major_formatter(FuncFormatter(one_decimal))
    plt.xlabel("Wind Speed [m/s]")
    plt.xlim(v0[0],v0[-1])
    plt.ylabel("Normalized flapwise displacement [-]")
    plot_name = "Static_displacement_flapwise"
    plot_saver(plot_name,save_plot)

    max_edge_no_stiff = np.max(edge_dis_no_stiff)
    max_edge_with_stiff = np.max(edge_dis_with_stiff)
    plt.plot(v0, edge_dis_no_stiff/max_edge_no_stiff, label='Without geometrical stiffening',  color=colors[1])
    plt.plot(v0,edge_dis_with_stiff/max_edge_with_stiff, label='With geometrical stiffening', color=colors[1],linestyle="--")
    plt.legend(loc='upper right')
    plt.grid()
    plt.gca().yaxis.set_major_formatter(FuncFormatter(one_decimal))
    plt.xlabel("Wind Speed [m/s]")
    plt.xlim(v0[0],v0[-1])
    plt.ylabel("Normalized edgewise displacement [-]")
    plot_name = "Static_displacement_edgewise"
    plot_saver(plot_name,save_plot)

    wind_speed_dynamic = 15 # m/s
    index_v0 = np.where(v0 == wind_speed_dynamic)[0][0]

def core_code(periodic_yes_no,blade_velocities_yes_no,geometric_stiffness_yes_no):
    if periodic_yes_no:
        add_period = True
        name_periodic = "periodic_15_ms"
    else:
        add_period = False
        name_periodic = "const_15_ms"
    if blade_velocities_yes_no:
        add_blade_vel = True
        name_periodic = "with_blade_vel"
    else:
        add_blade_vel = False
        name_periodic = "no_blade_vel"

    if geo_stiffness_yes_no:
        initial_conditions = np.array([flap_dis_with_stiff[index_v0],edge_dis_with_stiff[index_v0],flap_vel_with_stiff[index_v0],edge_vel_with_stiff[index_v0]])
        add_geo_stiff = True
        name_geo = "with_geo_stiff"
    else:
        initial_conditions = np.array([flap_dis_no_stiff[index_v0], edge_dis_no_stiff[index_v0], flap_vel_no_stiff[index_v0], edge_vel_no_stiff[index_v0]])
        add_geo_stiff = False
        name_geo = "no_geo_stiff"

    #Constant velocity of 15 m/s without taking into account the blade velocity
    res, Mn, FF, FE = structural_model.calculate_time_response_static_load(timestamps, initial_conditions, v0[index_v0],rotational_frequancies[index_v0],pitch_angles[index_v0],geometric_stiffness=add_geo_stiff)

    #Plotting Displacement
    plt.plot(res.t[int(start_plot * N_T):int(end_time * N_T)], res.y[0, int(start_plot * N_T):int(end_time * N_T)], label=f'{labels[0]}',color=colors[0])
    plt.legend(loc='upper right')
    plt.grid()
    adjust_x_y_ticks(start_plot,end_time)
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement [m]")
    plot_name = "Flapwise_disp_const_15ms_no_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)

    plt.plot(res.t[int(start_plot * N_T):int(end_time * N_T)], res.y[1, int(start_plot * N_T):int(end_time * N_T)], label=f'{labels[1]}',color=colors[1])
    plt.legend(loc='upper right')
    plt.grid()
    adjust_x_y_ticks(start_plot,end_time)
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement [m]")
    plot_name = "Edgewise_disp_const_15ms_no_blade_velo_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)

    #Constant velocity of 15 m/s with taking into account the blade velocity
    res, Mn_lst, FF_lst, FE_lst= structural_model.calculate_time_response_dynamic_load(timestamps,initial_conditions,v0[index_v0],rotational_frequancies[index_v0],pitch_angles[index_v0],periodic=False,blade_velocities=True,geometric_stiffness=add_geo_stiff)
    #Time Plot displacement
    plt.plot(timestamps[int(start_plot * N_T):int(end_time * N_T)], res[0, int(start_plot * N_T):int(end_time * N_T)],label=f'{labels[0]}', color=colors[0])
    plt.legend(loc='upper right')
    plt.grid()
    adjust_x_y_ticks(start_plot,end_time)
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement [m]")
    plot_name = "Flapwise_disp_const_15ms_with_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)

    plt.plot(timestamps[int(start_plot * N_T):int(end_time * N_T)], res[1, int(start_plot * N_T):int(end_time * N_T)],label=f'{labels[1]}', color=colors[1])
    plt.legend(loc='upper right')
    plt.grid()
    adjust_x_y_ticks(start_plot,end_time)
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement [m]")
    plot_name = "Edgewise_disp_const_15ms_with_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)

    #Time plot Forces and Moments
    plt.plot(timestamps[int(start_plot * N_T):int(end_time * N_T)-1],FF_lst[int(start_plot * N_T):int(end_time * N_T)],label=f'{labels[2]}', color=colors[2])
    plt.legend(loc='upper right')
    plt.grid()
    adjust_x_y_ticks(start_plot,end_time)
    plt.xlabel("Time [s]")
    plt.ylabel("Force at tip [N]")
    plot_name = "Flapwise_force_const_15ms_with_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)

    plt.plot(timestamps[int(start_plot * N_T):int(end_time * N_T)-1],FE_lst[int(start_plot * N_T):int(end_time * N_T)],label=f'{labels[3]}', color=colors[3])
    plt.legend(loc='upper right')
    plt.grid()
    adjust_x_y_ticks(start_plot,end_time)
    plt.xlabel("Time [s]")
    plt.ylabel("Force at tip [N]")
    plot_name = "Edgewise_force_const_15ms_with_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)

    plt.plot(timestamps[int(start_plot * N_T):int(end_time * N_T)-1],Mn_lst[int(start_plot * N_T):int(end_time * N_T)],label=f'{labels[4]}', color=colors[4])
    plt.legend(loc='upper right')
    plt.grid()
    adjust_x_y_ticks(start_plot, end_time)
    plt.xlabel("Time [s]")
    plt.ylabel("Moment [Nm m]")
    plot_name = "Root_bend_moment_const_15ms_with_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)

    #Frequency Plot
    for i in range(2):
        fft_values = np.fft.fft(res[i,int(start_plot * N_T):int(end_time * N_T)])
        magnitude = np.abs(fft_values)
        plt.plot(frequencies[0:N//2-1], 2*magnitude[0:N//2-1], label=f'{labels[i]}',color=colors[i])
    plt.yscale('log')
    plt.legend(loc='upper right')
    plt.grid()
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude [-] ")
    plot_name = "Frequency_plot_const_15ms_with_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)


    #Varying velocity around 15 m/s without taking into account the blade velocity
    res, Mn_lst, FF_lst, FE_lst = structural_model.calculate_time_response_dynamic_load(timestamps,initial_conditions,v0[index_v0],rotational_frequancies[index_v0],pitch_angles[index_v0],periodic=True,blade_velocities=False,geometric_stiffness=add_geo_stiff)
    #Time Plot displacement
    plt.plot(timestamps[int(start_plot * N_T):int(end_time * N_T)], res[0, int(start_plot * N_T):int(end_time * N_T)],label=f'{labels[0]}', color=colors[0])
    plt.legend(loc='upper right')
    plt.grid()
    adjust_x_y_ticks(start_plot, end_time)
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement [m]")
    plot_name = "Flapwise_disp_varying_15ms_no_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)

    plt.plot(timestamps[int(start_plot * N_T):int(end_time * N_T)], res[1, int(start_plot * N_T):int(end_time * N_T)],
             label=f'{labels[1]}', color=colors[1])
    plt.legend(loc='upper right')
    plt.grid()
    adjust_x_y_ticks(start_plot, end_time)
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement [m]")
    plot_name = "Edgewise_displ_varying_15ms_no_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)

    #Varying velocity around 15 m/s with taking into account the blade velocity
    res, Mn_lst, FF_lst, FE_lst = structural_model.calculate_time_response_dynamic_load(timestamps,initial_conditions,v0[index_v0],rotational_frequancies[index_v0],pitch_angles[index_v0],periodic=True,blade_velocities=True,geometric_stiffness=add_geo_stiff)
    #Time Plot Displacement
    plt.plot(timestamps[int(start_plot * N_T):int(end_time * N_T)], res[0, int(start_plot * N_T):int(end_time * N_T)],label=f'{labels[0]}', color=colors[0])
    plt.legend(loc='upper right')
    plt.grid()
    adjust_x_y_ticks(start_plot, end_time)
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement [m]")
    plot_name = "Flapwise_displ_varying_15ms_with_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)

    plt.plot(timestamps[int(start_plot * N_T):int(end_time * N_T)], res[1, int(start_plot * N_T):int(end_time * N_T)],label=f'{labels[1]}', color=colors[1])
    plt.legend(loc='upper right')
    plt.grid()
    adjust_x_y_ticks(start_plot, end_time)
    plt.xlabel("Time [s]")
    plt.ylabel("Tip displacement [m]")
    plot_name = "Edgewise_displ_varying_15ms_with_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)

    # Time plot Forces and Moments
    plt.plot(timestamps[int(start_plot * N_T):int(end_time * N_T)-1],FF_lst[int(start_plot * N_T):int(end_time * N_T)],label=f'{labels[2]}', color=colors[2])
    plt.legend(loc='upper right')
    plt.grid()
    adjust_x_y_ticks(start_plot, end_time)
    plt.xlabel("Time [s]")
    plt.ylabel("Force at tip [N]")
    plot_name = "Flapwise_force_varying_15ms_with_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)

    plt.plot(timestamps[int(start_plot * N_T):int(end_time * N_T)-1],FE_lst[int(start_plot * N_T):int(end_time * N_T)],label=f'{labels[3]}', color=colors[3])
    plt.legend(loc='upper right')
    plt.grid()
    adjust_x_y_ticks(start_plot, end_time)
    plt.xlabel("Time [s]")
    plt.ylabel("Force at tip [N]")
    plot_name = "Edgewise_force_varying_15ms_with_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)

    plt.plot(timestamps[int(start_plot * N_T):int(end_time * N_T)-1],Mn_lst[int(start_plot * N_T):int(end_time * N_T)],label=f'{labels[4]}', color=colors[4])
    plt.legend(loc='upper right')
    plt.grid()
    adjust_x_y_ticks(start_plot, end_time)
    plt.xlabel("Time [s]")
    plt.ylabel("Moment [Nm m]")
    plot_name = "Root_bend_moment_varying_15ms_with_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)

    #Frequency Plot
    for i in range(2):
        fft_values = np.fft.fft(res[i,int(start_plot * N_T):int(end_time * N_T)])
        magnitude = np.abs(fft_values)
        plt.plot(frequencies[0:N//2-1], 2*magnitude[0:N//2-1], label=f'{labels[i]}',color=colors[i])
    plt.yscale('log')
    plt.legend(loc='upper right')
    plt.grid()
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude [-] ")
    plot_name = "Frequency_plot_varying_15ms_with_blade_vel_{}".format(add_to_name)
    plot_saver(plot_name,save_plot)


# Generate all combinations of True and False for three variables
combinations_variables = list(itertools.product([True, False], repeat=3))

# Iterate over each combination and run the code
for combo in combinations_variables:
    run_code(*combo)

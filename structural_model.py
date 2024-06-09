from scipy import integrate
import numpy as np
import pandas as pd
from bem import BEM
from Shape_function import ShapeFunction

dof_settings = {
    'flap': {
        'EI': 'FlpStff\r',
        'damping_ratio': 0.477465/100,
        'coefficients': [0, 0, 0.0622, 1.7254, -3.2452, 4.7131, -2.2555]  # a + bx + cx**2 + dx**3 etc
    },
    'edge':
        {
            'EI': 'EdgStff\r',
            'damping_ratio': 0.477465/100,
            'coefficients': [0, 0, 0.3627, 2.5337, -3.5772, 2.2376, -0.6952]
        },
}


class StructuralModel:
    def __init__(self, filepath):
        structural_data = pd.read_csv(filepath, header=0)
        structural_data.columns = [col.split('\n')[0] for col in structural_data.columns]
        r = structural_data['Radius\r'].to_numpy()
        rho_A = structural_data['BMassDen\r'].to_numpy()
        self.R = structural_data['Radius\r'].iloc[-1]
        self.mass_matrix = np.eye(2)
        self.stiffness_matrix = np.eye(2)
        self.damping_matrix = np.eye(2)
        self.force_vector = []
        self.force_vector_list = []
        self.shape_functions = []
        for i, (dof, settings) in enumerate(dof_settings.items()):
            phi = ShapeFunction(coefficients=settings['coefficients'], R=self.R)
            EI = structural_data[settings['EI']].to_numpy()

            M = integrate.trapz(rho_A * (phi.f(r) ** 2), x=r)
            K = integrate.trapz(EI * (phi.d2f_dr2(r) ** 2), x=r)
            C = 2 * settings['damping_ratio'] * np.sqrt(K * M)

            self.shape_functions.append(phi)
            self.mass_matrix[i, i] = M
            self.stiffness_matrix[i, i] = K
            self.damping_matrix[i, i] = C

    @property
    def natural_frequencies(self):
        return np.sqrt(np.divide(self.stiffness_matrix, self.damping_matrix)).diagonal()

    def __set_force_vector(self, f_edge, f_flap, r):
        force_vect = []
        for phi, f in zip(self.shape_functions, [f_flap, f_edge]):
            force_vect.append(integrate.trapz(phi.f(r) * f, x=r))

        self.force_vector = np.array(force_vect)

    def calculate_time_response_static_load(self, timestamps,initial_conditions,V0,omega,pitch,Vf = 0, Ve = 0):
        r, Ff, Fe, Mn = BEM(V0,omega,pitch,Vf,Ve,self.shape_functions)
        self.__set_force_vector(f_flap=Ff, f_edge=Fe, r=r)
        for i, _ in enumerate(self.force_vector):
            def _ivp(t, y):
                x = y[:2]
                x_dot = y[2:]
                x_ddot = np.linalg.inv(self.mass_matrix) @ (self.force_vector - self.damping_matrix @ x_dot - self.stiffness_matrix @ x)
                return np.concatenate((x_dot, x_ddot))
        t_span = timestamps[0], timestamps[-1]
        res = integrate.solve_ivp(_ivp, t_span, initial_conditions, t_eval=timestamps)

        return res, Mn, self.force_vector[0], self.force_vector[1]

    def calculate_time_response_dynamic_load(self, timestamps,initial_conditions,V0,omega,pitch,periodic="False",blade_velocities="False"):
        input_conditions = initial_conditions
        res = np.array(input_conditions).reshape(-1, 1)
        v = V0
        Vf = 0
        Ve = 0
        Mn_lst = []
        FF_lst = []
        FE_lst = []
        for i in range(len(timestamps)-1):
            # print(i/len(timestamps)*100, "% done")
            t0 = timestamps[i]
            t1 = timestamps[i+1]
            thalf = (t0+t1)/2
            if periodic == "True":
                v = V0 + 0.5 * np.cos(1.267 * thalf) + 0.085 * np.cos(2.534 * thalf) + 0.015 * np.cos(3.801 * thalf)
            timestep = np.linspace(t0,t1,50)
            if blade_velocities == "True":
                Vf = input_conditions[2]
                Ve = input_conditions[3]
            res_step,Mn,FF,FE = self.calculate_time_response_static_load(timestep,input_conditions,v,omega,pitch,Vf,Ve)
            input_conditions = res_step.y[:, -1]
            res = np.concatenate((res, input_conditions.reshape(-1, 1)), axis=1)
            Mn_lst.append(Mn)
            FF_lst.append(FF)
            FE_lst.append(FE)
        return res, Mn_lst, FF_lst, FE_lst





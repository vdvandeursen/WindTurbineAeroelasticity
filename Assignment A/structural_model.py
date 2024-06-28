from scipy import integrate
import numpy as np
import pandas as pd
from bem import BEM
from shape_function import ShapeFunction

dof_settings = {
    'flap': {
        'EI': 'FlpStff',
        'damping_ratio': 0.477465/100,
        'coefficients': [0, 0, 0.0622, 1.7254, -3.2452, 4.7131, -2.2555]  # a + bx + cx**2 + dx**3 etc
    },
    'edge':
        {
            'EI': 'EdgStff',
            'damping_ratio': 0.477465/100,
            'coefficients': [0, 0, 0.3627, 2.5337, -3.5772, 2.2376, -0.6952]
        },
}




class StructuralModel:
    def __init__(self, filepath):
        structural_data = pd.read_csv(filepath, header=0)
        structural_data.columns = [col.split('\n')[0].replace('\r', '') for col in structural_data.columns]
        self.r = structural_data['Radius'].to_numpy()
        self.rho_A = structural_data['BMassDen'].to_numpy()
        self.R = structural_data['Radius'].iloc[-1]
        self.mass_matrix = np.eye(2)
        self.stiffness_matrix = np.eye(2)
        self.damping_matrix = np.eye(2)
        self.geometric_stiffness_matrix = np.eye(2)
        self.force_vector = []
        self.force_vector_list = []
        self.shape_functions = []
        for i, (dof, settings) in enumerate(dof_settings.items()):
            phi = ShapeFunction(coefficients=settings['coefficients'], R=self.R)
            EI = structural_data[settings['EI']].to_numpy()

            M = integrate.trapz(self.rho_A * (phi.f(self.r) ** 2), x=self.r)
            K = integrate.trapz(EI * (phi.d2f_dr2(self.r) ** 2), x=self.r)
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

    def __N(self, omega):
        result = []
        for i,r_min in enumerate(self.r):
            integral = integrate.trapz(self.r[i:]*omega**2*self.rho_A[i:],x=self.r[i:])
            result.append(integral)
        return np.array(result)

    def __set_geometric_stiffness_matrix(self,Omega):
        for i,phi in enumerate(self.shape_functions):
            K_geo = integrate.trapz(self.__N(Omega)*(phi.df_dr(self.r) ** 2), x=self.r)
            self.geometric_stiffness_matrix[i, i] = K_geo

    def calculate_time_response_static_load(self, timestamps,initial_conditions,V0,omega,pitch,Vf = 0, Ve = 0,geometric_stiffness=False):
        r, Ff, Fe, Mn = BEM(V0,omega,pitch,Vf,Ve,self.shape_functions)
        self.__set_force_vector(f_flap=Ff, f_edge=Fe, r=r)
        if geometric_stiffness == True:
            self.__set_geometric_stiffness_matrix(omega)
        else:
            self.geometric_stiffness_matrix=np.zeros([2,2])
        for i, _ in enumerate(self.force_vector):
            def _ivp(t, y):
                x = y[:2]
                x_dot = y[2:]
                x_ddot = np.linalg.inv(self.mass_matrix) @ (self.force_vector - self.damping_matrix @ x_dot - (self.stiffness_matrix+self.geometric_stiffness_matrix) @ x)
                return np.concatenate((x_dot, x_ddot))
        t_span = timestamps[0], timestamps[-1]
        res = integrate.solve_ivp(_ivp, t_span, initial_conditions, t_eval=timestamps)

        return res, Mn, self.force_vector[0], self.force_vector[1]

    def calculate_time_response_dynamic_load(self, timestamps,initial_conditions,V0,omega,pitch,periodic=False,blade_velocities=False,geometric_stiffness=False):
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
            if periodic == True:
                v = V0 + 0.5 * np.cos(1.267 * thalf) + 0.085 * np.cos(2.534 * thalf) + 0.015 * np.cos(3.801 * thalf)
            timestep = np.linspace(t0,t1,10)
            if blade_velocities == True:
                Vf = input_conditions[2]
                Ve = input_conditions[3]
            res_step,Mn,FF,FE = self.calculate_time_response_static_load(timestep,input_conditions,v,omega,pitch,Vf,Ve,geometric_stiffness)
            input_conditions = res_step.y[:, -1]
            res = np.concatenate((res, input_conditions.reshape(-1, 1)), axis=1)
            Mn_lst.append(Mn)
            FF_lst.append(FF)
            FE_lst.append(FE)
        return res, Mn_lst, FF_lst, FE_lst





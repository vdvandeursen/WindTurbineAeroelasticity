from scipy import integrate
import numpy as np
import pandas as pd
from bem import BEM

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

    def calculate_time_response_constant_velocity(self, timestamps, f_edge, f_flap, r):
        self.__set_force_vector(f_edge=f_edge, f_flap=f_flap, r=r)

        for i, _ in enumerate(self.force_vector):
            def _ivp(t, y):
                x = y[:2]
                x_dot = y[2:]
                x_ddot = np.linalg.inv(self.mass_matrix) @ (self.force_vector - self.damping_matrix @ x_dot - self.stiffness_matrix @ x)
                return np.concatenate((x_dot, x_ddot))

            initial_conditions = np.array([
                0, 0, 0, 0
            ])
            t_span = timestamps[0], timestamps[-1]
            res = integrate.solve_ivp(_ivp, t_span, initial_conditions, t_eval=timestamps)

            return res

    def calculate_time_response_varying_velocity(self, timestamps):
        self.force_vector_list = []
        for t in timestamps:
            v0 = 15 + 0.5 *np.cos(1.267*t)+ 0.085*np.cos(2.534*t)+0.015*np.cos(3.801*t)
            r, f_flap, f_edge, P = BEM(v0, 12.06/ 60 * 2 * np.pi, 10.45)
            self.__set_force_vector(f_edge=f_edge,f_flap=f_flap,r=r)
            self.force_vector_list.append(self.force_vector)

        for i, _ in enumerate(self.force_vector):
            print(i)

            def _ivp(t, y):
                x = y[:2]
                x_dot = y[2:]
                j = int(t/(timestamps[-1]-timestamps[-2]))
                x_ddot = np.linalg.inv(self.mass_matrix) @ (
                            self.force_vector_list[j] - self.damping_matrix @ x_dot - self.stiffness_matrix @ x)
                return np.concatenate((x_dot, x_ddot))

            initial_conditions = np.array([
                0, 0, 0, 0
            ])
            t_span = timestamps[0], timestamps[-1]
            res = integrate.solve_ivp(_ivp, t_span, initial_conditions, t_eval=timestamps)

            return res

class ShapeFunction:
    def __init__(self, coefficients: list, R):
        """ Polynomial of order x = a + bx + cx**2 + dx**3 + ... etc """
        self.coefficients = coefficients
        self.R = R

    def f(self, r):
        f = 0
        for i, c in enumerate(self.coefficients):
            f += c * (r/self.R) ** i

        return f

    # def df_dr(self, r):
    #     df_dr = 0
    #     for i, c in enumerate(self.coefficients):
    #         df_dr += i * c * r ** (i-1) / self.R ** i
    #
    #     return df_dr

    def d2f_dr2(self, r):
        d2f_dr2 = 0
        for i, c in enumerate(self.coefficients):
            d2f_dr2 += i * (i - 1) * c * r ** (i - 2) / self.R ** i

        return d2f_dr2

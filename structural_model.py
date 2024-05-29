from scipy import integrate
import numpy as np
import pandas as pd


dof_settings = {
    'flap': {
        'EI': 'FlpStff',
        'damping_ratio': 0.477465,
        'coefficients': [0, 0, 0.0622, 1.7254, -3.2452, 4.7131, -2.2555]  # a + bx + cx**2 + dx**3 etc
    },
    'edge':
        {
            'EI': 'EdgStff',
            'damping_ratio': 0.477465,
            'coefficients': [0, 0, 0.3627, 2.5337, -3.5772, 2.2376, -0.6952]
        },
}


class StructuralModel:
    def __init__(self, filepath):
        structural_data = pd.read_csv(filepath, header=0)
        structural_data.columns = [col.split('\n')[0] for col in structural_data.columns]

        r_R = structural_data['BlFract'].to_numpy()
        rho_A = structural_data['BMassDen'].to_numpy()

        self.mass_matrix = np.eye(2)
        self.stiffness_matrix = np.eye(2)
        self.damping_matrix = np.eye(2)
        self.force_vector = []
        self.shape_functions = []

        for i, (dof, settings) in enumerate(dof_settings.items()):
            phi = ShapeFunction(coefficients=settings['coefficients'])
            EI = structural_data[settings['EI']].to_numpy()

            M = integrate.trapz(rho_A * (phi.f(r_R) ** 2))
            K = integrate.trapz(EI * (phi.d2f_dr2(r_R) ** 2))
            C = 2 * settings['damping_ratio'] * np.sqrt(K * M)

            self.shape_functions.append(phi)
            self.mass_matrix[i, i] = M
            self.stiffness_matrix[i, i] = K
            self.damping_matrix[i, i] = C

    @property
    def natural_frequencies(self):
        return np.sqrt(np.divide(self.stiffness_matrix, self.damping_matrix)).diagonal()

    def __set_force_vector(self, f_edge, f_flap, r_R):
        force_vect = []
        for phi, f in zip(self.shape_functions, [f_flap, f_edge]):
            force_vect.append(integrate.trapz(phi.f(r_R) * f))

        self.force_vector = np.array(force_vect)

    def calculate_time_response(self, time_span, f_edge, f_flap, r_R):
        self.__set_force_vector(f_edge=f_edge, f_flap=f_flap, r_R=r_R)

        for i, _ in enumerate(self.force_vector):
            def _ivp(t, y):
                x = y[:2]
                x_dot = y[2:]

                x_ddot = np.linalg.inv(self.mass_matrix) @ (self.force_vector - self.damping_matrix @ x_dot - self.stiffness_matrix @ x)
                return np.concatenate((x_dot, x_ddot))

            initial_conditions = np.array([
                0, 0, 0, 0
            ])

            res = integrate.solve_ivp(_ivp, time_span, initial_conditions, method='Radau')
            return res




class ShapeFunction:
    def __init__(self, coefficients: list):
        """ Polynomial of order x = a + bx + cx**2 + dx**3 + ... etc """

        super().__init__()

        while len(coefficients) < 3:
            coefficients.append(0)

        self.coefficients = coefficients

    def f(self, r):
        f = 0
        for i, c in enumerate(self.coefficients):
            f += c * r ** i

        return f

    def df_dr(self, r):
        dr_dr = 0
        for i, c in enumerate(self.coefficients[1:]):
            dr_dr += c * r ** i

        return dr_dr

    def d2f_dr2(self, r):
        d2f_dr2 = 0
        for i, c in enumerate(self.coefficients[2:]):
            d2f_dr2 += c * r ** i

        return d2f_dr2

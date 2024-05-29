from scipy import integrate
import numpy as np
import pandas as pd

from shape_functions import Polynomial


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
        'coefficients':  [0, 0, 0.3627, 2.5337, -3.5772, 2.2376, -0.6952]
    },
}


if __name__ == '__main__':
    structural_data = pd.read_csv('structural_data.csv', header=0)
    structural_data.columns = [col.split('\n')[0] for col in structural_data.columns]

    r_R = structural_data['BlFract'].to_numpy()
    rho_A = structural_data['BMassDen'].to_numpy()

    M_matrix = np.eye(2)
    K_matrix = np.eye(2)
    C_matrix = np.eye(2)

    for i, (dof, settings) in enumerate(dof_settings.items()):
        phi_1 = Polynomial(coefficients=settings['coefficients'])
        EI = structural_data[settings['EI']].to_numpy()

        M = integrate.trapz(rho_A*(phi_1.f(r_R))**2)
        K = integrate.trapz(EI*(phi_1.d2f_dr2(r_R))**2)
        C = 2 * settings['damping_ratio'] * np.sqrt(K * M)

        omega_n = np.sqrt(K / M)

        M_matrix[i, i] = M
        K_matrix[i, i] = K
        C_matrix[i, i] = C


    print('break')
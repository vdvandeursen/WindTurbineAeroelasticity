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


    v0 = np.array([3   ,	4   ,	5   ,	6   ,	7   ,	8   ,	9   ,	10  ,	11  ,	11.4    ,	12  ,	13  ,	14  ,	15  ,	16  ,	17  ,	18  ,	19  ,	20  ,	21  ,	22  ,	23  ,	24  ,	25])
    RtSpeeds = np.array([5.7173  ,	6.9957  ,	7.470993355 ,	7.887136213 ,	8.396777409 ,	9.014936877 ,	10.14196013 ,	11.27405316 ,	11.85744186 ,	12.1    ,	12.10207641 ,	12.10166113 ,	12.10111296 ,	12.10069767 ,	12.10004983 ,	12.09983389 ,	12.09961794 ,	12.09928571 ,	12.09950166 ,	12.09960133 ,	12.09965116 ,	12.09975083 ,	12.09945183 ,	12.09956811 ])
    Pitchangles = np.array([0   ,	0   ,	0   ,	0   ,	0   ,	0   ,	0   ,	0   ,	0   ,	0   ,	3.83    ,	6.6 ,	8.7 ,	10.45   ,	12.06   ,	13.54   ,	14.92   ,	16.23   ,	17.47   ,	18.68   ,	19.85   ,	21.02   ,	22.12   ,	23.15])



    print('break')
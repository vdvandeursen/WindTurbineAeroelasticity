from scipy import integrate
import numpy as np

from shape_functions import Polynomial


if __name__ == '__main__':
    phi_1f = Polynomial(coefficients=[
        0, 0, 0.0622, 1.7254, -3.2452, 4.7131, -2.2555
    ])

    phi_1e = Polynomial(coefficients=[
        0, 0, 0.3627, 2.5337, -3.5772, 2.2376, -0.6952
    ])

    r_R = np.linspace(0, 1, 100)

    M_1f = integrate.trapz(r_R)

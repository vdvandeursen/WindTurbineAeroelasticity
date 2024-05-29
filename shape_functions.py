

class Polynomial:
    def __init__(self, coefficients: list):
        """ Polynomial of order x = a + bx + cx**2 + dx**3 + ... etc """

        super().__init__()

        while len(coefficients) < 3:
            coefficients.append(0)

        self.coefficients = coefficients

    def f(self, s):
        x = 0
        for i, c in enumerate(self.coefficients):
            x += c * s ** i

        return x

    def df_ds(self, s):
        df_ds = 0
        for i, c in enumerate(self.coefficients[1:]):
            df_ds += c * s ** i

        return df_ds

    def d2f_ds2(self, s):
        d2f_ds2 = 0
        for i, c in enumerate(self.coefficients[2:]):
            d2f_ds2 += c * s ** i

        return d2f_ds2



def phi_1f(r, R):
    ratio = r / R
    result = (
            0.0622 * ratio ** 2 +
            1.7254 * ratio ** 3 -
            3.2452 * ratio ** 4 +
            4.7131 * ratio ** 5 -
            2.2555 * ratio ** 6
    )
    return result

def phi_1e(r, R):
    ratio = r / R

    # Evaluate the polynomial
    result = (
            0.3627 * ratio ** 2 +
            2.5337 * ratio ** 3 -
            3.5772 * ratio ** 4 +
            2.2376 * ratio ** 5 -
            0.6952 * ratio ** 6
    )

    return result
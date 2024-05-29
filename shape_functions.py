

class Polynomial:
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

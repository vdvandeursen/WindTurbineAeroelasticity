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

    def df_dr(self, r):
        df_dr = 0
        for i, c in enumerate(self.coefficients):
            df_dr += i * c * r ** (i-1) / self.R ** i

        return df_dr

    def d2f_dr2(self, r):
        d2f_dr2 = 0
        for i, c in enumerate(self.coefficients):
            d2f_dr2 += i * (i - 1) * c * r ** (i - 2) / self.R ** i

        return d2f_dr2
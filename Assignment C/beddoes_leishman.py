import numpy as np


class DynamicStallModel:
    def __init__(self):
        self.dCn_dalpha = 2 * np.pi  # lift slope
        self.alpha0 = 0  # alpha for which normal load is zero

        self.c_normal = None
        self.c_normal_circ = None
        self.c_normal_noncirc = None
        self.c_normal_prime = None
        self.c_normal_f = None
        self.c_normal_v = None

        self.alpha_qs = None
        self.dalphaqs_dt = None
        self.alpha_equivalent = None
        self.s_array = None
        self.fprimeprime = None

    def run(self, time, chord, u_inf):
        # Todo: correct the alpha_qs for fluctuation in u_inf instead of heaving motion of blade. Formula below is probably not correct
        speed_blade_tangential = 0
        self.alpha_qs = np.arctan(u_inf / speed_blade_tangential)

        self.dalphaqs_dt = np.gradient(self.alpha_qs, time)  # calculate the time derivative of the quasi-steady alpha
        self.s_array = 2 * u_inf * time / chord  # define the array semi-chord time scale

        self.unsteady_flat_plate(time=time, chord=chord, u_inf=u_inf)
        self.non_linear_trailing_edge_separation(time=time)
        self.leading_edge_vortex_shedding(time=time)

        return self.c_normal + self.c_normal_f + self.c_normal_v

    def unsteady_flat_plate(self, time, chord, u_inf):
        def _duhamel_approx(Xi, Yi, delta_s, delta_alpha, order=2, A1=0.3, A2=0.7, b1=0.14, b2=0.53):
            """determining X and Y terms for recursive marching formula for approximation of Duhamel's integral"""

            # A1=0.165,A2=0.335,b1=0.0455,b2=0.3
            # determine the next values of X and Y, named Xip1 and Yip1
            if order == 1:
                Xip1 = Xi * np.exp(-b1 * delta_s) + A1 * delta_alpha
                Yip1 = Yi * np.exp(-b2 * delta_s) + A2 * delta_alpha
            elif order == 2:
                Xip1 = Xi * np.exp(-b1 * delta_s) + A1 * delta_alpha * np.exp(-b1 * delta_s / 2)
                Yip1 = Yi * np.exp(-b2 * delta_s) + A2 * delta_alpha * np.exp(-b2 * delta_s / 2)
            else:
                Xip1 = Xi * np.exp(-b1 * delta_s) + A1 * delta_alpha * (
                        (1 + 4 * np.exp(-b1 * delta_s / 2) + np.exp(-b1 * delta_s)) / 6)
                Yip1 = Yi * np.exp(-b2 * delta_s) + A2 * delta_alpha * (
                        (1 + 4 * np.exp(-b2 * delta_s / 2) + np.exp(-b2 * delta_s)) / 6)

            return Xip1, Yip1

        def _circulatory_normal_force():
            """define function for circulatory force, potential flow"""
            return self.dCn_dalpha * (self.alpha_equivalent - self.alpha0)

        def _deficiency_function(Dnoncirc_i, delta_dalpha_dt, delta_t, asound=343, kalpha=0.75):
            """deficiency function for non-circulatory normal force"""

            # a sound is the speed of sound
            TI = chord / asound
            Dnoncirc_ip1 = Dnoncirc_i * np.exp(-delta_t / (kalpha * TI)) + delta_dalpha_dt * np.exp(
                -delta_t / (2 * kalpha * TI))
            return Dnoncirc_ip1

        def _non_circulatory_normal_force():
            """non-circulatory normal force"""
            k_alpha = 0.75
            return 4 * k_alpha * chord / u_inf * (self.dalphaqs_dt - Dnoncirc)

        # define arrays for X,Y and alpha_equivalent
        Xarray = np.zeros(np.shape(time))
        Yarray = np.zeros(np.shape(time))

        # define the array of alpha_equivalent
        alpha_equivalent = np.zeros(np.shape(time))
        Dnoncirc = np.zeros(np.shape(time))

        alpha_equivalent[0] = self.alpha_qs[0]

        # march solution in time for alpha_E, and Dnoncirc, the deficiency function for non-circulatory loading
        for i, val in enumerate(time[:-1]):
            Xarray[i + 1], Yarray[i + 1] = _duhamel_approx(
                Xi=Xarray[i],
                Yi=Yarray[i],
                delta_s=self.s_array[i + 1] - self.s_array[i],
                delta_alpha=self.alpha_qs[i + 1] - self.alpha_qs[i]
            )
            Dnoncirc[i + 1] = _deficiency_function(
                Dnoncirc_i=Dnoncirc[i],
                delta_dalpha_dt=self.dalphaqs_dt[i + 1] - self.dalphaqs_dt[i],
                delta_t=time[i+1] - time[i],
            )

        self.alpha_equivalent = self.alpha_qs - Xarray - Yarray
        self.c_normal_circ = _circulatory_normal_force()
        self.c_normal_noncirc = _non_circulatory_normal_force()
        self.c_normal = self.c_normal_circ + self.c_normal_noncirc

    def non_linear_trailing_edge_separation(self, time):
        # TODO: rewrite formula to use radians?

        def _f_trailing_edge_separation_point(alpha, a1=7, a2=15, a3=21):
            # receives alpha in radians, converts to degrees
            alphadeg = alpha * 180 / np.pi
            if alphadeg <= a1:
                f = 1
            elif (alphadeg > a1) and (alphadeg <= a2):
                f = 1 - .8 * ((alphadeg - a1) / (a2 - a1))
            elif (alphadeg > a2) and (alphadeg < a3):
                f = .2 * (1 - ((alphadeg - a2) / (a3 - a2)) ** .3)
            else:
                f = 0
            return f

        # we will now do the time marching to solve for the pressure lag deficiency function
        def pressure_lag_deficiency(Dpress_i, delta_s, delta_CNpot, Tp=1.7):
            return Dpress_i * np.exp(-delta_s / Tp) + delta_CNpot * np.exp(-delta_s / 2 / Tp)

        Dpress = np.zeros(np.shape(time))
        for i, val in enumerate(time[:-1]):
            Dpress[i + 1] = pressure_lag_deficiency(Dpress[i], self.s_array[i + 1] - self.s_array[i],
                                                    self.c_normal[i + 1] - self.c_normal[i])

        # we now determine the normal force coefficient due to the pressure lag
        self.c_normal_prime = self.c_normal - Dpress

        # based on this c_normal_prime, we determine a new equivalent angle of attack to determine
        # the onset of trailing edge separation
        alpha_f = self.c_normal_prime / self.dCn_dalpha + self.alpha0

        # use this equivalent angle of attack alpha_f to determine a new trailing edge separation point effect f_prime
        fprime = []
        for a_f in alpha_f:
            fprime.append(_f_trailing_edge_separation_point(a_f))

        def boundary_layer_lag_deficiency(Dbl_i, delta_s, delta_fprime, Tf=3.0):
            return Dbl_i * np.exp(-delta_s / Tf) + delta_fprime * np.exp(-delta_s / 2 / Tf)

        # we will now do the time marching to solve for the boundary layer lag deficiency function
        Dbl = np.zeros(np.shape(time))
        for i, val in enumerate(time[:-1]):
            Dbl[i + 1] = boundary_layer_lag_deficiency(
                Dbl_i=Dbl[i],
                delta_s=self.s_array[i + 1] - self.s_array[i],
                delta_fprime=fprime[i + 1] - fprime[i]
            )

        # we now determine the new expression of fprimeprime due to the boundary layer lag
        self.fprimeprime = fprime - Dbl

        self.c_normal_f = self.dCn_dalpha * (
                (1 + np.sqrt(self.fprimeprime)) / 2
        ) ** 2 * (
                self.alpha_equivalent - self.alpha0
        ) + self.c_normal_noncirc

    def leading_edge_vortex_shedding(self, time):
        def vortime_function(vortime_i, delta_s, delta_alphaqs, Cnormal_prime, CN1=1.0093):
            if Cnormal_prime > CN1:
                vortime_ip1 = vortime_i + 0.45 * delta_s
            else:
                if (delta_alphaqs < 0 and vortime_i > 0):
                    vortime_ip1 = vortime_i + 0.45 * delta_s
                else:
                    vortime_ip1 = 0

            return vortime_ip1

        # we will now do the time marching to solve for the non-dimensional vortex-time parameter vortime
        vortime = np.zeros(np.shape(time))
        for i, val in enumerate(time[:-1]):
            vortime[i + 1] = vortime_function(
                vortime[i], self.s_array[i + 1] - self.s_array[i], self.dalphaqs_dt[i], self.c_normal_prime[i])

        # define the function for decay of the cumulative normal force due to the presence of the leading edge vortex
        def leading_edge_vortex_normal_force(c_normal_vortex_i, delta_s, delta_c_vortex, vor_time, TVL=11, TV=6):
            if 0.001 < vor_time < TVL:
                c_normal_vortex_ip1 = (
                        c_normal_vortex_i * np.exp(-delta_s / TV) +
                        delta_c_vortex * np.exp(-delta_s / 2 / TV)
                )
            else:
                c_normal_vortex_ip1 = c_normal_vortex_i * np.exp(-delta_s / TV)
            return c_normal_vortex_ip1

        # We will now solve for the cumulative normal force due to the leading edge vortex by marching in time.
        c_vortex = self.c_normal_circ * (1 - (((1 + np.sqrt(self.fprimeprime)) / 2) ** 2))
        c_normal_vortex = np.zeros(np.shape(time))
        c_normal_vortex[0] = c_vortex[0]

        for i, val in enumerate(time[:-1]):
            c_normal_vortex[i + 1] = leading_edge_vortex_normal_force(
                c_normal_vortex_i=c_normal_vortex[i],
                delta_s=self.s_array[i + 1] - self.s_array[i],
                delta_c_vortex=c_vortex[i + 1] - c_vortex[i],
                vor_time=vortime[i]
            )

        self.c_normal_v = c_normal_vortex

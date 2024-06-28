import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


class DynamicStallModel:
    def __init__(
            self,
            twist_func: callable,
            chord_func: callable,
            span_blade,
            ratio_blade_start,
            n_blades,
            airfoil
    ):
        # Solver input vectors
        self.r_R, self.dr = np.linspace(ratio_blade_start, 1, 100, retstep=True)

        # Rotor properties
        self.span_blade = span_blade
        self.n_blades = n_blades

        # Local blade properties
        self.twist = twist_func(self.r_R)
        self.chord = chord_func(self.r_R)

        self.dCn_dalpha = 2 * np.pi  # lift slope
        self.alpha0 = 0  # alpha for which normal load is zero

        # Lift drag polars
        col_names = ['angle_of_attack', 'lift_coefficient', 'drag_coefficient', 'moment_coefficient']
        self.airfoil = pd.read_csv(f'./{airfoil}.csv', names=col_names)

        # Solver output vectors
        self.c_lift = None
        self.c_lift_circ = None
        self.c_lift_noncirc = None
        self.c_lift_prime = None
        self.c_lift_f = None
        self.c_lift_v = None

        self.alpha_qs = None
        self.dalphaqs_dt = None
        self.alpha_equivalent = None
        self.s_array = None
        self.fprimeprime = None

    @property
    def c_lift_total(self):
        print(f'c_lift {self.c_lift[-5:]}')
        print(f'c_lift_f {self.c_lift_f[-5:]}')
        print(f'c_lift_v {self.c_lift_v[-5:]}')

        return self.c_lift + self.c_lift_f + self.c_lift_v

    def run(self, time_start, time_end, u_inf_func: callable, omega:float, pitch, max_iter=200, ctol=1e-2, num_time=50):
        """Returns the cn at each position along the blade span for 0<r/R<1, given rotor operation conditions

        :param time:
        :param u_inf_func: the local inflow speed
        :param omega: rotational frequency of blade

        """

        # Input vectors
        time = np.linspace(time_start, time_end, num_time)

        # Result vectors
        Fn = np.zeros(shape=(len(self.chord), num_time))
        Ft = np.zeros(shape=(len(self.chord), num_time))

        # Perform Beddoes Leishman for each anulus
        for i, (chord, twist) in enumerate(zip(self.chord, self.twist)):
            r = self.r_R[i] * self.span_blade

            # Initialize local flow properties
            induction_axial = 0.3
            induction_tangential = 0

            for n in range(max_iter):
                u_normal_local = (1 - induction_axial) * u_inf_func(time)
                u_tangential_local = (1 + induction_tangential) * r * omega

                self.alpha_qs = np.arctan(u_normal_local / u_tangential_local) - twist - pitch

                self.dalphaqs_dt = np.gradient(self.alpha_qs,time)  # calculate the time derivative of the quasi-steady alpha
                self.s_array = 2 * u_normal_local * time / chord  # define the array semi-chord time scale

                self.unsteady_flat_plate(time=time, chord=chord, u_inf=u_normal_local)
                self.non_linear_trailing_edge_separation(time=time)
                self.leading_edge_vortex_shedding(time=time)

                cl_total = self.c_lift_total  # local, for each time step

                # Perform BEM to calculate induction factors per timestep
                induction_axial_new, induction_tangential_new, Fn_local, Ft_local = self.blade_element_momentum(
                    index=i,
                    alpha=self.alpha_qs,
                    cl_total=cl_total,
                    u_normal_local=u_inf_func(time),
                    omega=omega,
                    pitch=pitch
                )

                if (np.abs(induction_axial_new - induction_axial) < ctol).all():
                    Ft[i, :] = Ft_local
                    Fn[i, :] = Fn_local
                    break  # solution converged
                else:
                    induction_axial = 0.9 * induction_axial + 0.1 * induction_axial_new
                    induction_tangential = 0.9 * induction_tangential + 0.1 * induction_tangential_new

                if n == max_iter - 1:
                    print(f'ERROR: Max iterations reached but solution has not converged for anulus at r={r:.2f} m.')
                    # return False

        return Ft, Fn

    def blade_element_momentum(self, index, alpha, cl_total, u_normal_local, omega, pitch):
        # Constants
        n_blades = self.n_blades  # number of blades
        span_blade = self.span_blade  # rotor radius
        hubrad = 1.5  # hub radius
        rou = 1.225  # density of air
        EPS = 0.001  # iterative precision tolerance
        MAX_ITER = 1000

        # Airfoil data
        alphas = self.airfoil['angle_of_attack'].to_numpy()
        drag_coefficients = self.airfoil['drag_coefficient'].to_numpy()

        interp_Cd = interp1d(alphas, drag_coefficients, fill_value="extrapolate")

        # 1d arrays, indexed ( t)
        a = np.zeros(len(alpha))
        a_prime = np.zeros(len(alpha))
        CN = np.zeros(len(alpha))
        CT = np.zeros(len(alpha))

        r = self.r_R[index] * span_blade
        twist = self.twist[index]
        chord = self.chord[index]

        Sigma = chord * n_blades / (2 * np.pi * r)  # Solidity

        for j, (alpha, u0, cl) in enumerate(zip(alpha, u_normal_local, cl_total)):  # For each timestep, determine induction factors

            phi = alpha + twist + pitch

            # phi = np.arctan(((1 - a) * u0) / ((1 + a_prime) * r * omega))

            # Find Cl and Cd
            lift_coefficient = cl
            drag_coefficient = interp_Cd(np.degrees(alpha))

            # Projection in and out of plane
            Cn = lift_coefficient * np.cos(phi) + drag_coefficient * np.sin(phi)
            Ct = lift_coefficient * np.sin(phi) - drag_coefficient * np.cos(phi)

            # Prandtl loss
            f_tiploss = n_blades / 2 * (span_blade - r) / (r * np.sin(phi))
            F_tiploss = (2 / np.pi) * np.arccos(np.exp(-f_tiploss))
            f_hubloss = n_blades / 2 * (r - hubrad) / (r * np.sin(phi))
            F_hubloss = (2 / np.pi) * np.arccos(np.exp(-f_hubloss))
            F = np.nan_to_num(F_tiploss * F_hubloss)

            # Determine a(t) and a'(t) at this location
            a[j] = 1 / (4 * F * (np.sin(phi)) ** 2 / (Sigma * Cn) + 1)
            a_prime[j] = 1 / (4 * F * np.sin(phi) * np.cos(phi) / (Sigma * Ct) - 1)
            CN[j] = Cn
            CT[j] = Ct

            # Glauert correction
            ac = 0.2
            if a[j] > ac:
                K = 4 * F * np.sin(phi) ** 2 / (Sigma * Cn)
                a[j] = 0.5 * (2 + K * (1 - 2 * ac) - np.sqrt((K * (1 - 2 * ac) + 2) ** 2 + 4 * (K * ac ** 2 - 1)))

            # Force in two directions and bending moment, as functions of time
        FN = (0.5 * rou * ((r * omega * (1 + a_prime)) ** 2 + (u_normal_local * (1 - a)) ** 2) * chord * CN * self.dr)
        FT = (0.5 * rou * ((r * omega * (1 + a_prime)) ** 2 + (u_normal_local * (1 - a)) ** 2) * chord * CT * self.dr)

        # TODO calculate loads
        # Rx = self.r_R * self.span_blade
        # Ff = np.cos(self.twist + pitch) * FN + np.sin(self.twist + pitch) * FT
        # Fe = np.cos(self.twist + pitch) * FT - np.sin(self.twist + pitch) * FN
        # MT = FT * Rx
        # MN = FN * Rx
        #
        # Mt = sum(MT)  # Rotor torque
        # Mn = sum(MN)  # Root moment
        # P = Mt * omega * 3 * 0.944  # Power calculation
        return a, a_prime, FN, FT

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
        self.c_lift_circ = _circulatory_normal_force()
        self.c_lift_noncirc = _non_circulatory_normal_force()
        self.c_lift = self.c_lift_circ + self.c_lift_noncirc

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
                                                    self.c_lift[i + 1] - self.c_lift[i])

        # we now determine the normal force coefficient due to the pressure lag
        self.c_lift_prime = self.c_lift - Dpress

        # based on this c_normal_prime, we determine a new equivalent angle of attack to determine
        # the onset of trailing edge separation
        alpha_f = self.c_lift_prime / self.dCn_dalpha + self.alpha0

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

        self.c_lift_f = np.nan_to_num(self.dCn_dalpha * (
                (1 + np.sqrt(self.fprimeprime)) / 2
        ) ** 2 * (
                self.alpha_equivalent - self.alpha0
        ) + self.c_lift_noncirc)

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
                vortime[i], self.s_array[i + 1] - self.s_array[i], self.dalphaqs_dt[i], self.c_lift_prime[i])

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
        c_vortex = self.c_lift_circ * (1 - (((1 + np.sqrt(self.fprimeprime)) / 2) ** 2))
        c_normal_vortex = np.zeros(np.shape(time))
        c_normal_vortex[0] = c_vortex[0]

        for i, val in enumerate(time[:-1]):
            c_normal_vortex[i + 1] = leading_edge_vortex_normal_force(
                c_normal_vortex_i=c_normal_vortex[i],
                delta_s=self.s_array[i + 1] - self.s_array[i],
                delta_c_vortex=c_vortex[i + 1] - c_vortex[i],
                vor_time=vortime[i]
            )

        self.c_lift_v = c_normal_vortex

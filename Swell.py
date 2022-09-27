import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import quad

g = 9.80665  # m/s2
alpha = 0.0081


class Swell:
    """A simple class to define the swell parameters, for the JONSWAP repartition, and the coefficient for the function
    of the wave repartition .
        :attributes
        ------------
        significant_wave_height: a float
            It is a float representing the maximum height of the waves (in meter)
        gamma: a float
            It is the gamma coefficient in the JONSWAP formula
        speed: a float
            It is the ship speed in the swell (in m/s)
        coeff_wave: a float
            It is the coefficient for the function of the angular repartition. The function is
            cosinus(angle)**(2*coeff_wave)
        depth: a float
            The depth of the water for the swell case (in meter)
        theta_mean: a float
            The direction of the swell in degrees. It represents the incident angle for which the wave will influence the
            most, the ship.
        """

    def __init__(self, significant_wave_height: float, gamma: float, speed: float, coeff_wave: int, depth: float,
                 theta_mean: float):
        self.significant_wave_height = significant_wave_height  # in meter
        self.gamma = gamma
        self.speed = speed  # in m/s
        self.coeff_wave = coeff_wave
        self.depth = depth  # in meter #the depth
        self.theta_mean = theta_mean  # in degree

    def JONSWAP(self, f: float):
        """This function permits to return the JONSWAP formula for the swell object and a frequency value.

        :argument
        -----------
        f: a float
            It is the value of the frequency (in Hertz) for which we want to compute the JONSWAP value

        :returns
        ---------
        JONSWAP: a float
            The value of the JONSWAP value for self and for the frequency f"""

        tau = self.tau(f)
        fp = self.fp()  # in Hertz
        gamma = self.gamma
        el1 = alpha * (g ** 2) * ((2 * np.pi) ** (-4)) * (f ** (-5))
        el2 = np.exp((-5 / 4) * ((f / fp) ** (-4)))
        el3 = gamma ** (np.exp(-(((f - fp) ** 2) / (2 * (tau ** 2) * (fp ** 2)))))
        # checked with EXCEL
        return el1 * el2 * el3

    def tau(self, f: float):
        """It is the tau function included in the JONSWAP computation.

        :argument
        ---------
        f: a float
            the frequency value of the wave we are interested in (a value in Hertz)

        :returns
        ----------
        val: 0.07 or 0.09
            If f<=fp returns 0.07 f>fp returns 0.09"""
        fp = self.fp()
        if f <= fp:
            return 0.07
        if f > fp:
            return 0.09

    def fp(self):
        """That function computes the peak frequency value, for the swell case.

        :argument
        ---------
        None

        :returns
        ---------
        fp: a float
            The frequency peak value (in Hertz)
        """
        Hs = self.significant_wave_height
        gamma = self.gamma
        # Checked with Excel
        return (1/2/np.pi) * np.sqrt(((1.555 + 0.2596 * gamma - 0.02231 * (gamma ** 2) + 0.001142 * (gamma ** 3)) * g * np.sqrt(alpha)) / Hs)

    def plot_JONSWAP(self, f_start: float, f_end: float, step: float):
        """That function permits to plot the graph of the JONSWAP, for the frequencies between f_start and f_end, with
        a step for the frequencies.

        :argument
        ----------
        f_start: a float
            The first frequency value for the plot, the frequency is a value in rad/s
        f_end: a float
            The last frequency value for the plot in rad/s
        step: a float
            The step between each frequency value in rad/s

        :returns
        ----------
        It prints the graph of the JONSWAP function"""
        les_f = np.arange(f_start, f_end, step)/2/np.pi
        les_y = []
        for f in les_f:
            y = self.JONSWAP(f)
            les_y.append(y)
        plt.plot(les_f, les_y)
        plt.xlabel("Frequencies in Hz")
        plt.ylabel("JONSWAP function")
        plt.show()

    def spreading_func(self, theta: float):
        """That functions permits to know the value of the spreading function a cos 2s function for a certain angle
        theta

        :argument
        ----------
        theta: a float
            the angle in degree
        :returns
        ------------
        spreading function: a float
            the value of the spreading function for the angle theta
        """
        s = self.coeff_wave
        gam1 = math.gamma(s + 1)
        gam2 = math.gamma(s + 0.5)
        el4 = (gam1 * np.cos((theta-self.theta_mean) * np.pi / 180 / 2) ** (2 * s)) / (2 * np.sqrt(np.pi) * gam2)
        return el4

    def spread_func_int(self,theta,step):
        """That functions permits to know the value of the spreading function a cos 2s function for a certain angle
        theta

        :argument
        ----------
        theta: a float
            the angle in degree
        :returns
        ------------
        spreading function: a float
            the value of the spreading function for the angle theta
        """
        s = self.coeff_wave
        gam1 = math.gamma(s + 1)
        gam2 = math.gamma(s + 0.5)
        int=quad(self.spreading_func,theta-step/2,theta+step/2)
        return int[0]*np.pi/180#*(gam1 * np.cos((theta) * np.pi / 180 / 2) ** (2 * s)) / (2 * np.sqrt(np.pi) * gam2)

    def wave_spectrum(self, theta: float, f: float):
        """That function for a certain wave with an angle theta in degree and a frequency f in Hertz, returns the
        value of the wave spectrum

        :argument
        ---------
        theta: a float
            The angle of the incident wave in degrees
        f: afloat
            The frequency of the incident wave in Hertz

        :returns
        ---------
        wave spectrum: a float
            The value represents the importance of the wave
        """
        return self.JONSWAP(f) * self.spreading_func(theta)

    def plot_spreading(self, theta1: float, theta2: float, step: float):
        """That function permits to print the spreading function in function of the angle.

        :argument
        ---------
        theta1: a float
            the first angle (in degree) of the plot
        theta2: a float
            the last angle (in degree) of the plot
        step: a float
            the gap between each point
        :returns
        ---------
        prints the graph of the spreading function"""
        les_theta = np.arange(theta1, theta2, step)
        les_y = []
        for theta in les_theta:
            y = self.spreading_func(theta)
            les_y.append(y)
        plt.plot(les_theta, les_y)
        plt.xlabel("Angles in degree")
        plt.ylabel("Spreading function")
        plt.show()

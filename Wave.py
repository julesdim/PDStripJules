import numpy as np
from numpy import inf

g = 9.80665  # in m/s2


class Wave:
    """A simple class to define a wave. That is useful because the loading uses different wave, with different
    frequencies and different angle.

    :attributes
    ------------
    wave_length: a float
        The length of a wave period (in m)
    angle: a float
        The angle measured between the x axis and the wave direction (in degree)
    speed: a float
        The speed of the boat (in m/s)
    k: a float
        The wave number m-1
    freq: a float
        This is the frequency of the wave (in Hertz)
    freq_e: a float
        This is the wave frequency encountered by the ship (in Hertz)
    depth: a float
        This is the depth of the water (in m)
    """

    def __init__(self, wave_length, angle, speed, depth):
        self.wave_length = wave_length
        self.angle = angle
        self.speed = speed
        self.k = 2 * np.pi / wave_length
        # check that coeff with Excel
        self.freq = np.sqrt(g * self.k * np.tanh(self.k * depth)) / (2 * np.pi)  # it is omega/2*pi, in s**(-1)
        self.freq_e = self.freq - self.k * speed * np.cos(np.pi * angle / 180)
        self.depth = depth

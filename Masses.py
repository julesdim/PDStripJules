import numpy as np
from scipy.integrate import quad


class Mass:
    """A simple class to define a weight by the weight, the beginning, the end, the center of gravity along all the axis
    the linear density at the beginning and at the end, it permits to compute different value of the weight if necessary

    :attributes
    ------------
    weight: a float
        the weight of the mass (in ton)
    x_start: a float
        x_coordinate of the mass start (in m)
    x_end: a float
        x_coordinate of the mass end (in m)
    x_coordinate_CoG: a float
        x coordinate of the center of gravity (in m)
    y_coordinate_CoG: a float
        y coordinate of the center of gravity (in m)
    z_coordinate_CoG: a float
        z coordinate of the center of gravity (in m)
    linear_density_start: a float
        linear density of the mass at x_start (in ton/m)
    linear_density_end: a float
        linear density of the mass at x_end (in ton/m)
    """

    def __init__(self, weight: float, x_start: float, x_end: float, x_coordinate_CoG: float, y_coordinate_CoG: float,
                 z_coordinate_CoG: float, linear_density_start: float, linear_density_end: float):
        self.weight = weight
        self.x_start = x_start
        self.x_end = x_end
        self.x_coordinate_CoG = x_coordinate_CoG
        self.y_coordinate_CoG = y_coordinate_CoG
        self.z_coordinate_CoG = z_coordinate_CoG
        self.linear_density_start = linear_density_start
        self.linear_density_end = linear_density_end

    def calcul_mass(self, x_start_section: float, x_end_section: float):
        """That function permits to the user to know what is the exact mass of the mass item between two coordinates,
        that is very useful for the input values of the ship.

        :parameter
        -----------
        x_start_section: a float
            It is the beginning of the section where we want to compute the mass (in m)
        x_end_section: a float
         It is the end of the section where the mass will be computed (in m)
        :returns
        ----------
        mass_section: a float
            The exact value of the mass between the coordinates (in ton)"""
        if x_start_section >= self.x_end or x_end_section <= self.x_start:
            return 0
        else:
            real_start = max([x_start_section, self.x_start])
            real_end = min([x_end_section, self.x_end])
            # we know the exact coordinates of the mass in the section
            if real_start < self.x_start:
                print("pb")
            if real_end > self.x_end:
                print("pb")
            if real_start!= self.x_start or real_end!=self.x_end:
                return quad(self.linear_density, real_start, real_end)[0]
            if real_start==self.x_start and real_end==self.x_end:
                return self.weight

    def calcul_xg_mass(self):
        """That function permits to the user to compute X_CoG in function of linear density when the program calcul the
        new formula.

        :argument
        ----------
        self: a mass object
            the current mass
        :returns
        ----------
        X_CoG: a float
            the current x_coordinate of the center of gravity (in m)"""
        x_end = self.x_end
        x_start = self.x_start
        return self.calcul_xg_coordinates(x_start, x_end)

    def calcul_linear_density_for_coordinates(self, x_start_part: float, x_end_part: float):
        """That function permits to know the value of the linear density at the start part and at the end part to
        compute the linear density.

        :argument
        ----------
        x_start_part: a float
            x coordinate of the section start (in m)
        x_end_part: a float
            x coordinate of the section end (in m)

        :returns
        -----------
        linear_density_at_x_start_part: a float
            the value for the x coordinate x_start_part (in ton/m)
        linear_density_at_x_end_part: a float
            the value for the x coordinate x_end_part (in ton/m)
        """
        x_end_of_the_mass = self.x_end
        x_start_of_the_mass = self.x_start
        linear_density_start = self.linear_density_start
        linear_density_end = self.linear_density_end
        coefficient_of_the_line = (linear_density_end - linear_density_start) / (x_end_of_the_mass -
                                                                                 x_start_of_the_mass)
        linear_density_at_x_start_part = linear_density_start + coefficient_of_the_line * (x_start_part -
                                                                                           x_start_of_the_mass)
        linear_density_at_x_end_part = linear_density_start + coefficient_of_the_line * (
                x_end_part - x_start_of_the_mass)
        return linear_density_at_x_start_part, linear_density_at_x_end_part

    def calcul_xg_not_the_mid(self, eps: float = 1e-12):
        """That functions permits to compute the new value of linear density, indeed at the beginning the linear density
        is defined equal at the beginning and at the end.

        :argument
        -----------
        eps: a float
            for the precision of the computation

        :returns
        ----------
        Nothing
        """
        weight = self.weight
        x_coordinate_CoG = self.x_coordinate_CoG
        x_end_mass = self.x_end
        x_start_mass = self.x_start
        linear_density_start = weight / (x_end_mass - x_start_mass)
        gap = x_end_mass - x_start_mass
        coefficient = 0
        i = 0
        while gap > eps and i < 50000:
            X_CoG_temp = self.calcul_xg_mass()
            ratio = abs(X_CoG_temp - x_coordinate_CoG) / (x_end_mass - x_start_mass)
            if X_CoG_temp > x_coordinate_CoG:
                coefficient -= ratio
            if X_CoG_temp < x_coordinate_CoG:
                coefficient += ratio
            self.linear_density_start = (1 - coefficient) * linear_density_start
            self.linear_density_end = (1 + coefficient) * linear_density_start
            gap = abs(X_CoG_temp - x_coordinate_CoG)
            i += 1
        self.linear_density_start = round(self.linear_density_start, 11)
        self.linear_density_end = round(self.linear_density_end, 11)

    def linear_density(self, x: float):
        """That function permits to know the exact value of the linear density for a certain coordinate.

        :argument
        ----------
        x: a float
            The coordinate where the linear density will be computed (in m)

        :returns
        ---------
        linear_density: a float
            The linear density at the coordinate x """
        x_start = self.x_start
        x_end = self.x_end
        if x > x_end or x < x_start:
            return 0
        elif x_start <= x <= x_end:
            linear_density = (self.linear_density_end - self.linear_density_start) * (x - x_start) / (x_end - x_start) + \
                             self.linear_density_start  # in ton/m
            return linear_density

    def linear_densityx(self, x: float):
        """That function returns the linear density, multiplied by the x coordinates. That is useful for the
        computation of the center of gravity.

        :argument
        ----------
        x: a float
            The x coordinate where the linear density will be computed (in m)

        :returns
        ----------
        linear_density*x: a float
            The linear density at x multiplied by x in ton.m/m"""

        return self.linear_density(x) * x

    def calcul_xg_coordinates(self, x_start: float, x_end: float):
        """That function permits to compute the x coordinate of the center of gravity for a mass between coordinates,
        in function of its linear density.

        :argument
        --------------
        x_start: a float
            The start of the section where xg will be computed (in m)
        x_end: a float
            The end of the section where the xg will be computed (in m)

        :returns
        ----------
        X_CoG: a float
            The x coordinate of the center of gravity (in m)"""

        int1, err1 = quad(self.linear_densityx, x_start, x_end)
        mass = self.calcul_mass(x_start, x_end)
        return int1 / mass

import numpy as np
import matplotlib.pyplot as plt
import Masses as mass


class Loading:
    """This Class is very simple because it's a list of mass, indeed the loading is nothing more

    :argument
    -----------
    masses: a list
        it's a list of mass object"""

    def __init__(self):
        self.masses = []

    def collection_of_mass(self, filename: str):
        """That functions reads a csv file and collects all the information, the weight in ton, the beginning of the weight,
        the end of the weight, the center of gravity with x coordinate, y and z, then it defines the weight per meter at the
        beginning of the weight repartition and the weight per meter at the end of the weight repartition (along the x axis), the
        code uses the localisation of the center of gravity to compute that values.

        :parameter
        -----------
        filename: a csv file
            every line are : x coordinate of the beginning, x coordinate of the end, x coordinate of the center of gravity
            y coordinate of CoG, z coordinate of the CoG

        :returns
        ----------
        list_of_masses: a loading file
            it is a list of weight, defined by the information mentioned above but with weight per
            meter computed at the beginning and at the end
         """
        file = open(filename, "rt")
        the_lines = file.readlines()
        total_mass = 0  # to check the total weight
        for line in the_lines:
            line_formatted = line.strip("\n").split(
                ";")  # we stop the line to the \n and we cut the information where there is a ";"
            current_weight = float(line_formatted[0])  # the first info is the object weight (in ton)
            total_mass += current_weight
            x_beginning = float(line_formatted[1])  # the second is the beginning (in m)
            x_end = float(line_formatted[2])  # the end (in m)
            x_coordinate_CoG = float(line_formatted[3])  # the exact center of gravity (in m)
            y_coordinate_CoG = float(line_formatted[4])  # the turning radius (in m)
            z_coordinate_CoG = float(line_formatted[5])  # the position along z axis of the center of gravity (in m)
            linear_density_x_beginning, linear_density_x_end = current_weight / (
                    x_end - x_beginning), current_weight / (x_end - x_beginning)  # in ton/m for the initialization
            if current_weight != 0:
                current_mass = mass.Mass(current_weight, x_beginning, x_end, x_coordinate_CoG, y_coordinate_CoG,
                                         z_coordinate_CoG, linear_density_x_beginning, linear_density_x_end)
                if x_coordinate_CoG != (x_beginning + x_end) / 2:
                    current_mass.calcul_xg_not_the_mid()  # we compute the good linear density of the start and the end
                    print(current_mass.x_coordinate_CoG)
                if x_coordinate_CoG == (x_beginning + x_end) / 2:
                    print(current_mass.x_coordinate_CoG)
                self.__append__(current_mass)
        print(total_mass, "TOTAL")

    def plot_loading(self, x_start: float, x_end: float):
        """That function permits to print the weight loading, for a list of weight, and with the
        boundaries of the ship.

        :argument
        -----------
        x_start: beginning of the ship part where the loading will be plotted (in m)

        x_end: the end of the ship part (in m)

        :returns
        -----------
        it prints a graph with the loading along the x-axis"""

        delt_x = 0.1  # value of the strip to calculate the ship loading
        list_x_coordinates = np.arange(x_start, x_end + delt_x, delt_x)  # coordinates of each strip
        total = len(list_x_coordinates)
        linear_density = [0]  # we initialize the weight loading for x=0, it's equal to 0
        for i in range(total - 1):
            x_inf = list_x_coordinates[i]  # in m
            x_up = list_x_coordinates[i + 1]  # in m
            # coordinates of the strip i
            mass_element = self.mass_calculation_for_coordinates(x_inf,
                                                                 x_up) / delt_x  # element of weight for the strip
            linear_density.append(mass_element)  # we add to the list
        plt.plot(list_x_coordinates, linear_density)
        plt.ylabel("Linear density in ton/m")
        plt.xlabel("Position along the x-axis in m")
        plt.title("Linear density in ton per meter along the x axis")
        plt.show()
        return

    def __append__(self, mass: float):
        """That functions append a new element to the loading.

        :argument
        -----------
        mass: a mass object
            the new element of the loading"""

        self.masses.append(mass)

    def mass_calculation_for_coordinates(self, x_start: float, x_end: float):
        """That function computes the value of the mass between two coordinates of a section.

        :argument
        -----------
        x_start: x coordinate of the start of the section (in m)
        x_end: x coordinate of the end of the section (in m)

        :returns
        -----------
        total_mass: the total mass of the section (in ton)"""

        number_of_masses = len(self.masses)
        total_mass = 0
        for i in range(number_of_masses):
            x_start_mass = self.masses[i].x_start
            x_end_mass = self.masses[i].x_end
            x_CoG_mass = self.masses[i].x_coordinate_CoG
            if x_start_mass < x_end and x_end_mass > x_start:
                real_start = np.max([x_start, x_start_mass])  # real beginning of the weight for the section
                real_end = np.min([x_end, x_end_mass])
                if real_start < x_start_mass:
                    print("pb")
                if real_end > x_end_mass:
                    print("pb")
                if x_CoG_mass != (x_start_mass + x_end_mass) / 2:
                    total_mass += self.masses[i].calcul_mass(real_start, real_end)
                if x_CoG_mass == (x_start_mass + x_end_mass) / 2:
                    # real end of the weight for the section, if the end of the weight is after the end of the section
                    total_mass += (real_end - real_start) * self.masses[i].linear_density_start  # in ton
        return total_mass

    def calcul_center_of_gravity_for_coordinates(self, x_start: float, x_end: float):
        """That function permits to know the center of gravity for a specific section.

        :argument
        ----------
        x_start: a float
            x coordinate of the section (in m)
        x_end: a float
            x coordinate of the section (in m)

        :returns
        -------------
        X_CoG: a float
            x coordinate of the center of gravity of the section (in m)
        Y_CoG: a float
            y coordinate of the center of gravity of the section (in m)
        Z_CoG: a float
            z coordinate of the center of gravity of the section (in m)"""

        total_mass = self.mass_calculation_for_coordinates(x_start, x_end)
        number_of_masses = len(self.masses)
        # initialization of the values
        X_CoG_glob = 0
        Y_CoG_glob = 0
        Z_CoG_glob = 0
        for i in range(number_of_masses):
            x_start_mass = self.masses[i].x_start  # beginning of the weight
            x_end_mass = self.masses[i].x_end  # end of the weight
            if x_start_mass < x_end and x_end_mass > x_start:
                real_start = np.max(
                    [x_start, x_start_mass])  # real beginning of the weight, if the weight begins before the frame
                real_end = np.min([x_end, x_end_mass])
                if real_start < x_start_mass:
                    print("pb")
                if real_end > x_end_mass:
                    print("pb")
                real_mass = self.masses[i].calcul_mass(real_start, real_end)
                # proportion of the weight situated between the section
                if real_start!=x_start_mass or real_end!=x_end_mass:
                    X_CoG = self.masses[i].calcul_xg_coordinates(real_start, real_end)
                if real_start==x_start_mass and real_end==x_end_mass:
                    X_CoG = self.masses[i].x_coordinate_CoG
                X_CoG_glob += real_mass * X_CoG
                Y_CoG_glob += real_mass * self.masses[i].y_coordinate_CoG
                Z_CoG_glob += real_mass * self.masses[i].z_coordinate_CoG
        return X_CoG_glob / total_mass, Y_CoG_glob / total_mass, Z_CoG_glob / total_mass

    def pdstrip_coordinates_conversion(self, midship: float):
        """That function converts every coordinate depending of the PD strip coordinate system.

        :argument
        ------------
        midship: a float
            the length of the midship, defined by the length between perpendiculars divide by 2 (in m)

        :returns
        -------------
        It changes the x coordinates by x-midship (the pd strip system the pd strip system places the origin in
        the middle of the boat)."""

        for mass in self.masses:
            mass.x_start = mass.x_start - midship
            mass.x_end = mass.x_end - midship
            mass.x_coordinate_CoG = mass.x_coordinate_CoG - midship

    def calcul_inertia_x(self, Y_Cog: float, Z_Cog: float, x_start: float, x_end: float):
        """ That functions permits to compute the square of the inertial radius relating to the axis through the center
        of gravity parallel to the y-axis for exact coordinates of the ship part. That will use each weight item in the
        ship.

        :parameter
        -----------
        Y_Cog: a float
            The y coordinate of the center of gravity (in m)
        Z_Cog: a float
            The z coordinate of the center of gravity (in m)
        x_start: a float
            The x coordinate of the beginning of the part where the user wants to compute the value (in m)
        x_end: a float
            The x coordinate of the end of the part where the user wants to compute the value (in m)

        :returns
        ----------
        rx2: a float
            The square of the inertial radius relating to the axis through the center of gravity parallel to the x-axis.
            (in m2)
        """

        sum = 0
        total_mass = self.mass_calculation_for_coordinates(x_start, x_end)
        if total_mass == 0:
            return 0
        n_loading = len(self.masses)
        for i in range(n_loading):
            mass = self.masses[i].calcul_mass(x_start, x_end)
            if mass != 0:
                y = self.masses[i].y_coordinate_CoG
                z = self.masses[i].z_coordinate_CoG
                sum += ((y - Y_Cog) ** 2 + (z - Z_Cog) ** 2) * mass
        return sum / total_mass

    def calcul_inertia_y(self, X_Cog: float, Z_Cog: float, x_start: float, x_end: float):
        """ That functions permits to compute the square of the inertial radius relating to the axis parallel to the
        y-axis for exact coordinates of the ship part. That will use each weight item in the ship.

                :parameter
                -----------
                X_Cog: a float
                    The x coordinate of the center of gravity (in m)
                Z_Cog: a float
                    The z coordinate of the center of gravity (in m)
                x_start: a float
                    The x coordinate of the beginning of the part where the user wants to compute the value (in m)
                x_end: a float
                    The x coordinate of the end of the part where the user wants to compute the value (in m)

                :returns
                ----------
                ry2: a float
                    The square of the inertial radius relating to the axis through the center of gravity parallel to the
                    y-axis. (in m2)
                """

        sum = 0
        total_mass = self.mass_calculation_for_coordinates(x_start, x_end)
        if total_mass == 0:
            return 0
        n_loading = len(self.masses)
        for i in range(n_loading):
            mass = self.masses[i].calcul_mass(x_start, x_end)
            if mass != 0:
                x = self.masses[i].calcul_xg_coordinates(x_start, x_end)
                z = self.masses[i].z_coordinate_CoG
                sum += ((x - X_Cog) ** 2 + (z - Z_Cog) ** 2) * (mass / total_mass)
        return sum

    def calcul_inertia_z(self, X_Cog: float, Y_Cog: float, x_start: float, x_end: float):
        """ That functions permits to compute the square of the inertial radius relating to the axis parallel to the
                z-axis for exact coordinates of the ship part. That will use each weight item in the ship.

                :parameter
                -----------
                X_Cog: a float
                    The x coordinate of the center of gravity (in m)
                Y_Cog: a float
                    The z coordinate of the center of gravity (in m)
                x_start: a float
                    The x coordinate of the beginning of the part where the user wants to compute the value (in m)
                x_end: a float
                    The x coordinate of the end of the part where the user wants to compute the value (in m)

                :returns
                ----------
                rz2: a float
                    The square of the inertial radius relating to the axis through the center of gravity parallel to the
                    z-axis. (in m2)
                        """

        sum = 0
        total_mass = self.mass_calculation_for_coordinates(x_start, x_end)
        if total_mass == 0:
            return 0
        n_loading = len(self.masses)
        for i in range(n_loading):
            mass = self.masses[i].calcul_mass(x_start, x_end)
            if mass != 0:
                x = self.masses[i].calcul_xg_coordinates(x_start, x_end)
                y = self.masses[i].y_coordinate_CoG
                sum += ((x - X_Cog) ** 2 + (y - Y_Cog) ** 2) * (mass / total_mass)
        return sum

    def calcul_inertia_xy(self, X_Cog: float, Y_Cog: float, x_start: float, x_end: float):
        """ That functions permits to compute the mass weighted average of (x-xG)(y-yG).

                :parameter
                -----------
                X_Cog: a float
                    The x coordinate of the center of gravity (in m)
                Y_Cog: a float
                    The y coordinate of the center of gravity (in m)
                x_start: a float
                    The x coordinate of the beginning of the part where the user wants to compute the value (in m)
                x_end: a float
                    The x coordinate of the end of the part where the user wants to compute the value (in m)

                :returns
                ----------
                (x-xg)(y-yg)_av: a float
                    The mass weighted average of (x-xG)(y-yG), used to determine the mixed moment of inertia,
                    Ixy=m*average((x-xG)(y-yG)) (in m2)
                """
        sum = 0
        total_mass = self.mass_calculation_for_coordinates(x_start, x_end)
        if total_mass == 0:
            return 0
        n_loading = len(self.masses)
        for i in range(n_loading):
            mass = self.masses[i].calcul_mass(x_start, x_end)
            if mass != 0:
                x = self.masses[i].calcul_xg_coordinates(x_start, x_end)
                y = self.masses[i].y_coordinate_CoG
                sum += ((x - X_Cog) * (y - Y_Cog)) * (mass / total_mass)
        return sum

    def calcul_inertia_yz(self, Y_Cog: float, Z_Cog: float, x_start: float, x_end: float):
        """ That functions permits to compute the mass weighted average of (y-yG)(z-zG).

                :parameter
                -----------
                Y_Cog: a float
                    The y coordinate of the center of gravity (in m)
                Z_Cog: a float
                    The z coordinate of the center of gravity (in m)
                x_start: a float
                    The x coordinate of the beginning of the part where the user wants to compute the value (in m)
                x_end: a float
                    The x coordinate of the end of the part where the user wants to compute the value (in m)

                :returns
                ----------
                (y-yg)(z-zg)_av: a float
                    The mass weighted average of (y-yG)(z-zG), used to determine the mixed moment of inertia,
                    Iyz=m*average((y-yG)(z-zG)) (in m2)
                """
        sum = 0
        total_mass = self.mass_calculation_for_coordinates(x_start, x_end)
        if total_mass == 0:
            return 0
        n_loading = len(self.masses)
        for i in range(n_loading):
            mass = self.masses[i].calcul_mass(x_start, x_end)
            if mass != 0:
                z = self.masses[i].z_coordinate_CoG
                y = self.masses[i].y_coordinate_CoG
                sum += ((y - Y_Cog) * (z - Z_Cog)) * (mass / total_mass)
        return sum

    def calcul_inertia_xz(self, X_Cog: float, Z_Cog: float, x_start: float, x_end: float):
        """ That functions permits to compute the mass weighted average of (x-xG)(z-zG).

                :parameter
                -----------
                X_Cog: a float
                    The x coordinate of the center of gravity (in m)
                Z_Cog: a float
                    The z coordinate of the center of gravity (in m)
                x_start: a float
                    The x coordinate of the beginning of the part where the user wants to compute the value (in m)
                x_end: a float
                    The x coordinate of the end of the part where the user wants to compute the value (in m)

                :returns
                ----------
                (x-xg)(z-zg)_av: a float
                    The mass weighted average of (x-xG)(z-zG), used to determine the mixed moment of inertia,
                    Ixz=m*average((x-xG)(z-zG)) (in m2)
                """

        sum = 0
        total_mass = self.mass_calculation_for_coordinates(x_start, x_end)
        if total_mass == 0:
            return 0
        n_loading = len(self.masses)
        for i in range(n_loading):
            mass = self.masses[i].calcul_mass(x_start, x_end)
            if mass != 0:
                x = self.masses[i].calcul_xg_coordinates(x_start, x_end)
                z = self.masses[i].z_coordinate_CoG
                sum += ((x - X_Cog) * (z - Z_Cog)) * (mass / total_mass)
        return sum

    def calcul_every_parameters(self, X_Cog: float, Y_Cog: float, Z_Cog: float, x_start: float, x_end: float):
        """That function permits to return every parameter needed by PD Strip for the inertia information.

        :parameter
        ----------
        X_Cog: a float
            The x coordinate of the center of gravity (in m)
        Y_Cog: a float
            The y coordinate of the center of gravity (in m)
        Z_Cog: a float
            The z coordinate of the center of gravity (in m)
        x_start: a float
            X coordinate of the beginning of the ship part where the information will be computed (in m)
        x_end: a float
            X coordinate of the end of the ship where the information will be computed (in m)

        :returns
        ----------
        rx2: a float
            The square of the inertial radius relating to the axis through the center of gravity parallel to the x-axis.
            (in m2)
        ry2: a float
                    The square of the inertial radius relating to the axis through the center of gravity parallel to the
                    y-axis. (in m2)
        rz2: a float
                    The square of the inertial radius relating to the axis through the center of gravity parallel to the
                    z-axis. (in m2)
        (x-xg)(y-yg)_av: a float
                    The mass weighted average of (x-xG)(y-yG), used to determine the mixed moment of inertia,
                    Ixy=m*average((x-xG)(y-yG)) (in m2)
        (y-yg)(z-zg)_av: a float
                    The mass weighted average of (y-yG)(z-zG), used to determine the mixed moment of inertia,
                    Iyz=m*average((y-yG)(z-zG)) (in m2)
        (x-xg)(z-zg)_av: a float
                    The mass weighted average of (x-xG)(z-zG), used to determine the mixed moment of inertia,
                    Ixz=m*average((x-xG)(z-zG)) (in m2)
        """
        return (
            self.calcul_inertia_x(Y_Cog, Z_Cog, x_start, x_end), self.calcul_inertia_y(X_Cog, Z_Cog, x_start, x_end),
            self.calcul_inertia_z(X_Cog, Y_Cog, x_start, x_end), self.calcul_inertia_xy(X_Cog, Y_Cog, x_start, x_end),
            self.calcul_inertia_yz(Y_Cog, Z_Cog, x_start, x_end), self.calcul_inertia_xz(X_Cog, Z_Cog, x_start, x_end))

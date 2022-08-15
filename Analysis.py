import csv
import Wave as Wv
import Hull as Hll
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import pandas as pd

class Analysis:
    """A simple class to define a file to analyse with the functions defined.

    :attributes
    ------------
    filename: a str
        The filename is the PD Strip result file which will be analysed
        Normally the filename ends with ".out.ok"
    """

    def __init__(self, filename: str):
        self.filename = filename

    def analysis_res(self, hull: Hll.Hull, distance_from_neutral_axis: float):
        """That function can make the strength analysis with the information of the hull. It returns at the end the
        forces and the moments along the x-axis.

        :argument
        ---------
        hull: a Hull.Hull object
            It is the hull object that represents the repartition of the hull
        distance_from_neutral_axis: a float
            distance between the baseline and the neutral axis for the computation of the bending moment (in m)
        :returns
        --------
        the_forces: an np.array
            That is a list with the strength in the ship at each x coordinate (in ton)
        the_moments: a np.array
            It is a list with the moments in the ship at each x coordinates of intersections (in ton.m)
        x_coordinates: a list
            It is a list with all the x_coordinates in the PD Strip coordinate system. That is the coordinates (in m)
            associated to the strength and the moments
        """
        x_coordinates = []
        filename = self.filename
        file = open(filename, "r", encoding="utf-8")
        the_lines = csv.reader(file)
        speed = hull.speed  # m/s
        significant_wave_height = hull.significant_wave_height  # meter
        gamma = hull.gamma
        coeff = hull.coeff_wave
        deep = hull.deep  # meter
        for line in the_lines:
            try:
                line_formatted = line[0].strip().split()
            except IndexError:
                line_formatted = ["error", "error", "error", "error", "error"]
                pass
            if len(line_formatted) == 6:
                if line_formatted[4] == "x=":
                    x_coordinates.append(float(line_formatted[5]))
        intersections = []
        for i in range(len(x_coordinates) - 1):
            intersections.append((x_coordinates[i] + x_coordinates[i + 1]) / 2)
        # saving the coordinates where the forces are computed
        the_forces = np.zeros(len(intersections))
        the_moments = np.zeros(len(intersections))
        example_wave_angle = 0
        example_wave_length = 0
        example_wave_speed = 0
        end_wave = 0
        file = open(filename, "r", encoding="utf-8")
        the_lines = csv.reader(file)
        bool_complete_vag = False
        count_line = 0
        count_force = 0
        count_moment = 0
        the_angles = []
        the_speed = []
        the_length = []
        count_wave = 0
        element_force_x = 0
        sum = 0
        for line in the_lines:
            count_line += 1
            try:
                line_formatted = line[0].strip().split()
            except IndexError:
                line_formatted = ["error", "error", "error", "error", "error"]
                pass
            if line_formatted[0] == "wave" and line_formatted[1] == "length":
                try:
                    example_wave_length = float(line_formatted[2])
                except:
                    example_wave_length = 100000
            if line_formatted[0] == "speed" and not bool_complete_vag:
                example_wave_speed = float(line_formatted[1])  # in meter per second
            if line_formatted[0] == "wave" and line_formatted[1] == "angle" and not bool_complete_vag:
                example_wave_angle = float(line_formatted[2])  # in degrees
            if line_formatted[0] == "wave" and line_formatted[2] == "frequency" and not bool_complete_vag:
                freq = float(line_formatted[3]) / 2 / np.pi  # in Hertz
            if example_wave_speed == speed and example_wave_length != 0 and not bool_complete_vag:
                bool_complete_vag = True
                coeff_vag = hull.wave_spectrum(example_wave_angle, freq)
            if bool_complete_vag:
                act_wave = Wv.Wave(example_wave_length, example_wave_angle, example_wave_speed, deep)
                wave_length = act_wave.wave_length
                angle = act_wave.angle
                if example_wave_speed == speed:
                    if example_wave_speed not in the_speed:
                        the_speed.append(example_wave_speed)
                    if example_wave_length not in the_length:
                        the_length.append(wave_length)
                    if example_wave_angle not in the_angles:
                        the_angles.append(angle)
                    if line_formatted[0] == "Force" and float(line_formatted[1]) == intersections[-1]:
                        end_wave = count_line + 1
                    if line_formatted[0] == "Force":
                        the_forces[count_force] += coeff_vag * float(
                            line_formatted[8])  # in ton
                        element_force_x = float(line_formatted[2])
                        count_force += 1
                    if line_formatted[0] == "Moment":
                        the_moments[count_moment] += coeff_vag * (
                                float(line_formatted[4]) - element_force_x * distance_from_neutral_axis)  # in ton
                        # We add a value with the forces in the x direction to compute the moment for the neutral axis
                        count_moment += 1
                    if count_line == end_wave:
                        sum += coeff_vag
                        bool_complete_vag = False  # we will compute a new wave so now the boolean is False
                        # We reinitialize every parameter
                        example_wave_length = 0
                        example_wave_angle = 0
                        example_wave_speed = None
                        element_force_x = 0
                        count_moment = 0
                        count_force = 0
                        count_wave += 1
                        # we re-initialize the wave values
        print(the_angles, the_speed, len(the_length), count_wave)  # Used to check the function
        print(sum, "sum", count_wave)  # same for the both speed the sum has to be the same, and the count wave too
        return the_forces, the_moments, intersections

    def computation_mx_BM_in_fct_mn_angle(self, speed: float, significant_wave_height: float, gamma: float,
                                          coeff_wave: float,
                                          water_depth: float, distance_from_neutral_axis: float):
        """That function can make the strength analysis with the information of the hull. It returns at the end the
                forces and the moments along the x-axis, we compute the maximum bending moments from -90° to 260°, with
                a step of 10°.

                :argument
                ---------
                speed: a float
                    It is the speed of the ship which is interesting for the computation of the maximum (in m/s)
                significant_wave_height: a float
                    It is the significant wave height of the hull, the Half of the maximum wave height (in m)
                gamma: a float
                    It is the gamma value for the JONSWAP wave distribution
                coeff_wave: a float
                    It is the coefficient for the spreading function cos 2s used in the hull object
                water_depth: a float
                    The depth for the hull condition (in m)
                distance_from_neutral_axis: a float
                    distance between the baseline and the neutral axis for the computation of the bending moment (in m)
                :returns
                --------
                Nothing it prints the maximum bending moment for the main mean direction, and the x coordinates where
                the bending moment is the greatest. It prints too, the graphs for the maximum bending moment for a
                certain mean direction and a certain speed.
                """
        # We just compute from 0 to 180 degrees because the values are symmetric from 180 to 360.
        angles = np.arange(35, 45, 1)  # degrees
        list_y1 = []
        list_y2 = []
        for angle in angles:
            print(angle, "act angle")
            hull = Hll.Hull(significant_wave_height, gamma, speed, coeff_wave, water_depth, angle)
            the_forces, the_moments, x_coordinates = self.analysis_res(hull, distance_from_neutral_axis)
            max_shear_forces = max(abs(max(the_forces)), abs(min(the_forces)))
            max_bending_moments = max(abs(max(the_moments)), abs(min(the_moments)))
            list_y1.append(max_shear_forces)
            list_y2.append(max_bending_moments)
        plt.plot(angles, list_y1)
        plt.title("Max shear forces in function of the mean direction")
        plt.show()
        plt.plot(angles, list_y2)
        plt.title("Max bending moments in function of the mean direction")
        plt.show()
        print(max(list_y2) * 9.80665, "max moment in kN.m")  # the result for an input file computed with that program
        # is in ton so we convert the result in kN
        i = list_y2.index(max(list_y2))
        print(angles[i], "angle of the maximum bending moment in degrees")
        hull = Hll.Hull(significant_wave_height, gamma, speed, coeff_wave, water_depth, angles[i])
        the_forces, the_moments, x_coordinates = self.analysis_res(hull, distance_from_neutral_axis)
        i = np.argmax(the_moments)
        print(x_coordinates[i], "coordinates of the maximum bending moment max")
        plt.plot(x_coordinates, the_forces)
        plt.title("Shear forces for the hull case the most limiting")
        plt.show()
        plt.plot(x_coordinates, the_moments)
        plt.title("Bending moments for the hull case the most limiting")
        plt.show()
        return

    def plot(self, hull: Hll.Hull, distance_from_neutral_axis):
        """That function permits to plot the graphs of the moments and the strengths along the x axis for the wave list
        and the hull associated. It plots the values with the units of the input file, in ton for the forces and in
        ton.m for the moments

        :argument
        ----------
        hull: a Hull.Hull object
            It is the hull object that represents the repartition of the wave list

        :returns
        --------
        forces: a list
            That is a list with the strength in the ship at each x coordinate (in ton)
        moments: a list
            It is a list with the moments in the ship at each x coordinates of intersections (in ton.m)
        list_x: a list
            It is a list with all the x_coordinates in the PD Strip coordinate system. That is the coordinates
            associated to the strength and the moments (in meter)
        """
        forces, moments, list_x = self.analysis_res(hull, distance_from_neutral_axis)
        plt.plot(list_x, forces)
        plt.title("Shear Forces in ton")
        plt.show()
        plt.plot(list_x, moments)
        plt.title("Bending Moments in ton.m")
        plt.show()
        return forces, moments, list_x

    def writing(self, hull: Hll.Hull, distance_from_neutral_axis):
        """That function permits to create a file named data_for_hull and all the values for every x coordinates of the
        forces and the moments will be saved. That writes directly the result computed with units of the file, in ton
        for the shear forces and in ton.m for the bending moments.

        :argument
        hull: Hull.Hull() object
            It's the specific hull, for the saved file

        :returns
        Nothing it creates a file with in a first column the x coordinates, the second one the forces, and the last one
        the moments for each x coordinates."""
        forces, moments, list_x = self.analysis_res(hull, distance_from_neutral_axis)
        n = len(forces)
        file_written = open("data_for_hull", "w")
        for i in range(n):
            file_written.write(str(list_x[i]) + " " + str(forces[i]) + " " + str(moments[i]) + "\n")

    def m_n(self, n: int, hull: Hll.Hull, distance_from_neutral_axis: float):
        g = 9.80665
        file = open(self.filename, "r", encoding="utf-8")
        the_lines = csv.reader(file)
        theta = []
        x_coordinates = []
        intersections = []
        prev_angle = False
        # calculation of the values of theta
        for line in the_lines:
            try:
                line_formatted = line[0].strip().split()
            except IndexError:
                line_formatted = ["error", "error", "error", "error", "error"]
                pass
            if line_formatted[0] == "Number" and line_formatted[1] == "of" and line_formatted[2] == "sections":
                nb_intersections = int(line_formatted[3]) - 1
            if line_formatted[0] == "Wave" and line_formatted[1] == "encounter":
                line_formatted.remove("Wave")
                line_formatted.remove("encounter")
                line_formatted.remove("angles")
                for angle in line_formatted:
                    val = float(angle)
                    theta.append(val)
                prev_angle = True
            if len(line_formatted) == 6:
                if line_formatted[4] == "x=":
                    x_coordinates.append(float(line_formatted[5]))
            if prev_angle and float(line_formatted[0]) != theta[0]:
                try:
                    float(line_formatted[0])
                    test = True
                except ValueError:
                    test = False
                if test:
                    for angle in line_formatted:
                        val = float(angle)
                        theta.append(val)
                prev_angle = False
        for i in range(len(x_coordinates) - 1):
            intersections.append(round((x_coordinates[i] + x_coordinates[i + 1]) / 2, 2))
        theta_f = theta.copy()
        for angle in theta:
            if (180 - angle) not in theta_f:
                theta_f.append(180 - angle)
        new = []
        for angle in theta_f:
            if angle >= 0:
                val = np.pi / 180 * angle
                new.append(val)
            else:
                val = (np.pi / 180) * (360 + angle)
                new.append(val)
        new.sort()
        theta_f.sort()
        if theta_f[0] == theta_f[-1] - 360:
            theta_f.remove(270)
        theta = new  # end of the calcul of the angles in radians
        m_int = {}
        print(intersections)
        for inter in intersections:
            m_n_mu = {}
            for angle in theta_f:
                file = open(self.filename, "r", encoding="utf-8")
                the_lines = csv.reader(file)
                forces = []
                count_mom = 0
                bool_complete_wv = False
                ex_wv_angle = None
                ex_wv_spd = None
                m_freq = {}
                for line in the_lines:
                    try:
                        line_formatted = line[0].strip().split()
                    except IndexError:
                        line_formatted = ["error", "error", "error", "error", "error"]
                        pass
                    if line_formatted[0] == "wave" and line_formatted[1] == "circ." and line_formatted[
                        2] == "frequency":
                        ex_wv_freq = float(line_formatted[3])
                        bool_complete_wv = False
                    if line_formatted[0] == "encounter" and line_formatted[1] == "frequency" and not bool_complete_wv:
                        ex_wv_freq_e = float(line_formatted[2])
                    if line_formatted[0] == "wave" and line_formatted[1] == "angle" and not bool_complete_wv:
                        ex_wv_angle = float(line_formatted[2])
                    if line_formatted[0] == "speed" and not bool_complete_wv:
                        ex_wv_spd = float(line_formatted[1])
                    if ex_wv_spd == hull.speed and ex_wv_angle == angle and not bool_complete_wv:
                        if line_formatted[0] == "Force":
                            ex_int = float(line_formatted[1])
                        if line_formatted[0] == "Moment" and ex_int == inter:
                            m_freq[ex_wv_freq] = float(line_formatted[4])
                            if count_mom == nb_intersections:
                                bool_complete_wv = True
                                ex_wv_angle = None
                                ex_wv_spd = None
                                ex_wv_freq = None
                freqs = list(m_freq.keys())
                vals = []
                for freq in freqs:
                    vals.append((abs(freq - ((freq) ** 2) * hull.speed * np.cos(angle) / g) ** n) * (
                            m_freq[freq] ** 2) * hull.JONSWAP(freq / 2 / np.pi))
                res = simpson(vals, freqs)
                if angle < 0:
                    mu = np.pi / 180 * (360 + angle)
                else:
                    mu = np.pi / 180 * angle
                m_n_mu[mu] = res
            the_angles = list(m_n_mu.keys())
            the_angles.sort()
            vals2 = []
            for mu2 in the_angles:
                vals2.append(m_n_mu[mu2] * hull.spreading_func(mu2 * 180 / np.pi))
            res2 = simpson(vals2, the_angles)
            m_int[inter] = res2
            print(inter)
        print(m_int)
        file2 = open("m_n_res", 'w', encoding="utf-8")
        for x, mn in m_int.items():
            file2.write(str(x) + " " + str(mn) + "\n")
        return max(m_int, key=m_int.get)

    def m_n_improved(self, n: int, hull: Hll.Hull, distance_from_neutral_axis: float):
        g = 9.80665
        file = open(self.filename, "r", encoding="utf-8")
        the_lines = csv.reader(file)
        theta = []
        x_coordinates = []
        intersections = []
        prev_angle = False
        ex_wv_angle = None
        ex_wv_spd = None
        ex_wv_freq = None
        # calculation of the values of theta
        for line in the_lines:
            try:
                line_formatted = line[0].strip().split()
            except IndexError:
                line_formatted = ["error", "error", "error", "error", "error"]
                pass
            if line_formatted[0] == "Number" and line_formatted[1] == "of" and line_formatted[2] == "sections":
                nb_intersections = int(line_formatted[3]) - 1
            if line_formatted[0] == "Wave" and line_formatted[1] == "encounter":
                line_formatted.remove("Wave")
                line_formatted.remove("encounter")
                line_formatted.remove("angles")
                for angle in line_formatted:
                    val = float(angle)
                    theta.append(val)
                prev_angle = True
            if len(line_formatted) == 6:
                if line_formatted[4] == "x=":
                    x_coordinates.append(float(line_formatted[5]))
            if prev_angle and float(line_formatted[0]) != theta[0]:
                try:
                    float(line_formatted[0])
                    test = True
                except ValueError:
                    test = False
                if test:
                    for angle in line_formatted:
                        val = float(angle)
                        theta.append(val)
                prev_angle = False
        for i in range(len(x_coordinates) - 1):
            intersections.append(round((x_coordinates[i] + x_coordinates[i + 1]) / 2, 2))
        theta_f = theta.copy()
        for angle in theta:
            if (180 - angle) not in theta_f:
                theta_f.append(180 - angle)
        new = []
        for angle in theta_f:
            if angle >= 0:
                val = np.pi / 180 * angle
                new.append(val)
            else:
                val = (np.pi / 180) * (360 + angle)
                new.append(val)
        new.sort()
        theta_f.sort()
        if theta_f[0] == theta_f[-1] - 360:
            theta_f.remove(270)
        theta = new  # end of the calcul of the angles in radians
        m_int = {}
        count_mom = 0
        for inter in intersections:
            m_int[inter] = {}
            for theta3 in theta_f:
                m_int[inter][theta3] = {}
        file = open(self.filename, "r", encoding="utf-8")
        the_lines = csv.reader(file)
        for line in the_lines:
            try:
                line_formatted = line[0].strip().split()
            except IndexError:
                line_formatted = ["error", "error", "error", "error", "error"]
                pass
            if line_formatted[0] == "wave" and line_formatted[1] == "circ." and line_formatted[
                2] == "frequency":
                ex_wv_freq = float(line_formatted[3])
                bool_complete_wv = False
            if line_formatted[0] == "encounter" and line_formatted[1] == "frequency" and not bool_complete_wv:
                ex_wv_freq_e = float(line_formatted[2])
            if line_formatted[0] == "wave" and line_formatted[1] == "angle" and not bool_complete_wv:
                ex_wv_angle = float(line_formatted[2])
            if line_formatted[0] == "speed" and not bool_complete_wv:
                ex_wv_spd = float(line_formatted[1])
            if ex_wv_spd == hull.speed and not bool_complete_wv:
                if line_formatted[0] == "Force":
                    ex_int = float(line_formatted[1])
                if line_formatted[0] == "Moment":
                    m_int[ex_int][ex_wv_angle][(ex_wv_freq, ex_wv_freq_e)] = float(line_formatted[4])*9.80665
                    count_mom += 1
                    if count_mom == nb_intersections:
                        bool_complete_wv = True
                        ex_wv_angle = None
                        ex_wv_spd = None
                        ex_wv_freq = None
                        count_mom = 0
        res_func_int = {}
        for x, m_angle in m_int.items():
            res2 = {}
            res_func_angle = {}
            for theta2, m_freq in m_angle.items():
                freqs = list(m_freq.keys())
                freqs = sorted(freqs, key=lambda x: x[0])
                vals = {}
                list_freq1 = []
                for freq in freqs:
                    freq_th = freq[0]
                    freq_e = freq[1]
                    list_freq1.append(freq_th)
                    vals[freq_th]=((abs(freq_e) ** n) * (m_int[x][theta2][freq] ** 2) * hull.JONSWAP(abs(freq_e) / 2 / np.pi))
                list_freq1.sort()
                list_int1=[]
                for f1 in list_freq1:
                    list_int1.append(vals[f1])
                res = simpson(list_int1, list_freq1)
                if pd.isna(res):
                    res=0
                res_func_angle[theta2] = res
            angles = list(m_angle.keys())
            for angle in angles:
                res2[angle] = (res_func_angle[angle] * hull.spread_func_int(angle,10))
            les_res = {}
            for theta4, res4 in res2.items():
                if theta4 < 0:
                    angle_fin = (theta4 + 360) * np.pi / 180
                else:
                    angle_fin = theta4 * np.pi / 180
                les_res[angle_fin] = res4
            mu = list(les_res.keys())
            mu.sort()
            les_y = []
            for val1 in mu:
                les_y.append(les_res[val1])
            res3 = simpson(les_y, mu)
            res_func_int[x] = res3
        les_x = list(res_func_int.keys())
        les_x.sort()
        file2 = open("mn_opt_res", "w", encoding="utf-8")
        for x in les_x:
            file2.write(str(x) + " " + str(res_func_int[x]) + "\n")
        max_key = max(res_func_int, key=res_func_int.get)
        print(max_key, res_func_int[max_key])
        return res_func_int

    def max_BM_func_dir(self, significant_wav_height, gamma, speed, coeff_wave, deep, distance_from_neutral_axis, D,
                        alpha):
        the_angles = np.arange(0, 190, 45)
        file = open("res_BM_func_dir"+str(speed), "w", encoding="utf-8")
        for angle in the_angles:
            hull = Hll.Hull(significant_wav_height, gamma, speed, coeff_wave, deep, angle)
            m0 = self.m_n_improved(0, hull, distance_from_neutral_axis)
            m2 = self.m_n_improved(2, hull, distance_from_neutral_axis)
            inter = list(m0.keys())
            max_BM = {}
            file.write(str(angle) + "\n")
            for x in inter:
                max_BM[x] = (np.sqrt(2 * np.log(D / 2 / np.pi / alpha * np.sqrt(m2[x] / m0[x]))) * np.sqrt(
                    m0[x]))
                file.write(
                    str(x) + " " + str(np.sqrt(2 * np.log(D / 2 / np.pi / alpha * np.sqrt(m2[x] / m0[x]))) * np.sqrt(
                        m0[x])) + "\n")
            max_key = max(max_BM, key=max_BM.get)
            print(max_key, max_BM[max_key])
            print(angle)
            plt.plot(inter, max_BM.values(), label=str(angle))
            plt.legend()
        plt.show()
        file.close()

    def max_BM_func_spd(self, significant_wav_height, gamma, list_speed, coeff_wave, deep, distance_from_neutral_axis,
                        angle, D, alpha):
        file = open("res_BM_func_spd", "w", encoding="utf-8")
        for spd in list_speed:
            hull = Hll.Hull(significant_wav_height, gamma, spd, coeff_wave, deep, angle)
            m0 = self.m_n_improved(0, hull, distance_from_neutral_axis)
            m2 = self.m_n_improved(2, hull, distance_from_neutral_axis)
            inter = list(m0.keys())
            max_BM = {}
            file.write(str(spd) + "\n")
            for x in inter:
                max_BM[x] = (np.sqrt(2 * np.log(D / 2 / np.pi / alpha * np.sqrt(m2[x] / m0[x]))) * np.sqrt(
                    m0[x]))
                file.write(
                    str(x) + " " + str(np.sqrt(2 * np.log(D / 2 / np.pi / alpha * np.sqrt(m2[x] / m0[x]))) * np.sqrt(
                        m0[x])) + "\n")
            max_key = max(max_BM, key=max_BM.get)
            print(max_key, max_BM[max_key])
            print(spd)
            plt.plot(inter, max_BM.values(), label=str(spd))
            plt.legend()
        plt.show()
        file.close()

    def m_n_func_wv_height(self, n, list_significant_wav_height, gamma, speed, coeff_wave, deep,
                           distance_from_neutral_axis, angle):
        for wv_height in list_significant_wav_height:
            hull = Hll.Hull(wv_height, gamma, speed, coeff_wave, deep, angle)
            res = self.m_n_improved(n, hull, distance_from_neutral_axis)
            plt.plot(list(res.keys()), list(res.values()), label=str(wv_height))
            plt.legend()
        plt.show()

    def print_filewritting_max_BM(self, hull, dist, D, alpha):
        m0 = self.m_n_improved(0, hull, dist)
        m2 = self.m_n_improved(2, hull, dist)
        inter = list(m0.keys())
        max_BM = {}
        file = open("max_BM", "w", encoding="utf-8")
        for x in inter:
            max_BM[x] = (np.sqrt(2 * np.log(D / 2 / np.pi / alpha * np.sqrt(m2[x] / m0[x]))) * np.sqrt(m0[x]))
            file.write(str(x) + " " + str(
                np.sqrt(2 * np.log(D / 2 / np.pi / alpha * np.sqrt(m2[x] / m0[x]))) * np.sqrt(m0[x])) + "\n")
        file.close()
        max_key = max(max_BM, key=max_BM.get)
        print(max_key, max_BM[max_key])
        plt.plot(inter, max_BM.values())
        plt.show()

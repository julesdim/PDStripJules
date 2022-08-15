import unittest
import numpy as np
import Loading as ld
import Masses as ms
import Shape as sh
import Frames as fr
import csv
import Hull as Hll
import Wave as Wv
import Analysis as Anls


def arrondir(val, nb):
    ret = []
    for i in range(len(val)):
        x = round(val[i], nb)
        ret.append(x)
    return tuple(ret)

def combinaisons():
    les_angles = np.arange(-90, 270, 10)
    les_freq = np.arange(0.1, 2.58, 0.08)
    poss = []
    for angle in les_angles:
        for freq in les_freq:
            poss.append((angle, round(freq, 3)))
    return poss

class Test(unittest.TestCase):
    def test_calcul_mass(self):
        mass = ms.Mass(50, 0, 50, 25, 0, 1, 1, 1)
        mass_test = 25
        mass_th = mass.calcul_mass(0, 25)
        self.assertEqual(mass_th, mass_test)
        mass_th = mass.calcul_mass(0, 50)
        mass_test = 50
        self.assertEqual(mass_th, mass_test)
        mass_th = mass.calcul_mass(-50, 25)
        mass_test = 25
        self.assertEqual(mass_th, mass_test)
        mass_th = mass.calcul_mass(25, 75)
        mass_test = 25
        self.assertEqual(mass_th, mass_test)
        mass_th = mass.calcul_mass(50, 100)
        mass_test = 0
        self.assertEqual(mass_th, mass_test)

    def test_calcul_xg_mass(self):
        mass = ms.Mass(50, 0, 50, 25, 0, 1, 1, 1)
        xg_th = mass.calcul_xg_mass()
        self.assertEqual(25, xg_th)
        mass = ms.Mass(50, 0, 25, 10, 0, 1, 2, 2)
        mass.calcul_xg_not_the_mid()
        xg_th = mass.calcul_xg_mass()
        self.assertTrue(mass.x_coordinate_CoG - 0.000001 < xg_th < mass.x_coordinate_CoG + 0.000001)

    def test_calcul_linear_density_coordinates(self):
        mass = ms.Mass(50, 0, 50, 25, 0, 1, 1, 1)
        begin, end = mass.calcul_linear_density_for_coordinates(0, 50)
        self.assertEqual(1, begin)
        self.assertEqual(1, end)

    def test_calcul_linear_density(self):
        mass = ms.Mass(50, 0, 50, 45, 0, 1, 1, 1)
        begin, end = mass.calcul_linear_density_for_coordinates(0, 50)
        mass_th = (end + begin) * (mass.x_end - mass.x_start) / 2
        self.assertEqual(mass.weight, mass_th)

    def test_mass_calc_for_coord(self):
        mass = ms.Mass(50, 0, 50, 25, 0, 1, 1, 1)
        load = ld.Loading()
        load.__append__(mass)
        mass_th = load.mass_calculation_for_coordinates(0, 50)
        mass_test = 50
        self.assertEqual(mass_th, mass_test)
        mass = ms.Mass(50, 0, 50, 30, 0, 1, 1, 1)
        mass.calcul_xg_not_the_mid()
        mass_th = load.mass_calculation_for_coordinates(0, 50)
        self.assertEqual(mass_th, mass_test)

    def test_CoG_coord(self):
        mass = ms.Mass(50, 0, 50, 25, 0, 1, 1, 1)
        load = ld.Loading()
        load.__append__(mass)
        xg_test = (10, 0, 1)
        xg_th = load.calcul_center_of_gravity_for_coordinates(0, 20)
        print(xg_th, xg_test)
        self.assertEqual(xg_th, xg_test)

    def test_pdstrip_conversion(self):
        mass = ms.Mass(50, 0, 50, 25, 0, 1, 1, 1)
        load = ld.Loading()
        load.__append__(mass)
        load.pdstrip_coordinates_conversion(25)
        self.assertEqual(load.masses[0].x_coordinate_CoG, 0)

    def test_inertia(self):
        mass = ms.Mass(50, 0, 50, 25, 0, 1, 1, 1)
        load = ld.Loading()
        load.__append__(mass)
        test = load.calcul_every_parameters(25, 0, 1, 0, 50)
        les_x = np.linspace(0, 50, 10)
        for i in range(len(les_x) - 1):
            start = les_x[i]
            end = les_x[i + 1]
            X_CoG, Y_CoG, Z_CoG = load.calcul_center_of_gravity_for_coordinates(start, end)
            test = load.calcul_every_parameters(X_CoG, 0, 1, start, end)
            self.assertEqual((0, 0, 0, 0, 0, 0), arrondir(test, 20))
        self.assertEqual((0, 0, 0, 0, 0, 0), arrondir(test, 20))

    def test_frames(self):
        frames = fr.Frames(20)
        frames.__append__((0, 0))
        self.assertEqual(frames.x_coordinate, 20)
        self.assertEqual(frames.coordinates[0], (0, 0))

    def test_shape(self):
        frame1 = fr.Frames(20)
        frame1.__append__((0, 0))
        frame2 = fr.Frames(25)
        frame2.__append__((1, 0))
        form = sh.Form()
        form.__append__(frame1)
        form.__append__(frame2)
        CoG = form.center_of_gravity_no_mass_for_coordinates(10, 30)
        self.assertEqual(CoG, (20, 0.5, 0))

    def test_every_waves(self):
        filename = "pdstrip.out.ok"
        speed = 0
        significant_wave_height = 1.2
        gamma = 1.8
        coeff = 5
        deep=17.56
        list_poss_vague = combinaisons()
        x_coordinates = []
        file = open(filename, "r", encoding="utf-8")
        the_lines = csv.reader(file)
        hull = Hll.Hull(significant_wave_height, gamma, speed, coeff, 17.56, 0)
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
        the_forces = np.zeros(len(x_coordinates))
        the_moments = np.zeros(len(x_coordinates))
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
            if line_formatted[0] == "speed":
                example_wave_speed = float(line_formatted[1])
            if line_formatted[0] == "wave" and line_formatted[1] == "angle":
                example_wave_angle = float(line_formatted[2])
            if line_formatted[0] == "wave" and line_formatted[2] == "frequency" and not bool_complete_vag:
                freq = float(line_formatted[3])  # in Hertz
            if example_wave_speed == speed and example_wave_length != 0:
                bool_complete_vag = True
            if bool_complete_vag:
                act_wave = Wv.Wave(example_wave_length, example_wave_angle, example_wave_speed, deep)
                wave_length = act_wave.wave_length
                angle = act_wave.angle
                freq=act_wave.freq*2*np.pi
                coeff_repartition_angul = hull.spreading_coeff(act_wave.angle,act_wave.freq)
                if example_wave_speed == speed:
                    if example_wave_speed not in the_speed:
                        the_speed.append(example_wave_speed)
                    if example_wave_length not in the_length:
                        sum += hull.JONSWAP(freq) * coeff_repartition_angul
                        the_length.append(example_wave_length)
                    if example_wave_angle not in the_angles:
                        the_angles.append(example_wave_angle)
                        sum += hull.JONSWAP(freq) * coeff_repartition_angul
                    if line_formatted[0] == "Force" and float(line_formatted[1]) == intersections[-1]:
                        end_wave = count_line + 1
                    if line_formatted[0] == "Force":
                        the_forces[count_force] += coeff_repartition_angul * hull.JONSWAP(freq) * float(
                            line_formatted[8])
                        element_force_x = float(line_formatted[2])
                        count_force += 1
                    if line_formatted[0] == "Moment":
                        the_moments[count_moment] += coeff_repartition_angul * hull.JONSWAP(freq) * (
                                float(line_formatted[4]) - element_force_x * 1.821)
                        count_moment += 1
                    if count_line == end_wave:
                        comb = (example_wave_angle, round(freq,3))
                        list_poss_vague.remove(comb)
                        bool_complete_vag = False
                        example_wave_length = 0
                        example_wave_angle = 0
                        example_wave_speed = None
                        element_force_x = 0
                        count_moment = 0
                        count_force = 0
                        count_wave += 1
        file.close()
        self.assertEqual(list_poss_vague, [])

    def test_analysis_res(self):
        x_coordinates = []
        filename = "pdstrip.out.ok"
        file = open(filename, "r", encoding="utf-8")
        the_lines = csv.reader(file)
        hull = Hll.Hull(0.6, 1.8, 0, 5, 7.56,0)
        speed = hull.speed  # m/s
        significant_wave_height = hull.significant_wave_height  # meter
        gamma = hull.gamma
        coeff = hull.coeff_wave
        deep = hull.deep  # meter
        distance_from_neutral_axis=1.821
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
                coeff_vag = hull.spreading_coeff(example_wave_angle, freq)
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
                        self.assertEqual(count_moment,count_force)
                        self.assertEqual(count_moment,len(the_moments))
                        # We reinitialize every parameter
                        example_wave_length = 0
                        example_wave_angle = 0
                        example_wave_speed = None
                        element_force_x = 0
                        count_moment = 0
                        count_force = 0
                        count_wave += 1
        self.assertEqual(count_wave,len(the_length)*len(the_angles))

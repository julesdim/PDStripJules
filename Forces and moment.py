import matplotlib.pyplot as plt
import csv
import numpy as np


def graph_file_for_one_wave(filename: str, wave_length: float, wave_angle: float,
                            wave_speed: float, text: str, boolean_print, bool_initialisation, text_fonction,
                            Lpp: float):
    """That function plots 2 different graphs. The first is the shear forces along the x-axis and the second is the
    bending moment. It can be the real, the immaginary or the absolute value of that force.

    :argument
    ---------
    filename: a text object
        that's the name of the pdstrip results file, it ends by ".out.ok"
    wave_length: a float
        that's the value of the length of the wave that we want to print the results (in m)
    wave_angle: a float
        that's the value of the angle of the wave that we want to print the results (in degree)
    wave_speed: a float
        that's the value of the speed of the wave that we want to print the results (in m/s)
    text: a text object
        3 possibilities, "real" to plot the real part
        "imaginary" to plot the imaginary part
        "absolute" to plot the absolute value"
    boolean_print: a boolean
        if that boolean is True, the figure will be printed
        if not the figure will be saved and will be printed when the boolean will be True
    bool_initialisation: a boolean
        it permits to know if that is the first time we use the function to initialise the plot and the figures in the
        subplots
    text_function: a text object
        2 possibilities "speed" or "wave length" to plot the graph in function of that characteristics
    Lpp: a float
        length between perpendiculars (in m)
    :returns
    --------
    It plots 2 graphs for the wave selected and for the part selected (real, imaginary or absolute), it prints the
    shear forces and the bending moments.
    """
    if text == "real":
        constant = 0
    elif text == "imaginary":
        constant = 1
    elif text == "absolute":
        constant = 2
    file = open(filename, "r", encoding="utf-8")
    file_written = open("data_results.csv", "w")
    the_lines = csv.reader(file)
    line_counter = 0
    list_x_coordinates = []
    forces_along_z = []
    moment_along_y = []
    number_of_section = 0
    example_wave_length = 0
    example_wave_angle = 0
    example_wave_speed = 0
    for line in the_lines:
        try:
            line_formatted = line[0].strip().split()
        except IndexError:
            pass
        if len(line_formatted) == 0:
            line_formatted = ["error", "error", "error"]
        if line_formatted[0] == "wave" and line_formatted[2] == "frequency":
            example_wave_frequency = float(line_formatted[3])
        if line_formatted[0] == "wave" and line_formatted[1] == "length":
            try:
                example_wave_length = float(line_formatted[2])
                # if example_wave_length == wave_length:
                # print(example_wave_frequency)
            except:
                example_wave_length = 100000
        if line_formatted[0] == "wave" and line_formatted[1] == "angle":
            example_wave_angle = float(line_formatted[2])
        if line_formatted[0] == "speed":
            example_wave_speed = float(line_formatted[1])
        if line_formatted[0] == "Number":
            number_of_section = float(line_formatted[3])
        if example_wave_length == wave_length and example_wave_angle == wave_angle and example_wave_speed == wave_speed:
            if line_formatted[0] == "Force":
                x_coordinate = float(line_formatted[1]) + Lpp / 2
                element_force_z = float(line_formatted[8 + constant])
                forces_along_z.append(element_force_z)
                if x_coordinate not in list_x_coordinates:
                    list_x_coordinates.append(x_coordinate)
            if line_counter == (number_of_section - 1):
                break
            if line_formatted[0] == "Moment":
                line_counter += 1
                element_moment_y = float(line_formatted[4 + constant])
                moment_along_y.append(element_moment_y)
    forces_along_z = np.array(forces_along_z)
    all_graph = [forces_along_z, moment_along_y]
    n_all_graph = len(all_graph)
    for i in range(len(list_x_coordinates)):
        file_written.write(str(list_x_coordinates[i]) + " ")
        for graph in all_graph:
            file_written.write(str(graph[i]) + " ")
        file_written.write("\n")
    if bool_initialisation:
        fig = plt.figure()
        subplot1 = fig.add_subplot(1, 2, 1)
        subplot2 = fig.add_subplot(1, 2, 2)
        global list_subplot
        list_subplot = [subplot1, subplot2]
    if text_fonction == "speed":
        legend = str(wave_speed)
    if text_fonction == "wave length":
        legend = str(wave_length)
    if text_fonction == "wave angle":
        legend = str(wave_angle)
    for i in range(n_all_graph):
        list_subplot[i].plot(list_x_coordinates, all_graph[i], label=legend)
        plt.grid()
        if i == 0:
            list_subplot[i].set_title("Forces along z axis")
            list_subplot[i].legend()
        if i == 1:
            list_subplot[i].set_title("Bending moment")
            list_subplot[i].legend()
    if boolean_print:
        plt.show()
    return


def graph_many_wave_length(filename: str, first_wave_length: float, last_wave_length: float, speed: float, angle: float,
                           text: str, Lpp: float):
    """That function plots 2 graphs with the first the shear forces for different wave length and the second one the
    bending moment for the same wave lengths.

    :argument
    -----------
    filename: a str
        that is the name of the result file where the results for all the wave length is saved
    first_wave_length: a float
        that is the first wave length that we want plot. We will print all the graphs for the wave lengths greater than
        this value (in m)
    last_wave_length: a float
        that is the last wave length that we want plot. We will print all the graphs for the wave length lower than this
        value (in m)
    speed: a float
        that is the speed of the ship for the results plotted (in m/s)
    angle: a float
        that is the angle of the incident waves, for every graph (in degree)
    text: 3 possibilities, "real", "imaginary" or "absolute", to plot the value wanted by the user of the shear forces
    or the bending moment.
    Lpp: a float
        length between perpendiculars (in m)"""
    list_wave_length = []
    bool_initialisation = True
    file = open(filename, "r", encoding="utf-8")
    the_lines = csv.reader(file)
    boolean_print = False
    example_wave_length = 0
    text_fonction = "wave length"
    for line in the_lines:
        try:
            line_formatted = line[0].strip().split()
        except IndexError:
            pass
        if len(line_formatted) == 0:
            line_formatted = ["error", "error", "error"]
        if line_formatted[0] == "wave" and line_formatted[1] == "length":
            example_wave_length = float(line_formatted[2])
        if first_wave_length < example_wave_length < last_wave_length and example_wave_length not in list_wave_length:
            list_wave_length.append(example_wave_length)
    if len(list_wave_length) == 0:
        return
    for i in range(len(list_wave_length)):
        wave_length = list_wave_length[i]
        if i == 1:
            bool_initialisation = False
        if i == len(list_wave_length) - 1:
            boolean_print = True
        graph_file_for_one_wave(filename, wave_length, angle, speed, text, boolean_print, bool_initialisation,
                                text_fonction, Lpp)
    return

def graph_many_wave_angle(filename: str, speed: float, wave_length,
                           text: str, Lpp: float):
    """That function plots 2 graphs with the first the shear forces for different wave length and the second one the
    bending moment for the same wave lengths.

    :argument
    -----------
    filename: a str
        that is the name of the result file where the results for all the wave length is saved
    first_wave_length: a float
        that is the first wave length that we want plot. We will print all the graphs for the wave lengths greater than
        this value (in m)
    last_wave_length: a float
        that is the last wave length that we want plot. We will print all the graphs for the wave length lower than this
        value (in m)
    speed: a float
        that is the speed of the ship for the results plotted (in m/s)
    angle: a float
        that is the angle of the incident waves, for every graph (in degree)
    text: 3 possibilities, "real", "imaginary" or "absolute", to plot the value wanted by the user of the shear forces
    or the bending moment.
    Lpp: a float
        length between perpendiculars (in m)"""
    list_wave_length = []
    bool_initialisation = True
    file = open(filename, "r", encoding="utf-8")
    the_lines = csv.reader(file)
    boolean_print = False
    example_wave_length = 0
    text_fonction = "wave angle"
    list_wave_angle=np.arange(0,190,30)
    if len(list_wave_angle) == 0:
        return
    for i in range(len(list_wave_angle)):
        wave_angle = list_wave_angle[i]
        if i == 1:
            bool_initialisation = False
        if i == len(list_wave_length) - 1:
            boolean_print = True
        graph_file_for_one_wave(filename, wave_length, wave_angle, speed, text, boolean_print, bool_initialisation,
                                text_fonction, Lpp)
    return

def graph_many_speed(filename: str, first_speed: float, last_speed: float, wave_length: float, angle: float, text: str,
                     Lpp: float):
    """That function plots 2 graphs with the first the shear forces for different speeds and the second one the
    bending moment for the same speeds.

    :argument
    -----------
    filename: a str
        that is the name of the result file where the results for all the wave length is saved
    first_speed: a float
        that is the first speed that we want plot. We will print all the graphs for the speeds greater than
        this value (in m/s)
    last_speed: a float
        that is the last speed that we want plot. We will print all the graphs for the speeds lower than this
        value (in m/s)
    wave_length: a float
        that is the wave length of the incident wave of the graphs plotted (in m)
    angle: a float
        that is the angle of the incident waves, for every graph (in degree)
    text: a str
        3 possibilities, "real", "imaginary" or "absolute", to plot the value wanted by the user of the shear forces
        or the bending moment.
    Lpp: a float
        length between perpendiculars (in m)"""
    list_speed = []
    bool_initialisation = True
    file = open(filename, "r", encoding="utf-8")
    the_lines = csv.reader(file)
    boolean_print = False
    example_wave_length = 0
    text_fonction = "speed"
    for line in the_lines:
        try:
            line_formatted = line[0].strip().split()
        except IndexError:
            pass
        if len(line_formatted) == 0:
            line_formatted = ["error", "error", "error"]
        if line_formatted[0] == "speed":
            example_wave_speed = float(line_formatted[1])
            if first_speed <= example_wave_speed <= last_speed and example_wave_speed not in list_speed:
                list_speed.append(example_wave_speed)
    if len(list_speed) == 0:
        return
    for i in range(len(list_speed)):
        speed = list_speed[i]
        if i == len(list_speed) - 1:
            boolean_print = True
        graph_file_for_one_wave(filename, wave_length, angle, speed, text, boolean_print, bool_initialisation,
                                text_fonction, Lpp)
        bool_initialisation = False
    return


def printing_the_extreme_point(filename: str, text: str, Lpp: float):
    """That functions permits to return every extreme point for the values of the shear forces and the bending, for real
    imaginary or absolute value.
    moments.

    :argument
    -----------
    filename: a str
        the name of the file where the results are saved
    text: a str
        3 possibilities "real", "imaginary" or "absolute", to print the part of the extreme point the user is
        interested in
    Lpp: a float
        The length between perpendiculars of the ship (in m)
    :returns
    ---------
    extreme_point_shear_forces: a float
        the value of the greatest shear forces in the file (in ton.m)
    wave_length_extremesf: a float
        the value of the wave length corresponding to the extreme point (in m)
    speed_extremesf: a float
        the value of the speed corresponding to the extreme point (in m/s)
    angle_extremesf: a float
        the value of the angle corresponding to the extreme point (in degree)
    x_coordinate_extremesf: a float
        the x_coordinate, corresponding to the extreme point of shear forces (in m)
    extreme_point_bending_moments: a float
        the value of the greatest bending moment in the file, that value is the one computed from the point of the
        neutral axis (in ton.m)
    wave_length_extremebm: a float
        the value of the wave length corresponding to the extreme point of bending moment (in m)
    speed_extremebm: a float
        the value of speed corresponding to the extreme point of bending moment (in m)
    angle_extremebm: a float
        the value of the angle corresponding to the extreme point of bending moment (in degree)
    x_coordinate_extremebm: a float
        the x_coordinate, corresponding to the extreme point of bending moment (in m)
    """
    if text == "real":
        constant = 0
    elif text == "imaginary":
        constant = 1
    elif text == "absolute":
        constant = 2
    file = open(filename, "r", encoding="utf-8")
    extreme_point_shear_forces = 0
    extreme_point_bending_moments = 0
    wave_length_extremesf = 0
    wave_length_extremebm = 0
    speed_extremebm = 0
    speed_extremesf = 0
    angle_extremesf = 0
    angle_extremebm = 0
    the_lines = csv.reader(file)
    line_counter = 0
    number_of_section = 0
    example_wave_frequency = 0
    example_wave_length = 0
    example_wave_angle = 0
    example_wave_speed = 0
    for line in the_lines:
        try:
            line_formatted = line[0].strip().split()
        except IndexError:
            pass
        if len(line_formatted) == 0:
            line_formatted = ["error", "error", "error"]
        if line_formatted[0] == "wave" and line_formatted[2] == "frequency":
            example_wave_frequency = float(line_formatted[3])
        if line_formatted[0] == "wave" and line_formatted[1] == "length":
            try:
                example_wave_length = float(line_formatted[2])
                # if example_wave_length == wave_length:
                # print(example_wave_frequency)
            except:
                example_wave_length = 100000
        if line_formatted[0] == "wave" and line_formatted[1] == "angle":
            example_wave_angle = float(line_formatted[2])
        if line_formatted[0] == "speed":
            example_wave_speed = float(line_formatted[1])
        if line_formatted[0] == "Number":
            number_of_section = float(line_formatted[3])
        if line_formatted[0] == "Force":
            x_coordinate = float(line_formatted[1]) + Lpp / 2
            element_force_x = float(line_formatted[2 + constant])
            element_force_z = float(line_formatted[8 + constant])
            if abs(element_force_z) > extreme_point_shear_forces:
                extreme_point_shear_forces = abs(element_force_z)
                wave_length_extremesf = example_wave_length
                speed_extremesf = example_wave_speed
                angle_extremesf = example_wave_angle
                x_coordinate_extremesf = x_coordinate
        if line_formatted[0] == "Moment":
            line_counter += 1
            element_moment_y = float(line_formatted[4 + constant]) - element_force_x * 1.821
            if abs(element_moment_y) > extreme_point_bending_moments:
                extreme_point_bending_moments = abs(element_moment_y)
                wave_length_extremebm = example_wave_length
                speed_extremebm = example_wave_speed
                angle_extremebm = example_wave_angle
                x_coordinate_extremebm = x_coordinate
    return (extreme_point_shear_forces, wave_length_extremesf, speed_extremesf, angle_extremesf, x_coordinate_extremesf,
            extreme_point_bending_moments, wave_length_extremebm, speed_extremebm, angle_extremebm,
            x_coordinate_extremebm)

graph_many_wave_angle("pdstrip.out.ok",0,540.3,"absolute",135)
# graph_file_for_one_wave("pdstrip.out.ok", 183.17, 0, 0, "real", True, True, "speed")
# graph_many_wave_length("pdstrip.out.ok", 100, 7000, 0, 0, "real", 135)
# graph_many_wave_length("pdstrip.out.ok", 100, 7000, 4.17, 0, "real", 135)
# graph_many_speed("pdstrip.out.ok", 0, 5, 540.3, 0, "real", 135)
# graph_many_speed("pdstrip.out.ok", 0, 5, 299.31, 0, "real", 135)
# graph_many_speed("pdstrip.out.ok", 0, 5, 206.27, 0, "real", 135)
# extreme_point = printing_the_extreme_point("pdstrip.out.ok", "real", 135)
# frequence de resonnance de 0.269 et frequence de 0.218 (pour 246.47 de frequence de resonance) et de 0.316 pour 171.16
# graph_file_for_one_wave("pdstrip.out.ok", 104.72, 180, 0, "real", True, True, "speed")
# graph_file_for_one_wave("pdstrip.out.ok", 385.11, 0, 0, "real", True, True, "speed")
# print(extreme_point)


import numpy as np
import matplotlib.pyplot as plt
import Frames as fr
import Loading as ld
import csv


class Form:
    """That a class to define the form of the ship, it's a list of frame

    :argument
    ---------
    shape: a list
        it's a list of frame for x-axis"""

    def __init__(self):
        self.shape = []

    def __append__(self, new_frame: fr.Frames):
        """That function appends a new element to the existing list

        :argument
        -----------
        new_frame: a frame object
            It's a new frame defined by an x coordinates and different points for that x. coordinates are in meter

        :returns
        -----------
        A new element is added to the existing shape"""
        self.shape.append(new_frame)

    def collection_of_coordinates(self, filename: str):
        """That functions read a Pias file of coordinates, and it reads every coordinate to return a shape object

        :parameter
        -----------
        filename: a text variable corresponding to the pias file of coordinates

        :returns
        ----------
        form: A Form object corresponding to a list of frame, defined by a x coordinate and coordinates of point for that
        frame
        """
        file = open(filename, "r", encoding="utf-8")
        the_lines = csv.reader(file)
        line_counter = 0
        start_line_of_the_frame = 0
        end_line_of_the_frame = 1000  # no importance just for an initialization
        boolean_for_beginning_of_the_frame = True
        # we can know if we are at the beginning of the frame to get the x coordinate of the frame
        for line in the_lines:
            line_formatted = line[0].strip().split()  # formatting the line
            if line_counter == 0:
                total_number_of_frame = float(line_formatted[0])  # to know the number of frame
            if line_counter != 0 and len(line_formatted) == 1 and boolean_for_beginning_of_the_frame:
                x_coordinate = (float(line_formatted[0]))
                # if there is just one coordinate, it is the position along x-axis, if beg_frame==True
                start_line_of_the_frame = line_counter
                current_frame = fr.Frames(x_coordinate)
                boolean_for_beginning_of_the_frame = False  # As we just passed the start of the frame, next is false
            if line_counter != 0 and len(line_formatted) == 3:
                y_coordinate = (float(line_formatted[0]))
                z_coordinate = (float(line_formatted[1]))
                current_frame.__append__((y_coordinate, z_coordinate))  # we append to the list
                if y_coordinate != 0:
                    current_frame.__append__((-y_coordinate, z_coordinate))
                # We suppose that the ship is symmetric
                boolean_for_beginning_of_the_frame = True
                # the next time len(new)==1 it will be the beginning of a new frame
            if len(line_formatted) == 1 and line_counter == start_line_of_the_frame + 1:
                number_points_of_the_current_frame = int(line_formatted[0])
                end_line_of_the_frame = line_counter + number_points_of_the_current_frame
            if line_counter == end_line_of_the_frame:
                self.__append__(current_frame)
            line_counter += 1
        file.close()

    def center_of_gravity_no_mass_for_coordinates(self, x_start: float, x_end: float):
        """That function allows the user to know the center of gravity if there's no mass between the x_start and the
        x_end of the section. It returns the center of volume instead of the center of gravity, because we suppose there
        is no mass between start and en.

        :argument
        ---------
        x_start: a float
            the x coordinate at the start of the section (in m)
        x_end: a float
            the x coordinate at the end of the section (in m)

        :returns
        -----------
        x_CoG: a float
            it corresponds to the center of the volume of the x-axis (in m)
        y_CoG: a float
            it corresponds to the center of the volume of the y-axis (in m)
        z_CoG: a float
            it corresponds to the center of the volume of the z-axis (in m)"""
        list_y_coordinates = []
        list_z_coordinates = []
        x_CoG = (x_start + x_end) / 2
        for frame in self.shape:
            if x_end >= frame.x_coordinate >= x_start:
                for coordinate in frame.coordinates:
                    list_y_coordinates.append(coordinate[0])
                    list_z_coordinates.append(coordinate[1])
        return x_CoG, np.mean(list_y_coordinates), np.mean(list_z_coordinates)

    def correction_of_coordinates_for_up(self):
        """That function corrects the coordinates if it misses a point at the maximum z to get the good Z_CoG, in the
        case of no mass.

        :argument
        ----------
        self: a form object

        :returns
        -----------
        self: a form object
            the point needed are added to the current form """
        n = len(self.shape)
        for j in range(n):
            list_z_coordinates_of_frame = []  # all the z of the frame
            list_y_coordinates_of_frame = []
            coordinates = self.shape[j].coordinates
            number_coordinates = len(coordinates)
            for i in range(number_coordinates):
                list_z_coordinates_of_frame.append(coordinates[i][1])
                list_y_coordinates_of_frame.append(coordinates[i][0])
            maximum_z = max(list_z_coordinates_of_frame)  # we save the max of z in pd strip coordinates, that means
            # z to the ground
            maximum_y = max(list_y_coordinates_of_frame)
            for coordinate in coordinates:
                nb_point_coord = 0
                for coordinate2 in coordinates:
                    if coordinate2[0] == coordinate[0]:
                        nb_point_coord += 1
                if coordinate[1] < maximum_z and abs(coordinate[0]) < maximum_y and nb_point_coord < 2:
                    self.shape[j].__append__((coordinate[0], maximum_z))
        return self

    def conversion_coordinate_to_pdstrip(self, midship: float):
        """That function changes the x coordinate into the PD strip system of coordinate.

        :argument
        ---------
        midship: a float
            that's the middle of the ship, defined by the length between perpendicular divided by 2 (in m)

        :returns
        ---------
        self: a form object
            that's the same frames but with different x_coordinates, x-midship for every frame"""
        for frame in self.shape:
            frame.x_coordinate = frame.x_coordinate - midship
        return self

    def x_coordinates(self):
        """That functions permits to get each x coordinate of the frame to check if there's no missing frame

        :argument
        ---------
        self: a form object
            the current form

        :returns
        ---------
        x_coordinate: a list
            the list of x coordinate of the frames
        """
        x_coordinates = []
        for frame in self.shape:
            x_coordinates.append(frame.x)
        return x_coordinates

    def plotting(self):
        """That function plots the current form in 3D.

        :argument
        ----------
        self: a form object
            the current form

        :returns
        ------------
            It prints the current form in 3D"""
        list_x_coordinates = []
        list_y_coordinates = []
        list_z_coordinates = []
        for frame in self.shape:
            x = frame.x_coordinate
            for coordinate in frame.coordinates:
                y = coordinate[0]
                z = coordinate[1]
                list_x_coordinates.append(x)
                list_y_coordinates.append(y)
                list_z_coordinates.append(z)
        list_z_coordinates = np.array(list_z_coordinates)
        list_y_coordinates = np.array(list_y_coordinates)
        list_x_coordinates = np.array(list_x_coordinates)
        fig = plt.figure()
        ax = fig.gca(projection="3d")
        ax.scatter(list_x_coordinates, list_y_coordinates, list_z_coordinates, label="courbe", marker='d')
        ax.set_title('Hull shape')
        plt.tight_layout()
        plt.show()

    def plot_one_frame(self, x: float):
        """That function plots one frame of the form for an x written by the user.

        :argument
        ------------
        x: a float
            corresponds to an x coordinate that is part of the frames (in m)

        :returns
        ----------
        It plots a graph in 2D for one frame
        """
        list_x_coordinates = []
        list_y_coordinates = []
        for frame in self.shape:
            if frame.x_coordinate == x:
                for coordinate in frame.coordinates:
                    list_x_coordinates.append(coordinate[0])
                    list_y_coordinates.append(coordinate[1])
        plt.plot(list_x_coordinates, list_y_coordinates)
        plt.show()

    def checking(self):
        """That functions checks if there are not 2 same points, with the same coordinates

        :argument
        ---------
        self: a form object
            the current form

        :returns
        -----------
        It prints an error message if there are two same coordinates
        """
        for frame in self.shape:
            for i in range(len(frame.coordinates)):
                for j in range(len(frame.coordinates)):
                    if frame.coordinates[i][0] == frame.coordinates[j][0] and \
                            frame.coordinates[i][1] == frame.coordinates[j][1] and i != j:
                        if frame.coordinates[i][0] == frame.coordinates[j][0]:
                            print("pb y")
                        if frame.coordinates[i][1] == frame.coordinates[j][0]:
                            print("pb x")

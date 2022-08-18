import csv
import Masses as mass
import Shape as shape
import Frames as fr
import Loading as ld
import matplotlib.pyplot as plt
import Wave as wv
import Analysis as Anls
import Hull as Hll
import numpy as np
from numpy import inf


def conversion_for_pdstrip_xaxis(x: float, midship: float):
    """That functions converts a x coordinate with an origin at the PPAR into a coordinate with an origin at
     middle ship

     :parameter
     ----------
     x: a float
        the coordinate to convert (in meter)
     midship: a float
        the localisation of midship from the PPAR, corresponds to Lpp/2 (in meter)

    :returns
    ----------
    x: a float
        the coordinates from the middle ship (in meter)
     """
    return x - midship


def Writing_of_the_PDstrip_input_file(masses_filename: str, coordinates_filename: str, Lpp: float):
    """That function writes a file titled data_pdstrip_prev.csv with all the input data from the frames that needs the
    PDstrip program, that function needs the loading informations, the shape of the ship, and the length between
    perpendiculars.

    :parameter
    -----------
    masses_filename: a text variable
        it corresponds to a csv file
        every line are : x coordinate of the beginning, x coordinate of the end, x coordinate of the center of gravity
        y coordinate of CoG, z coordinate of the CoG
    coordinates_filename:a text variable
        it corresponds to the pias file of coordinates
    Lpp: a number
        it is the length between perpendiculars (in m)

    :returns
    ----------
    nothing, it writes a file with the data needed named "data_pdstrip.csv" """
    hull_form = shape.Form()
    hull_form.collection_of_coordinates(coordinates_filename)
    hull_form.conversion_coordinate_to_pdstrip(Lpp / 2)
    hull_form.correction_of_coordinates_for_up()
    hull_form.checking()
    hull_form.plotting()
    weight_loading = ld.Loading()
    weight_loading.collection_of_mass(masses_filename)
    weight_loading.pdstrip_coordinates_conversion(Lpp / 2)
    weight_loading.plot_loading(-Lpp / 2, Lpp / 2)
    file = open("data_pdstrip.csv", "w")  # writing of the info in the file "data_pdstrip.csv"
    # for every section we have the backward and the forward
    back_section_total_ship = conversion_for_pdstrip_xaxis(0, Lpp / 2)
    front_section = conversion_for_pdstrip_xaxis(Lpp, Lpp / 2)
    weight_of_the_current_part = weight_loading.mass_calculation_for_coordinates(back_section_total_ship, front_section)
    x_coordinate_CoG_glob, y_coordinate_CoG_glob, z_coordinate_CoG_glob = \
        weight_loading.calcul_center_of_gravity_for_coordinates(back_section_total_ship, front_section)
    radius_of_inertia_x_square, radius_of_inertia_y_square, radius_of_inertia_z_square, xy, yz, xz = \
        weight_loading.calcul_every_parameters(x_coordinate_CoG_glob, y_coordinate_CoG_glob, z_coordinate_CoG_glob,
                                               back_section_total_ship, front_section)
    data = [weight_of_the_current_part, x_coordinate_CoG_glob, y_coordinate_CoG_glob, z_coordinate_CoG_glob,
            radius_of_inertia_x_square, radius_of_inertia_y_square, radius_of_inertia_z_square,
            xy, yz, xz]
    print(x_coordinate_CoG_glob, y_coordinate_CoG_glob, z_coordinate_CoG_glob)
    for input_value in data:
        # we write every input for the section
        file.write(str(input_value) + " ")
    file.write("\n")
    n = len(hull_form.shape)
    for i in range(len(hull_form.shape) - 1):
        print(i)
        back_section = (hull_form.shape[i].x_coordinate + hull_form.shape[i + 1].x_coordinate) / 2
        front_section = conversion_for_pdstrip_xaxis(Lpp, Lpp / 2)
        weight_of_the_current_part = weight_loading.mass_calculation_for_coordinates(back_section, front_section)
        try:
            x_coordinate_CoG, y_coordinate_CoG, z_coordinate_CoG = \
                weight_loading.calcul_center_of_gravity_for_coordinates(back_section, front_section)
        except ZeroDivisionError:
            x_coordinate_CoG, y_coordinate_CoG, z_coordinate_CoG = \
                hull_form.center_of_gravity_no_mass_for_coordinates(back_section, front_section)
        radius_of_inertia_x_square, radius_of_inertia_y_square, radius_of_inertia_z_square, xy, yz, xz = \
            weight_loading.calcul_every_parameters(x_coordinate_CoG, y_coordinate_CoG, z_coordinate_CoG, back_section,
                                                   front_section)
        data = [weight_of_the_current_part, x_coordinate_CoG, y_coordinate_CoG, z_coordinate_CoG,
                radius_of_inertia_x_square, radius_of_inertia_y_square, radius_of_inertia_z_square,
                xy, yz, xz]
        for input_value in data:
            # we write every input for the section
            file.write(str(input_value) + " ")
        file.write("\n")
    file.close()  # total weight is checked
    return


masses1 = "masses1.csv"
shape1 = "barge_standaard_pias_text_file.txt"
masses2 = "test2.csv"
# masses2 = "North route.csv"
shape2 = "oural_standaard_pias_text_file.txt"
# Writing_of_the_PDstrip_input_file(masses2, shape2, 135)
D = 11000 / 4.17 * 100
pdstrip = Anls.Analysis("pdstrip.out.ok")
hull_oural = Hll.Hull(0.6, 1.8, 4.17, 5, 7.56, 0)
# hull_oural.plot_spreading(-90,270,10)
# hull_oural.plot_JONSWAP(0.1,2.5,0.08)
# res=pdstrip.m_n_improved(2, hull_oural, 1.821)
# pdstrip.print_filewritting_max_BM(hull_oural,1.821,D,0.01)
# list_wv_height=np.arange(0.5,3,0.5)
# pdstrip.max_BM_SF_func_dir(0.6,1.8,0,5,7.56,1.821,D,0.01,"SF")
# pdstrip.max_BM_SF_func_dir(0.6, 1.8, 4.17, 5, 7.56, 1.821, D, 0.01,"SF")
pdstrip.max_BM_SF_func_spd(0.6,1.8,[0,4.17],5,7.56,1.821,135,D,0.01,"SF")


# hull_oural = Hll.Hull(0.61, 1.8, 0, 5, 7.56, 230)
# pdstrip.computation_mx_BM_in_fct_mn_angle(0, 0.61, 1.8, 5, 17.56, 1.821)
# pdstrip.computation_mx_BM_in_fct_mn_angle(4.17, 0.61, 1.8, 5, 7.56, 1.821)
# hull_max0= Hll.Hull(0.61, 1.8, 4.17, 5, 7.56, 58)
# pdstrip.writing(hull_max0, 1.821)
# hull_max4=Hll.Hull(0.6, 1.8, 4.17, 5, 7.56,128)
# pdstrip.writing(hull_max4,1.821)

# hull_oural.plot_spreading(0,360,10)

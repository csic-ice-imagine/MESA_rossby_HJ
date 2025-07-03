import numpy as np
import mesa_reader
import os
import matplotlib
from matplotlib import rc

compile = False

column_depth = 200

mesh_delta_coeffs = 0.25 # This should be the same as the create and core inlists
time_delta_coeffs = 0.25 # This should be the same as the create and core inlists

fluxes = np.array([1.6576e+10,8.3787e+09,4.2353e+09,2.1408e+09,1.0822e+09,
                   5.4701e+08,2.7650e+08,1.3977e+08,7.0648e+07,3.5711e+07])
names = ["P_0_316_F_1_65e10_T_2848",
         "P_0_528_F_8_38e9_T_2401",
         "P_0_880_F_4_24e9_T_2025",
         "P_1_47_F_2_14e9_T_1707",
         "P_2_45_F_1_08e9_T_1440",
         "P_4_08_F_5_47e8_T_1214",
         "P_6_81_F_2_77e8_T_1023",
         "P_11_4_F_1_40e8_T_863",
         "P_19_0_F_7_06e7_T_728",
         "P_31_6_F_3_57e7_T_614"]

fluxes = np.array([9.1246e+09,
                    8.3787e+09,7.6938e+09,7.0648e+09,6.4873e+09,5.9570e+09,5.4701e+09,5.0229e+09,4.6123e+09,
                    4.2353e+09,3.8890e+09,3.5711e+09,3.2792e+09,3.0111e+09,2.7650e+09,2.5390e+09,2.3314e+09,
                    2.1408e+09,1.9658e+09,1.8051e+09,1.6576e+09,1.5221e+09,1.3977e+09,1.2834e+09,1.1785e+09,
                    1.0822e+09,9.9369e+08,9.1246e+08,8.3787e+08,7.6938e+08,7.0648e+08,6.4873e+08,5.9570e+08,
                    5.4701e+08,5.0229e+08,4.6123e+08,4.2353e+08,3.8890e+08,3.5711e+08,3.2792e+08,3.0111e+08,
                    2.7650e+08,2.5390e+08,2.3314e+08,2.1408e+08,1.9658e+08,1.8051e+08,1.6576e+08,1.5221e+08,
                    1.3977e+08,1.2834e+08,1.1785e+08,1.0822e+08,9.9369e+07,9.1246e+07,8.3787e+07,7.6938e+07,
                    7.0648e+07,6.4873e+07,5.9570e+07,5.4701e+07,5.0229e+07,4.6123e+07,4.2353e+07,3.8890e+07,
                    3.5711e+07])

names = ["P_0_495_F_9_12e9_T_2453",
        "P_0_527_F_8_38e9_T_2401","P_0_562_F_7_69e9_T_2351","P_0_599_F_7_06e9_T_2301","P_0_639_F_6_49e9_T_2252",
        "P_0_681_F_5_96e9_T_2205","P_0_726_F_5_47e9_T_2158","P_0_774_F_5_02e9_T_2113","P_0_825_F_4_61e9_T_2068",
        "P_0_880_F_4_24e9_T_2025","P_0_938_F_3_89e9_T_1982","P_1_000_F_3_57e9_T_1940","P_1_066_F_3_28e9_T_1899",
        "P_1_136_F_3_01e9_T_1859","P_1_212_F_2_76e9_T_1820","P_1_292_F_2_54e9_T_1782","P_1_377_F_2_33e9_T_1744",
        "P_1_468_F_2_14e9_T_1707","P_1_565_F_1_97e9_T_1671","P_1_668_F_1_81e9_T_1636","P_1_778_F_1_66e9_T_1601",
        "P_1_896_F_1_52e9_T_1568","P_2_021_F_1_40e9_T_1535","P_2_154_F_1_28e9_T_1502","P_2_297_F_1_18e9_T_1471",
        "P_2_448_F_1_08e9_T_1440","P_2_610_F_9_94e8_T_1409","P_2_783_F_9_12e8_T_1379","P_2_966_F_8_38e8_T_1350",
        "P_3_162_F_7_69e8_T_1322","P_3_371_F_7_06e8_T_1294","P_3_594_F_6_49e8_T_1267","P_3_831_F_5_96e8_T_1240",
        "P_4_084_F_5_47e8_T_1214","P_4_354_F_5_02e8_T_1188","P_4_642_F_4_61e8_T_1163","P_4_948_F_4_24e8_T_1139",
        "P_5_275_F_3_89e8_T_1115","P_5_623_F_3_57e8_T_1091","P_5_995_F_3_28e8_T_1068","P_6_391_F_3_01e8_T_1046",
        "P_6_813_F_2_76e8_T_1023","P_7_263_F_2_54e8_T_1002","P_7_743_F_2_33e8_T_981","P_8_254_F_2_14e8_T_960",
        "P_8_799_F_1_97e8_T_940","P_9_380_F_1_81e8_T_920","P_10_000_F_1_66e8_T_901","P_10_661_F_1_52e8_T_882",
        "P_11_365_F_1_40e8_T_863","P_12_115_F_1_28e8_T_845","P_12_915_F_1_18e8_T_827","P_13_769_F_1_08e8_T_809",
        "P_14_678_F_9_94e7_T_792","P_15_647_F_9_12e7_T_776","P_16_681_F_8_38e7_T_759","P_17_783_F_7_69e7_T_743",
        "P_18_957_F_7_06e7_T_728","P_20_209_F_6_49e7_T_712","P_21_544_F_5_96e7_T_697","P_22_967_F_5_47e7_T_683",
        "P_24_484_F_5_02e7_T_668","P_26_102_F_4_61e7_T_654","P_27_826_F_4_24e7_T_640","P_29_663_F_3_89e7_T_627",
        "P_31_623_F_3_57e7_T_614"]

fluxes = fluxes[::-1]
names = names[::-1]

if compile:
    os.system("./clean")
    os.system("./mk")

for index,flux in enumerate(fluxes):

    file_read = open('inlist_relax', "r")
    lines = file_read.readlines()
    index_line=0
    for line in lines:
        if '      relax_to_this_irrad_flux = ' in line:
            lines[index_line] = "      relax_to_this_irrad_flux = " + "{:1.4e}".format(flux) + " ! in erg/cm^2/s\n"
        if '      irrad_col_depth = ' in line:
            lines[index_line] = "      irrad_col_depth = " + "{:1.4e}".format(column_depth) + " ! in erg/cm^2/s\n"
        if '    mesh_delta_coeff = ' in line:
            lines[index_line] = "    mesh_delta_coeff = " + "{:1.4e}".format(mesh_delta_coeffs) + " ! space resolution \n"
        if '    time_delta_coeff = ' in line:
            lines[index_line] = "    time_delta_coeff = " + "{:1.4e}".format(time_delta_coeffs) + " ! time resolution \n"
        index_line+=1
    file_read.close()
    with open('inlist_relax', 'w') as file: 
        file.writelines(lines)

    file_read = open('inlist_evolve', "r")
    lines = file_read.readlines()
    index_line=0
    for line in lines:
        if '    irradiation_flux = ' in line:
            lines[index_line] = "    irradiation_flux = " + "{:1.4e}".format(flux) + " ! in erg/cm^2/s\n"
        if '    column_depth_for_irradiation = ' in line:
            lines[index_line] = "	column_depth_for_irradiation = " + "{:1.4e}".format(column_depth) + " ! in erg/cm^2/s\n"
        if '    mesh_delta_coeff = ' in line:
            lines[index_line] = "    mesh_delta_coeff = " + "{:1.4e}".format(mesh_delta_coeffs) + " ! space resolution \n"
        if '    time_delta_coeff = ' in line:
            lines[index_line] = "    time_delta_coeff = " + "{:1.4e}".format(time_delta_coeffs) + " ! time resolution \n"
        index_line+=1
    file_read.close()
    with open('inlist_evolve', 'w') as file: 
        file.writelines(lines)

    file_read = open('inlist_pgstar', "r")
    lines = file_read.readlines()
    index_line=0
    for line in lines:
        if '    Grid1_file_dir = ' in line:
            lines[index_line] = "    Grid1_file_dir = 'png_5RJ_" + names[index] + "'\n"
        index_line+=1
    file_read.close()
    with open('inlist_pgstar', 'w') as file:
        file.writelines(lines) 
    
    file_read = open('rn', "r")
    lines = file_read.readlines()
    index_line=0
    for line in lines:
        if 'do_one inlist_evolve_header' in line:
            lines[index_line] = "do_one inlist_evolve_header planet_evolve_1.0_MJ_10.0_ME_5.0_RJ.mod LOGS_5RJ_" + names[index] + "\n"
        index_line+=1
    file_read.close()
    with open('rn', 'w') as file: 
        file.writelines(lines)

    os.system("./rn")


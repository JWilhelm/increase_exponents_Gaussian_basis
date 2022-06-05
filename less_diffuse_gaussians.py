#/bin/bahs
import sys
import numpy as np
import os
import re

orig_basis_set_file_name = "/dss/dsshome1/07/di76zil/30_CP2K_ham_k_from_Gamma_and_make_ADMM_in_Sigx_work/cp2k/data/BASIS_MOLOPT"
new_basis_set_file_name  = "./BASIS_SET_MoS_MOLOPT_robust"
basis_sets_to_modify     = ["DZVP-MOLOPT-SR-GTH"]
elements                 = ['S', 'Mo']
desired_min_exponent     = 0.13
min_expansion_exponent   = 0.13
max_expansion_exponent   = 0.15
delta_exponent           = 0.02

#-----------------------------------------------
# from here, no parameters need to be modified |
#-----------------------------------------------

def get_inv_mat_S(min_expansion_exponent, max_expansion_exponent, delta_exponent):
    small_delta = 1.0e-6
    list_exponents = np.arange(min_expansion_exponent, max_expansion_exponent+small_delta, delta_exponent)
    num_exponents = np.size(list_exponents)
    mat_S = np.zeros([num_exponents, num_exponents])
    for i in range(num_exponents):
        for j in range(num_exponents):
            mat_S[i,j] = 1.0/np.sqrt(list_exponents[i] + list_exponents[j])

    inv_mat_S = np.linalg.inv(mat_S)

    print("inv_mat_S")
    print(inv_mat_S)

    identity = np.matmul(inv_mat_S, mat_S)

    print("inv_mat_S * mat_S")
    print(identity)

    return inv_mat_S, list_exponents, num_exponents

basis_set_file  = open(orig_basis_set_file_name, 'r')

inv_mat_S, list_exponents, num_exponents = \
        get_inv_mat_S(min_expansion_exponent, max_expansion_exponent, delta_exponent)

basis_set_found = False
basis_set_lines_so_far = []
basis_set_lines_new = []
basis_set_lines_new_int = []
unchanged_exponent_lines = []
num_new_gaussians = 0

for line in basis_set_file:
    split_line = line.split(" ")
    # remove empty entries originating from white spaces
    split_line = list(filter(None, split_line))

    # we have found next basis set, no numbers any more
    if split_line[0].isalpha() and basis_set_found:

        # print everything
        with open(new_basis_set_file_name, 'a') as the_file:
            the_file.write(header_line)
            the_file.write(basis_set_lines_so_far[0])

            split_third_line = basis_set_lines_so_far[1].split(" ")
            split_third_line = list(filter(None, split_third_line))

            split_third_line[3] = str(int(split_third_line[3]) + num_new_gaussians)

            third_line_string = ' '.join(map(str, split_third_line)) 

            the_file.write(third_line_string)

            for count, unchanged_line in enumerate(basis_set_lines_so_far):
                if count > 1: 
                    the_file.write(unchanged_line)

            for new_line in basis_set_lines_new_int:
                new_line_to_write=("      "+ ' '.join(map(str, new_line)) + ' \n')

                the_file.write(new_line_to_write)

            the_file.write("\n")

        basis_set_found = False
        basis_set_lines_so_far = []
        basis_set_lines_new = []
        basis_set_lines_new_int = []
        unchanged_exponent_lines = []
        num_new_gaussians = 0
        first_replacement = False

    # modify the basis set
    if basis_set_found: 
        exponent = float(split_line[0])

        print(exponent)
        if exponent < desired_min_exponent:

            if num_new_gaussians == 0:
               num_new_gaussians = num_exponents-1
               first_replacement = True
            else:
               num_new_gaussians = num_new_gaussians-1
               first_replacement = False

            print(split_line[1:])

            old_expansion_coeff = np.array([float(x) for x in split_line[1:] ])

            vec_P = 1.0/np.sqrt(np.ones(num_exponents)*exponent + list_exponents)

            expansion_coeff = inv_mat_S.dot(vec_P)

            for count, exponent_new in enumerate(list_exponents):
                if first_replacement:
                    basis_set_line_new = np.zeros(np.size(old_expansion_coeff)+1)
                    basis_set_line_new[0] = exponent_new
                    basis_set_line_new[1:] = old_expansion_coeff[:]*expansion_coeff[count]
                    basis_set_lines_new_int.append(basis_set_line_new)

                else:
                    basis_set_lines_new_int[count][1:] = basis_set_lines_new_int[count][1:] + old_expansion_coeff[:]*expansion_coeff[count]

        elif exponent > desired_min_exponent and split_line[0].isdigit:

            basis_set_lines_so_far.append(line)

    # check whether a new basis set needs to be modified
    if split_line[0] in elements:
        if split_line[1] in basis_sets_to_modify:
            basis_set_found = True

            header_line = split_line[0]+" "+split_line[1]+'_robust\n'


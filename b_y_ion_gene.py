"""
Calculate the theoritical b and y ions given one peptide sequence
"""
import numpy as np


def b_y_ion_gen(peptide_seq_charge):
    """
    generate a list of tuple with b ion and y ion
    :param peptide_seq_charge: PEPTIDEcharge
    :return:
    """
    from aa_mass_table import aa_mass_table,h_oh_mass_dict
    peptide_seq = peptide_seq_charge[:-1]
    ion_list = []
    rev_peptide_seq = peptide_seq[::-1]

    for i in range(len(peptide_seq)):
        b_ion = h_oh_mass_dict['H']  # N-terminal adding H
        y_ion = h_oh_mass_dict['OH']  # C-terminal adding OH
        pep = peptide_seq[:i+1]
        rev_pep = rev_peptide_seq[:i+1]

        for each_pep in pep:
            b_ion += aa_mass_table[each_pep]

        for each_rev in rev_pep:
            y_ion += aa_mass_table[each_rev]
        ion_list.append((b_ion,y_ion))

    return ion_list


def b_y_ion_bins_gen(b_y_ion_list, ppm=50):
    """
    get b/y ions into bins
    :param b_y_ion_list: returned by last function
    :param ppm: default is 50
    :return:
    """
    # add ppm to each b/y ions
    b_y_ion_tolerance_list = [(each*(1-ppm/1000000),each*(1+ppm/1000000))
                              for each_tp in b_y_ion_list
                              for each in each_tp]

    # sorted the bin as monotonic
    b_y_ion_monotonic_bins = sorted([each
                                     for each_tp in b_y_ion_tolerance_list
                                     for each in each_tp])

    return b_y_ion_monotonic_bins


def dump_mass_into_ion_bins(mass_array, ion_bins):
    """
    dump b and y ions of a peptide into ion bins
    :param mass_array:
    :param ion_bins:
    :return: index of each dumped mass
    """
    bin_index = np.digitize(mass_array,ion_bins)

    return bin_index


def vector_gen(int_array,bin_index,b_y_ion_monotonic_bins):
    """
    get a vector for cosine similarity computation
    :param int_array:
    :param bin_index: returned by function dump_mass_into_ion_bins
    :param b_y_ion_monotonic_bins: returned by function b_y_ion_bins_gen
    :return: a vector with same size as b/y ion list
    """
    v1 = np.zeros(len(b_y_ion_monotonic_bins)+1)
    for ind, val in enumerate(bin_index):
        # if val%2 == 0:
        #     continue
        # else:
        #print(ind, val)
        v1[val] += int_array[ind]
    #print(v1)
    v1 = [v1[i] for i in range(1, len(v1), 2)]  # get ppm tolerance
    return v1


def cosine_sim_calculating(v1, v2):
    """
    calculate the cosine similarity beweeen two b/y ions binned vectors
    :param v1:
    :param v2:
    :return:
    """
    from scipy import spatial
    return 1-spatial.distance.cosine(v1,v2)


def single_usage(peptide_seq, predicted_mass_array, predicted_int_array, precursor_mass_array, precursor_int_array,ppm:int):
    """
    given a peptide, get cosine sim between predicted and actual ms2 after binned into b/y ions
    :param peptide_seq: peptide sequence, string
    :param predicted_mass_array: 1D array, from prosit prediction for the peptide
    :param predicted_int_array: 1D array, from prosit prediction for the peptide
    :param precursor_mass_array: 1D array, from MS2 spectrum
    :param precursor_int_array: 1D array, from MS2 spectrum
    :return: cosine similarity for b/y ions binned ms2 between predicted and precursor ms2 spectrum
    """
    ion_list = b_y_ion_gen(peptide_seq)  # get theoretical ion mass
    b_y_ion_bins = b_y_ion_bins_gen(ion_list,ppm=ppm) # get ion mass bin edges

    predicted_bin_index = dump_mass_into_ion_bins(predicted_mass_array,b_y_ion_bins)  # get bin index for each mass
    precusor_bin_index = dump_mass_into_ion_bins(precursor_mass_array,b_y_ion_bins)  #

    v_predicted = vector_gen(predicted_int_array,predicted_bin_index,b_y_ion_bins) # get binned intensity array
    v_precusor = vector_gen(precursor_int_array,precusor_bin_index,b_y_ion_bins) # get binned intensity array

    return cosine_sim_calculating(v_predicted,v_precusor)


if __name__ == '__main__':
    import numpy as np
    peptide ='NVIFLNK'

    b_y_ion_list = b_y_ion_gen(peptide)

    b_y_ion_clean_list = [(each*(1-100/1000000),each*(1+100/1000000)) for each_tp in b_y_ion_list for each in each_tp ]  # 50 ppm

    b_y_ion_monotonic_bins = sorted([each for each_tp in b_y_ion_clean_list for each in each_tp])
    print (b_y_ion_monotonic_bins)
    print (len(b_y_ion_monotonic_bins))
    test_array = np.array([147.11281, 115.0502, 261.15573, 214.11862, 374.2398, 327.20267, 521.3082, 261.15775, 474.2711, 237.63919, 634.3923, 317.69977, 733.4607])
    int_array = [0.11995926, 0.00087668106, 0.64565444, 0.48343208, 0.45107266, 0.08600334, 0.7407538, 0.009013158, 0.008180862, 0.00033641688, 1.0, 0.0143297, 0.042357314]
    # print (int_array)
    bin_index = np.digitize(test_array,b_y_ion_monotonic_bins)
    print (bin_index)
    v1 = np.zeros(len(b_y_ion_monotonic_bins)+1)
    for ind,val in enumerate(bin_index):
        # if val%2 == 0:
        #     continue
        # else:
        print (ind,val)
        v1[val] += int_array[ind]
    print (v1)
    v1 = [v1[i] for i in range(1,len(v1),2)]
    print (v1)



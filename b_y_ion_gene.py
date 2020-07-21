"""
Calculate the theoritical b and y ions given one peptide sequence
"""
import numpy as np


def b_y_ion_gen(peptide_seq):
    """
    generate a list of tuple with b ion and y ion
    :param peptide_seq:
    :return:
    """
    from aa_mass_table import aa_mass_table,h_oh_mass_dict
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


def vector_gen(int_array,bin_index,b_y_ion_list):
    """
    get a vector for cosine similarity computation
    :param int_array:
    :param bin_index: returned by function dump_mass_into_ion_bins
    :return: a vector with same size as b/y ion list
    """
    v1 = np.zeros(len(b_y_ion_list))
    for ind,val in enumerate(bin_index):
        if val%2 == 0:  #
            continue
        else:
            v1[val-1] += int_array[ind]
    return v1

if __name__=='__main__':
    import numpy as np
    peptide ='GPGTYDTSDPSIYK'

    b_y_ion_list = b_y_ion_gen(peptide)

    b_y_ion_clean_list = [(each*(1-50/1000000),each*(1+50/1000000)) for each_tp in b_y_ion_list for each in each_tp ]

    b_y_ion_monotonic_bins = sorted([each for each_tp in b_y_ion_clean_list for each in each_tp])
    print (b_y_ion_monotonic_bins)
    test_array = np.array([57.9,58.02333,58.0289,69.9, 147.10001,155.0791,200.08,212.1005,310.1777])
    int_array = np.random.rand(9)
    print (int_array)
    bin_index = np.digitize(test_array,b_y_ion_monotonic_bins)
    print (bin_index)
    v1 = np.zeros(len(b_y_ion_list))
    for ind,val in enumerate(bin_index):
        if val%2 == 0:
            continue
        else:
            v1[val-1] += int_array[ind]
    print (v1)




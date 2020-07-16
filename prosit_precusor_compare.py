from msp_reader import msp_reader, prosit_csv_output
from ms2_reader import target_peptide_file_spec_getter, target_ms2_info_reader
from psm_reader import peptide_file_spec_dict_of_dict
from collections import defaultdict
from glob import glob
import numpy as np

pep_file = ''  # tsv or dta file
csv_out = ''
peptide_list = []  # target peptide list
prosit_csv = prosit_csv_output(peptide_list,csv_out,pep_file) # output csv file subject to prosit prediction


msp_file_path = ''
prosit_info_dict = msp_reader(msp_file_path)  # info dictionary with predicted ms2 m/z intensity array

# prosit dictionary key and value change
new_prosit_info_dict = defaultdict(list)  # make pep_seq as key, (m/z array, int array) as value
for each in prosit_info_dict:
    pep_seq, mz_array, int_array = prosit_info_dict[each][0], prosit_info_dict[each][-2], prosit_info_dict[each][-1]
    new_prosit_info_dict[pep_seq].append((mz_array,int_array))

# target peptides file and spectra info
psm_path = ''
target_pep_file_spec_dict_of_dict = {each:peptide_file_spec_dict_of_dict(psm_path)[each] for each in peptide_list}

# compare similarity between prosit prediction and acquired precusor ms2
ms2_dict_of_dict = {}
cos_sim_score_dict = {}
for each_pep in new_prosit_info_dict:
    for each_tuple in new_prosit_info_dict[each_pep]:
        predicted_mz_array, predicted_int_array = each_tuple
        v1 = np.zeros(280000)  # m/z from 200-3000, 100 bins, make a int vector to compare similarity.
        v2 = np.zeros(280000)
        for each_pred_mz, each_pred_int in zip(predicted_mz_array,predicted_int_array):

        file_spec_list_dict = target_pep_file_spec_dict_of_dict[each_pep]  # file_spec_list dict {file:[spec1, spec2],}
        for each_file in file_spec_list_dict:
            for each_spec in file_spec_list_dict[each_file]:
                precusor_mz_array, precusor_int_array = ms2_dict_of_dict[each_file][each_spec][-2:]

def ms2_info_dict_generator(psm_tsv_path, target_pep_list, ms2_path, pickle_saved_path = None):
    """
    write ms2 info of target peptides in a dictionary of dictionary, {file:{spec1:(),spec2(),...}}
    :param psm_tsv_path:
    :param target_pep_list:
    :param ms2_path:
    :param pickle_saved_path:
    :return:
    """
    ms2_list = glob(ms2_path + '*_clean.ms2')

    if '_clean.ms2' not in ms2_list[0]:
        raise ValueError('ms2 file not clean, use ms2_cleaner.py to clean it')

    file_spec_no_list_dict = target_peptide_file_spec_getter(target_pep_list, psm_tsv_path)
    print(file_spec_no_list_dict)

    ms2_dict_of_dict = {
        each: target_ms2_info_reader(each, file_spec_no_list_dict[each.split('\\')[-1].split('_clean')[0]]) for each in
        ms2_list}

    if pickle_saved_path:
        import pickle as ppp
        ppp.dump(ms2_dict_of_dict,open(pickle_saved_path,'wb'))

    return ms2_dict_of_dict


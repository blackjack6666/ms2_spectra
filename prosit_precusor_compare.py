from msp_reader import msp_reader, prosit_csv_output
from ms2_reader import target_peptide_file_spec_getter, target_ms2_info_reader
from psm_reader import peptide_file_spec_dict_of_dict
from collections import defaultdict
from glob import glob
import numpy as np
from scipy import spatial
import pickle as ppp

peptide_list = ppp.load(open('C:/Users/gao lab computer/PycharmProjects/extend_different_species/PXD001723_ext_pep_list.p','rb'))  # target peptide list
peptide_list = [each for each in peptide_list if len(each)<=30]
"""
pep_file = 'D:/data/ext_evo_pj/gb_ext_search_7_11_PXD001723/peptide.tsv'  # tsv or dta file
csv_out = 'D:/data/ext_evo_pj/gb_ext_search_7_11_PXD001723/PXD001723_ext_prosit.csv'
prosit_csv = prosit_csv_output(peptide_list,csv_out,pep_file) # output csv file subject to prosit prediction
"""

msp_file_path = 'D:/data/ext_evo_pj/gb_ext_search_7_11_PXD001723/myPrositLib.msp'
prosit_info_dict = msp_reader(msp_file_path)  # info dictionary with predicted ms2 m/z intensity array

# prosit dictionary key and value change
new_prosit_info_dict = defaultdict(list)  # make pep_seq as key, (m/z array, int array) as value
for each in prosit_info_dict:
    pep_seq, mz_array, int_array = prosit_info_dict[each][0], prosit_info_dict[each][-2], prosit_info_dict[each][-1]
    new_prosit_info_dict[pep_seq].append((mz_array,int_array))

print (new_prosit_info_dict)

# target peptides file and spectra info
psm_path = 'D:/data/ext_evo_pj/gb_ext_search_7_11_PXD001723/psm.tsv'
pep_file_spec_dict_of_dict = peptide_file_spec_dict_of_dict(psm_path)
print ('1')
target_pep_file_spec_dict_of_dict = {each:pep_file_spec_dict_of_dict[each] for each in peptide_list}
print (target_pep_file_spec_dict_of_dict)

# compare similarity between prosit prediction and acquired precusor ms2
ms2_dict_of_dict = ppp.load(open('D:/data/ext_evo_pj/gb_ext_search_7_11_PXD001723/PXD001723_ms2_dict_of_dict_7_13.p','rb'))
cos_sim_score_dict = defaultdict(list)
for each_pep in new_prosit_info_dict:
    print (each_pep)

    # file_spec_list dict {file:[spec1, spec2],}
    file_spec_list_dict = target_pep_file_spec_dict_of_dict[each_pep]

    # get m/z, int array for peptide seq from prosit
    for each_tuple in new_prosit_info_dict[each_pep]:
        predicted_mz_array, predicted_int_array = each_tuple
        v1 = np.zeros(290000)  # m/z from 200-3000, 100 bins, make a int vector to compare similarity.

        for each_pred_mz, each_pred_int in zip(predicted_mz_array,predicted_int_array):
            v1[int((each_pred_mz-200)*100)] += each_pred_int

        for each_file in file_spec_list_dict:
            for each_spec in file_spec_list_dict[each_file]:
                print (each_file,each_spec)
                precusor_mz_array, precusor_int_array = ms2_dict_of_dict['F:/XS/c_elegans/PXD001723'+'\\'+each_file+'_clean.ms2'][each_spec][-2:]
                max_int = max(precusor_int_array)
                # normalize the intensity array to be compatible with prosit result
                precusor_int_array = [float(each)/max_int for each in precusor_int_array]
                v2 = np.zeros(290000)

                for each_mz, each_int in zip(precusor_mz_array,precusor_int_array):
                    v2[int((each_mz-200)*100)] += each_int
                cos_sim = 1-spatial.distance.cosine(v1,v2)
                print (cos_sim)
                cos_sim_score_dict[each_pep].append((each_file+'_'+str(each_spec),cos_sim))
print (cos_sim_score_dict)

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


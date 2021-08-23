"""
usage: 1. use prosit_csv_input to generate a csv file for your interested peptides subject to prosit ms2 spectra predict
       2. follow the steps at the end of this script
"""


from msp_reader import msp_reader, prosit_csv_output
from ms2_reader import target_peptide_file_spec_getter, target_ms2_info_reader
from psm_reader import peptide_file_spec_dict_of_dict
from collections import defaultdict
from glob import glob
import numpy as np
from scipy import spatial
from predfull import mgf_file_reader
import pickle as ppp


def prosit_csv_input(pep_file,csv_out,peptide_list):
    """
    write peptide seq and charge info into a csv file subject to prosit predict
    :param pep_file: could be tsv or dta file
    :param csv_out:
    :param peptide_list:
    :return: csv file for prosit prediction
    """
    prosit_csv = prosit_csv_output(peptide_list,csv_out,pep_file) # output csv file subject to prosit prediction
    return prosit_csv


def msp_info_dict_gen(msp_file_path, file='msp'):
    """
    get info from prosit output msp file into dictionary
    :param msp_file_path:
    :return: dictionary with peptide as key and list of m/z array and int array tuple as value
    """
    if file=='msp':
        prosit_info_dict = msp_reader(msp_file_path)  # info dictionary with predicted ms2 m/z intensity array
    elif file=='mgf':
        prosit_info_dict = mgf_file_reader(msp_file_path)
    # prosit dictionary key and value change
    new_prosit_info_dict = defaultdict(list)  # make pep_seq as key, (m/z array, int array) as value
    for each in prosit_info_dict:
        pep_seq, charge, mz_array, int_array = prosit_info_dict[each][0], prosit_info_dict[each][1], prosit_info_dict[each][-2], prosit_info_dict[each][-1]
        new_prosit_info_dict[pep_seq+str(charge)].append((mz_array,int_array))
    return new_prosit_info_dict


def target_pep_files_spectra_gen(peptide_list,psm_path):
    """
    get a dictionary of target peptides with file and spectras
    :param peptide_list:
    :param psm_path:
    :return: dict of dict, {peptide+charge:{file:[spectra_#1, spectra_#2]}}
    """

    pep_file_spec_dict_of_dict = peptide_file_spec_dict_of_dict(psm_path)

    target_pep_file_spec_dict_of_dict = {each:pep_file_spec_dict_of_dict[each] for each in peptide_list}
    return target_pep_file_spec_dict_of_dict


def cosine_similarity_compare(new_prosit_info_dict, target_pep_file_spec_dict_of_dict, ms2_dict_of_dict):
    """
    compare cosine similarity between prosit prediciton and acquired precursor
    :param new_prosit_info_dict:
    :param target_pep_file_spec_dict_of_dict:
    :param ms2_dict_of_dict: returned by function ms2_info_dict_generator, a dictionary with file name as key, info_dict as value
    :return:
    """

    cos_sim_score_dict = defaultdict(list)
    for each_pep in new_prosit_info_dict:
        print (each_pep)

        # file_spec_list dict {file:[spec1, spec2],}
        file_spec_list_dict = target_pep_file_spec_dict_of_dict[each_pep]

        # get mass, int array for peptide seq from prosit
        for each_tuple in new_prosit_info_dict[each_pep]:
            predicted_mz_array, predicted_int_array = each_tuple
            v1 = np.zeros(290000)  # m/z from 200-3000, 100 bins, make a int vector to compare similarity.

            for each_pred_mz, each_pred_int in zip(predicted_mz_array,predicted_int_array):
                v1[int((each_pred_mz-200)*100)] += each_pred_int

            for each_file in file_spec_list_dict:
                for each_spec in file_spec_list_dict[each_file]:
                    print (each_file,each_spec)
                    precusor_mz_array, precusor_int_array = ms2_dict_of_dict['C:/uic/lab/mankin/ms2_files/api_ms2/api05'+'\\'+each_file+'_clean.ms2'][each_spec][-2:]
                    max_int = max(precusor_int_array)
                    # normalize the intensity array to be compatible with prosit result
                    precusor_int_array = [float(each)/max_int for each in precusor_int_array]
                    v2 = np.zeros(290000)

                    for each_mz, each_int in zip(precusor_mz_array,precusor_int_array):
                        v2[int((each_mz-200)*100)] += each_int
                    cos_sim = 1-spatial.distance.cosine(v1,v2)
                    print (cos_sim)
                    cos_sim_score_dict[each_pep].append((each_file+'_'+str(each_spec),cos_sim))
    #print (cos_sim_score_dict)
    return cos_sim_score_dict


def b_y_ion_cos_sim_compare(new_prosit_info_dict, target_pep_file_spec_dict_of_dict, ms2_dict_of_dict,ppm:int):
    """
    calculate cosine sim score in batch mode
    :param new_prosit_info_dict:
    :param target_pep_file_spec_dict_of_dict:
    :param ms2_dict_of_dict:
    :return: cosine sim dictionary with pep seq as key, a list of tuple as value
    """
    from b_y_ion_gene import single_usage
    cos_sim_score_dict = defaultdict(list)
    for each_pep in new_prosit_info_dict:
        print(each_pep)

        # file_spec_list dict {file:[spec1, spec2],}
        file_spec_list_dict = target_pep_file_spec_dict_of_dict[each_pep]

        # get mass, int array for peptide seq from prosit
        for each_tuple in new_prosit_info_dict[each_pep]:
            predicted_mz_array, predicted_int_array = each_tuple

            for each_file in file_spec_list_dict:
                for each_spec in file_spec_list_dict[each_file]:
                    print(each_file, each_spec)
                    precusor_mz_array, precusor_int_array = \
                    ms2_dict_of_dict['C:/uic/lab/mankin/ms2_files/api_ms2/api05' + '\\' + each_file + '_clean.ms2'][each_spec][-2:]
                    max_int = max(precusor_int_array)

                    # normalize the intensity array to be compatible with prosit result
                    precusor_int_array = [float(each) / max_int for each in precusor_int_array]
                    # get b_y_ion binned cosine similarity score
                    b_y_binned_cos_score = single_usage(each_pep,predicted_mz_array,predicted_int_array,precusor_mz_array,precusor_int_array,ppm)

                    cos_sim_score_dict[each_pep].append((each_file + '_' + str(each_spec), b_y_binned_cos_score))
    # print (cos_sim_score_dict)
    return cos_sim_score_dict

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


if __name__=='__main__':
    import b_y_ion_gene
    import matplotlib.pyplot as plt
    from predfull import mgf_file
    import pickle as ppp
    msp_file_path = 'D:/uic/lab/code/PredFull-master/test_prediction.mgf'
    msp_info_dict = msp_info_dict_gen(msp_file_path, file='mgf')

    # print('predicted',[i for i in zip(*msp_info_dict['ESTIDETTRYGPI'][0])])
    # # b/y compare
    # ion_list = b_y_ion_gene.b_y_ion_gen('ESTIDETTRYGPI')
    # b_y_bins = b_y_ion_gene.b_y_ion_bins_gen(ion_list,ppm=50)
    # print ('bins',b_y_bins)
    # mass_array,int_array = msp_info_dict['ESTIDETTRYGPI'][0]
    # print ('predicted mass',mass_array)
    # print ('predicted int',int_array)
    # bin_index = b_y_ion_gene.dump_mass_into_ion_bins(mass_array,b_y_bins)
    # print (bin_index)
    # v_predicted = b_y_ion_gene.vector_gen(int_array,bin_index,b_y_bins)
    # print ('predicted vector',v_predicted)


    # peptide_list = ppp.load(
    #     open('C:/Users/gao lab computer/PycharmProjects/extend_different_species/PXD001723_ext_pep_list.p',
    #          'rb'))  # target peptide list
    peptide_list = ['AEHLVFWNGGR2','AEHLVFWNGGR3','VPVTDESPATR2','WKNPTPSYSK2']
    peptide_list = [each for each in peptide_list if len(each) <= 30] # peptides longer than 30aa are not compatible with prosit
    print (len(peptide_list))

    ms2_dict_of_dict = ppp.load(open('D:/uic/lab/mankin/ms2_files/api_ms2/api05/api05_ms2_dict_of_dict.p','rb'))

    #print ('precursor',[i for i in zip(*ms2_dict_of_dict['F:/XS/c_elegans/PXD001364'+'\\20091013_Velos3_DiWa_SA_Celegans_Hsf1-Day12_Offgel07_clean.ms2'][7473][5:7])])

    # prec_mass_array, prec_int_array = ms2_dict_of_dict['F:/XS/c_elegans/PXD001364'+'\\20091013_Velos3_DiWa_SA_Celegans_Hsf1-Day12_Offgel01_clean.ms2'][10159][5:7]
    # max_prec_int = max(prec_int_array)
    # prec_int_array = [float(each)/max_prec_int for each in prec_int_array]
    # print ('precusor', [i for i in zip(*(prec_mass_array,prec_int_array))])
    # prec_bin_index = b_y_ion_gene.dump_mass_into_ion_bins(prec_mass_array,b_y_bins)
    # v_precusor = b_y_ion_gene.vector_gen(prec_int_array,prec_bin_index,b_y_bins)
    # print ('precusor vector',v_precusor)
    # print (1-spatial.distance.cosine(v_predicted,v_precusor))
    # print (b_y_ion_gene.single_usage('ESTIDETTRYGPI',mass_array,int_array,prec_mass_array,prec_int_array,ppm=50))

    psm_path = 'D:/uic/lab/mankin/20200302_3_2_db_search/api05/psm.tsv'  # one psm file could include multiple file search result
    target_pep_file_spec_dict_of_dict = target_pep_files_spectra_gen(peptide_list,psm_path)

    b_y_ion_binned_cos_sim_dict = b_y_ion_cos_sim_compare(msp_info_dict,target_pep_file_spec_dict_of_dict,ms2_dict_of_dict,50)
    print(len(b_y_ion_binned_cos_sim_dict), b_y_ion_binned_cos_sim_dict)

    # caculate avearge cosine score for each pep
    cosine_average_dict = {}
    for each_pep in b_y_ion_binned_cos_sim_dict:  # each_pep: peptide_seq + charge
        total = 0
        for each_tp in b_y_ion_binned_cos_sim_dict[each_pep]:
            total+=each_tp[1]
        aver = float(total)/len(b_y_ion_binned_cos_sim_dict[each_pep])
        cosine_average_dict[each_pep]=aver
    print (cosine_average_dict)
    # plt.hist([v for v in cosine_average_dict.values()])
    # ppp.dump(cosine_average_dict,open('PXD001723_ext_pep_cosScore_dict.p','wb'))
    # print (target_pep_file_spec_dict_of_dict)
    # cosine_score_dict = cosine_similarity_compare(msp_info_dict,target_pep_file_spec_dict_of_dict,ms2_dict_of_dict)
    # print (len(cosine_score_dict))
    #
    # high_conf_ext_pep_list = [(pep,each_tu[0],each_tu[1]) for pep in cosine_score_dict for each_tu in cosine_score_dict[pep] if each_tu[1]>0.5]
    # print (high_conf_ext_pep_list)

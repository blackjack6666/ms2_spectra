"""
plot ms spectrum from ms2 files
"""
import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from pyteomics import mgf
import pickle as ppp
from psm_reader import peptide_file_spectra_generator, dta_file_spectra_generator
from collections import defaultdict
from glob import glob
import numpy as np
import pandas as pd
import os

'''
# Read the spectrum from an MGF file using Pyteomics.
spectrum_dict = mgf.get_spectrum(
    'spectra.mgf',
    'mzspec:PXD004732:01650b_BC2-TUM_first_pool_53_01_01-3xHCD-1h-R2:scan:'
    '41840:WNQLQAFWGTGK/2')
identifier = spectrum_dict['params']['title']
precursor_mz = spectrum_dict['params']['pepmass'][0]
precursor_charge = spectrum_dict['params']['charge'][0]
mz = spectrum_dict['m/z array']
intensity = spectrum_dict['intensity array']
retention_time = float(spectrum_dict['params']['rtinseconds'])
peptide = 'WNQLQAFWGTGK'
'''

def ms2_visulizer(ms2_info_dict:dict,spectra_no:int, peptide:str,saved_file_path:str,input_file:str, peptide_with_mod):
    """
    plot the ms2 peak give spectrum number and peptide sequence
    -----
    :param ms2_info_dict: spectra_no as key, (mz, ret_time, charge1,charge2, mass_array,intensity arrary) as value
    :param spectra_no: spectrum number
    :param peptide: the peptide to be visulized
    :return: ms2 peak plot
    """

    m_over_z,ret_time,charge1,charge2,pre_cursor_mass,m_array,int_array = ms2_info_dict[spectra_no]

    m_array = np.array(m_array, dtype=float)
    int_array = np.array(int_array,dtype=float)
    # add one value before and after m_array and int_array to make the range same

    m_array = np.append(m_array,2000)
    m_array = np.insert(m_array,0,1)
    int_array = np.append(int_array,np.max(int_array)*0.006)
    int_array = np.insert(int_array,0,np.max(int_array)*0.006)
    #print(m_array, int_array)
    # Create the MS/MS spectrum.
    test_spec = sus.MsmsSpectrum('identifier', m_over_z, int(charge1), m_array, int_array)
    pep_seq_1 = peptide.replace('M(ox)', 'M[+15.9949]')
    test_spec.annotate_proforma(pep_seq_1, fragment_tol_mass=50, fragment_tol_mode="ppm", ion_types="aby")

    fig, ax = plt.subplots(figsize=(12, 8))

    sup.spectrum(test_spec, grid=False, ax=ax)

    # Process the MS/MS spectrum.
    # fragment_tol_mass = 50
    # fragment_tol_mode = 'ppm'
    # spectrum = (spectrum.set_mz_range(min_mz=50, max_mz=2005)
    #             .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
    #             .filter_intensity(min_intensity=0.0005, max_num_peaks=8000)
    #             .scale_intensity()
    #             .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
    #                                     ion_types='by'))
    #
    # # Plot the MS/MS spectrum.
    # fig, ax = plt.subplots(figsize=(12, 6))
    # sup.spectrum(spectrum,grid=False,ax=ax)
    plt.title(input_file+' '+peptide_with_mod+' spec: '+str(spectra_no)+' ret time: '+str(ret_time)+ ' charge: '+str(charge1)+
              ' precursor mass: '+str(pre_cursor_mass))
    plt.savefig(saved_file_path+peptide+'_'+input_file+'_'+str(spectra_no)+'.png')
    #plt.show()
    plt.close()


def ms2_visualize_target(ms2_dict_of_dict, pep_file_spec_no_dict_of_dict,target_psm_list, save_folder):
    """
    plot ms2 spectrum
    :param ms2_dict_of_dict: from ms2_reader.py {file:{spec_no:()}}
    :param pep_file_spec_no_dict_of_dict: {peptidecharge:{file_name:[spec_no1, spec_no2]},peptide:{}}
    :param save_folder: folder to output pngs
    :param target_psm_list: target PSMs list to plot. [psm1charge,psm2charge...]
    :return:
    """
    for pepcharge in target_psm_list:
        print (pepcharge)
        pep_seq, charge = ''.join(pepcharge[:-1]), pepcharge[-1]
        for file in pep_file_spec_no_dict_of_dict[pepcharge]:
            for spec_no in pep_file_spec_no_dict_of_dict[pepcharge][file]:
                if not os.path.exists(save_folder + pep_seq +'_'+ charge + '_' + file + '_' + str(spec_no) + '.png'):
                    m_over_z, ret_time, charge1, charge2, pre_cursor_mass, m_array, int_array = ms2_dict_of_dict[file][spec_no]
                    m_array = np.array(m_array, dtype=float)
                    int_array = np.array(int_array, dtype=float)
                    # add one value before and after m_array and int_array to make the range same

                    m_array = np.append(m_array, 2000)
                    m_array = np.insert(m_array, 0, 1)
                    int_array = np.append(int_array, np.max(int_array) * 0.006)
                    int_array = np.insert(int_array, 0, np.max(int_array) * 0.006)
                    test_spec = sus.MsmsSpectrum('identifier', m_over_z, int(charge1), m_array, int_array)
                    pep_seq_1 = pep_seq.replace('M(ox)', 'M[+15.9949]')
                    # spectrum utilis
                    test_spec.annotate_proforma(pep_seq_1, fragment_tol_mass=50, fragment_tol_mode="ppm", ion_types="aby")
                    fig, ax = plt.subplots(figsize=(12, 8))
                    sup.spectrum(test_spec, grid=False, ax=ax)
                    plt.title(file + ' ' + pep_seq + ' spec: ' + str(spec_no) + ' ret time: ' + str(
                        ret_time) + ' charge: ' + charge +
                              ' precursor mass: ' + str(pre_cursor_mass))
                    plt.savefig(save_folder + pep_seq +'_'+ charge + '_' + file + '_' + str(spec_no) + '.png')
                    plt.close()


if __name__=='__main__':
    # control_pep_list = ['SHPQFEKAARLMSAAA']
    # ms2_info_dict_of_dict = ppp.load(open('F:/alanine_tailing/2022_03_07/SHPQFEKAARLMSAAA_psms.p','rb'))
    #
    # # print ([key for key in ms2_info_dict_of_dict])
    # print (ms2_info_dict_of_dict)
    #
    # psm_path = 'F:/alanine_tailing/search/open_search/chymo_open_search/Tarpt_HS_chymo/psm.tsv'
    # peptides_info_dict = peptide_file_spectra_generator(psm_path)
    # #
    # # dta_path = 'C:/uic/lab/mankin/dta_results/dta_242_4_29/api4/'
    # # peptides_info_dict = dta_file_spectra_generator(dta_path)
    #
    #
    #
    # total_len = 0
    # except_peptides = defaultdict(list)
    # for each_pep in list(set(control_pep_list)):
    #     print (each_pep)
    #     #peptide_seq=peptides_info_dict[each_pep][0][2]
    #     file_spctra_dict = defaultdict(list)
    #     pep_spectra_dict = defaultdict(list)
    #     for each in peptides_info_dict[each_pep.split('_')[0]]:
    #         file_spctra_dict[each[0]].append(each[1])
    #         pep_spectra_dict[each[0]+str(each[1])]=each[2]
    #
    #     for each_file in file_spctra_dict:
    #         print (each_file)
    #         ms2_info_dict = ms2_info_dict_of_dict['F:/alanine_tailing/2022_03_07'+'\\'+each_file+'_clean.ms2']
    #         for each_spectra in file_spctra_dict[each_file]:
    #             print (each_spectra)
    #             peptide_seq=pep_spectra_dict[each_file+str(each_spectra)]  # peptide with mod
    #             try:
    #                 ms2_visulizer(ms2_info_dict,each_spectra,each_pep.split('_')[0],'F:/alanine_tailing/prosit/'
    #                             ,'_'.join(each_file.split('_')[-2:]), each_pep)
    #             except ValueError:
    #                 except_peptides[each_pep].append((each_file,each_spectra))
    #     total_len+=1
    #     print (total_len)
    #
    # print (f'peptides that are not processed: {except_peptides}')


    ## colon dataset plotting
    ms2_dict_of_dict = ppp.load(open('F:/Colon/ms2/ms2_dict_of_dict_target.p','rb'))
    pep_file_spec_dict_of_dict = ppp.load(open('F:/Colon/prosit/pep_file_spec_dict_of_dict.p','rb'))
    cos_sim_df = pd.read_csv('F:/Colon/prosit/cos_sim.csv').sort_values(by=['cos similarity'], ascending=False).iloc[
                 :100, :]
    psm_charge_list = [psm + str(charge) for psm, charge in zip(cos_sim_df['PSM'], cos_sim_df['charge'])]
    output_folder = 'F:/Colon/prosit/top100_real_spectrum/'
    ms2_visualize_target(ms2_dict_of_dict,pep_file_spec_dict_of_dict,psm_charge_list,output_folder)
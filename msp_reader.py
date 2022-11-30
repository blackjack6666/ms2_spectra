"""
read msp file from prosit for ms2 visualization
"""
import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from tsv_reader import peptide_charge_getter, peptide_counting, dta_charge_reader
import os

def prosit_csv_output(extended_pep_list,output_path,peptide_file_path):
    """
    -----
    output a csv file to prosit peptide prediction
    -----
    :param extended_pep_list:
    :param pep_tsv_path:
    :param peptide_file_path:could be peptide.tsv or .dta
    :return:
    """
    import pandas as pd
    if peptide_file_path.endswith('.tsv'):
        peptide_charge_dict = peptide_charge_getter(peptide_file_path)

    else:
        peptide_charge_dict = dta_charge_reader(peptide_file_path)
    peptide_charge_list = [[each, 30, charge] for each in extended_pep_list for charge in
                           peptide_charge_dict[each]]
    #print (peptide_charge_list)
    df = pd.DataFrame(peptide_charge_list,
                      columns=['modified_sequence', 'collision_energy', 'precursor_charge'])
    df = df.set_index('modified_sequence')
    df.to_csv(output_path)
    return df


def msp_reader(msp_file_path):
    """
    -----
    read peptide sequence, m/z, charge, mass and intensity array from msp file predicted by Prosit
    :param msp_file_path:
    :return: a dictionary
    """
    import numpy as np
    info_dict = {}
    with open(msp_file_path, 'r') as f_open:
        f_split = f_open.read().split('Name: ')
        for each in f_split[1:]:
            each_split = each.split('\n')
            info_line = each_split[2]
            print (info_line)
            pep_seq = info_line.split('ModString=')[1].split('//')[0]
            charge = int(each_split[0].split('/')[1])

            if 'Oxidation@M' in info_line:  # decide oxidation index
                pep_seq_moxi = ''
                oxi_M = info_line.split('Oxidation@M')[1:]
                if len(oxi_M) == 1:
                    oxi_m_idx = [int(oxi_M[0].split('/')[0])] if '; ' not in oxi_M[0] else [int(oxi_M[0].split('; ')[0].split('/')[0])]
                else:
                    oxi_m_idx = [int(chunk.split('; ')[0]) for chunk in oxi_M[:-1]]
                    if 'Carbamidomethyl@' not in oxi_M[-1]:
                        oxi_m_idx.append(int(oxi_M[-1].split('/')[0]))
                    else:
                        oxi_m_idx.append(int(oxi_M[-1].split('/')[0].split('; ')[0]))
                # add M(ox) to original sequence
                for i, j in enumerate(pep_seq):
                    if i in oxi_m_idx:
                        pep_seq_moxi += 'M(ox)'
                    else:
                        pep_seq_moxi += j
            else:
                pep_seq_moxi = pep_seq
            print (pep_seq_moxi)
            m_weight = float(each_split[1].split('MW: ')[1])
            m_over_z = m_weight
            mass_array = [float(line.split('\t')[0]) for line in each_split[4:-1]]
            int_array = [float(line.split('\t')[1]) for line in each_split[4:-1]]
            info_dict[pep_seq_moxi+str(charge)] = (pep_seq_moxi,charge,m_weight,m_over_z,mass_array,int_array)
    return info_dict


def ms2_visulizer(msp_info_dict:dict, saved_file_path:str, input_file:str, target_list=None):
    """
    plot the ms2 peak give spectrum number and peptide sequence (this function is not working anymore due to upgrade
    for spectrum utils, use ms2_visulizer2 instead)
    -----
    :param msp_info_dict: spectra_no as key, (mz, ret_time, charge1,charge2, mass_array,intensity arrary) as value
    :return: ms2 peak plot
    """
    import numpy as np
    if target_list == None:
        target_list = msp_info_dict


    for each in target_list:

        pep_seq, charge1, pre_cursor_mass, m_over_z, m_array, int_array = msp_info_dict[each]


        m_array = np.array(m_array, dtype=float)
        int_array = np.array(int_array,dtype=float)
        # add one value before and after m_array and int_array to make the range same

        m_array = np.append(m_array,2000)
        m_array = np.insert(m_array,0,1)
        int_array = np.append(int_array,np.max(int_array)*0.006)
        int_array = np.insert(int_array,0,np.max(int_array)*0.006)
        #print(m_array, int_array)
        # Create the MS/MS spectrum.
        spectrum = sus.MsmsSpectrum(
            'identifier', m_over_z, int(charge1), m_array, int_array,
             )

    # Process the MS/MS spectrum.
        fragment_tol_mass = 50
        fragment_tol_mode = 'ppm'
        try:
            spectrum = (spectrum.set_mz_range(min_mz=None, max_mz=None)
                        .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
                        .filter_intensity(min_intensity=0.0005, max_num_peaks=8000)
                        .scale_intensity()
                        .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
                                        ion_types='by'))
        except ValueError:
            print (each)
    # Plot the MS/MS spectrum.
        fig, ax = plt.subplots(figsize=(12, 6))
        sup.spectrum(spectrum,grid=False,ax=ax)
        plt.title(str(pep_seq)+' '+input_file+' '+ ' charge: '+str(charge1)+
                ' precursor mass: '+str(pre_cursor_mass*charge1))
        plt.savefig(saved_file_path+pep_seq+'_'+str(charge1)+'_'+input_file+'.png')
        #plt.show()
        plt.close()


def ms2_visulizer2(msp_info_dict:dict, saved_file_path:str, input_file:str, target_list=None):
    """
    plot the ms2 spectram from prosit predicted spectrum
    -----
    :param msp_info_dict: spectra_no as key, (mz, ret_time, charge1,charge2, mass_array,intensity arrary) as value
    :return: ms2 peak plot
    """
    import numpy as np
    if target_list == None:
        target_list = msp_info_dict

    for each in target_list:

        pep_seq, charge1, precursor_moz, m_over_z, m_array, int_array = msp_info_dict[each]
        if not os.path.exists(saved_file_path+pep_seq+'_'+str(charge1)+'_'+input_file+'.png'):
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
            pep_seq_1 = pep_seq.replace('M(ox)','M[+15.9949]')
            test_spec.annotate_proforma(pep_seq_1, fragment_tol_mass=50, fragment_tol_mode="ppm", ion_types="aby")

            fig, ax = plt.subplots(figsize=(12, 8))

            sup.spectrum(test_spec, grid=False, ax=ax)

        # Process the MS/MS spectrum.
        #     fragment_tol_mass = 50
        #     fragment_tol_mode = 'ppm'
        #     try:
        #         spectrum = (spectrum.set_mz_range(min_mz=None, max_mz=None)
        #                     .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
        #                     .filter_intensity(min_intensity=0.0005, max_num_peaks=8000)
        #                     .scale_intensity()
        #                     .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
        #                                     ion_types='by'))
        #     except ValueError:
        #         print (each)
        # Plot the MS/MS spectrum.
        #     fig, ax = plt.subplots(figsize=(12, 6))
        #     sup.spectrum(spectrum,grid=False,ax=ax)
            plt.title(str(pep_seq)+' '+input_file+' '+ ' charge: '+str(charge1)+
                    ' precursor mass: '+str(precursor_moz*charge1))
            plt.savefig(saved_file_path+pep_seq+'_'+str(charge1)+'_'+input_file+'.png')
            #plt.show()
            plt.close()


if __name__=='__main__':
    import pickle as ppp
    import pandas as pd
    # peptsv = 'F:/alanine_tailing/search/open_search/chymo_open_search/Tarpt_HS_chymo/peptide.tsv'
    #pep_list = peptide_counting(tsv_path)
    # peptide = ['TSYSEFLSQLANQYASCLKGDG']
    #dta_path = 'C:/uic/lab/mankin/dta_results/dta_242_20aa_normal_fs/api05/'
    # pep_list = ppp.load(open('C:/Users/gao lab computer/PycharmProjects/extend_different_species/PXD001364_ext_pep_list.p','rb'))
    # pep_list = [each for each in pep_list if len(each)<=30]
    # pep_list = ['SHPQFEKAARLMSAAA']
    # prosit_csv_output(pep_list,'F:/alanine_tailing/SHPQFEKAARLMSAAA.csv',peptsv)


    # info_dict = msp_reader('F:/alanine_tailing/prosit/SHPQFEKAARLMSAAA.msp')
    #
    #
    # ms2_visulizer(info_dict,'F:/alanine_tailing/prosit/','prosit_predict')

    # Colon datasets from Cornell
    info_dict = msp_reader('F:/Colon/prosit/myPrositLib.msp')
    cos_sim_df = pd.read_csv('F:/Colon/prosit/cos_sim.csv').sort_values(by=['cos similarity'],ascending=False).iloc[100:200,:]
    psm_charge_list = [psm+str(charge) for psm,charge in zip(cos_sim_df['PSM'],cos_sim_df['charge'])]
    print (psm_charge_list)
    #
    ms2_visulizer2(info_dict, 'F:/Colon/prosit/top100_200_prosit_spectrum/', 'prosit_predict',target_list=psm_charge_list)
    # from b_y_ion_gene import b_y_ion_gen
    # print (b_y_ion_gen('LM(ox)ITAM(ox)RPK3'))

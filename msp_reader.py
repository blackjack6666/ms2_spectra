"""
read msp file from prosit for ms2 visualization
"""
import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from tsv_reader import peptide_charge_getter, peptide_counting, dta_charge_reader


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
            key = each_split[0]
            pep_seq,charge = key.split('/')[0], int(key.split('/')[1])
            m_weight = float(each_split[1].split('MW: ')[1])
            m_over_z = m_weight/charge
            mass_array = [float(line.split('\t')[0]) for line in each_split[4:-1]]
            int_array = [float(line.split('\t')[1]) for line in each_split[4:-1]]
            info_dict[key] = (pep_seq,charge,m_weight,m_over_z,mass_array,int_array)
    return info_dict


def ms2_visulizer(msp_info_dict:dict, saved_file_path:str, input_file:str):
    """
    plot the ms2 peak give spectrum number and peptide sequence
    -----
    :param msp_info_dict: spectra_no as key, (mz, ret_time, charge1,charge2, mass_array,intensity arrary) as value
    :return: ms2 peak plot
    """
    import numpy as np
    for each in msp_info_dict:

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
            'identifier', m_over_z, charge1, m_array, int_array,
             peptide=pep_seq)

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
                ' precursor mass: '+str(pre_cursor_mass))
        plt.savefig(saved_file_path+pep_seq+'_'+str(charge1)+'_'+input_file+'.png')
        #plt.show()
        plt.close()


if __name__=='__main__':
    """
    tsv_path = 'C:/uic/lab/Irina/search_3_22_sorf/T/peptide.tsv'
    pep_list = peptide_counting(tsv_path)
    peptide = ['TSYSEFLSQLANQYASCLKGDG']
    dta_path = 'C:/uic/lab/mankin/dta_results/dta_242_20aa_normal_fs/api05/'
    prosit_csv_output(peptide,'C:/uic/lab/mankin/dta_results/new.csv',dta_path)
    """

    info_dict = msp_reader('C:/uic/lab/mankin/dta_results/myPrositLib.msp')


    ms2_visulizer(info_dict,'C:/uic/lab/mankin/dta_results/','api05')

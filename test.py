import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from pyteomics import mgf

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
    import numpy as np
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
    spectrum = sus.MsmsSpectrum(
        'identifier', m_over_z, charge1, m_array, int_array,
        retention_time=ret_time, peptide=peptide)

    # Process the MS/MS spectrum.
    fragment_tol_mass = 50
    fragment_tol_mode = 'ppm'
    spectrum = (spectrum.set_mz_range(min_mz=50, max_mz=2005)
                .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
                .filter_intensity(min_intensity=0.0005, max_num_peaks=8000)
                .scale_intensity()
                .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
                                        ion_types='by'))

    # Plot the MS/MS spectrum.
    fig, ax = plt.subplots(figsize=(12, 6))
    sup.spectrum(spectrum,grid=False,ax=ax)
    plt.title(input_file+' '+peptide_with_mod+' spec: '+str(spectra_no)+' ret time: '+str(ret_time)+ ' charge: '+str(charge1)+
              ' precursor mass: '+str(pre_cursor_mass))
    plt.savefig(saved_file_path+peptide+'_'+input_file+'_'+str(spectra_no)+'.png')
    #plt.show()
    plt.close()



import pickle as ppp
from psm_reader import peptide_file_spectra_generator, dta_file_spectra_generator
from collections import defaultdict
from glob import glob


control_pep_list = ppp.load(open('C:/uic/lab/Irina/sorf_dta_updated_db/sorf_pep_list_4_23.p','rb'))
ms2_info_dict_of_dict = ppp.load(open('C:/uic/lab/Irina/2020-03-12/sorf_file_ms2_info_dict_4_23.p','rb'))

print ([key for key in ms2_info_dict_of_dict])


#psm_path = 'D:/data/Mankin/search_result/20200129_merged_gln_tyr/ctrl/psm.tsv'
#peptides_info_dict = peptide_file_spectra_generator(psm_path)

dta_path = 'C:/uic/lab/Irina/sorf_dta_updated_db/'
peptides_info_dict = dta_file_spectra_generator(dta_path)



total_len = 0
except_peptides = defaultdict(list)
for each_pep in list(set(control_pep_list)):
    print (each_pep)
    #peptide_seq=peptides_info_dict[each_pep][0][2]
    file_spctra_dict = defaultdict(list)
    pep_spectra_dict = defaultdict(list)
    for each in peptides_info_dict[each_pep]:
        file_spctra_dict[each[0]].append(each[1])
        pep_spectra_dict[each[0]+str(each[1])]=each[2]

    for each_file in file_spctra_dict:
        print (each_file)
        ms2_info_dict = ms2_info_dict_of_dict['C:/uic/lab/Irina/2020-03-12'+'\\'+each_file+'_clean.ms2']
        for each_spectra in file_spctra_dict[each_file]:
            print (each_spectra)
            peptide_seq=pep_spectra_dict[each_file+str(each_spectra)]  # peptide with mod
            try:
                ms2_visulizer(ms2_info_dict,each_spectra,each_pep,'C:/uic/lab/Irina/sorf_dta_updated_db/sorf_ms2_spec/'
                            ,'_'.join(each_file.split('_')[-2:]), peptide_seq)
            except ValueError:
                except_peptides[each_pep].append((each_file,each_spectra))
    total_len+=1
    print (total_len)

print (except_peptides)

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

def ms2_visulizer(ms2_info_dict:dict,spectra_no:int, peptide:str,saved_file_path:str,input_file:str):
    """
    plot the ms2 peak give spectrum number and peptide sequence
    -----
    :param ms2_info_dict: spectra_no as key, (mz, ret_time, charge1,charge2, mass_array,intensity arrary) as value
    :param spectra_no: spectrum number
    :param peptide: the peptide to be visulized
    :return: ms2 peak plot
    """

    m_over_z,ret_time,charge1,charge2,m_array,int_array = ms2_info_dict[spectra_no]

    # Create the MS/MS spectrum.
    spectrum = sus.MsmsSpectrum(
        'identifier', m_over_z, charge1, m_array, int_array,
        retention_time=ret_time, peptide=peptide)

    # Process the MS/MS spectrum.
    fragment_tol_mass = 10
    fragment_tol_mode = 'ppm'
    spectrum = (spectrum.set_mz_range(min_mz=100, max_mz=1400)
                .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
                .filter_intensity(min_intensity=0.05, max_num_peaks=50)
                .scale_intensity('root')
                .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
                                        ion_types='aby'))

    # Plot the MS/MS spectrum.
    fig, ax = plt.subplots(figsize=(12, 6))
    sup.spectrum(spectrum, ax=ax)
    plt.title(peptide+' spectrum number: '+str(spectra_no)+' retention time: '+str(ret_time))
    plt.savefig(saved_file_path+peptide+'_'+input_file+'_'+str(spectra_no)+'.png')
    #plt.show()
    plt.close()


from ms2_reader import  ms2_info_reader
from glob import glob
from psm_reader import dta_info_reader
from  parameters import candidate_peptides
ms2_path = glob('D:/data/Mankin/Shura_Ribo_2020/2020-01-09/'+'*_clean.ms2')
print (ms2_path)
dta_path = glob('D:/data/Mankin/search_result/dta_result2/'+'*.dta')
print (dta_path)

'''
for ms2,dta in zip(ms2_path,dta_path):
    ms2_info_dict = ms2_info_reader(ms2)
    peptide_spectrum_dict = dta_info_reader(dta)
    for each_pep in candidate_peptides:
        for each_spectrum in peptide_spectrum_dict[each_pep]:
            ms2_visulizer(ms2_info_dict,each_spectrum,each_pep,
                          'D:/data/Mankin/ms2_spectrum/'+ms2.split('\\')[1].split('.')[0]+'/')
'''
import pickle as ppp
from psm_reader import peptide_file_spectra_generator
from collections import defaultdict
from glob import glob
peptide_list = ppp.load(open('frac_peptide_list_of_list.p', 'rb'))
control_pep_list = ppp.load(open('extend_pep_in_control.p','rb'))
psm_path = 'D:/data/Mankin/search_result/20200129_merged_gln_tyr/ctrl/psm.tsv'
peptides_info_dict = peptide_file_spectra_generator(psm_path)



ms2_info_dict_of_dict = ppp.load(open('frac_ms2_dict.p','rb'))


total_len = 0
except_peptides = defaultdict(list)
for each_pep in list(set(control_pep_list)):
    print (each_pep)
    
    file_spctra_dict = defaultdict(list)
    for each in peptides_info_dict[each_pep]:
        file_spctra_dict[each[0]].append(each[1])
    for each_file in file_spctra_dict:
        print (each_file)
        ms2_info_dict = ms2_info_dict_of_dict['D:/data/Mankin/Shura_Ribo_2020/2020_01_24_ms2'+'\\'+each_file+'_clean.ms2']
        for each_spectra in file_spctra_dict[each_file]:
            try:
                ms2_visulizer(ms2_info_dict,each_spectra,each_pep,'D:/data/Mankin/Shura_Ribo_2020/2020_01_24_ms2/ctrl/'
                              ,'_'.join(each_file.split('_')[-2:]))
            except ValueError:
                except_peptides[each_pep].append((each_file,each_spectra))
    total_len+=1
    print (total_len)

print (except_peptides)
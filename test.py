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

def ms2_visulizer(ms2_info_dict:dict,spectra_no:int, peptide:str,saved_file_path:str):
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
    plt.savefig(saved_file_path+peptide+str(spectra_no)+'.png')
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

for ms2,dta in zip(ms2_path,dta_path):
    ms2_info_dict = ms2_info_reader(ms2)
    peptide_spectrum_dict = dta_info_reader(dta)
    for each_pep in candidate_peptides:
        for each_spectrum in peptide_spectrum_dict[each_pep]:
            ms2_visulizer(ms2_info_dict,each_spectrum,each_pep,
                          'D:/data/Mankin/ms2_spectrum/'+ms2.split('\\')[1].split('.')[0]+'/')


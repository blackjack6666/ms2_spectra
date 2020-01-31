def peptide_spectra_dict(psm_tsv:str):
    """
    get peptide-spectra numbers information
    :param psm_tsv:
    :return: a dictionary with peptide as key, matched spectra number as value
    """
    from collections import defaultdict
    info_dict = defaultdict(list)
    with open(psm_tsv,'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            pep_seq = line_split[1]
            spectra_number = line_split[0].split('.')[-2]
            info_dict[pep_seq].append(spectra_number)
    return info_dict

def spectra_info_generator(psm_tsv):
    """
    read psm line by line, file_spectra as key, the rest as value
    :param psm_path:
    :return:
    """
    info_dict = {}
    with open(psm_tsv,'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            file_spectra_number = line_split[0]
            info = '\t'.join(line_split[1:])
            info_dict[file_spectra_number] = info
    return info_dict

def peptide_spectra_dict2(psm_tsv:str):
    """

    :param psm_tsv:
    :return: a dictionary with peptide as key, matched file_spectra number as value
    """
    from collections import defaultdict
    info_dict = defaultdict(list)
    with open(psm_tsv, 'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            pep_seq = line_split[1]
            spectra_number = line_split[0]  # file_name_spectra_number
            info_dict[pep_seq].append(spectra_number)
    return info_dict
def peptide_file_spectra_generator(psm_tsv:str):
    """
    get peptide-file-spectra-infomation, for multiple file search together
    -----
    :param psm_tsv:
    :return: a dictionary of dictionary with peptide as key, file and matched spectra number as value
    """
    from collections import defaultdict
    info_dict = defaultdict(list)

    with open(psm_tsv, 'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            pep_seq = line_split[1]
            file_name = line_split[0].split('.')[0]
            spectra_number = line_split[0].split('.')[-2]
            info_dict[pep_seq].append((file_name,spectra_number))

    return info_dict
def dta_info_reader(dta_file):

    info_list = [] # list of tuple
    from collections import defaultdict

    with open(dta_file, 'r') as file_open:
        for i in range(29):
            next(file_open)
        Reverse_start = 0
        for line in file_open:
            line_split = line.split('\t')
            if line.startswith('Reverse_'):
                Reverse_start = 1
            elif line.startswith('sp'):
                Reverse_start = 0
                uniprot_id = line_split[0].split('|')[1]
                seq_coverage = line_split[3]
            elif len(line_split) == 15 and Reverse_start == 0:
                pep_seq = line_split[-1].split('.')[1]
                length = len(pep_seq)
                unique = line_split[0] == '*'
                file_spect_no = line_split[1]
                spect_no = int(line_split[1].split('.')[-2])
                xcorr = line_split[2]
                DeltCN = line_split[3]
                mass_obs = line_split[5]
                mass_calc = line_split[6]
                intensity = line_split[8]
                SpR = line_split[9]
                prob = line_split[10]
                charge = line_split[1].split('.')[-1]
                info_list.append((unique, pep_seq, length, charge, file_spect_no, spect_no,
                                    xcorr, DeltCN, seq_coverage, mass_obs, mass_calc,
                                    SpR, prob, uniprot_id, intensity))
    info_dict = defaultdict(list)
    for each in info_list:
        info_dict[each[1]].append(each[5])
    return info_dict

def psm_compare(psm_path1,psm_path2,target_peptide):
    """
    given target peptide, check if the spectrum in psm path1 also is identified into something else
    -----
    :param psm_path1: psm path in which the target peptide has spectrum
    :param psm_path2: psm path to be compared
    :param target_peptide:
    :return:
    """
    spectrum_list = peptide_spectra_dict2(psm_path1)[target_peptide]
    print ('%s matched spectrum: %s' %(target_peptide,spectrum_list))
    info_dict2 = spectra_info_generator(psm_path2)
    same_spectrum = False
    same_spectrum_count = 0
    for each_spectrum in spectrum_list:
        if each_spectrum in info_dict2:
            print (each_spectrum,info_dict2[each_spectrum])
            same_spectrum = True
            same_spectrum_count+=1
    if not same_spectrum:
        print('Total spectrum = %i, No same spectrum found' % len(spectrum_list))
    else:
        print('Total spectrum = %i, found %i same spectrum' % (len(spectrum_list),same_spectrum_count))
    return None

psm_path1 = 'D:/data/Mankin/search_result/20200129_merged_gln_tyr/api1/psm.tsv'
psm_path2 = 'D:/data/Mankin/search_result/20200129_noextension_databs/api1/psm.tsv'
target_peptide = 'KPTLPVAAWEIR'
psm_compare(psm_path1,psm_path2,target_peptide)
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
            spectra_number = int(line_split[0].split('.')[-2])
            info_dict[pep_seq].append(spectra_number)
    return info_dict

def peptide_file_spectra_dict(psm_tsv:str, fragpipe=15.0):
    """
    get peptide-file-spectra dictionary
    :param psm_tsv:
    :return:
    """
    from collections import defaultdict
    info_dict = defaultdict(list)
    with open(psm_tsv, 'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            pep_seq = line_split[2] if fragpipe == 15.0 else line_split[1]
            spectra_number = int(line_split[0].split('.')[-2])
            file_name = line_split[0].split('.')[0]
            info_dict[pep_seq].append((file_name,spectra_number))
    return info_dict


def peptide_file_spec_dict_of_dict(psm_tsv:str):
    """
    {pep1:{file1:[spec1,spec2], file2:[spec1,spec2]}, pep2:{}}
    :param psm_tsv:
    :return:
    """
    from collections import defaultdict
    info_dict = defaultdict(list)

    with open(psm_tsv, 'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            pep_seq = line_split[1]
            charge = line_split[0].split('.')[-1]
            spectra_number = int(line_split[0].split('.')[-2])
            file_name = line_split[0].split('.')[0]
            info_dict[pep_seq+charge].append((file_name,spectra_number))

    new_info_dict = {}  # dictionary structure change
    for each in info_dict:
        file_spec_list_dict = defaultdict(list)
        for each_tuple in info_dict[each]:
            file_spec_list_dict[each_tuple[0]].append(each_tuple[1])
        new_info_dict[each] = file_spec_list_dict

    return new_info_dict


def  spectra_info_generator(psm_tsv):
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


def peptide_file_spectra_generator(psm_tsv:str,fragpipe=15.0):
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
            pep_seq = line_split[2] if fragpipe == 15.0 else line_split[1]
            file_name = line_split[0].split('.')[0]
            spectra_number = int(line_split[0].split('.')[-2])
            pep_seq_mod = pep_seq if line_split[3] == '' else line_split[3]
            info_dict[pep_seq].append((file_name,spectra_number,pep_seq_mod))

    return info_dict


def dta_file_spectra_generator(file_location_path:str):
    """
    -----
    get peptide-file-spectra-infomation from dta files
    -----
    :param dta_file_location:
    :return:
    """
    from glob import glob
    from collections import defaultdict

    dta_files = []

    # single dta read
    if isinstance(file_location_path, str) and file_location_path.endswith('.dta'):
        dta_files = [file_location_path]
        print('reading single dta file')

    # dta folder read
    elif isinstance(file_location_path, str) and not file_location_path.endswith('.dta'):
        dta_files = glob(file_location_path + '*.dta')
        print('reading single folder')

    # multiple dta files read or multiple dta folder read
    elif isinstance(file_location_path, list):
        for each_dta in file_location_path:
            # when parameter is dta file

            if each_dta.endswith('.dta'):
                dta_files.append(each_dta)

            # when parameter is a folder path
            else:
                dta_files += glob(each_dta + '*.dta')
    else:
        raise ValueError('parameter should be string folder path or dta file path or list of dta files or folder paths')

    # exclude wash and hela files
    clean_dta_files = []
    for each in dta_files:
        wash_hela = 0
        for word in ['wash', 'Wash', 'WASH', 'Hela', 'hela', 'HELA']:
            if word in each:
                wash_hela += 1
                break
        if wash_hela == 0:
            clean_dta_files.append(each)

    print(clean_dta_files)

    # read info.

    info_dict = defaultdict(list)
    for dta_file in clean_dta_files:
        with open(dta_file, 'r') as file_open:
            for i in range(29):
                next(file_open)
            Reverse_start = 0
            for line in file_open:
                line_split = line.split('\t')
                if line.startswith('Reverse_') or line.startswith('Rev_'):
                    Reverse_start = 1
                elif line.startswith('sp') or line.startswith('tr'):
                    Reverse_start = 0

                elif len(line_split) == 15 and Reverse_start == 0:
                    pep_seq = line_split[-1].split('.')[1]
                    file_name = line_split[1].split('.')[0]
                    spectra_number = int(line_split[1].split('.')[-2])
                    pep_seq_mod = pep_seq
                    info_dict[pep_seq].append((file_name,spectra_number,pep_seq_mod))

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
    return same_spectrum_count


if __name__== "__main__":
    psm_path = 'D:/data/ext_evo_pj/gb_ref_search_7_7_PXD001723/psm.tsv'
    print (peptide_file_spec_dict_of_dict(psm_path))
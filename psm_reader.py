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

print(dta_info_reader('D:/data/Mankin/search_result/dta_result/XS_Shura_Ribo_Api_4.dta')['KPTLPVAAWEIR'])
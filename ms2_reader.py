def ms2_info_reader(ms2_path:str):
    """
    read a ms2 file and get it's info.
    :param ms2_path: a ms2 path string
    :return: a info. dictionary
    """

    info_dict = {}
    with open(ms2_path,'r') as f:
        f_read = f.read().rstrip('\r\n')
        # split file
        f_split = f_read.split('\nS')[1:]
        for each in f_split:
            each_split = each.split('\n')
            spectra_no = int(each_split[0].split('\t')[1])
            print (spectra_no)
            m_over_z = float(each_split[0].split('\t')[-1])
            ret_time = float(each_split[1].split('\t')[-1])
            charge1 = int(each_split[9].split('\t')[1])
            pre_cursor_mass = float(each_split[9].split('\t')[2])
            charge2 = int(each_split[10].split('\t')[1]) if each_split[10].startswith('Z') else None
            m_array = [float(line.split(' ')[0]) for line in each_split[10:] if line[0].isdigit()]
            int_array = [float(line.split(' ')[1]) for line in each_split[10:] if line[0].isdigit()]
            info_dict[spectra_no] = (m_over_z,ret_time,charge1,charge2,pre_cursor_mass,m_array,int_array)
    return info_dict

def target_ms2_info_reader(ms2_path:str, spec_list:list):
    """
    -----
    only read ms2 info of the ones in para spec_list
    -----
    :param ms2_path: a ms2 file path string
    :param spec_list: target peptide list
    :return: a info dictionary, only includes the spec in para spec_list
    """
    info_dict = {}
    with open(ms2_path,'r') as f:
        f_read = f.read().rstrip('\r\n')
        # split file
        f_split = f_read.split('\nS')[1:]
        for each in f_split:
            each_split = each.split('\n')
            spectra_no = int(each_split[0].split('\t')[1])
            if spectra_no in spec_list:
                spec = True
            else:
                spec = False
            if spec:

                print (spectra_no)
                m_over_z = float(each_split[0].split('\t')[-1])
                ret_time = float(each_split[1].split('\t')[-1])
                charge1 = int(each_split[9].split('\t')[1])
                pre_cursor_mass = float(each_split[9].split('\t')[2])
                charge2 = int(each_split[10].split('\t')[1]) if each_split[10].startswith('Z') else None
                m_array = [float(line.split(' ')[0]) for line in each_split[10:] if line[0].isdigit()]
                int_array = [float(line.split(' ')[1]) for line in each_split[10:] if line[0].isdigit()]
                info_dict[spectra_no] = (m_over_z,ret_time,charge1,charge2,pre_cursor_mass,m_array,int_array)
                spec_list.remove(spectra_no)
    return info_dict

def target_peptide_file_spec_getter(target_pep_list:list,psm_tsv_path:str):
    """
    -----
    get the file and spectra number info only for target peptide, make ms2 info reader faster and save space
    -----
    :param target_pep_list: a list of target peptides
    :param psm_tsv_path: psm_tsv result path
    :return: a dictionary with file name as key, a list of target spec numbers as value
    """
    from psm_reader import peptide_file_spectra_dict
    from collections import defaultdict
    pep_info_dict = peptide_file_spectra_dict(psm_tsv_path)
    target_pep_info_dict = {each:pep_info_dict[each] for each in target_pep_list}
    file_spec_no_list_dict = defaultdict(list)
    for each in target_pep_info_dict:
        for each_tuple in target_pep_info_dict[each]:
            file_spec_no_list_dict[each_tuple[0]].append(each_tuple[1])

    return file_spec_no_list_dict

if __name__=='__main__':
    from glob import glob
    import pickle as ppp

    psm_tsv = 'F:/alanine_tailing/search/open_search/chymo_open_search/Tarpt_HS_chymo/psm.tsv'


    # pep_list = ppp.load(open('PXD001723_ext_pep_list.p', 'rb'))
    pep_list = ['SHPQFEKAARLMSAAA']
    # print (pep_list)
    # ms2_path = 'F:/XS/c_elegans/PXD001723/'
    #
    # ms2_list = glob(ms2_path + '*_clean.ms2')
    ms2_list = glob('F:/alanine_tailing/2022_03_07/*Chymo_clean.ms2')
    file_spec_no_list_dict = target_peptide_file_spec_getter(pep_list, psm_tsv)
    print(file_spec_no_list_dict)

    ms2_dict_of_dict = {
        each: target_ms2_info_reader(each, file_spec_no_list_dict[each.split('\\')[-1].split('_clean')[0]]) for each in
        ms2_list}

    ppp.dump(ms2_dict_of_dict,
             open('F:/alanine_tailing/2022_03_07/SHPQFEKAARLMSAAA_psms.p', 'wb'))

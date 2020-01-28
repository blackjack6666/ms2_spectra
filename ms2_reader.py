def ms2_info_reader(ms2_path:str):
    """
    read a ms2 file and get it's info.
    :param ms2_path: a ms2 path string
    :return: a info. dictionary
    """

    info_dict = {}
    with open(ms2_path,'r') as f:
        f_read = f.read()
        # split file
        f_split = f_read.split('\nS')[1:]
        for each in f_split:
            each_split = each.split('\n')
            spectra_no = int(each_split[0].split('\t')[1])
            print (spectra_no)
            m_over_z = float(each_split[0].split('\t')[-1])
            ret_time = float(each_split[1].split('\t')[-1])
            charge1 = int(each_split[9].split('\t')[1])
            charge2 = int(each_split[10].split('\t')[1]) if each_split[10].startswith('Z') else None
            m_array = [float(line.split(' ')[0]) for line in each_split[10:] if line[0].isdigit()]
            int_array = [float(line.split(' ')[1]) for line in each_split[10:] if line[0].isdigit()]
            info_dict[spectra_no] = (m_over_z,ret_time,charge1,charge2,m_array,int_array)
    return info_dict


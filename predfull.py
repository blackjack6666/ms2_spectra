from tsv_reader import peptide_charge_getter, peptide_counting, dta_charge_reader
from msp_reader import ms2_visulizer

def prosit_csv_output(pep_list,output_path,peptide_file_path):
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
    peptide_charge_list = [[each, charge, 'HCD', 30] for each in pep_list for charge in
                           peptide_charge_dict[each]]
    #print (peptide_charge_list)
    df = pd.DataFrame(peptide_charge_list,
                      columns=['Peptide', 'Charge', 'Type', 'NCE'])
    df = df.set_index('Peptide')
    df.to_csv(output_path,sep = '\t')
    return df


def mgf_file_reader(mgf_file):
    """
    read mgf file that contains peptide info and predicted m/z and int
    :param mgf_file: returned by predfull.py
    :return:
    """
    info_dict = {}
    with open(mgf_file, 'r') as f_open:
        file_split = f_open.read().split('END IONS\n\n')

        for each in file_split[:-1]:
            each_split = each.split('\n')
            pep_seq, charge, m_over_z = each_split[2][8:], int(each_split[3][7]), float(each_split[4][8:])

            precusor_mass = m_over_z*charge
            mass_array = [float(line.split(' ')[0]) for line in each_split[5:-1]]
            int_array = [float(line.split(' ')[1]) for line in each_split[5:-1]]
            info_dict[pep_seq+str(charge)] = pep_seq,charge,precusor_mass,m_over_z,mass_array,int_array
    return info_dict

# peptsv = 'D:/uic/lab/mankin/20200302_3_2_db_search/api05/peptide.tsv'
#
# peptide_list = ['AEHLVFWNGGR','VPVTDESPATR','WKNPTPSYSK']
# prosit_csv_output(peptide_list,'D:/uic/lab/mankin/predfull/test.tsv', peptsv)

mgf_file = 'D:/uic/lab/code/PredFull-master/test_prediction.mgf'
info_dict = mgf_file_reader(mgf_file)
ms2_visulizer(info_dict,'D:/uic/lab/mankin/predfull/', 'predfull_predict')
import re
from collections import defaultdict
import pandas as pd

def protein_tsv_reader(protein_tsv_file):
    protein_list = []
    with open(protein_tsv_file, 'r') as file_open:
        next(file_open)
        for line in file_open:
            line_split = line.split("\t")
            protein_ID = line_split[3]
            protein_list.append(protein_ID)
    return protein_list

def protein_phospho_counting(protein_csv_file):
    phospho_protein_dict = {}
    with open(protein_csv_file, 'r') as file_open:
        next(file_open)
        pattern = re.compile('\d+\w{1}\(79\.9663\)')  # look for phosphorylation
        for line in file_open:
            protein_ID = line.split('\t')[3]
            regex = re.findall(pattern, line)
            for ele in regex:
                if ele != '':
                    phospho_protein_dict[protein_ID] = regex
    return phospho_protein_dict


def info_getter(protein_tsv_file):
    info_dict = {}
    with open(protein_tsv_file, 'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            protein_ID = line_split[3]
            length = int(line_split[6])
            coverage = float(line_split[7])
            total_spec_count = int(line_split[-9])
            total_intensity = line_split[-5]
            info_dict[protein_ID] = (length,coverage,total_spec_count,total_intensity)
    return info_dict


def info_getter1(protein_tsv_file):
    df = pd.read_csv(protein_tsv_file, delimiter='\t', header=0)
    return df

def peptide_counting(peptide_tsv_file):
    peptide_list = []
    with open(peptide_tsv_file, 'r') as file_open:
        next(file_open)
        for line in file_open:
            peptide_seq = line.split("\t")[0]
            peptide_list.append(peptide_seq)
    return peptide_list


def peptide_phospho_reader(peptide_tsv_file, mod=79.9663): # 79.9663 is the delta mass of phosphorylation on STY
    pep_phos_dict = defaultdict()
    with open(peptide_tsv_file) as file_open:
        for i in range(1):
            next(file_open)
        for line in file_open:
            pep_seq = line.split('\t')[0]

            pattern = re.compile('\d+\w{1}\('+str(mod)+'\)')
            regex = re.findall(pattern, line)
            for ele in regex:
                if ele != '':
                    pep_phos_dict[pep_seq]=regex
    return pep_phos_dict

def venn_diagram_gen(dictionary, title=''): # parameter could be a dictionary of proteins or peptides from different samples {'sample1': [], 'sample2': []}
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2, venn3

    value_list_of_sets = [set(l) for l in dictionary.values()]
    sample_name_list = [n for n in dictionary.keys()]
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    if len(dictionary) == 2:  # two samples for venn diagram
        venn2(value_list_of_sets, set_labels=sample_name_list)

    elif len(dictionary) == 3:  # 3 samples for venn diagram
        venn3(value_list_of_sets, set_labels=sample_name_list)

    else:
        print ('Error: only 2 or 3 comparison for venn diagram are accepted in this script.')

    plt.title(title)
    plt.show()

def psm_reader(psm_path):
    pep_spec_count_dict = defaultdict(int)
    ret_pep_dict = {}
    with open(psm_path, 'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            pep_seq = line_split[1]
            retention_time = float(line_split[4])/60  # in minute
            pep_spec_count_dict[pep_seq]+=1
            ret_pep_dict[retention_time] = pep_seq
    return pep_spec_count_dict, ret_pep_dict


def peptide_charge_getter(peptide_tsv_file:str):
    """
    -----
    get charge for each peptide sequence
    -----
    :param tsv_path:
    :return:
    """
    peptide_charge_dict = defaultdict(list)
    with open(peptide_tsv_file, 'r') as file_open:
        next(file_open)
        for line in file_open:
            line_split = line.split("\t")
            peptide_seq = line_split[0]
            charge = line_split[2]
            if ',' in charge:
                peptide_charge_dict[peptide_seq].append(int(charge.split(',')[0]))
                peptide_charge_dict[peptide_seq].append(int(charge.split(',')[1]))
            else:
                peptide_charge_dict[peptide_seq].append(int(charge))
    return peptide_charge_dict


def dta_charge_reader(file_location_path):
    """
    -----
    read peptide seq, PSM, protein ID from dta files
    -----
    :param file_location_path: could be str folder path or dta file path or list of dta files or folder path
    :return: peptide list, PSM dictionary showing number of PSM, protein ID list
    """
    from glob import glob
    from collections import defaultdict

    dta_files = []

    # single dta read
    if isinstance(file_location_path,str) and file_location_path.endswith('.dta'):
        dta_files = [file_location_path]
        print ('reading single dta file')

    # dta folder read
    elif isinstance(file_location_path,str) and not file_location_path.endswith('.dta'):
        dta_files = glob(file_location_path + '*.dta')
        print ('reading single folder')

    # multiple dta files read or multiple dta folder read
    elif isinstance(file_location_path,list):
        for each_dta in file_location_path:
            # when parameter is dta file

            if each_dta.endswith('.dta'):
               dta_files.append(each_dta)

            # when parameter is a folder path
            else:
                dta_files += glob(each_dta+'*.dta')
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

    print (clean_dta_files)

    # read info.

    seq_charge_dict = defaultdict(list)
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
                    charge = int(line_split[1].split('.')[-1])
                    seq_charge_dict[pep_seq].append(charge)

    seq_charge_dict = {each:list(set(seq_charge_dict[each])) for each in seq_charge_dict}
    return seq_charge_dict
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn3, venn3_circles
    import pandas as pd
    '''
    SC_1_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_1/protein.tsv'
    SC_1_protein_set = (protein_tsv_reader(SC_1_tsv_path))

    SC_2_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_2_3/SC_2/protein.tsv'
    SC_2_protein_set = (protein_tsv_reader(SC_2_tsv_path))
    SC_3_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_2_3/SC_3/protein.tsv'
    SC_3_protein_set = (protein_tsv_reader(SC_3_tsv_path))
    SC_combined_set = set(SC_1_protein_set+SC_2_protein_set+SC_3_protein_set)
    print (len(SC_combined_set))
    '''

    SC_1_peptide_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_phos_2nd_search/SC_phos_1/peptide.tsv'
    SC_2_peptide_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_phos_2nd_search/SC_phos_2/peptide.tsv'
    SC_3_peptide_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_phos_2nd_search/SC_phos_3/peptide.tsv'

    SC_1_peptide_list = pd.read_csv(SC_1_peptide_tsv_path, delimiter="\t")['Peptide'].values
    SC_2_peptide_list = pd.read_csv(SC_2_peptide_tsv_path, delimiter="\t")['Peptide'].values
    SC_3_peptide_list = pd.read_csv(SC_3_peptide_tsv_path, delimiter="\t")['Peptide'].values
    #print (SC_1_peptide_list)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    venn3([set(SC_1_peptide_list), set(SC_2_peptide_list), set(SC_3_peptide_list)], set_labels=('SC_1_peptide', 'SC_2_peptide', 'SC_3_peptide'))
    plt.title("Venn diagram for peptides identified in SC after phospho-peptide enrichment")
    plt.show()

    '''
    fig, ax = plt.subplots(1,1,figsize=(10,10))
    venn3([SC_1_protein_set, SC_2_protein_set, SC_3_protein_set],set_labels=('SC_1', 'SC_2', 'SC_3'))
    c = venn3_circles([SC_1_protein_set, SC_2_protein_set, SC_3_protein_set],linewidth=1)
    c[0].set_color('magenta')
    c[0].set_edgecolor('none')
    c[0].set_alpha(0.4)
    plt.title('Venn diagram for proteins identified in spinal cord without phosphopeptide enrichment')
    plt.show()
    '''
    '''
    SC_phos_1_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_phos_2nd_search/SC_phos_1/protein.tsv'
    SC_phos_2_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_phos_2nd_search/SC_phos_2/protein.tsv'
    SC_phos_3_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_phos_2nd_search/SC_phos_3/protein.tsv'

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    venn3([set(protein_tsv_reader(SC_phos_1_tsv_path)), set(protein_tsv_reader(SC_phos_2_tsv_path)),
           set(protein_tsv_reader(SC_phos_3_tsv_path))], set_labels=('SC_1', "SC_2", "SC_3"))
    plt.title('Venn diagram for proteins identified in spinal cord after phospho-peptide enrichment')
    plt.show()
    '''

    '''
    SC_1_peptide_tsv_path = 'D:/data/phospho_wang/9_17_2019_search_result/SC_1/peptide.tsv'
    pep_phos_dict = peptide_phospho_reader(SC_1_peptide_tsv_path)
    print (len(pep_phos_dict))
    '''
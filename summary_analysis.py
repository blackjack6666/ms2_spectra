import pandas as pd
import numpy as np
from collections import defaultdict
def fasta_reverse_generator(ID_descrip_dict, protein_dict, fasta_file_out):
    print ('reverse_algorithm = null, or other protease')
    # read protein sequence into dic


    # read description into dic


    # write id and reverse sequence into fasta_file_out
    with open(fasta_file_out, 'w') as file_open:
        for id in protein_dict:

            rev_seq = protein_dict[id][::-1]


            block = np.arange(0, len(rev_seq), 60)
            block = np.append(block, len(rev_seq))

            # write forward
            file_open.write('>sp'+ '|'+id+'|'+ID_descrip_dict[id]+'\n')
            for i in range(len(block)-1):
                file_open.write(protein_dict[id][block[i]:block[i+1]]+'\n')

            # write reverse
            file_open.write('>Rev_sp'+ '|'+id+'|'+ID_descrip_dict[id]+'\n')
            for i in range(len(block)-1):
                file_open.write(rev_seq[block[i]:block[i+1]]+'\n')
    return fasta_file_out

"""
from protein_coverage import read_fasta_into_dict, read_description_into_dict
fasta_path = 'D:/data/Mankin/fasta/gln_tyr_reverse_2_13.fasta'

protein_dict = read_fasta_into_dict(fasta_path)[0]
for each in protein_dict:
    if not protein_dict[each].startswith('M'):
        #new_seq = 'M'+''.join(list(protein_dict[each])[1:])
        #protein_dict[each] = new_seq
        print(each)
#id_description_dict = read_description_into_dict(fasta_path)
#fasta_reverse_generator(id_description_dict,protein_dict,'D:/data/Mankin/fasta/gln_tyr_reverse_2_13.fasta')
"""
from tsv_reader import venn_diagram_gen
import matplotlib.pyplot as plt
path = 'D:/data/Mankin/3_9_images/apis_control_extended_psm_summary_3_9.xlsx'

df = pd.read_excel(path)
df_api05 = df[df['file_name'].str.contains("api05")]
df_api1 = df[df['file_name'].str.contains("api1")]
df_api4 = df[df['file_name'].str.contains("api4")]
df_ctrl = df[df['file_name'].str.contains("ctrl")]
df_apis = df[df['file_name'].str.contains("api")]




"""
uniport_api05,uniport_api1,uniport_api4,uniport_control = \
    df_api05['uniprot_id'].tolist(), df_api1['uniprot_id'].tolist(),df_api4['uniprot_id'].tolist(),df_ctrl['uniprot_id'].tolist()


uniprot_apis = uniport_api05+uniport_api1+uniport_api4

venn_diagram_gen({'Api05':uniport_api05, 'Api1':uniport_api1, 'Api4':uniport_api4},
                 title='Proteins identified by extended peptides in all Apis')
"""


df_apis_normal_count, df_apis_forward_count,df_apis_backward_count = \
    df_api4[df_api4['protein id'].str.contains('normal')].shape[0], \
    df_api4[df_api4['protein id'].str.contains('forward')].shape[0], \
    df_api4[df_api4['protein id'].str.contains('backward')].shape[0]
labels = 'in frame', '+1 frameshift', '-1 frameshift'
colors = ['gold', 'yellowgreen', 'lightcoral']
explode = (0.1,0,0)
patches, texts,autotexts= plt.pie([df_apis_normal_count,df_apis_forward_count,df_apis_backward_count], explode=explode,
                                  colors=colors,autopct='%1.1f%%',shadow=True, startangle=140)

plt.legend(patches, labels, loc="best")
plt.axis('equal')
plt.title('Ratio of modes of extended PSMs in Api4')
plt.tight_layout()

plt.show()

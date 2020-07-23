import pickle as ppp
import matplotlib.pyplot as plt
from tsv_reader import psm_reader, pep_probability_getter
from collections import defaultdict
import pandas as pd
import seaborn as sns

PXD001364_psm_path = 'D:/data/ext_evo_pj/c_elegans/gb_ext_search_7_11_PXD001364/psm.tsv'
PXD001723_psm_path = 'D:/data/ext_evo_pj/c_elegans/gb_ext_search_7_11_PXD001723/psm.tsv'
PXD001364_pep_path = 'D:/data/ext_evo_pj/c_elegans/gb_ext_search_7_11_PXD001364/peptide.tsv'
PXD001723_pep_path = 'D:/data/ext_evo_pj/c_elegans/gb_ext_search_7_11_PXD001723/peptide.tsv'

PXD001364_prob_dict = pep_probability_getter(PXD001364_pep_path)
PXD001723_prob_dict = pep_probability_getter(PXD001723_pep_path)
PXD001364_prob_dict.update(PXD001723_prob_dict)

PXD001364_psm_dict = psm_reader(PXD001364_psm_path)[0]
PXD001723_psm_dict = psm_reader(PXD001723_psm_path)[0]
total_psm_dict = defaultdict(int)
for each in PXD001364_psm_dict:
    total_psm_dict[each]+=PXD001364_psm_dict[each]
for each in PXD001723_psm_dict:
    total_psm_dict[each]+=PXD001723_psm_dict[each]

PXD001364_ext_dict = ppp.load(open('PXD001364_ext_pep_cosScore_dict.p', 'rb'))
PXD001723_ext_dict = ppp.load(open('PXD001723_ext_pep_cosScore_dict.p', 'rb'))
PXD001364_ext_dict.update(PXD001723_ext_dict)


PXD001364_ext_pep_list = [key for key in PXD001364_ext_dict]
PXD001723_ext_pep_list = [key for key in PXD001723_ext_dict]
total_ext_pep_list = list(set(PXD001723_ext_pep_list+PXD001364_ext_pep_list))

info_list = [[each, PXD001364_prob_dict[each],PXD001364_ext_dict[each],total_psm_dict[each]] for each in total_ext_pep_list]
df = pd.DataFrame(info_list,columns=['peptide_seq', 'probability', 'cosine_score', 'PSM'])
print (df.shape)
print (df)
fig, ax = plt.subplots(figsize=(10, 5))
sns.scatterplot(x='cosine_score', y='PSM', hue='probability',data=df,ax=ax)
plt.title('ext_pep_c_elegans_cos_PSM_prob')
plt.show()


# overlap = [each for each in PXD001364_ext_pep_list if each in PXD001723_ext_pep_list]
#
# for each in overlap:
#     print (each, PXD001723_ext_dict[each], PXD001364_ext_dict[each])
#
# all_cosScore = [PXD001364_ext_dict[each] for each in PXD001364_ext_dict] + [PXD001723_ext_dict[each] for each in PXD001723_ext_dict]
#
# plt.hist(all_cosScore,100)
# #plt.xlim(0,1)
# plt.xlabel('cosScore')
# plt.ylabel('Frequency')
# plt.title('Histogram of ext pep cosScore from C_elegans')
# plt.show()

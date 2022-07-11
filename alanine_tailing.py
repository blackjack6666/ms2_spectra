from tsv_reader import fasta_reader, read_description, combined_proteintsv_map


def add_alanine(fasta_in, fasta_out, num_alanine=5):
    """
    add a defined number of A to every sequence in the fasta file
    :param fasta_in:
    :param fasta_out:
    :param num_alanine:
    :return:
    """
    protein_dict = fasta_reader(fasta_in)
    print (f'{len(protein_dict)} proteins in fasta in')
    descript_dict = read_description(fasta_in)

    with open(fasta_out,'w',newline='\n') as f_out:
        for prot in protein_dict:
            seq = protein_dict[prot]
            block = range(0, len(seq) + 60, 60)
            f_out.write('>'+ descript_dict[prot][0]+'|'+prot+'|'+descript_dict[prot][1]+'\n')
            for i in range(len(block)-1):
                f_out.write(seq[block[i]:block[i+1]]+'\n')
            for i in range(num_alanine):
                new_id = prot+'_'+str(i+1)+'A'
                new_seq = seq+'A'*(i+1)
                new_block = range(0, len(new_seq) + 60, 60)
                f_out.write('>' + descript_dict[prot][0] + '|' + new_id + '|' + descript_dict[prot][1] +' '+ str(i+1)+'Alanines'+'\n')
                for i in range(len(new_block) - 1):
                    f_out.write(new_seq[new_block[i]:new_block[i + 1]] + '\n')
    return f_out


if __name__ == '__main__':
    # fastain = 'F:/alanine_tailing/20220709/gfp_truncation.txt'
    # fastaout = 'F:/alanine_tailing/20220709/gfp_truncation_alanines.fasta'
    # add_alanine(fastain,fastaout)
    combined_prot_tsv = 'F:/alanine_tailing/20220709/ecoli_db_search/combined_protein.tsv'
    info_dict = combined_proteintsv_map(combined_prot_tsv)
    for each in info_dict:
        print (each, len(info_dict[each]))

def lorikeet_view_gen(pep_seq,
                      scan_num,
                      charge,
                      precursorMZ,
                      file_name,
                      mz_int_2d_array,
                      html_out,
                      var_mods=None,
                      ntermmod=0,
                      cterm_mod=None,
                      staticmod=None,
                      ):
    """
    html out must be put at the same directory with js scripts needed for specview
    :param pep_seq:
    :param scan_num:
    :param charge:
    :param precursorMZ:
    :param file_name:
    :param mz_int_2d_array:
    :param var_mods: ex. varMods[0] = {index: 14, modMass: 16.0, aminoAcid: 'M'};
    :param ntermmod:
    :param html_out:
    :return:
    """

    with open(html_out, 'w') as f_write:
        with open('lorikeet_template.html', 'r') as f_open:
            f_str = f_open.read()
            head = f_str.split('$("#lorikeet").specview(')[0]
            tail = '];\n</script>\n\n</body>\n\n</html>'
            info_dict = {'sequence': 'sequence',
                         'scanNum': scan_num,
                         'charge': charge,
                         'precursorMz': precursorMZ,
                         'fileName': file_name,
                         'peaks': 'peaks'}
            if var_mods:
                info_dict['variableMods'] = 'varMods'
            if ntermmod:
                info_dict['ntermMod'] = 'ntermMod'
            if cterm_mod:
                info_dict['//ctermMod'] = 'ctermMod'
            if staticmod:
                info_dict['//staticMods'] = 'staticMods'

            info_dict_str = str(info_dict).replace("'sequence'",'sequence').replace("'ntermMod'","ntermMod").replace("'ctermMod'","ctermMod").replace("'peaks'","peaks").replace("'varMods'","varMods").replace("'staticMods'","staticMods")

            middle = ');\n\n});\n\n\n'
            f_write.write(head+
                          '$("#lorikeet").specview('+
                          info_dict_str+
                          middle+
                          'var sequence = "'+
                          str(pep_seq)+
                          '";\n')
            if var_mods:
                f_write.write('var varMods = [];\n// modification index = 14; modification mass = 16.0; modified residue = "M"\n'+
                          str(var_mods)+';\n')

            f_write.write('// mass to be added to the N-terminus\n'+
                          'var ntermMod = '+str(ntermmod)+';\n\n'+
                          '// peaks in the scan: [m/z, intensity] pairs.\n'+
                          'var peaks = [')


            for each_mz_int_pair in mz_int_2d_array:
                f_write.write(str(each_mz_int_pair)+',\n')
            f_write.write(tail)
    return html_out

seq,charge,scan_num,precusor_mz,filename, = 'KFSLIFMLMLK',2,1110,685.900562804995,'test'
mz_2d_array = [[129.10223,1],[138.58896,0.00109026],[147.11281,0.21086179],[182.10498,0.000305613],[260.19687,0.20461063],[276.17065,0.16967659],[363.20267,0.04228261],[368.72324,0.002124215],
               [391.23737,0.14057337],[434.2435,0.017290255],[476.28674,0.062178917],[490.78552,0.020696316],[504.3214,0.1092823],[504.80286,7.49E-05],[556.3058,0.002918057],[589.3708,0.10352384],[635.3619,0.107186384],
               [736.4392,0.068093434],[782.4303,0.3423374],[867.47974,0.11775447],[895.5144,0.26976636],[980.5638,0.19554809],[1008.59845,0.12986062],[1095.6305,0.7183024],[1111.6042,0.20025733],[1224.6884,0.0836296],
               [1242.6989,0.4526999]]

lorikeet_view_gen(seq,scan_num,charge,precusor_mz,filename,mz_2d_array,html_out='C:/Users/gao lab computer/Downloads/UWPR-Lorikeet-d2a4a49/html/test2.html')
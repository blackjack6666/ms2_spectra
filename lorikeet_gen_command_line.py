import os
import argparse
import sys


class ProgramError(Exception):
    """An :py:class:`Exception` raised in case of a problem.
    :param msg: the message to print to the user before exiting.
    :type msg: string
    """
    def __init__(self, msg):
        """Construction of the :py:class:`ProgramError` class.
        :param msg: the message to print to the user
        :type msg: string
        """
        self.message = str(msg)

    def __str__(self):
        """Creates a string representation of the message."""
        return self.message


def ms2_csv_reader(mz_int_csv_path):
    """
    -----
    Read a csv file that includes m/z and int pairs for one spectra, mz and int separated by \t,
    each pair separated by \n
    -----
    csv example: mz1    int1
                 mz2    int2
                 mz3    int3
    :param mz_int_csv_path:
    :return:
    """
    with open(mz_int_csv_path,'r') as file_open:
        try:
            return [[float(line.split("\t")[0]),float(line.split("\t")[1])] for line in file_open]
        except:
            message = "some error, check mz_int_csv file"
            raise ProgramError(message)


def lorikeet_view_gen(args):
    """
    function to plot, html out must be put at the same directory with js scripts needed for specview
    :param pep_seq:
    :param scan_num:
    :param charge:
    :param precursorMZ:
    :param file_name: original ms2 file name, just the file name, not absolute path
    :param mz_int_2d_array:
    :param var_mods: ex. varMods[0] = {index: 14, modMass: 16.0, aminoAcid: 'M'};
    :param ntermmod:
    :param html_out:
    :return:
    """

    with open(args.html_out, 'w') as f_write:
        try:
            with open('lorikeet_template.html', 'r') as f_open:
                f_str = f_open.read()
                head = f_str.split('$("#lorikeet").specview(')[0]
                tail = '];\n</script>\n\n</body>\n\n</html>'
                info_dict = {'sequence': 'sequence',
                             'scanNum': args.scan_number,
                             'charge': args.charge,
                             'precursorMz': args.precusor_mz,
                             'fileName': args.file_name,
                             'peaks': 'peaks'}
                if args.var_mod:
                    info_dict['variableMods'] = 'varMods'
                if args.ntermmod:
                    info_dict['ntermMod'] = 'ntermMod'
                if args.ctermmod:
                    info_dict['//ctermMod'] = 'ctermMod'
                if args.staticmod:
                    info_dict['//staticMods'] = 'staticMods'

                info_dict_str = str(info_dict).replace("'sequence'", 'sequence').replace("'ntermMod'",
                                                                                         "ntermMod").replace(
                    "'ctermMod'", "ctermMod").replace("'peaks'", "peaks").replace("'varMods'", "varMods").replace(
                    "'staticMods'", "staticMods")

                middle = ');\n\n});\n\n\n'
                f_write.write(head +
                              '$("#lorikeet").specview(' +
                              info_dict_str +
                              middle +
                              'var sequence = "' +
                              str(args.peptide_seq) +
                              '";\n')
                if args.var_mod:
                    f_write.write(
                        'var varMods = [];\n// modification index = 14; modification mass = 16.0; modified residue = "M"\n' +
                        str(args.var_mod) + ';\n')

                f_write.write('// mass to be added to the N-terminus\n' +
                              'var ntermMod = ' + str(args.ntermmod) + ';\n\n' +
                              '// peaks in the scan: [m/z, intensity] pairs.\n' +
                              'var peaks = [')

                if not args.mz_int_csv_path:
                    mess = "No mz_int csv file, create one"
                    raise ProgramError(mess)

                mz_int_2d_array = ms2_csv_reader(args.mz_int_csv_path)
                #print (mz_int_2d_array)
                for each_mz_int_pair in mz_int_2d_array:
                    f_write.write(str(each_mz_int_pair) + ',\n')
                f_write.write(tail)
                print ("Plotting done.")
        except FileNotFoundError:
            sys.exit(
                'FileNotFoundError: lorikeet_template.html was not put in the same directory as this script, please do so')

    return args.html_out


def parse_args():
    """
    parse the command line arguments
    ============================  =======  ====================================
         Options                   Type                   Description
    ============================  =======  ====================================
    ``--peptide_seq``             string    the peptide sequence you want to plot
    ``--scan_number``             int       the ms2 scan number for peptide seq
    ``--charge``                  int       the charge of peptide
    ``--precusor_mz``             float     m/z value for peptide
    ``--file_name``               string    ms2 file name
    ``--mz_int_csv_path``         file_path a csv path of m/z and intensity pairs,
                                            m/z and intensity separated by \t,
                                            pair separated by \n, could copy and
                                            paste from *.ms2
    ``--html_out``                file_path the path to output a html file plotting
                                            ms2 spectrum using lorikeet
    ``--var_mod``                 dict      ex. varMods[0] =
                                            {index: 14, modMass: 16.0, aminoAcid: 'M'};
    ``--ntermmod``                float     n terminal mods
    ``--ctermmod``                float     c terminal mods
    ``--staticmod``               dict      still working on it
    ============================  =======  ====================================

    :return:
    """
    # Creating the parser object
    parser = argparse.ArgumentParser(
        description="This script produces ms2 spectrum based on lorikeet html script "
    )

    # Peptide sequence
    parser.add_argument(
        "--peptide_seq", type=str, default="PEPTIDE",
        help="The aa sequence of a given peptide "
        "[Default: %(default)s]."
    )

    # scan number
    parser.add_argument(
        "--scan_number", type=int,
        help="The ms2 scan number of what you want to plot"
    )

    # charge
    parser.add_argument(
        "--charge", type=int, default=2,
        help="the charge of peptide"
    )

    # precusor mz value
    parser.add_argument(
        "--precusor_mz", type=float,
        help="the peptide m/z"
    )

    # ms2 file name
    parser.add_argument(
        "--file_name", type=str, default="ms2 file name",
        help="the ms2 file name, not essential"
    )

    # mz int csv path
    parser.add_argument(
        "--mz_int_csv_path", type=str,
        help="the csv having ms2 info (mz and int pairs) from .ms2"
    )

    # html output
    parser.add_argument(
        "--html_out", type=str, default='ms2_lorikeet.html',
        help="path to output ms2 plotting html"
    )

    # var mod
    parser.add_argument(
        "--var_mod", type=dict, default=None,
        help="variable modification"
    )

    # Ntermial modification
    parser.add_argument(
        "--ntermmod", type=float, default=0.0,
        help="n terminal modification"
    )

    # Cterminal modification
    parser.add_argument(
        "--ctermmod", type=float,default=0.0,
        help="c terminal modification"
    )

    # static modificaiton
    parser.add_argument(
        "--staticmod", type=dict, default=None,
        help="static modification"
    )

    return parser.parse_args()


def check_args(args):
    """
    check some args
    :param args:
    :return:
    """

    if args.mz_int_csv_path and not os.path.isfile(args.mz_int_csv_path):
        message = "m/z and intensity file not found"
        raise ProgramError(message)


def main():

    try:
        args = parse_args()
        print ("checking arguments...")
        check_args(args)
        print ("plotting...")
        lorikeet_view_gen(args)

        print ("If the output html is empty, the most possible reason is html_out path is not same as the js script for "
               "specview, for more infomation, please contact Xinhao Shao, Gao lab at xshao8@uic.edu")

    except KeyboardInterrupt:
        print ('Cancelled by keyboardinterruption')
        sys.exit(0)

if __name__=="__main__":
    main()
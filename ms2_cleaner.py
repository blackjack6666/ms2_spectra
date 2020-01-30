def ms_cleaner(ms2_file_name):
    with open(ms2_file_name) as file_read:
        with open(ms2_file_name.replace('.ms2', '_clean.ms2'), 'w') as file_write:
            for line in file_read:
                if line[0].isdigit() and len(line.split(' ')) == 4:
                    file_write.write(' '.join(line.split(' ')[:2]) + '\n')
                else:
                    file_write.write(line)


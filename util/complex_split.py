def split_complex_intelligently(complex_path='',binder_path='',receptor_path='',chain_info_file = ""):
    def is_type(line,c_t_lst):
        try:
            #return line.split()[4].upper() in 'HLIMJN'
            return line[21] in c_t_lst
        except:
            return False

    with open(chain_info_file,"r") as f:
        r_chn = f.readline().split()[1].split(",")
        b_chn = f.readline().split()[1].split(",")
    with open(receptor_path,'w') as rec_f:
        with open(binder_path,'w') as bnd_f:
            with open(complex_path) as f:
                for line in f.readlines():
                    if start_with_atom(line):
                        if is_type(line,r_chn):
                            rec_f.write(line[:60]+'\n')
                        elif is_type(line,b_chn):                            
                            bnd_f.write(line[:60]+'\n')


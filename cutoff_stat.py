from collections import Counter,defaultdict
l_count = 1
a_count_list = []
env_count_dict = defaultdict(list)
with open("flu.txt","r") as f:
    for l in f.readlines():
        if ( l_count + 2 ) % 4 == 0:
            print l_count
            print l
            a_count_list += eval(l)
        elif ( l_count + 1 ) % 4 == 0:
            tmp_d = eval(l)
            for t,l_ in tmp_d.items():
                env_count_dict[t] += l_
            print env_count_dict
        l_count += 1
stat = Counter(a_count_list)
p_stat = Counter(env_count_dict["p_list"])
c_stat = Counter(env_count_dict["c_list"])
a_stat = Counter(env_count_dict["a_list"])
h_stat = Counter(env_count_dict["h_list"])
print stat,p_stat ,c_stat,a_stat,h_stat

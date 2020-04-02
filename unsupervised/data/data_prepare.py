
def make_geno(hap1, hap2):
    hap = hap1 + '~' + hap2
    hap = ('~').join(sorted(hap.split('~')))

    return hap


f_out_hap = open('new_data_hap_test.txt', 'w')
f_out_geno = open('new_data_geno_test.txt', 'w')
dict_ids_genos = {}
dict_ids_haps = {}
#with open('out.pullimp.2019-01-14.txt') as f:
with open('test') as f:
    for line in f:
        splited_line = line.strip().split(',')
        id = splited_line[1]

        hap_line = splited_line[1] + ',' + splited_line[2] + ';' + splited_line[3] + ',' + splited_line[5] + ';' + \
                   splited_line[6] + ',' + splited_line[4] + ',' + splited_line[7] + '\n'
        geno = splited_line[1] + ',' + make_geno(splited_line[2], splited_line[5]) + ',' + splited_line[3] + ',' + \
               splited_line[6] + ',' + splited_line[4] + ',' + splited_line[7] + '\n'
        prob = float(splited_line[4]) * float(splited_line[7])
        if splited_line[2] != splited_line[5]:
            prob = prob * 2
        if not id in dict_ids_genos:
            dict_ids_genos[id] = {geno:prob}
            dict_ids_haps[id] = {hap_line:prob}
        else:
            dict_ids_genos[id][geno] =  prob
            dict_ids_haps[id][hap_line] =  prob


for key,val in dict_ids_haps.items():
    sorted_dict = sorted(val.items(), key=lambda kv: kv[1], reverse=True)
    f_out_hap.write(sorted_dict[0][0])

for key,val in dict_ids_genos.items():
    sorted_dict = sorted(val.items(), key=lambda kv: kv[1], reverse=True)
    f_out_geno.write(sorted_dict[0][0])

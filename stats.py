import matplotlib.pyplot as plt
import numpy as np
import json

dict_races = {}
dict_race_count_mr_sr = {}#list for each race - 0idx - count of SR, 1 idx - count of MR

count_mr = 0

with open('params_file.json') as f:
    params = json.load(f)

f_data = params.get("data_file", 'new_data_geno.txt')

with open(f_data) as in_f:
    for line in in_f:
        line = line.strip().split(',')
        race1 = line[2]
        race2 = line[3]

        races = [race1, race2 ]
        if race1 != race2:
            if race1 in dict_race_count_mr_sr:
                dict_race_count_mr_sr[race1][1] +=1
            else:
                dict_race_count_mr_sr[race1] = [0,1]

            if race2 in dict_race_count_mr_sr:
                dict_race_count_mr_sr[race2][1] += 1
            else:
                dict_race_count_mr_sr[race2] = [0, 1]

            count_mr +=1
        else:
            if race1 in dict_race_count_mr_sr:
                dict_race_count_mr_sr[race1][0] += 1
                #dict_race_count_mr_sr[race1] += 1
            else:
                dict_race_count_mr_sr[race1] = [1, 0]
                #dict_race_count_mr_sr[race1] = 1
        races = ('_').join(sorted(races))

        dict_races[races] = dict_races.get(races, 0) + 1



print(count_mr)
print(len(dict_races))

sorted_dict_races = sorted(dict_races.items(), key=lambda kv: kv[1], reverse=True)
sorted_dict = sorted(dict_race_count_mr_sr.items(), key=lambda kv: kv[1], reverse=True)
sum = 0

x, y = [], []


for key, val in sorted_dict:
    value = val[0]+val[1]
    x.append(key)
    y.append(value)

    print(key + '-  ' + str(val))
    sum += val[0]

plt.figure()
plt.bar(x, y)
plt.xticks(np.arange(len(sorted_dict)), x, rotation='vertical',fontsize=7)
plt.yticks(fontsize=7)
plt.ylabel('Number In Category')
plt.xlabel('Races')
plt.tight_layout()
plt.show()
print(sum)

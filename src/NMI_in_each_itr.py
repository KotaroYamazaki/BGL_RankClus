from sklearn.metrics.cluster import normalized_mutual_info_score
import csv

rr = []
rp = []
for v in csv.reader(open("results/cluster_result_rankclus.csv", "r")):
    rr.append([int(elm) for elm in v[1:]])
for v in csv.reader(open("results/cluster_result_proposal.csv", "r")):
    rp.append([int(elm) for elm in v[1:]])

nmi_list = []
f = open('results/nmi_in_each_iteration.csv', 'a')

for i in range(len(rr)):
    nmi_list.append(normalized_mutual_info_score(rr[i], rp[i]))
    f.write("{}".format(str(nmi_list[i])))
    if i != len(rr) - 1:
    	f.write(",")
    print(nmi_list[i])

f.write("\n")
f.close()

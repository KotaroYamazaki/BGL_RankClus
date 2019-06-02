from sklearn.metrics.cluster import normalized_mutual_info_score
import csv
result = [[ int(elm) for elm in v] for v in csv.reader(open("results/result.csv", "r"))]
correct = [[ int(elm) for elm in v] for v in csv.reader(open("results/correct.csv", "r"))]
NMI = normalized_mutual_info_score(correct[0],result[0])
f = open('results/result_time_compare.csv','a')
f.write(",{}".format(str(NMI)))
f.write("\n")
f.close
print(NMI)

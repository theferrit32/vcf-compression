import json
import matplotlib.pyplot as plt
from functools import cmp_to_key
filename = 'genotype-frequencies-chr22.json'
with open(filename) as fin:
    j = json.load(fin)

x = [0]

plots = []
legend = []
items = list(j.items())
def gt_compare(gt1, gt2):
    gt1a1 = gt1[0].split('|')[0]
    gt1a2 = gt1[0].split('|')[1]
    gt1val1 = int(gt1a1)
    gt1val2 = int(gt1a2)

    gt2a1 = gt2[0].split('|')[0]
    gt2a2 = gt2[0].split('|')[1]
    gt2val1 = int(gt2a1)
    gt2val2 = int(gt2a2)

    return (gt1val1 + gt1val2) - (gt2val1 + gt2val2)

items = sorted(items, key=cmp_to_key(gt_compare))
# stacked bar
# for i in range(len(items)):
#     gt, count = items[i]
#     #print(gt, count)
#     if i == 0:
#         bottom = None
#     else:
#         bottom = items[i-1][1]
#     p = plt.bar([0], [count], bottom=None)
#     plots.append(p)
#     legend.append(gt)

# print('legend: ' + str(legend))
# plt.legend(legend)
# plt.show()

total = sum([i[1] for i in items])

x = [i[0] for i in items]
y = [i[1] for i in items]
plt.bar(x, y)
#plt.yticks([i for i in range(0, max(y), int(max(y)/100))])
plt.xticks(x, x, rotation='vertical')
#plt.yscale('log')
plt.tight_layout()
plt.grid()
plt.xlabel('Allele Genotype Values')
plt.ylabel('Occurences in Chromosome 22')
plot_margin = 0.25

plt.margins(x=0.005, y=0.005, tight=True)
#lt.savefig(filename + '.png', bbox_inches='tight')
plt.show()

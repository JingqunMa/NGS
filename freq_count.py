import sys
from collections import defaultdict

fqFile = open(sys.argv[1])
freqs = defaultdict(int)
for line in fqFile:
    seq = line.strip()
    for i, c in enumerate(seq):
        freqs[(i,c)] += 1

fqFile.close()

maxLen = max([x[0] for x in freqs])

for pos in range(maxLen):
    data = [(x,c) for x,c in freqs.items() if x[0] == pos]
    total = sum([x[1] for x in data])
    for x,c in data:
        print('\t'.join(map(str,[pos,x[1],c,total,float(c) / total])))

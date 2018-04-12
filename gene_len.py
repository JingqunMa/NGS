import sys, os
from collections import defaultdict

gtfFile=open("hg38_plus.gtf", "r")
# gene_len_file=open("gene_len_from_hg38", "w")
# gene_len_file.write("gene_id\tgene_name\tgene_length")
print("gene_id"+"\t"+"gene_name"+"\t"+"gene_length")

for line in gtfFile:
	if line[0]!="#":
		_, _, feature, start,end, _,_, _, info = line.strip().split('\t')
		info = dict([(x.split('"')[0].strip(),x.split('"')[1].strip()) for x in info.split(';')[:-1]])
		if feature=="gene":
			gene_id=info['gene_id']
			gene_name=info['gene_name']
			gene_length=int(end)-int(start)
			# gene_len_file.write("\t".join(gene_id, gene_name, gene_length))
			print(gene_id+"\t"+gene_name+"\t"+str(gene_length))

gtfFile.close()
# gene_len_file.close()
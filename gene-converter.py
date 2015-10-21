

probemap = {}
with open('data/probe2gene.map') as pmap:
		for line in pmap:
			l = line.replace('"', '').strip().split()
			probe = l[0]
			genes = l[1:]
			probemap[probe] = genes

genematrix = {}
header = ""
with open('data/sample_probe_profile.matrix') as pmtrx:
	header = pmtrx.readline()
	for line in pmtrx:
		arr = line.replace('"', '').strip().split()
		genes = probemap[ arr[-1]]
		newrow = arr[:-1]
		for gene in genes:
			if gene in genematrix:
				(weight,gene,oldrow) = genematrix[gene]
				weight = int(weight)
				newrow = [str(float(n) + float(o)*weight) for n, o in zip(newrow, oldrow)]
				genematrix[gene] = (str(weight+1),gene,newrow)
			else:
				genematrix[gene] = (1,gene,newrow)

with open('data/sample_gene_profile.matrix', 'w') as gmtrx:
	gmtrx.write(header)
	gmtrx.write('\n')
	for (weight,gene,row) in genematrix.values():
		# print("  ".join(row) + " " + gene)
		gmtrx.write("  ".join(row) + " " + gene)
		gmtrx.write('\n')

#!/usr/env/python3

import sys
import os

hgnc_path = sys.argv[1]
in_bed_path = sys.argv[2]
basename = os.path.splitext(os.path.basename(in_bed_path))[0]

hgnc_file = open(hgnc_path, 'r') #use the protein-coding only filtered hgnc file here so there is no header

hgnc_gene_list = []
for line in hgnc_file:
    hgnc_line_list = line.strip().split('\t')
    #print(len(hgnc_line_list))
    #print(hgnc_line_list)
    hgnc_gene_name = hgnc_line_list[1]
    hgnc_gene_list.append(hgnc_gene_name)
hgnc_file.close()

in_bed_file = open(in_bed_path, 'r')
bed_header = in_bed_file.readline()
bed_gene_list = []
bed_gene_line_list = []

for line in in_bed_file:
    bed_line_list = line.strip().split('\t')
    bed_gene_name = bed_line_list[3]
    bed_gene_list.append(bed_gene_name)
    bed_gene_line_list.append(bed_line_list)
in_bed_file.close()

out_bed_file = open(basename + '.protein-coding_hgnc_filtered.txt', 'w')
ensembl_100_missing_file = open(basename + '.ensembl_100_missing_from_hgnc.txt', 'w')
ensembl_100_missing_file.write(bed_header)
for i, bed_name in enumerate(bed_gene_list):
    if bed_name in hgnc_gene_list:
        #bed_gene_line_list[i][1] = "chr" + bed_gene_line_list[i][1]
        for j, item in enumerate(bed_gene_line_list[i]):
            if j == len(bed_gene_line_list[i])-1:
                bed_gene_line_list[i][j] = item+'\n'
            else:
                bed_gene_line_list[i][j] = item+'\t'
        out_bed_file.writelines(bed_gene_line_list[i])
    else:
        for j, item in enumerate(bed_gene_line_list[i]):
            if j == len(bed_gene_line_list[i])-1:
                bed_gene_line_list[i][j] = item+'\n'
            else:
                bed_gene_line_list[i][j] = item+'\t'
        ensembl_100_missing_file.writelines(bed_gene_line_list[i])
ensembl_100_missing_file.close()
out_bed_file.close()

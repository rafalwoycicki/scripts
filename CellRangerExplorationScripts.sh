###counting number of reads per cells from the 10x Genomics Cell Ranger output matrix

zcat matrix.mtx.gz | grep -v '%' | awk '$3>0' | sed 's/ /\t/g' | tail -n +2 | cut -f2,3 | awk '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }' | sed 's/ /\t/' | sort -k2,2n > no_reads_per_cells.txt

cat no_reads_per_cells.txt | cut -f2 | sort -k1,1n | uniq -c | sed 's/^     / /' | sed 's/^    / /' | sed 's/^   / /' | sed 's/^  / /' | sed 's/^ //' | sed 's/ /\t/' > no_reads_per_cells1.dist

cat *no_reads_per_cells1.dist | awk '$2<2000&&$1<200' > no_reads_per_cells1.4graph

#gplot.pl -style points -using 2:1 -xlabel no_of_reads -ylabel frequency_of_cells -title no_of_cells_with_number_of_reads_L5 -type ps -outfile "no_of_cells_with_number_of_reads.ps" *no_reads_per_cells1.dist

#gplot.pl -style points -using 2:1 -xlabel no_of_reads -ylabel frequency_of_cells -title no_of_cells_with_number_of_reads_L5 -type ps -outfile "no_of_cells_with_number_of_reads_1.ps" *no_reads_per_cells1.4graph 


###counting number of reads per genes from the 10x Genomics Cell Ranger output matrix

zcat matrix.mtx.gz | grep -v '%' | awk '$3>0' | sed 's/ /\t/g' | tail -n +2 | cut -f1,3 | awk '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }' | sed 's/ /\t/' | sort -k2,2n > no_reads_per_genes.txt

cat no_reads_per_genes.txt | cut -f2 | sort -k1,1n | uniq -c | sed 's/^     / /' | sed 's/^    / /' | sed 's/^   / /' | sed 's/^  / /' | sed 's/^ //' | sed 's/ /\t/' > no_reads_per_genes1.dist

cat *no_reads_per_genes1.dist | awk '$2<200' > no_reads_per_genes1.4graph

#gplot.pl -style points -using 2:1 -xlabel with_no_of_reads -ylabel frequency_of_genes -title no_of_genes_with_number_of_reads_L5 -type ps -outfile "no_of_genes_with_number_of_reads.ps" *no_reads_per_genes1.dist

#gplot.pl -style points -using 2:1 -xlabel with_no_of_reads -ylabel frequency_of_genes -title no_of_genes_with_number_of_reads_L5 -type ps -outfile "no_of_genes_with_number_of_reads_1.ps" *no_reads_per_genes1.4graph


###counting number of cells with genes from the 10x Genomics Cell Ranger output matrix

zcat matrix.mtx.gz | grep -v '%' | awk '$3>0' | sed 's/ /\t/g' | tail -n +2 | cut -f1 | sort -k1,1 | uniq -c | sed 's/^     / /' | sed 's/^    / /' | sed 's/^   / /' | sed 's/^  / /' | sed 's/^ //' | sed 's/ /\t/' | sort -k1,1n > no_cells_with_genes.txt

cat no_cells_with_genes.txt | cut -f1 | sort -k1,1n | uniq -c | sed 's/^     / /' | sed 's/^    / /' | sed 's/^   / /' | sed 's/^  / /' | sed 's/^ //' | sed 's/ /\t/' > no_cells_with_genes1.dist

cat *no_cells_with_genes1.dist | awk '$2<200' > no_cells_with_genes1.4graph

#gplot.pl -style points -using 2:1 -xlabel no_of_cells_with_gene -ylabel frequency_of_genes -title no_of_genes_in_number_of_cells_L5 -type ps -outfile "no_of_genes_in_number_of_cells.ps" *no_cells_with_genes1.dist

#gplot.pl -style points -using 2:1 -xlabel no_of_cells_with_gene -ylabel frequency_of_genes -title no_of_genes_in_number_of_cells_L5 -type ps -outfile "no_of_genes_in_number_of_cells_1.ps" *no_cells_with_genes1.4graph


###counting number of genes per cells from the 10x Genomics Cell Ranger output matrix

zcat matrix.mtx.gz | grep -v '%' | awk '$3>0' | sed 's/ /\t/g' | tail -n +2 | cut -f2 | sort -k1,1 | uniq -c | sed 's/^     / /' | sed 's/^    / /' | sed 's/^   / /' | sed 's/^  / /' | sed 's/^ //' | sed 's/ /\t/' | sort -k1,1n > no_genes_per_cells.txt

cat no_genes_per_cells.txt | cut -f1 | sort -k1,1n | uniq -c | sed 's/^     / /' | sed 's/^    / /' | sed 's/^   / /' | sed 's/^  / /' | sed 's/^ //' | sed 's/ /\t/' > no_genes_per_cells1.dist

cat *no_genes_per_cells1.dist | awk '$2<2000&&$1<200' > no_genes_per_cells1.4graph

#gplot.pl -style points -using 2:1 -ylabel frequency_of_cells -xlabel no_gene_per_cell -title no_of_cells_with_number_of_genes_L5 -type ps -outfile "no_of_cells_with_number_of_genes.ps" *no_genes_per_cells1.dist

#gplot.pl -style points -using 2:1 -ylabel frequency_of_cells -xlabel no_gene_per_cell -title no_of_cells_with_number_of_genes_L5 -type ps -outfile "no_of_cells_with_number_of_genes_1.ps" *no_genes_per_cells1.4graph #OK





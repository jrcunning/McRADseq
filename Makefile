#output/topten_contigs_snps.txt: output/topten_contigs.csv
#	filter_fasta.py -f data/contigs/filtered_CSR_.92.fasta -o output/topten_contigs.fasta -s output/topten_contigs.csv
#	R --vanilla < R/get_topten_SNPs.R

output/BC_NC_summary_plus.txt: output/BC_NC_summary.txt
	R --vanilla < R/get_annot.R

output/BC_NC_summary.txt: output/BC_NC_contigs2nt_blastout_tabular.txt output/BC_NC_cohits2sp_blastout_tabular.txt output/BC_NC_cohits2tr_blastout_tabular.txt output/BC_NC_contigs2bac_blastout_tabular.txt output/BC_NC_contigs2vir_blastout_tabular.txt R/summ_data.R
	R --vanilla < R/summ_data.R

output/BC_NC_cohits2tr_blastout_tabular.txt: output/BC_NC_cohits.faa
	blastp -db /Volumes/CoralReefFutures/ref/uniprot_trembl/uniprot_trembl.fasta \
	-query output/BC_NC_cohits.faa -num_threads 4 \
	-outfmt 11 -out output/BC_NC_cohits2tr_blastout_archive.txt
	blast_formatter -archive output/BC_NC_cohits2tr_blastout_archive.txt \
	-outfmt 6 -num_alignments 1 -out output/BC_NC_cohits2tr_blastout_tabular.txt
	blast_formatter -archive output/BC_NC_cohits2tr_blastout_archive.txt \
	-outfmt 0 -num_alignments 1 -out output/BC_NC_cohits2tr_blastout_pairwise.txt

output/BC_NC_cohits2sp_blastout_tabular.txt: output/BC_NC_cohits.faa
	blastp -db data/ref/uniprot_sprot.fasta \
	-query output/BC_NC_cohits.faa -num_threads 4 \
	-outfmt 11 -out output/BC_NC_cohits2sp_blastout_archive.txt
	blast_formatter -archive output/BC_NC_cohits2sp_blastout_archive.txt \
	-outfmt 6 -num_alignments 1 -out output/BC_NC_cohits2sp_blastout_tabular.txt
	blast_formatter -archive output/BC_NC_cohits2sp_blastout_archive.txt \
	-outfmt 0 -num_alignments 1 -out output/BC_NC_cohits2sp_blastout_pairwise.txt

output/BC_NC_cohits.faa: output/BC_NC_contigs2co_blastout_tabular.txt
	sort -u -k1,1 output/BC_NC_contigs2co_blastout_tabular.txt | cut -f2 > output/BC_NC_contigs2co_besthitnames.txt
	filter_fasta.py -f data/ref/all_coral.faa -o output/BC_NC_cohits.faa -s output/BC_NC_contigs2co_besthitnames.txt
	rm output/BC_NC_contigs2co_besthitnames.txt

output/BC_NC_contigs2co_blastout_tabular.txt: output/BC_NC_contigs.fasta
	blastx -db data/ref/all_coral.faa \
	-query output/BC_NC_contigs.fasta -num_threads 4 \
	-evalue 1 -outfmt 11 \
	-out output/BC_NC_contigs2co_blastout_archive.txt
	blast_formatter -archive output/BC_NC_contigs2co_blastout_archive.txt \
	-outfmt 6 -num_alignments 1 -out output/BC_NC_contigs2co_blastout_tabular.txt
	blast_formatter -archive output/BC_NC_contigs2co_blastout_archive.txt \
	-outfmt 0 -num_alignments 1 -out output/BC_NC_contigs2co_blastout_pairwise.txt

output/BC_NC_contigs2nt_blastout_tabular.txt: output/BC_NC_contigs.fasta
	blastn -db /Volumes/CoralReefFutures/ref/ncbi_nt/nt \
	-query output/BC_NC_contigs.fasta -num_threads 4 \
	-evalue 1 -word_size 11 -gapopen 5 -gapextend 2 -outfmt 11 \
	-out output/BC_NC_contigs2nt_blastout_archive.txt
	blast_formatter -archive output/BC_NC_contigs2nt_blastout_archive.txt \
	-outfmt 6 -num_alignments 1 -out output/BC_NC_contigs2nt_blastout_tabular.txt
	blast_formatter -archive output/BC_NC_contigs2nt_blastout_archive.txt \
	-outfmt 0 -num_alignments 1 -out output/BC_NC_contigs2nt_blastout_pairwise.txt

output/BC_NC_contigs2bac_blastout_tabular.txt: output/BC_NC_contigs.fasta
	blastn -db /Volumes/CoralReefFutures/ref/bacteria/all_complete_Gb_bac.fasta \
	-query output/BC_NC_contigs.fasta -num_threads 4 \
	-evalue 1 -word_size 11 -gapopen 5 -gapextend 2 -outfmt 11 \
	-out output/BC_NC_contigs2bac_blastout_archive.txt
	blast_formatter -archive output/BC_NC_contigs2bac_blastout_archive.txt \
	-outfmt 6 -num_alignments 1 -out output/BC_NC_contigs2bac_blastout_tabular.txt
	blast_formatter -archive output/BC_NC_contigs2bac_blastout_archive.txt \
	-outfmt 0 -num_alignments 1 -out output/BC_NC_contigs2bac_blastout_pairwise.txt

output/BC_NC_contigs2vir_blastout_tabular.txt: output/BC_NC_contigs.fasta
	blastn -db /Volumes/CoralReefFutures/ref/viruses/phantome_db.fasta \
	-query output/BC_NC_contigs.fasta -num_threads 4 \
	-evalue 1 -word_size 11 -gapopen 5 -gapextend 2 -outfmt 11 \
	-out output/BC_NC_contigs2vir_blastout_archive.txt
	blast_formatter -archive output/BC_NC_contigs2vir_blastout_archive.txt \
	-outfmt 6 -num_alignments 1 -out output/BC_NC_contigs2vir_blastout_tabular.txt
	blast_formatter -archive output/BC_NC_contigs2vir_blastout_archive.txt \
	-outfmt 0 -num_alignments 1 -out output/BC_NC_contigs2vir_blastout_pairwise.txt

output/BC_NC_contigs.fasta: R/filter_SNPs.R
	R --vanilla < R/filter_SNPS.R


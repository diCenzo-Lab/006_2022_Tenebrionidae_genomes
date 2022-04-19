## Repeat the masking

perl ../switch_contig_names.pl superworm_cleaned.fasta > superworm_cleaned_renamed.fasta
sdust -w 64 -t 28 superworm_cleaned_renamed.fasta > superworm_cleaned_renamed_sdust.bed
RepeatMasker -pa 7 -species Insecta -xsmall superworm_cleaned_renamed.fasta
bedtools maskfasta -fi superworm_cleaned_renamed.fasta.masked -bed superworm_cleaned_renamed_sdust.bed -fo superworm_fixed_masked.fasta -soft
sed -i 's/contig_4179\b/contig_4179 [location=mitochondrion]/' superworm_fixed_masked.fasta
perl ../switch_scaffold_names.pl superworm_fixed_masked.fasta > superworm_fixed_masked_renamed.fasta

##Annotation

# Extract nuclear genome
perl ../extract_nuclear_genome.pl superworm_fixed_masked_renamed.fasta > superworm_nuclear_genome.fasta

# Trim RNA-seq reads
mkdir processed_reads/
bbduk.sh in=/datadisk1/Sequencing_data/Illumina_data/20220307-Insect-RNA-pilot/NS.1812.002.NEBNext_dual_i7_A12---NEBNext_dual_i5_A12.Superworm_test_RNA_R1.fastq.gz in2=/datadisk1/Sequencing_data/Illumina_data/20220307-Insect-RNA-pilot/NS.1812.002.NEBNext_dual_i7_A12---NEBNext_dual_i5_A12.Superworm_test_RNA_R2.fastq.gz ref=adapters,artifacts,phix,lambda out=processed_reads/forward.bbduk.fastq.gz out2=processed_reads/reverse.bbduk.fastq.gz
java -jar /datadisk1/Bioinformatics_programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 28 processed_reads/forward.bbduk.fastq.gz processed_reads/reverse.bbduk.fastq.gz processed_reads/forward.bbduk.trimmed.fastq.gz processed_reads/forward.bbduk.unpaired.fastq.gz processed_reads/reverse.bbduk.trimmed.fastq.gz processed_reads/reverse.bbduk.unpaired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Map RNA-seq data
mkdir star_genome/
STAR --runMode genomeGenerate --runThreadN 28 --genomeDir star_genome --genomeFastaFiles superworm_nuclear_genome.fasta --genomeSAindexNbases 13
STAR --runMode alignReads --runThreadN 28 --genomeDir star_genome --readFilesIn processed_reads/forward.bbduk.trimmed.fastq.gz,processed_reads/reverse.bbduk.trimmed.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --alignIntronMin 20 --alignIntronMax 1000000 --outFilterType BySJout --outFileNamePrefix star_aligned_pass_1_ --limitBAMsortRAM 20000000000
STAR --runMode alignReads --runThreadN 28 --genomeDir star_genome --readFilesIn processed_reads/forward.bbduk.trimmed.fastq.gz,processed_reads/reverse.bbduk.trimmed.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --alignIntronMin 20 --alignIntronMax 1000000 --outFilterType BySJout --outFileNamePrefix star_aligned_pass_2_ --limitBAMsortRAM 20000000000 --sjdbFileChrStartEnd star_aligned_pass_1_SJ.out.tab

# Annotate using RNA-seq data
braker.pl --species superworm_rna --genome superworm_fixed_masked_renamed.fasta --bam star_aligned_pass_2_Aligned.sortedByCoord.out.bam --softmasking --cores 28
mv braker braker_rna

# Annotate using proteins
cat ../beetle_reference_proteins.fasta ../superworm_busco_complete_proteins.fasta > protein_reference_set.fasta
sed -i 's/\*//g' protein_reference_set.fasta
braker.pl --species superworm_prot --genome superworm_fixed_masked_renamed.fasta --prot_seq protein_reference_set.fasta --softmasking --cores 28
mv braker braker_protein

# Combine outputs and prep for submission
cp ../tsebra_configuration.cfg .
tsebra.py -g braker_rna/augustus.hints.gtf,braker_protein/augustus.hints.gtf -c tsebra_configuration.cfg -e braker_rna/hintsfile.gff,braker_protein/hintsfile.gff -o braker_combined.gtf
rename_gtf.py --gtf braker_combined.gtf --prefix Zmor --translation_tab translation.tab --out braker_combined_renamed.gtf
gt gtf_to_gff3 <(grep -P "\tCDS\t|\texon\t" braker_combined_renamed.gtf ) > braker_combined_renamed.gff
gt gff3 -sort -tidy -fixregionboundaries -force -o temp.gff braker_combined_renamed.gff
mv temp.gff braker_combined_renamed.gff
mkdir final_annotation/
cp superworm_fixed_masked_renamed.fasta final_annotation/superworm_masked.fasta
cp braker_combined_renamed.gtf final_annotation/
cp braker_combined_renamed.gff final_annotation/
cp ../Zmor-template.sbt final_annotation/

# Fix annotation based on table2asn results
cd final_annotation/
perl removeContained.pl braker_combined_renamed.gff > test_removeContained.gff
grep 'CDS' test_removeContained.gff | cut -f2,2 -d ';' | cut -f2,2 -d'=' | sort -u > CDS_transcripts.txt
grep 'ID=mRNA' test_removeContained.gff | cut -f3,3 -d ';' | cut -f2,2 -d'=' | sort -u > mRNA_transcripts.txt
grep -v -f CDS_transcripts.txt mRNA_transcripts.txt > mRNA_missing_CDS.txt
grep -v -f mRNA_missing_CDS.txt test_removeContained.gff > test_removeContained_2.gff
tac test_removeContained_2.gff > test_removeContained_3.gff
perl removeContained.pl test_removeContained_3.gff > test_removeContained_4.gff
tac test_removeContained_4.gff > test_removeContained_5.gff
grep 'CDS' test_removeContained_5.gff | cut -f2,2 -d ';' | cut -f2,2 -d'=' | sort -u > CDS_transcripts.txt
grep 'ID=mRNA' test_removeContained_5.gff | cut -f3,3 -d ';' | cut -f2,2 -d'=' | sort -u > mRNA_transcripts.txt
grep -v -f CDS_transcripts.txt mRNA_transcripts.txt > mRNA_missing_CDS.txt
grep -v -f mRNA_missing_CDS.txt test_removeContained_5.gff > test_removeContained_6.gff
gt gff3 -sort -tidy -fixregionboundaries -force -o test_removeContained_7.gff test_removeContained_6.gff
grep -v ',' contained_locus.txt | grep -v 'Zmor_021870' | grep -v 'Zmor_022573' | grep -v 'Zmor_027130' | grep -v 'Zmor_003754' | grep -v 'Zmor_003794' | grep -v 'Zmor_004558' | grep -v 'Zmor_022443' | grep -v 'Zmor_004049' | grep -v 'Zmor_004163' | grep -v 'Zmor_004167' | grep -v 'Zmor_012104' | cut -f4,4 > contained_locus_toRemove.txt # contained_locus.txt taken from previous table2asn output
echo Zmor_023085 >> contained_locus_toRemove.txt
echo Zmor_027187 >> contained_locus_toRemove.txt
grep 'ID=gene' test_removeContained_7.gff | sed 's/_/\t/' | sed 's/ID=gene//' | sed 's/;/\t/' | sort -k2,2 -k10,10 -g | cut -f2,2 -d'=' > geneList.txt
grep 'locus-tag' Zophobas_morio.sqn | cut -f2,2 -d'"' > locusList.txt
paste geneList.txt locusList.txt > gene_to_locus.txt
grep -f contained_locus_toRemove.txt gene_to_locus.txt | cut -f1,1 > contained_genes_toRemove.txt
grep -v -w -f contained_genes_toRemove.txt test_removeContained_7.gff > test_removeContained_8.gff
gt gff3 -sort -tidy -fixregionboundaries -force -o test_removeContained_9.gff test_removeContained_8.gff
sed -i 's/\ ->\ /\t/' internal_stop.txt # internal_stop.txt comes from discrepancy report when submitted to NCBI
cut -f2,2 internal_stop.txt | sed 's/\[gnl|Zmor|cds.//' | sed 's/\]//' | sed 's/\[lcl|//' > internal_to_remove.txt
grep -v -w -f internal_to_remove.txt test_removeContained_9.gff > test_removeContained_10.gff
gt gff3 -sort -tidy -fixregionboundaries -force -o test_removeContained_11.gff test_removeContained_10.gff
grep -v -w 'mRNA21853' test_removeContained_11.gff > test_removeContained_12.gff
gt gff3 -sort -tidy -fixregionboundaries -force -o test_removeContained_13.gff test_removeContained_12.gff

# Run table2asn
table2asn -M n -J -c w -euk -t Zmor-template.sbt -gaps-min 10 -gaps-unknown 100 -l align-genus -j "[organism=Zophobas morio]" -i superworm_masked.fasta -f test_removeContained_13.gff -o Zophobas_morio.sqn -Z -locus-tag-prefix Zmor -V vb

# Get some stats
gt stat -addintrons -force -genelengthdistri -o superworm_gene_length_distribution.txt test_removeContained_13.gff
gt stat -addintrons -force -exonlengthdistri -o superworm_exon_length_distribution.txt test_removeContained_13.gff
gt stat -addintrons -force -exonnumberdistri -o superworm_exon_number_distribution.txt test_removeContained_13.gff
gt stat -addintrons -force -intronlengthdistri -o superworm_intron_length_distribution.txt test_removeContained_13.gff
gt stat -addintrons -force -cdslengthdistri -o superworm_cds_length_distribution.txt test_removeContained_13.gff

# Get proteins
cp ../../extractFAA.pl .
perl extractFAA.pl Zophobas_morio.gbf > Zophobas_morio.faa
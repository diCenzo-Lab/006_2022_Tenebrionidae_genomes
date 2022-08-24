## Find and fix mitochondrion
#Zophobas atratus downloaded from NCBI as Zata_mito.fasta

# Get files
mkdir annotation/
cd annotation/
cp ../00_final_assemblies/superworm.masurca-purged.flye-purged.quickmerge.purged.fasta.gz scaffolded_assembly.fasta.gz
cp ../03_merge_and_purge/16_final_assemblies/superworm.masurca-purged.flye-purged.quickmerge.purged.fasta.gz original_assembly.fasta.gz
gunzip *.fasta.gz

# Make blast database of both assemblies
makeblastdb -in scaffolded_assembly.fasta -out scaffolded_assembly -title scaffolded_assembly -dbtype 'nucl'
makeblastdb -in original_assembly.fasta -out original_assembly -title original_assembly -dbtype 'nucl'

# Run blast to find mitochondrion
blastn -query Zata_mito.fasta -db scaffolded_assembly -out scaffolded_blastresults.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send bitscore evalue sstrand'
blastn -query Zata_mito.fasta -db original_assembly -out original_blastresults.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send bitscore evalue sstrand'
head -15702 original_assembly.fasta | tail -2 > original_assembly_mito.fasta
blastn -query original_assembly_mito.fasta -db scaffolded_assembly -out compared_blastresults.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send bitscore evalue sstrand'

# Extract contaminated contig and remove from assembly
grep -A 1 'JAKCWB010000011_1_RagTag' scaffolded_assembly.fasta > mito_contaminated_contig.fasta
grep -v -f mito_contaminated_contig.fasta scaffolded_assembly.fasta > mito_removed_assembly.fasta

# Fix contaminated contig and add back to assembly
awk '{ print substr($1,27250,100000000) }' mito_contaminated_contig.fasta > temp.fasta
tail -1 temp.fasta > temp2.fasta
echo '>JAKCWB010000011_1_RagTag' > temp3.fasta
cat temp3.fasta temp2.fasta > mito_cleaned_contig.fasta
rm temp.fasta
rm temp2.fasta
rm temp3.fasta
cat mito_removed_assembly.fasta mito_cleaned_contig.fasta > cleaned_assembly_without_mitochondrion.fasta

# Extract mitochondrion, circularize, and add back to assembly
awk '{ print substr($1,1,27149) }' mito_contaminated_contig.fasta > mitochondrion_full_sequence.fasta
awk '{ print substr($1,700,15494) }' mitochondrion_full_sequence.fasta > temp.fasta
tail -1 temp.fasta > temp2.fasta
echo '>mitochondrion' > temp3.fasta
cat temp3.fasta temp2.fasta > circularized_mitochondrion.fasta
rm temp.fasta
rm temp2.fasta
rm temp3.fasta
cat cleaned_assembly_without_mitochondrion.fasta circularized_mitochondrion.fasta > superworm_fixed_assembly.fasta

# Confirm all looks good
makeblastdb -in superworm_fixed_assembly.fasta -out superworm_fixed_assembly -title superworm_fixed_assembly -dbtype 'nucl'
blastn -query Zata_mito.fasta -db superworm_fixed_assembly -out fixed_blastresults.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send bitscore evalue sstrand'

## Switch contig names

perl switch_contig_names.pl superworm_fixed_assembly.fasta > superworm_renamed_contigs.fasta

## Masking

sdust -w 64 -t 28 superworm_renamed_contigs.fasta > superworm_renamed_sdust.bed
RepeatMasker -pa 7 -species Insecta -xsmall superworm_renamed_contigs.fasta
bedtools maskfasta -fi superworm_renamed_contigs.fasta.masked -bed superworm_renamed_sdust.bed -fo superworm_masked.fasta -soft
sed -i 's/contig_4662\b/contig_4662 [location=mitochondrion]/' superworm_masked.fasta

## Fix contamination

# Remove mitochondrion and contaminated contigs
perl ../modifyFasta.pl superworm_masked.fasta > superworm_masked.modified.fasta
grep -w -v 'contig_2113' superworm_masked.modified.fasta | grep -w -v 'contig_2942' | grep -w -v 'contig_4043' | grep -w -v 'contig_4280' | grep -w -v 'contig_4327' | grep -w -v 'contig_4384' | grep -w -v 'contig_4472' | grep -w -v 'contig_4654' > superworm_masked.modified.reduced.fasta
grep 'mitochondrion' superworm_masked.modified.reduced.fasta > superworm_mitochondrion_final.fasta
sed -i 's/\t/\n/' superworm_mitochondrion_final.fasta
grep -v 'mitochondrion' superworm_masked.modified.reduced.fasta > superworm_masked.modified.reduced.noMito.fasta
grep -w -v 'contig_1' superworm_masked.modified.reduced.noMito.fasta | grep -w -v 'contig_19' | grep -w -v 'contig_2914' | grep -w -v 'contig_4521' | grep -w -v 'contig_5' | grep -w -v 'contig_963' > superworm_masked.modified.reduced.noMito.noContamScaffold.fasta
sed -i 's/\t/\n/' superworm_masked.modified.reduced.noMito.noContamScaffold.fasta
cp superworm_masked.modified.reduced.noMito.noContamScaffold.fasta superworm_masked.fixed.fasta

# Split contaminated contigs and add back the uncontaminated parts
grep -w 'contig_1' superworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,15404031)}' temp.fasta > temp_a.fasta
awk '{print substr($0,15404062,55083720)}' temp.fasta > temp_b.fasta
echo '>contig_1a' >> superworm_masked.fixed.fasta
cat superworm_masked.fixed.fasta temp_a.fasta > superworm_masked.fixed.fasta.tmp
echo '>contig_1b' >> superworm_masked.fixed.fasta.tmp
cat superworm_masked.fixed.fasta.tmp temp_b.fasta > superworm_masked.fixed.fasta

grep -w 'contig_19' superworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,4452)}' temp.fasta > temp_a.fasta
awk '{print substr($0,4479,66824)}' temp.fasta > temp_b.fasta
echo '>contig_19a' >> superworm_masked.fixed.fasta
cat superworm_masked.fixed.fasta temp_a.fasta > superworm_masked.fixed.fasta.tmp
echo '>contig_19b' >> superworm_masked.fixed.fasta.tmp
cat superworm_masked.fixed.fasta.tmp temp_b.fasta > superworm_masked.fixed.fasta

grep -w 'contig_2914' superworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,358)}' temp.fasta > temp_a.fasta
awk '{print substr($0,675,1877)}' temp.fasta > temp_b.fasta
echo '>contig_2914a' >> superworm_masked.fixed.fasta
cat superworm_masked.fixed.fasta temp_a.fasta > superworm_masked.fixed.fasta.tmp
echo '>contig_2914b' >> superworm_masked.fixed.fasta.tmp
cat superworm_masked.fixed.fasta.tmp temp_b.fasta > superworm_masked.fixed.fasta

grep -w 'contig_4521' superworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,454)}' temp.fasta > temp_a.fasta
echo '>contig_4521a' >> superworm_masked.fixed.fasta
cat superworm_masked.fixed.fasta temp_a.fasta > superworm_masked.fixed.fasta.tmp
mv superworm_masked.fixed.fasta.tmp superworm_masked.fixed.fasta

grep -w 'contig_963' superworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,420,6330)}' temp.fasta > temp_a.fasta
echo '>contig_963a' >> superworm_masked.fixed.fasta
cat superworm_masked.fixed.fasta temp_a.fasta > superworm_masked.fixed.fasta.tmp
mv superworm_masked.fixed.fasta.tmp superworm_masked.fixed.fasta

grep -w 'contig_5' superworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,957599)}' temp.fasta > temp_a.fasta
awk '{print substr($0,958311,31805947)}' temp.fasta > temp_b.fasta
echo '>contig_5a' >> superworm_masked.fixed.fasta
cat superworm_masked.fixed.fasta temp_a.fasta > superworm_masked.fixed.fasta.tmp
echo '>contig_5b' >> superworm_masked.fixed.fasta.tmp
cat superworm_masked.fixed.fasta.tmp temp_b.fasta > superworm_masked.fixed.fasta

# Add back in mitochondrion
cat superworm_masked.fixed.fasta superworm_mitochondrion_final.fasta > superworm_masked.fixed.fasta.tmp
mv superworm_masked.fixed.fasta.tmp superworm_masked.fixed.fasta

## Rescaffold the nuclear genome

# Unmask
sed -i 's/a/A/g' superworm_masked.fixed.fasta
sed -i 's/t/T/g' superworm_masked.fixed.fasta
sed -i 's/g/G/g' superworm_masked.fixed.fasta
sed -i 's/c/C/g' superworm_masked.fixed.fasta
sed -i 's/n/N/g' superworm_masked.fixed.fasta
sed -i 's/miToChoNdrioN/mitochondrion/' superworm_masked.fixed.fasta

# Extract nuclear genome
perl ../extract_nuclear_genome.pl superworm_masked.fixed.fasta > superworm_nuclear_genome.fasta
tail -2 superworm_masked.fixed.fasta > superworm_mito_genome.fasta

# Scaffold
cp ../17_final_scaffolds/GCA_022388445.1_Zophobas_atratus_KSU_v.1_genomic.fna.gz .
gunzip GCA_022388445.1_Zophobas_atratus_KSU_v.1_genomic.fna.gz
ragtag.py scaffold -t 24 -u -o ragtag_output GCA_022388445.1_Zophobas_atratus_KSU_v.1_genomic.fna superworm_nuclear_genome.fasta

# Add back in the mitochondrion
cat ragtag_output/ragtag.scaffold.fasta superworm_mito_genome.fasta > superworm_cleaned.fasta

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
cp ../../removeContained.pl .
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

# Add back in accidentally removed exons
grep '##gff-version 3' test_removeContained_9.gff > test_removeContained_10.gff
grep '##sequence-region' test_removeContained_9.gff >> test_removeContained_10.gff
grep 'ID=gene' test_removeContained_9.gff | cut -f2,2 -d';' > gene_ids_to_keep.txt
grep 'ID=mRNA' test_removeContained_9.gff | cut -f3,3 -d';' > transcript_ids_to_keep.txt
grep 'ID=mRNA' braker_combined_renamed.gff | cut -f3,3 -d';' > all_original_transcript_ids.txt
diff all_original_transcript_ids.txt transcript_ids_to_keep.txt | grep '<' | sed 's/<\ //' > transcript_ids_to_remove.txt
grep -w -f 'gene_ids_to_keep.txt' braker_combined_renamed.gff | grep -w -v -f 'transcript_ids_to_remove.txt' >> test_removeContained_10.gff
gt gff3 -sort -tidy -fixregionboundaries -force -o test_removeContained_11.gff test_removeContained_10.gff

# Remove genes fully contained within another gene on the same strand
cp ../../removeContainedGenes.pl .
perl removeContainedGenes.pl test_removeContained_11.gff | cut -f2,2 -d';' > contained_genes.txt
grep -v -w -f 'contained_genes.txt' test_removeContained_11.gff > test_removeContained_12.gff
gt gff3 -sort -tidy -fixregionboundaries -force -o test_removeContained_13.gff test_removeContained_12.gff

# Remove contained complete coding regions on the same strand
tac test_removeContained_13.gff > test_removeContained_14.gff
perl removeContained.pl test_removeContained_14.gff > test_removeContained_15.gff
tac test_removeContained_15.gff > test_removeContained_16.gff
grep 'CDS' test_removeContained_16.gff | cut -f2,2 -d ';' | cut -f2,2 -d'=' | sort -u > CDS_transcripts.txt
grep 'ID=mRNA' test_removeContained_16.gff | cut -f3,3 -d ';' | cut -f2,2 -d'=' | sort -u > mRNA_transcripts.txt
grep -v -f CDS_transcripts.txt mRNA_transcripts.txt > mRNA_missing_CDS.txt
grep -v -f mRNA_missing_CDS.txt test_removeContained_16.gff > test_removeContained_17.gff
gt gff3 -sort -tidy -fixregionboundaries -force -o test_removeContained_18.gff test_removeContained_17.gff
grep '##gff-version 3' test_removeContained_18.gff > test_removeContained_19.gff
grep '##sequence-region' test_removeContained_18.gff >> test_removeContained_19.gff
grep 'ID=gene' test_removeContained_18.gff | cut -f2,2 -d';' > gene_ids_to_keep.txt
grep 'ID=mRNA' test_removeContained_18.gff | cut -f3,3 -d';' > transcript_ids_to_keep.txt
grep 'ID=mRNA' braker_combined_renamed.gff | cut -f3,3 -d';' > all_original_transcript_ids.txt
diff all_original_transcript_ids.txt transcript_ids_to_keep.txt | grep '<' | sed 's/<\ //' > transcript_ids_to_remove.txt
grep -w -f 'gene_ids_to_keep.txt' braker_combined_renamed.gff | grep -w -v -f 'transcript_ids_to_remove.txt' >> test_removeContained_19.gff
gt gff3 -sort -tidy -fixregionboundaries -force -o test_removeContained_20.gff test_removeContained_19.gff

# Run table2asn
table2asn -M n -J -c w -euk -t Zmor-template.sbt -gaps-min 10 -gaps-unknown 100 -l align-genus -j "[organism=Zophobas morio]" -i superworm_masked.fasta -f test_removeContained_20.gff -o Zophobas_morio.sqn -Z -locus-tag-prefix Zmor -V vb -verbose

# Get some stats
gt stat -addintrons -force -genelengthdistri -o superworm_gene_length_distribution.txt test_removeContained_20.gff
gt stat -addintrons -force -exonlengthdistri -o superworm_exon_length_distribution.txt test_removeContained_20.gff
gt stat -addintrons -force -exonnumberdistri -o superworm_exon_number_distribution.txt test_removeContained_20.gff
gt stat -addintrons -force -intronlengthdistri -o superworm_intron_length_distribution.txt test_removeContained_20.gff
gt stat -addintrons -force -cdslengthdistri -o superworm_cds_length_distribution.txt test_removeContained_20.gff

# Get proteins
cp ../../extractFAA.pl .
perl extractFAA.pl Zophobas_morio.gbf > Zophobas_morio.faa

## Find and fix mitochondrion
#Zophobas atratus mitochondrion (Genbank accession MK140669.1) downloaded from NCBI as Zata_mito.fasta

# Get files
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

##Annotation

# Extract nuclear genome
perl extract_nuclear_genome.pl superworm_masked.fasta > superworm_nuclear_genome.fasta

# Trim RNA-seq reads
mkdir processed_reads/
bbduk.sh in=/datadisk1/Sequencing_data/Illumina_data/20220307-Insect-RNA-pilot/NS.1812.002.NEBNext_dual_i7_A12---NEBNext_dual_i5_A12.Superworm_test_RNA_R1.fastq.gz in2=/datadisk1/Sequencing_data/Illumina_data/20220307-Insect-RNA-pilot/NS.1812.002.NEBNext_dual_i7_A12---NEBNext_dual_i5_A12.Superworm_test_RNA_R2.fastq.gz ref=adapters,artifacts,phix,lambda out=processed_reads/forward.bbduk.fastq.gz out2=processed_reads/reverse.bbduk.fastq.gz
java -jar /datadisk1/Bioinformatics_programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 28 processed_reads/forward.bbduk.fastq.gz processed_reads/reverse.bbduk.fastq.gz processed_reads/forward.bbduk.trimmed.fastq.gz processed_reads/forward.bbduk.unpaired.fastq.gz processed_reads/reverse.bbduk.trimmed.fastq.gz processed_reads/reverse.bbduk.unpaired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Map RNA-seq data
mkdir star_genome/
STAR --runMode genomeGenerate --runThreadN 28 --genomeDir star_genome --genomeFastaFiles superworm_nuclear_genome.fasta --genomeSAindexNbases 13
STAR --runMode alignReads --runThreadN 28 --genomeDir star_genome --readFilesIn processed_reads/forward.bbduk.trimmed.fastq.gz,processed_reads/reverse.bbduk.trimmed.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --alignIntronMin 20 --alignIntronMax 1000000 --outFilterType BySJout --outFileNamePrefix star_aligned_pass_1_ --limitBAMsortRAM 11000000000
STAR --runMode alignReads --runThreadN 28 --genomeDir star_genome --readFilesIn processed_reads/forward.bbduk.trimmed.fastq.gz,processed_reads/reverse.bbduk.trimmed.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --alignIntronMin 20 --alignIntronMax 1000000 --outFilterType BySJout --outFileNamePrefix star_aligned_pass_2_ --limitBAMsortRAM 13000000000 --sjdbFileChrStartEnd star_aligned_pass_1_SJ.out.tab

# Annotate using RNA-seq data
braker.pl --species superworm_rna --genome superworm_masked.fasta --bam star_aligned_pass_2_Aligned.sortedByCoord.out.bam --softmasking --cores 28
mv braker braker_rna

# Annotate using proteins
# The protein fasta files of the T. castaneum genome annotation (ENA accession PRJNA12540) and a previously published T. molitor genome annotation (ENA accession PRJEB44755) were downloaded and concatenated to form the beetle_reference_proteins.fasta
# superworm_busco_complete_proteins.fasta contains the single-copy and multi-copy complete Endopterygota BUSCO genes identified in the assembly to be annotated
cat beetle_reference_proteins.fasta superworm_busco_complete_proteins.fasta > protein_reference_set.fasta
braker.pl --species superworm_prot --genome superworm_masked.fasta --prot_seq protein_reference_set.fasta --softmasking --cores 28
mv braker braker_protein

# Combine outputs and prep for submission
tsebra.py -g braker_rna/augustus.hints.gtf,braker_protein/augustus.hints.gtf -c tsebra_configuration.cfg -e braker_rna/hintsfile.gff,braker_protein/hintsfile.gff -o braker_combined.gtf
rename_gtf.py --gtf braker_combined.gtf --prefix Zmor --translation_tab translation.tab --out braker_combined_renamed.gtf
gt gtf_to_gff3 <(grep -P "\tCDS\t|\texon\t" braker_combined_renamed.gtf ) > braker_combined_renamed.gff
mkdir final_annotation/
cp superworm_masked.fasta final_annotation/
cp braker_combined_renamed.gtf final_annotation/
cp braker_combined_renamed.gff final_annotation/
cp Zmor-template.sbt final_annotation/
cd final_annotation/
table2asn -M n -J -c w -euk -t Zmor-template.sbt -gaps-min 10 -gaps-unknown 100 -l align-genus -j "[organism=Zophobas morio]" -i superworm_masked.fasta -f braker_combined_renamed.gff -o Zophoas_morio.sqn -Z
cd ..

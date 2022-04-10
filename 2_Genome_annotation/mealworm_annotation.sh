##To find mitochondria
#Tenebrio molitor mitochondrion (GenBank accession KF418153.1) was downloaded from NCBI as Tmol_mito.fasta

# Get files
cp ../00_final_assemblies/mealworm.flye.assembly.purged.fasta.gz scaffolded_assembly.fasta.gz
gunzip *.fasta.gz

# Make blast database of both assemblies
makeblastdb -in scaffolded_assembly.fasta -out scaffolded_assembly -title scaffolded_assembly -dbtype 'nucl'

# Run blast to find mitochondrion
blastn -query Tmol_mito.fasta -db scaffolded_assembly -out scaffolded_blastresults.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send bitscore evalue sstrand'

# Extract mitochondrion and add back to end of assembly
grep -A 1 'contig_17307_1_RagTag' scaffolded_assembly.fasta > mito_contig.fasta
grep -v -f mito_contig.fasta scaffolded_assembly.fasta > mito_removed_assembly.fasta
tail -1 mito_contig.fasta > temp2.fasta
echo '>mitochondrion' > temp3.fasta
cat temp3.fasta temp2.fasta > circularized_mitochondrion.fasta
rm temp.fasta
rm temp2.fasta
rm temp3.fasta
cat mito_removed_assembly.fasta circularized_mitochondrion.fasta > mealworm_fixed_assembly.fasta

# Confirm all looks good
makeblastdb -in mealworm_fixed_assembly.fasta -out mealworm_fixed_assembly -title mealworm_fixed_assembly -dbtype 'nucl'
blastn -query Tmol_mito.fasta -db mealworm_fixed_assembly -out fixed_blastresults.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send bitscore evalue sstrand'

## Switch contig names
perl switch_contig_names.pl mealworm_fixed_assembly.fasta > mealworm_renamed_contigs.fasta

## Masking
sdust -w 64 -t 28 mealworm_renamed_contigs.fasta > mealworm_renamed_sdust.bed
RepeatMasker -pa 7 -species Insecta -xsmall mealworm_renamed_contigs.fasta
bedtools maskfasta -fi mealworm_renamed_contigs.fasta.masked -bed mealworm_renamed_sdust.bed -fo mealworm_masked.fasta -soft
sed -i 's/contig_2142\b/contig_2142 [location=mitochondrion]/' mealworm_masked.fasta

##Annotation

# Extract nuclear genome
perl extract_nuclear_genome.pl mealworm_masked.fasta > mealworm_nuclear_genome.fasta

# Trim RNA-seq reads
mkdir processed_reads/
bbduk.sh in=/datadisk1/Sequencing_data/Illumina_data/20220307-Insect-RNA-pilot/NS.1812.002.NEBNext_dual_i7_G11---NEBNext_dual_i5_G11.Mealworm_test_RNA_R1.fastq.gz in2=/datadisk1/Sequencing_data/Illumina_data/20220307-Insect-RNA-pilot/NS.1812.002.NEBNext_dual_i7_G11---NEBNext_dual_i5_G11.Mealworm_test_RNA_R2.fastq.gz ref=adapters,artifacts,phix,lambda out=processed_reads/forward.bbduk.fastq.gz out2=processed_reads/reverse.bbduk.fastq.gz
java -jar /datadisk1/Bioinformatics_programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 28 processed_reads/forward.bbduk.fastq.gz processed_reads/reverse.bbduk.fastq.gz processed_reads/forward.bbduk.trimmed.fastq.gz processed_reads/forward.bbduk.unpaired.fastq.gz processed_reads/reverse.bbduk.trimmed.fastq.gz processed_reads/reverse.bbduk.unpaired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Map RNA-seq data
mkdir star_genome/
STAR --runMode genomeGenerate --runThreadN 28 --genomeDir star_genome --genomeFastaFiles mealworm_nuclear_genome.fasta --genomeSAindexNbases 12
STAR --runMode alignReads --runThreadN 28 --genomeDir star_genome --readFilesIn processed_reads/forward.bbduk.trimmed.fastq.gz,processed_reads/reverse.bbduk.trimmed.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --alignIntronMin 20 --alignIntronMax 1000000 --outFilterType BySJout --outFileNamePrefix star_aligned_pass_1_ --limitBAMsortRAM 11000000000
STAR --runMode alignReads --runThreadN 28 --genomeDir star_genome --readFilesIn processed_reads/forward.bbduk.trimmed.fastq.gz,processed_reads/reverse.bbduk.trimmed.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --alignIntronMin 20 --alignIntronMax 1000000 --outFilterType BySJout --outFileNamePrefix star_aligned_pass_2_ --limitBAMsortRAM 13000000000 --sjdbFileChrStartEnd star_aligned_pass_1_SJ.out.tab

# Annotate using RNA-seq data
braker.pl --species mealworm_rna --genome mealworm_masked.fasta --bam star_aligned_pass_2_Aligned.sortedByCoord.out.bam --softmasking --cores 28
mv braker braker_rna

# Annotate using proteins
# The protein fasta files of the T. castaneum genome annotation (ENA accession PRJNA12540) and a previously published T. molitor genome annotation (ENA accession PRJEB44755) were downloaded and concatenated to form the beetle_reference_proteins.fasta
# superworm_busco_complete_proteins.fasta contains the single-copy and multi-copy complete Endopterygota BUSCO genes identified in the assembly to be annotated
cat beetle_reference_proteins.fasta mealworm_busco_complete_proteins.fasta > protein_reference_set.fasta
braker.pl --species mealworm_prot --genome mealworm_masked.fasta --prot_seq protein_reference_set.fasta --softmasking --cores 28
mv braker braker_protein

# Combine outputs and prep for submission
tsebra.py -g braker_rna/augustus.hints.gtf,braker_protein/augustus.hints.gtf -c tsebra_configuration.cfg -e braker_rna/hintsfile.gff,braker_protein/hintsfile.gff -o braker_combined.gtf
rename_gtf.py --gtf braker_combined.gtf --prefix Tmol --translation_tab translation.tab --out braker_combined_renamed.gtf
gt gtf_to_gff3 <(grep -P "\tCDS\t|\texon\t" braker_combined_renamed.gtf ) > braker_combined_renamed.gff
mkdir final_annotation/
cp mealworm_masked.fasta final_annotation/
cp braker_combined_renamed.gtf final_annotation/
cp braker_combined_renamed.gff final_annotation/
cp Tmol-template.sbt final_annotation/
cd final_annotation/
table2asn -M n -J -c w -euk -t Tmol-template.sbt -gaps-min 10 -gaps-unknown 100 -l align-genus -j "[organism=Tenebrio molitor]" -i mealworm_masked.fasta -f braker_combined_renamed.gff -o Tenebrio_molitor.sqn -Z
cd ..

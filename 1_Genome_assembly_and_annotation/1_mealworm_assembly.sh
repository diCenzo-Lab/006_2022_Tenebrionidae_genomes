## Trim the illumina reads

# Trim using platanus
~/scratch/platanus_trim -t 64 sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq

## Perform the flye assembly and polishing

# Assembly with flye
mkdir genome_assembly/
mkdir genome_assembly/mealworm_assembly/
mkdir genome_assembly/mealworm_assembly/01_flye/
cd genome_assembly/mealworm_assembly/01_flye/
mkdir 01_flye_assembly/
flye --nano-raw /datadisk1/Nanopore_data/20210827-Mealworm-genome/Basecalled_data/passed.gz --outdir 01_flye_assembly/ --threads 8
pigz -r -p 24 01_flye_assembly/
gunzip 01_flye_assembly/assembly.fasta.gz

# Polishing with racon and nanopore data
mkdir 02_racon_correction/
minimap2 -x map-ont -d 02_racon_correction/index.mni 01_flye_assembly/assembly.fasta
minimap2 -t 16 -a 02_racon_correction/index.mni /datadisk1/Nanopore_data/20210827-Mealworm-genome/Basecalled_data/passed.gz > 02_racon_correction/mapped.sam
racon -m 8 -x -6 -g -8 -w 500 -t 14 /datadisk1/Nanopore_data/20210827-Mealworm-genome/Basecalled_data/passed.gz 02_racon_correction/mapped.sam 01_flye_assembly/assembly.fasta > 02_racon_correction/assembly_racon.fasta
pigz -r -p 24 02_racon_correction/
gunzip 02_racon_correction/assembly_racon.fasta.gz

# Polishing with medaka and nanopore data
mkdir 03_medaka_correction/
medaka_consensus -i /datadisk1/Nanopore_data/20210827-Mealworm-genome/Basecalled_data/passed.gz -d 02_racon_correction/assembly_racon.fasta -o 03_medaka_correction -t 16 -b 50 -m r941_min_hac_g507
pigz -r -p 24 03_medaka_correction/
gunzip 03_medaka_correction/consensus.fasta.gz

# Polish with pilon and illumina data
mkdir 04_pilon_polishing/
cd ../superworm/
bowtie2-build --threads 64 03_medaka_correction/consensus.fasta 04_pilon_polishing/medaka_consensus
bowtie2 -p 64 -x 04_pilon_polishing/medaka_consensus -1 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq.trimmed -2 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq.trimmed -S 04_pilon_polishing/illumina_mapped.sam 2> 04_pilon_polishing/illumina_mapped.txt
samtools sort -@ 64 04_pilon_polishing/illumina_mapped.sam > 04_pilon_polishing/illumina_mapped_sorted.sam
samtools view -S -b -@ 64 04_pilon_polishing/illumina_mapped_sorted.sam > 04_pilon_polishing/illumina_mapped_sorted.bam
samtools index -@ 64 04_pilon_polishing/illumina_mapped_sorted.bam
bioawk -c fastx '{print $name}' 03_medaka_correction/consensus.fasta > 04_pilon_polishing/nms
mkdir -p 04_pilon_polishing/split
samtools faidx 03_medaka_correction/consensus.fasta
for CHR in $(cat 04_pilon_polishing/nms);do
  samtools faidx 03_medaka_correction/consensus.fasta $CHR > 04_pilon_polishing/split/$CHR\.fa;
done
mkdir 04_pilon_polishing/segments
for NUM in 1 2 3 4 5 6 7 8 9;do
  cat 04_pilon_polishing/split/contig_${NUM}*.fa > 04_pilon_polishing/segments/${NUM}.fa
done
mkdir 04_pilon_polishing/polished_segments/
ls 04_pilon_polishing/segments/*fa > 04_pilon_polishing/toprocess
JAVA_TOOL_OPTIONS="-Xmx450G -Xss2560k"
for CHRFILE in $(cat 04_pilon_polishing/toprocess);do
  echo $CHRFILE
  CHRNAME=$(basename $CHRFILE | cut -f 1 -d '.')
  CHROUTDIR=04_pilon_polishing/polished_segments/$CHRNAME
  java -Xmx450G -jar $EBROOTPILON/pilon-1.24.jar --threads 12 --diploid --genome $CHRFILE --frags 04_pilon_polishing/illumina_mapped_sorted.bam --output $CHRNAME --outdir $CHROUTDIR --changes --vcf > 04_pilon_polishing/assembly_pilon.txt
done
cat 04_pilon_polishing/polished_segments/*/*.fasta > 04_pilon_polishing/assembly_pilon.fasta
pigz -r -p 64 04_pilon_polishing/
gunzip 04_pilon_polishing/assembly_pilon.fasta.gz

# Polish with hapo-g and illumina data
mkdir 05_hapog_polishing/
bowtie2-build --threads 64 04_pilon_polishing/assembly_pilon.fasta 05_hapog_polishing/assembly_pilon
bowtie2 -p 64 -x 05_hapog_polishing/assembly_pilon -1 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq.trimmed -2 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq.trimmed -S 05_hapog_polishing/illumina_mapped.sam 2> 05_hapog_polishing/illumina_mapped.txt
samtools sort -@ 64 05_hapog_polishing/illumina_mapped.sam > 05_hapog_polishing/illumina_mapped_sorted.sam
samtools view -S -b -@ 64 05_hapog_polishing/illumina_mapped_sorted.sam > 05_hapog_polishing/illumina_mapped_sorted.bam
samtools index -@ 64 05_hapog_polishing/illumina_mapped_sorted.bam
hapog.py --genome 05_hapog_polishing/assembly_pilon.fasta -b 05_hapog_polishing/illumina_mapped_sorted.bam -o 05_hapog_polishing/output -t 24 -u &
pigz -r -p 64 05_hapog_polishing/
gunzip 05_hapog_polishing/output/hapog_results/hapog.fasta.gz

# Finalize assembly
cp 05_hapog_polishing/output/hapog_results/hapog.fasta 00_final_assembly/polished_genome.fasta
sed -i 's/_pilon_polished//' 00_final_assembly/polished_genome.fasta
sed -i 's/_pilon//' 00_final_assembly/polished_genome.fasta
cd ..

## Perform the masurca assembly and polishing

# Assembly with masurca
mkdir 02_masurca/
cd 02_masurca/
mkdir 01_masurca_assembly/
cd 01_masurca_assembly/
masurca -t 64 -i ../../../../sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq,../../../../sequencing_data/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq -r ../../../../sequencing_data/mealworm/passed.fastq
# Masurca froze while creating mr.41.17.15.0.02.txt.tmp. The job was therefore killed and restarted, which made a backup of the mr.41.17.15.0.02.txt.tmp file and tried to continue making the rest of the file from where it left off.
masurca -t 64 -i ../../../../sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq,../../../../sequencing_data/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq -r ../../../../sequencing_data/mealworm/passed.fastq
# Masurca again quickly froze and so the job was again killed. The mr.41.17.15.0.02.txt.tmp file and the backup of this file were combined, and reads in the longest_reads.25x.fa file that were not present in the combined mr.41.17.15.0.02.txt.tmp file were removed from the the longest_reads.25x.fa file. The command was run once again.
masurca -t 64 -i ../../../../sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq,../../../../sequencing_data/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq -r ../../../../sequencing_data/mealworm/passed.fastq
cd ..
pigz -r -p 64 01_masurca_assembly/
gunzip 01_masurca_assembly/CA.mr.41.17.15.0.02/primary.genome.scf.fasta.gz

# Polish with pilon and illumina data
mkdir 02_pilon_polishing/
bowtie2-build  --threads 64 01_masurca_assembly/CA.mr.41.17.15.0.02/primary.genome.scf.fasta 02_pilon_polishing/masurca_assembly
bowtie2 -p 64 -x 02_pilon_polishing/masurca_assembly -1 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq.trimmed -2 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq.trimmed -S 02_pilon_polishing/illumina_mapped.sam 2> 02_pilon_polishing/illumina_mapped.txt
samtools sort -@ 64 02_pilon_polishing/illumina_mapped.sam > 02_pilon_polishing/illumina_mapped_sorted.sam
samtools view -S -b -@ 64 02_pilon_polishing/illumina_mapped_sorted.sam > 02_pilon_polishing/illumina_mapped_sorted.bam
samtools index -@ 64 02_pilon_polishing/illumina_mapped_sorted.bam
sh /home/gd38/scratch/insect_genomes/genome_assembly/flye/superworm/run_pilon.sh
bioawk -c fastx '{print $name}' 01_masurca_assembly/CA.mr.41.17.15.0.02/primary.genome.scf.fasta > 02_pilon_polishing/nms
mkdir -p 02_pilon_polishing/split
samtools faidx 01_masurca_assembly/CA.mr.41.17.15.0.02/consensus.fasta
for CHR in $(cat 02_pilon_polishing/nms);do
  samtools faidx 01_masurca_assembly/CA.mr.41.17.15.0.02/primary.genome.scf.fasta $CHR > 02_pilon_polishing/split/$CHR\.fa;
done
mkdir 02_pilon_polishing/segments
for NUM in 5 6 7;do
  cat 02_pilon_polishing/split/scf71800000${NUM}*.fa > 02_pilon_polishing/segments/${NUM}.fa
done
mkdir 02_pilon_polishing/polished_segments/
ls 02_pilon_polishing/segments/*fa > 02_pilon_polishing/toprocess
JAVA_TOOL_OPTIONS="-Xmx450G -Xss2560k"
for CHRFILE in $(cat 02_pilon_polishing/toprocess);do
  echo $CHRFILE
  CHRNAME=$(basename $CHRFILE | cut -f 1 -d '.')
  CHROUTDIR=02_pilon_polishing/polished_segments/$CHRNAME
  java -Xmx450G -jar $EBROOTPILON/pilon-1.24.jar --threads 12 --diploid --genome $CHRFILE --frags 02_pilon_polishing/illumina_mapped_sorted.bam --output $CHRNAME --outdir $CHROUTDIR --changes --vcf > 02_pilon_polishing/assembly_pilon.txt
done
cat 02_pilon_polishing/polished_segments/*/*.fasta > 02_pilon_polishing/assembly_pilon.fasta
pigz -r -p 64 02_pilon_polishing/
gunzip 02_pilon_polishing/assembly_pilon.fasta.gz

# Polish with hapo-g and illumina data
mkdir 05_hapog_polishing/
bowtie2-build --threads 64 02_pilon_polishing/assembly_pilon.fasta 03_hapog_polishing/assembly_pilon
bowtie2 -p 64 -x 03_hapog_polishing/assembly_pilon -1 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq.trimmed -2 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq.trimmed -S 03_hapog_polishing/illumina_mapped.sam 2> 03_hapog_polishing/illumina_mapped.txt
samtools sort -@ 64 03_hapog_polishing/illumina_mapped.sam > 03_hapog_polishing/illumina_mapped_sorted.sam
samtools view -S -b -@ 64 03_hapog_polishing/illumina_mapped_sorted.sam > 03_hapog_polishing/illumina_mapped_sorted.bam
samtools index -@ 64 03_hapog_polishing/illumina_mapped_sorted.bam
hapog.py --genome 03_hapog_polishing/assembly_pilon.fasta -b 03_hapog_polishing/illumina_mapped_sorted.bam -o 03_hapog_polishing/output -t 24 -u &
pigz -r -p 64 03_hapog_polishing/
gunzip 03_hapog_polishing/output/hapog_results/hapog.fasta.gz

# Finalize assembly
cp 03_hapog_polishing/output/hapog_results/hapog.fasta 00_final_assembly/polished_genome.fasta
sed -i 's/_pilon_polished//' 00_final_assembly/polished_genome.fasta
sed -i 's/_pilon//' 00_final_assembly/polished_genome.fasta
cd ..

## Merge genomes and purge duplicates

# Prepare folders and input genomes
mkdir 03_merge_and_purge/
cd 03_merge_and_purge/
mkdir 00_input_genomes/
mkdir 01_purged_flye/
mkdir 02_purged_masurca/
mkdir 03_quickmerge_original/
mkdir 04_quickmerge_purged/
mkdir 05_ragtag_original/
mkdir 07_purged_quickmerge_original_1/
mkdir 08_purged_quickmerge_original_2/
mkdir 09_purged_quickmerge_purged_1/
mkdir 10_purged_quickmerge_purged_2/
mkdir 11_purged_ragtag_original_1/
mkdir 12_purged_ragtag_original_2/
mkdir 13_purged_ragtag_purged_1/
mkdir 14_purged_ragtag_purged_2/
mkdir 16_final_assemblies/
mkdir 17_final_scaffolds/
cp ../01_flye/00_final_assembly/polished_genome.fasta 00_input_genomes/mealworm.flye.assembly.fasta
cp ../02_masurca/00_final_assembly/polished_genome.fasta 00_input_genomes/mealworm.masurca.assembly.fasta

# Purge flye assemblies
cd 01_purged_flye/
bowtie2-build  --threads 64 ../00_input_genomes/mealworm.flye.assembly.fasta genome_assembly
bowtie2 -p 64 -x genome_assembly -1 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq.trimmed -2 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq.trimmed -S illumina_mapped.sam 2> illumina_mapped.txt
samtools view -S -b -@ 64 illumina_mapped.sam > illumina_mapped.bam
split_fa ../00_input_genomes/mealworm.flye.assembly.fasta > assembly.split
minimap2 -xasm5 -DP assembly.split assembly.split | gzip -c - > assembly.split.self.paf.gz
ngscstat illumina_mapped.bam
calcuts TX.stat > cutoffs 2> calcults.log
purge_dups -2 -T cutoffs -c TX.base.cov assembly.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed ../00_input_genomes/mealworm.flye.assembly.fasta
mv purged.fa mealworm.flye.assembly.purged.fasta
pigz -p 64 *
cp mealworm.flye.assembly.purged.fasta.gz ../00_input_genomes/mealworm.flye.assembly.purged.fasta.gz
gunzip ../00_input_genomes/mealworm.flye.assembly.purged.fasta.gz
cd ..

# Purge masurca assemblies
bowtie2-build  --threads 64 ../00_input_genomes/mealworm.masurca.assembly.fasta genome_assembly
bowtie2 -p 64 -x genome_assembly -1 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq.trimmed -2 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq.trimmed -S illumina_mapped.sam 2> illumina_mapped.txt
samtools view -S -b -@ 64 illumina_mapped.sam > illumina_mapped.bam
split_fa ../00_input_genomes/mealworm.masurca.assembly.fasta > assembly.split
minimap2 -xasm5 -DP assembly.split assembly.split | gzip -c - > assembly.split.self.paf.gz
ngscstat illumina_mapped.bam
calcuts TX.stat > cutoffs 2> calcults.log
purge_dups -2 -T cutoffs -c TX.base.cov assembly.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed ../00_input_genomes/mealworm.masurca.assembly.fasta
mv purged.fa mealworm.masurca.assembly.purged.fasta
pigz -p 64 *
cp mealworm.masurca.assembly.purged.fasta.gz ../00_input_genomes/mealworm.masurca.assembly.purged.fasta.gz
gunzip ../00_input_genomes/mealworm.masurca.assembly.purged.fasta.gz
cd ..

# Merge assemblies using quickmerge
cd 03_quickmerge_original/
merge_wrapper.py -l 220000 -ml 10000 -t 16 -v --prefix mealworm.flye.masurca.quickmerge ../00_input_genomes/mealworm.flye.assembly.fasta ../00_input_genomes/mealworm.masurca.assembly.fasta
merge_wrapper.py -l 190000 -ml 10000 -t 16 -v --prefix mealworm.masurca.flye.quickmerge ../00_input_genomes/mealworm.masurca.assembly.fasta ../00_input_genomes/mealworm.flye.assembly.fasta
cd ../04_quickmerge_purged/
merge_wrapper.py -l 240000 -ml 10000 -t 16 -v --prefix mealworm.flye-purged.masurca-purged.quickmerge ../00_input_genomes/mealworm.flye.assembly.purged.fasta ../00_input_genomes/mealworm.masurca.assembly.purged.fasta
merge_wrapper.py -l 200000 -ml 10000 -t 16 -v --prefix mealworm.masurca-purged.flye-purged.quickmerge ../00_input_genomes/mealworm.masurca.assembly.purged.fasta ../00_input_genomes/mealworm.flye.assembly.purged.fasta
cd ..
pigz -r -p 24 03_quickmerge_original/
pigz -r -p 24 04_quickmerge_purged/

# Patch assemblies using RagTag
cd 05_ragtag_original
ragtag.py patch --nucmer-params '--maxmatch -l 100 -c 500 -t 16' -o mealworm.flye.masurca.ragtag ../00_input_genomes/mealworm.flye.assembly.fasta ../00_input_genomes/mealworm.masurca.assembly.fasta
ragtag.py patch --nucmer-params '--maxmatch -l 100 -c 500 -t 16' -o mealworm.masurca.flye.ragtag ../00_input_genomes/mealworm.masurca.assembly.fasta ../00_input_genomes/mealworm.flye.assembly.fasta
cd ../06_ragtag_purged/
ragtag.py patch --nucmer-params '--maxmatch -l 100 -c 500 -t 16' -o mealworm.flye-purged.masurca-purged.ragtag ../00_input_genomes/mealworm.flye.assembly.purged.fasta ../00_input_genomes/mealworm.masurca.assembly.purged.fasta
ragtag.py patch --nucmer-params '--maxmatch -l 100 -c 500 -t 16' -o mealworm.masurca-purged.flye-purged.ragtag ../00_input_genomes/mealworm.masurca.assembly.purged.fasta ../00_input_genomes/mealworm.flye.assembly.purged.fasta
cd ..
pigz -r -p 24 05_ragtag_original/
pigz -r -p 24 06_ragtag_purged/

# Copy genomes
cp 00_input_genomes/*.fasta 16_final_assemblies/
cp 03_quickmerge_original/merged_mealworm.flye.masurca.quickmerge.fasta.gz 16_final_assemblies/mealworm.flye.masurca.quickmerge.fasta.gz
cp 03_quickmerge_original/merged_mealworm.masurca.flye.quickmerge.fasta.gz 16_final_assemblies/mealworm.masurca.flye.quickmerge.fasta.gz
cp 04_quickmerge_purged/merged_mealworm.flye-purged.masurca-purged.quickmerge.fasta.gz 16_final_assemblies/mealworm.flye-purged.masurca-purged.quickmerge.fasta.gz
cp 04_quickmerge_purged/merged_mealworm.masurca-purged.flye-purged.quickmerge.fasta.gz 16_final_assemblies/mealworm.masurca-purged.flye-purged.quickmerge.fasta.gz
cp 05_ragtag_original/mealworm.flye.masurca.ragtag/ragtag.patch.fasta.gz 16_final_assemblies/mealworm.flye.masurca.ragtag.fasta.gz
cp 05_ragtag_original/mealworm.masurca.flye.ragtag/ragtag.patch.fasta.gz 16_final_assemblies/mealworm.masurca.flye.ragtag.fasta.gz
cp 06_ragtag_purged/mealworm.flye-purged.masurca-purged.ragtag/ragtag.patch.fasta.gz 16_final_assemblies/mealworm.flye-purged.masurca-purged.ragtag.fasta.gz
cp 06_ragtag_purged/mealworm.masurca-purged.flye-purged.ragtag/ragtag.patch.fasta.gz 16_final_assemblies/mealworm.masurca-purged.flye-purged.ragtag.fasta.gz
gunzip 16_final_assemblies/*

# Purge merged assemblies 1
cd 07_purged_quickmerge_original_1/
bowtie2-build  --threads 64 ../16_final_assemblies/mealworm.flye.masurca.quickmerge.fasta genome_assembly
bowtie2 -p 64 -x genome_assembly -1 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq.trimmed -2 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq.trimmed -S illumina_mapped.sam 2> illumina_mapped.txt
samtools view -S -b -@ 64 illumina_mapped.sam > illumina_mapped.bam
split_fa ../16_final_assemblies/mealworm.flye.masurca.quickmerge.fasta > assembly.split
minimap2 -xasm5 -DP assembly.split assembly.split | gzip -c - > assembly.split.self.paf.gz
ngscstat illumina_mapped.bam
calcuts TX.stat > cutoffs 2> calcults.log
purge_dups -2 -T cutoffs -c TX.base.cov assembly.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed ../16_final_assemblies/mealworm.flye.masurca.quickmerge.fasta
mv purged.fa mealworm.flye.masurca.quickmerge.purged.fasta
pigz -p 64 *

# Purge merged assemblies 2
cd ../08_purged_quickmerge_original_2/
bowtie2-build  --threads 64 ../16_final_assemblies/mealworm.masurca.flye.quickmerge.fasta genome_assembly
bowtie2 -p 64 -x genome_assembly -1 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq.trimmed -2 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq.trimmed -S illumina_mapped.sam 2> illumina_mapped.txt
samtools view -S -b -@ 64 illumina_mapped.sam > illumina_mapped.bam
split_fa ../16_final_assemblies/mealworm.masurca.flye.quickmerge.fasta > assembly.split
minimap2 -xasm5 -DP assembly.split assembly.split | gzip -c - > assembly.split.self.paf.gz
ngscstat illumina_mapped.bam
calcuts TX.stat > cutoffs 2> calcults.log
purge_dups -2 -T cutoffs -c TX.base.cov assembly.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed ../16_final_assemblies/mealworm.masurca.flye.quickmerge.fasta
mv purged.fa mealworm.masurca.flye.quickmerge.purged.fasta
pigz -p 64 *

# Purge merged assemblies 3
cd ../09_purged_quickmerge_purged_1/
bowtie2-build  --threads 64 ../16_final_assemblies/mealworm.flye-purged.masurca-purged.quickmerge.fasta genome_assembly
bowtie2 -p 64 -x genome_assembly -1 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq.trimmed -2 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq.trimmed -S illumina_mapped.sam 2> illumina_mapped.txt
samtools view -S -b -@ 64 illumina_mapped.sam > illumina_mapped.bam
split_fa ../16_final_assemblies/mealworm.flye-purged.masurca-purged.quickmerge.fasta > assembly.split
minimap2 -xasm5 -DP assembly.split assembly.split | gzip -c - > assembly.split.self.paf.gz
ngscstat illumina_mapped.bam
calcuts TX.stat > cutoffs 2> calcults.log
purge_dups -2 -T cutoffs -c TX.base.cov assembly.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed ../16_final_assemblies/mealworm.flye-purged.masurca-purged.quickmerge.fasta
mv purged.fa mealworm.flye-purged.masurca-purged.quickmerge.purged.fasta
pigz -p 64 *

# Purge merged assemblies 4
cd ../10_purged_quickmerge_purged_2/
bowtie2-build  --threads 64 ../16_final_assemblies/mealworm.masurca-purged.flye-purged.quickmerge.fasta genome_assembly
bowtie2 -p 64 -x genome_assembly -1 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq.trimmed -2 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq.trimmed -S illumina_mapped.sam 2> illumina_mapped.txt
samtools view -S -b -@ 64 illumina_mapped.sam > illumina_mapped.bam
split_fa ../16_final_assemblies/mealworm.masurca-purged.flye-purged.quickmerge.fasta > assembly.split
minimap2 -xasm5 -DP assembly.split assembly.split | gzip -c - > assembly.split.self.paf.gz
ngscstat illumina_mapped.bam
calcuts TX.stat > cutoffs 2> calcults.log
purge_dups -2 -T cutoffs -c TX.base.cov assembly.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed ../16_final_assemblies/mealworm.masurca-purged.flye-purged.quickmerge.fasta
mv purged.fa mealworm.masurca-purged.flye-purged.quickmerge.purged.fasta
pigz -p 64 *

# Purge merged assemblies 5
cd ../11_purged_ragtag_original_1/
bowtie2-build  --threads 64 ../16_final_assemblies/mealworm.flye.masurca.ragtag.fasta genome_assembly
bowtie2 -p 64 -x genome_assembly -1 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq.trimmed -2 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq.trimmed -S illumina_mapped.sam 2> illumina_mapped.txt
samtools view -S -b -@ 64 illumina_mapped.sam > illumina_mapped.bam
split_fa ../16_final_assemblies/mealworm.flye.masurca.ragtag.fasta > assembly.split
minimap2 -xasm5 -DP assembly.split assembly.split | gzip -c - > assembly.split.self.paf.gz
ngscstat illumina_mapped.bam
calcuts TX.stat > cutoffs 2> calcults.log
purge_dups -2 -T cutoffs -c TX.base.cov assembly.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed ../16_final_assemblies/mealworm.flye.masurca.ragtag.fasta
mv purged.fa mealworm.flye.masurca.ragtag.purged.fasta
pigz -p 64 *

# Purge merged assemblies 6
cd ../12_purged_ragtag_original_2/
bowtie2-build  --threads 64 ../16_final_assemblies/mealworm.masurca.flye.ragtag.fasta genome_assembly
bowtie2 -p 64 -x genome_assembly -1 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq.trimmed -2 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq.trimmed -S illumina_mapped.sam 2> illumina_mapped.txt
samtools view -S -b -@ 64 illumina_mapped.sam > illumina_mapped.bam
split_fa ../16_final_assemblies/mealworm.masurca.flye.ragtag.fasta > assembly.split
minimap2 -xasm5 -DP assembly.split assembly.split | gzip -c - > assembly.split.self.paf.gz
ngscstat illumina_mapped.bam
calcuts TX.stat > cutoffs 2> calcults.log
purge_dups -2 -T cutoffs -c TX.base.cov assembly.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed ../16_final_assemblies/mealworm.masurca.flye.ragtag.fasta
mv purged.fa mealworm.masurca.flye.ragtag.purged.fasta
pigz -p 64 *

# Purge merged assemblies 7
cd ../13_purged_ragtag_purged_1/
bowtie2-build  --threads 64 ../16_final_assemblies/mealworm.flye-purged.masurca-purged.ragtag.fasta genome_assembly
bowtie2 -p 64 -x genome_assembly -1 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq.trimmed -2 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq.trimmed -S illumina_mapped.sam 2> illumina_mapped.txt
samtools view -S -b -@ 64 illumina_mapped.sam > illumina_mapped.bam
split_fa ../16_final_assemblies/mealworm.flye-purged.masurca-purged.ragtag.fasta > assembly.split
minimap2 -xasm5 -DP assembly.split assembly.split | gzip -c - > assembly.split.self.paf.gz
ngscstat illumina_mapped.bam
calcuts TX.stat > cutoffs 2> calcults.log
purge_dups -2 -T cutoffs -c TX.base.cov assembly.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed ../16_final_assemblies/mealworm.flye-purged.masurca-purged.ragtag.fasta
mv purged.fa mealworm.flye-purged.masurca-purged.ragtag.purged.fasta
pigz -p 64 *

# Purge merged assemblies 8
cd ../14_purged_ragtag_purged_2/
bowtie2-build  --threads 64 ../16_final_assemblies/mealworm.masurca-purged.flye-purged.ragtag.fasta genome_assembly
bowtie2 -p 64 -x genome_assembly -1 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R1.fastq.trimmed -2 ~/scratch/insect_genomes/sequencing_data/mealworm/NS.1706.004.IDT_i7_71---IDT_i5_71.Mealworm_R2.fastq.trimmed -S illumina_mapped.sam 2> illumina_mapped.txt
samtools view -S -b -@ 64 illumina_mapped.sam > illumina_mapped.bam
split_fa ../16_final_assemblies/mealworm.masurca-purged.flye-purged.ragtag.fasta > assembly.split
minimap2 -xasm5 -DP assembly.split assembly.split | gzip -c - > assembly.split.self.paf.gz
ngscstat illumina_mapped.bam
calcuts TX.stat > cutoffs 2> calcults.log
purge_dups -2 -T cutoffs -c TX.base.cov assembly.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed ../16_final_assemblies/mealworm.masurca-purged.flye-purged.ragtag.fasta
mv purged.fa mealworm.masurca-purged.flye-purged.ragtag.purged.fasta
pigz -p 64 *
cd ..

# Copy the assemblies
cp 07_*/*.fasta.gz 16_final_assemblies/
cp 08_*/*.fasta.gz 16_final_assemblies/
cp 09_*/*.fasta.gz 16_final_assemblies/
cp 10_*/*.fasta.gz 16_final_assemblies/
cp 11_*/*.fasta.gz 16_final_assemblies/
cp 12_*/*.fasta.gz 16_final_assemblies/
cp 13_*/*.fasta.gz 16_final_assemblies/
cp 14_*/*.fasta.gz 16_final_assemblies/
gunzip 16_final_assemblies/*

## Scaffold assemblies

# Scaffold using RagTag and the reference CAJRHG02.fasta
cd 17_final_scaffolds/
cp ~/Downloads/CAJRHG02.fasta.gz .
gunzip CAJRHG02.fasta.gz
ragtag.py scaffold -t 24 -u -o mealworm.flye.assembly.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.flye.assembly.fasta
ragtag.py scaffold -t 24 -u -o mealworm.flye.assembly.purged.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.flye.assembly.purged.fasta
ragtag.py scaffold -t 24 -u -o mealworm.flye.masurca.quickmerge.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.flye.masurca.quickmerge.fasta
ragtag.py scaffold -t 24 -u -o mealworm.flye.masurca.quickmerge.purged.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.flye.masurca.quickmerge.purged.fasta
ragtag.py scaffold -t 24 -u -o mealworm.flye.masurca.ragtag.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.flye.masurca.ragtag.fasta
ragtag.py scaffold -t 24 -u -o mealworm.flye.masurca.ragtag.purged.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.flye.masurca.ragtag.purged.fasta
ragtag.py scaffold -t 24 -u -o mealworm.flye-purged.masurca-purged.quickmerge.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.flye-purged.masurca-purged.quickmerge.fasta
ragtag.py scaffold -t 24 -u -o mealworm.flye-purged.masurca-purged.quickmerge.purged.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.flye-purged.masurca-purged.quickmerge.purged.fasta
ragtag.py scaffold -t 24 -u -o mealworm.flye-purged.masurca-purged.ragtag.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.flye-purged.masurca-purged.ragtag.fasta
ragtag.py scaffold -t 24 -u -o mealworm.flye-purged.masurca-purged.ragtag.purged.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.flye-purged.masurca-purged.ragtag.purged.fasta
ragtag.py scaffold -t 24 -u -o mealworm.masurca.assembly.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.masurca.assembly.fasta
ragtag.py scaffold -t 24 -u -o mealworm.masurca.assembly.purged.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.masurca.assembly.purged.fasta
ragtag.py scaffold -t 24 -u -o mealworm.masurca.flye.quickmerge.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.masurca.flye.quickmerge.fasta
ragtag.py scaffold -t 24 -u -o mealworm.masurca.flye.quickmerge.purged.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.masurca.flye.quickmerge.purged.fasta
ragtag.py scaffold -t 24 -u -o mealworm.masurca.flye.ragtag.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.masurca.flye.ragtag.fasta
ragtag.py scaffold -t 24 -u -o mealworm.masurca.flye.ragtag.purged.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.masurca.flye.ragtag.purged.fasta
ragtag.py scaffold -t 24 -u -o mealworm.masurca-purged.flye-purged.quickmerge.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.masurca-purged.flye-purged.quickmerge.fasta
ragtag.py scaffold -t 24 -u -o mealworm.masurca-purged.flye-purged.quickmerge.purged.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.masurca-purged.flye-purged.quickmerge.purged.fasta
ragtag.py scaffold -t 24 -u -o mealworm.masurca-purged.flye-purged.ragtag.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.masurca-purged.flye-purged.ragtag.fasta
ragtag.py scaffold -t 24 -u -o mealworm.masurca-purged.flye-purged.ragtag.purged.fasta CAJRHG02.fasta ../16_final_assemblies/mealworm.masurca-purged.flye-purged.ragtag.purged.fasta
cd ..
pigz -r -p 24 17_final_scaffolds/
cd ../

## Get final assemblies

mkdir 00_final_assemblies
cp 03_merge_and_purge/17_final_scaffolds/mealworm.flye.assembly.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.flye.assembly.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.flye.assembly.purged.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.flye.assembly.purged.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.flye.masurca.quickmerge.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.flye.masurca.quickmerge.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.flye.masurca.quickmerge.purged.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.flye.masurca.quickmerge.purged.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.flye.masurca.ragtag.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.flye.masurca.ragtag.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.flye.masurca.ragtag.purged.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.flye.masurca.ragtag.purged.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.flye-purged.masurca-purged.quickmerge.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.flye-purged.masurca-purged.quickmerge.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.flye-purged.masurca-purged.quickmerge.purged.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.flye-purged.masurca-purged.quickmerge.purged.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.flye-purged.masurca-purged.ragtag.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.flye-purged.masurca-purged.ragtag.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.flye-purged.masurca-purged.ragtag.purged.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.flye-purged.masurca-purged.ragtag.purged.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.masurca.assembly.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.masurca.assembly.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.masurca.assembly.purged.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.masurca.assembly.purged.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.masurca.flye.quickmerge.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.masurca.flye.quickmerge.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.masurca.flye.quickmerge.purged.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.masurca.flye.quickmerge.purged.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.masurca.flye.ragtag.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.masurca.flye.ragtag.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.masurca.flye.ragtag.purged.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.masurca.flye.ragtag.purged.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.masurca-purged.flye-purged.quickmerge.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.masurca-purged.flye-purged.quickmerge.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.masurca-purged.flye-purged.quickmerge.purged.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.masurca-purged.flye-purged.quickmerge.purged.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.masurca-purged.flye-purged.ragtag.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.masurca-purged.flye-purged.ragtag.fasta.gz
cp 03_merge_and_purge/17_final_scaffolds/mealworm.masurca-purged.flye-purged.ragtag.purged.fasta/ragtag.scaffold.fasta.gz 00_final_assemblies/mealworm.masurca-purged.flye-purged.ragtag.purged.fasta.gz
gunzip 00_final_assemblies/*
sed -i 's/|/_/g' 00_final_assemblies/*
sed -i 's/\./_/g' 00_final_assemblies/*

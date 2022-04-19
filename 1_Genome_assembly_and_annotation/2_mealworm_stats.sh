##Genome quality assessment
stats.sh mealworm.flye.assembly.fasta
stats.sh mealworm.flye.assembly.purged.fasta
stats.sh mealworm.flye.masuca.quickmerge.fasta
stats.sh mealworm.flye.masuca.quickmerge.purged.fasta
stats.sh mealworm.flye.masurca.ragtag.fasta
stats.sh mealworm.flye.masurca.ragtag.purged.fasta
stats.sh mealworm.flye-purged.masurca-purged.quickmerge.fasta
stats.sh mealworm.flye-purged.masurca-purged.quickmerge.purged.fasta
stats.sh mealworm.flye-purged.masurca-purged.ragtag.fasta
stats.sh mealworm.flye-purged.masurca-purged.ragtag.purged.fasta
stats.sh mealworm.masurca.assembly.fasta
stats.sh mealworm.masurca.assembly.purged.fasta
stats.sh mealworm.masurca.flye.quickmerge.fasta
stats.sh mealworm.masurca.flye.quickmerge.purged.fasta
stats.sh mealworm.masurca.flye.ragtag.fasta
stats.sh mealworm.masurca.flye.ragtag.purged.fasta
stats.sh mealworm.masurca-purged.flye-purged.quickmerge.fasta
stats.sh mealworm.masurca-purged.flye-purged.quickmerge.purged.fasta
stats.sh mealworm.masurca-purged.flye-purged.ragtag.fasta
stats.sh mealworm.masurca-purged.flye-purged.ragtag.purged.fasta

## BUSCO analysis
#Commands to load dataset and unzip the files
wget https://busco-data.ezlab.org/v5/data/lineages/endopterygota_odb10.2020-09-10.tar.gz
wget https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2020-09-10.tar.gz
mkdir busco_downloads/
mkdir busco_downloads/lineages/
mv endopterygota_odb10.2020-09-10.tar.gz busco_downloads/lineages/
mv eukaryota_odb10.2020-09-10.tar.gz busco_downloads/lineages/
cd busco_downloads/lineages/
tar -xvf endopterygota_odb10.2020-09-10.tar.gz
tar -xvf eukaryota_odb10.2020-09-10.tar.gz
cd ../..

#To run BUSCO
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye.assembly.fasta --out mealworm.flye.assembly.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye.assembly.fasta --out mealworm.flye.assembly.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye.assembly.purged.fasta --out mealworm.flye.assembly.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye.assembly.purged.fasta --out mealworm.flye.assembly.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye.masurca.quickmerge.fasta --out mealworm.flye.masurca.quickmerge.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye.masurca.quickmerge.fasta --out mealworm.flye.masurca.quickmerge.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye.masurca.quickmerge.purged.fasta --out mealworm.flye.masurca.quickmerge.purged.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye.masurca.quickmerge.purged.fasta --out mealworm.flye.masurca.quickmerge.purged.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/ mealworm.flye.masurca.ragtag.fasta --out mealworm.flye.masurca.ragtag.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/ mealworm.flye.masurca.ragtag.fasta --out mealworm.flye.masurca.ragtag.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye.masurca.ragtag.purged.fasta --out mealworm.flye.masurca.ragtag.purged.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/ mealworm.flye.masurca.ragtag.purged.fasta --out mealworm.flye.masurca.ragtag.purged.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye-purged.masurca-purged.quickmerge.fasta --out mealworm.flye-purged.masurca-purged.quickmerge.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/ mealworm.flye-purged.masurca-purged.quickmerge.fasta --out mealworm.flye-purged.masurca-purged.quickmerge.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye-purged.masurca-purged.quickmerge.purged.fasta --out mealworm.flye-purged.masurca-purged.quickmerge.purged.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye-purged.masurca-purged.quickmerge.purged.fasta --out mealworm.flye-purged.masurca-purged.quickmerge.purged.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye-purged.masurca-purged.ragtag.fasta --out mealworm.flye-purged.masurca-purged.ragtag.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye-purged.masurca-purged.ragtag.fasta --out mealworm.flye-purged.masurca-purged.ragtag.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye-purged.masurca-purged.ragtag.purged.fasta --out mealworm.flye-purged.masurca-purged.ragtag.purged.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.flye-purged.masurca-purged.ragtag.purged.fasta --out mealworm.flye-purged.masurca-purged.ragtag.purged.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca.assembly.fasta --out mealworm.masurca.assembly.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca.assembly.fasta --out mealworm.masurca.assembly.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca.assembly.purged.fasta --out mealworm.masurca.assembly.purged.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca.assembly.purged.fasta --out mealworm.masurca.assembly.purged.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca.flye.quickmerge.fasta --out mealworm.masurca.flye.quickmerge.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca.flye.quickmerge.fasta --out mealworm.masurca.flye.quickmerge.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca.flye.quickmerge.purged.fasta --out mealworm.masurca.flye.quickmerge.purged.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca.flye.quickmerge.purged.fasta --out mealworm.masurca.flye.quickmerge.purged.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca.flye.ragtag.fasta --out mealworm.masurca.flye.ragtag.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca.flye.ragtag.fasta --out mealworm.masurca.flye.ragtag.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca.flye.ragtag.purged.fasta --out mealworm.masurca.flye.ragtag.purged.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca.flye.ragtag.purged.fasta --out mealworm.masurca.flye.ragtag.purged.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca-purged.flye-purged.quickmerge.fasta --out mealworm.masurca-purged.flye-purged.quickmerge.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca-purged.flye-purged.quickmerge.fasta --out mealworm.masurca-purged.flye-purged.quickmerge.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca-purged.flye-purged.quickmerge.purged.fasta --out mealworm.masurca-purged.flye-purged.quickmerge.purged.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca-purged.flye-purged.quickmerge.purged.fasta --out mealworm.masurca-purged.flye-purged.quickmerge.purged.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca-purged.flye-purged.ragtag.fasta --out mealworm.masurca-purged.flye-purged.ragtag.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca-purged.flye-purged.ragtag.fasta --out mealworm.masurca-purged.flye-purged.ragtag.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca-purged.flye-purged.ragtag.purged.fasta --out mealworm.masurca-purged.flye-purged.ragtag.purged.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 64 -f
busco --offline --in assembled_genomes/mealworm_assemblies/mealworm.masurca-purged.flye-purged.ragtag.purged.fasta --out mealworm.masurca-purged.flye-purged.ragtag.purged.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 64 -f

##Jellyfish - for Genomescope
jellyfish count -C -m 21 -s 3000000000 -t 64 mealworm_illumina.fastq -o reads.jf
jellyfish histo -t 10 reads.jf > reads.histo

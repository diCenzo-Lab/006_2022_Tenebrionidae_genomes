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

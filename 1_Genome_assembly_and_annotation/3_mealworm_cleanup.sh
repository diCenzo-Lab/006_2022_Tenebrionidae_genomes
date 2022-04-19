##To find mitochondria
#Tenebrio molitor downloaded from NCBI as Tmol_mito.fasta

# Get files
mkdir annotation/
cd annotation/
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

perl ../switch_contig_names.pl mealworm_fixed_assembly.fasta > mealworm_renamed_contigs.fasta

## Masking

sdust -w 64 -t 28 mealworm_renamed_contigs.fasta > mealworm_renamed_sdust.bed
RepeatMasker -pa 7 -species Insecta -xsmall mealworm_renamed_contigs.fasta
bedtools maskfasta -fi mealworm_renamed_contigs.fasta.masked -bed mealworm_renamed_sdust.bed -fo mealworm_masked.fasta -soft
sed -i 's/contig_2142\b/contig_2142 [location=mitochondrion]/' mealworm_masked.fasta

## Fix contamination

# Remove mitochondrion and contaminated contigs
perl ../modifyFasta.pl mealworm_masked.fasta > mealworm_masked.modified.fasta
grep -w -v 'contig_685' mealworm_masked.modified.fasta > mealworm_masked.modified.reduced.fasta
grep 'mitochondrion' mealworm_masked.modified.reduced.fasta > mealworm_mitochondrion_final.fasta
sed -i 's/\t/\n/' mealworm_mitochondrion_final.fasta
grep -v 'mitochondrion' mealworm_masked.modified.reduced.fasta > mealworm_masked.modified.reduced.noMito.fasta
grep -w -v 'contig_1' mealworm_masked.modified.reduced.noMito.fasta | grep -w -v 'contig_11' | grep -w -v 'contig_114' | grep -w -v 'contig_14' | grep -w -v 'contig_15' | grep -w -v 'contig_1897' | grep -w -v 'contig_2' | grep -w -v 'contig_2050' | grep -w -v 'contig_227' | grep -w -v 'contig_331' | grep -w -v 'contig_376' | grep -w -v 'contig_421' | grep -w -v 'contig_474' | grep -w -v 'contig_5' | grep -w -v 'contig_56' | grep -w -v 'contig_6' | grep -w -v 'contig_642' | grep -w -v 'contig_7' | grep -w -v 'contig_8' | grep -w -v 'contig_859' | grep -w -v 'contig_9' > mealworm_masked.modified.reduced.noMito.noContamScaffold.fasta
sed -i 's/\t/\n/' mealworm_masked.modified.reduced.noMito.noContamScaffold.fasta
cp mealworm_masked.modified.reduced.noMito.noContamScaffold.fasta mealworm_masked.fixed.fasta

# Split contaminated contigs and add back the uncontaminated parts
grep -w 'contig_1' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,22440986)}' temp.fasta > temp_a.fasta
awk '{print substr($0,22441019,32026178)}' temp.fasta > temp_b.fasta
echo '>contig_1a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_1b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_11' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,6018224)}' temp.fasta > temp_a.fasta
awk '{print substr($0,6021436,9729241)}' temp.fasta > temp_b.fasta
echo '>contig_11a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_11b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_114' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,3486,4965)}' temp.fasta > temp_a.fasta
echo '>contig_114a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
mv mealworm_masked.fixed.fasta.tmp mealworm_masked.fixed.fasta

grep -w 'contig_2050' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,9416)}' temp.fasta > temp_a.fasta
awk '{print substr($0,9446,9761)}' temp.fasta > temp_b.fasta
echo '>contig_2050a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_2050b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_15' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,3638717)}' temp.fasta > temp_a.fasta
awk '{print substr($0,3639654,5916947)}' temp.fasta > temp_b.fasta
echo '>contig_15a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_15b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_1897' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,34,4169)}' temp.fasta > temp_a.fasta
echo '>contig_1897a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
mv mealworm_masked.fixed.fasta.tmp mealworm_masked.fixed.fasta

grep -w 'contig_227' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,4491)}' temp.fasta > temp_a.fasta
echo '>contig_227a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
mv mealworm_masked.fixed.fasta.tmp mealworm_masked.fixed.fasta

grep -w 'contig_331' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,17422)}' temp.fasta > temp_a.fasta
echo '>contig_331a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
mv mealworm_masked.fixed.fasta.tmp mealworm_masked.fixed.fasta

grep -w 'contig_376' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,1081)}' temp.fasta > temp_a.fasta
awk '{print substr($0,1108,1833)}' temp.fasta > temp_b.fasta
echo '>contig_376a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_376b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_421' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,1327)}' temp.fasta > temp_a.fasta
awk '{print substr($0,1358,5633)}' temp.fasta > temp_b.fasta
echo '>contig_421a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_421b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_474' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,15613)}' temp.fasta > temp_a.fasta
awk '{print substr($0,15644,16693)}' temp.fasta > temp_b.fasta
echo '>contig_474a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_474b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_5' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,10407455)}' temp.fasta > temp_a.fasta
awk '{print substr($0,10407488,21544811)}' temp.fasta > temp_b.fasta
echo '>contig_5a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_5b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_56' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,258)}' temp.fasta > temp_a.fasta
awk '{print substr($0,288,8427)}' temp.fasta > temp_b.fasta
echo '>contig_56a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_56b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_6' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,1712607)}' temp.fasta > temp_a.fasta
awk '{print substr($0,1712642,20827433)}' temp.fasta > temp_b.fasta
echo '>contig_6a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_6b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_642' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,392)}' temp.fasta > temp_a.fasta
awk '{print substr($0,419,2229)}' temp.fasta > temp_b.fasta
echo '>contig_642a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_642b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_7' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,1456861)}' temp.fasta > temp_a.fasta
awk '{print substr($0,1460148,16874873)}' temp.fasta > temp_b.fasta
echo '>contig_7a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_7b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_859' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,619)}' temp.fasta > temp_a.fasta
awk '{print substr($0,657,4639)}' temp.fasta > temp_b.fasta
echo '>contig_859a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_859b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_9' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,11219748)}' temp.fasta > temp_a.fasta
awk '{print substr($0,11219775,12982946)}' temp.fasta > temp_b.fasta
echo '>contig_9a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_9b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_14' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,2912305)}' temp.fasta > temp_a.fasta
awk '{print substr($0,2912333,130557)}' temp.fasta > temp_b.fasta
awk '{print substr($0,3046337,4980929)}' temp.fasta > temp_c.fasta
echo '>contig_14a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_14b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta
echo '>contig_14c' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_c.fasta > mealworm_masked.fixed.fasta.tmp
mv mealworm_masked.fixed.fasta.tmp mealworm_masked.fixed.fasta

grep -w 'contig_2' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,7467433)}' temp.fasta > temp_a.fasta
awk '{print substr($0,7467469,1183714)}' temp.fasta > temp_b.fasta
awk '{print substr($0,8651218,2139476)}' temp.fasta > temp_c.fasta
awk '{print substr($0,10792608,3061560)}' temp.fasta > temp_d.fasta
awk '{print substr($0,13856498,5092858)}' temp.fasta > temp_e.fasta
awk '{print substr($0,18950830,24300819)}' temp.fasta > temp_f.fasta
echo '>contig_2a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_2b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta
echo '>contig_2c' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_c.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_2d' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_d.fasta > mealworm_masked.fixed.fasta
echo '>contig_2e' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_e.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_2f' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_f.fasta > mealworm_masked.fixed.fasta

grep -w 'contig_8' mealworm_masked.modified.reduced.noMito.fasta | sed 's/\t/\n/' | tail -1 > temp.fasta
awk '{print substr($0,1,12899099)}' temp.fasta > temp_a.fasta
awk '{print substr($0,12899133,773256)}' temp.fasta > temp_b.fasta
awk '{print substr($0,13672422,15733334)}' temp.fasta > temp_c.fasta
echo '>contig_8a' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_a.fasta > mealworm_masked.fixed.fasta.tmp
echo '>contig_8b' >> mealworm_masked.fixed.fasta.tmp
cat mealworm_masked.fixed.fasta.tmp temp_b.fasta > mealworm_masked.fixed.fasta
echo '>contig_8c' >> mealworm_masked.fixed.fasta
cat mealworm_masked.fixed.fasta temp_c.fasta > mealworm_masked.fixed.fasta.tmp
mv mealworm_masked.fixed.fasta.tmp mealworm_masked.fixed.fasta

# Add back in mitochondrion
cat mealworm_masked.fixed.fasta mealworm_mitochondrion_final.fasta > mealworm_masked.fixed.fasta.tmp
mv mealworm_masked.fixed.fasta.tmp mealworm_masked.fixed.fasta

## Rescaffold the nuclear genome

# Unmask
sed -i 's/a/A/g' mealworm_masked.fixed.fasta
sed -i 's/t/T/g' mealworm_masked.fixed.fasta
sed -i 's/g/G/g' mealworm_masked.fixed.fasta
sed -i 's/c/C/g' mealworm_masked.fixed.fasta
sed -i 's/n/N/g' mealworm_masked.fixed.fasta
sed -i 's/miToChoNdrioN/mitochondrion/' mealworm_masked.fixed.fasta

# Extract nuclear genome
perl ../extract_nuclear_genome.pl mealworm_masked.fixed.fasta > mealworm_nuclear_genome.fasta
tail -2 mealworm_masked.fixed.fasta > mealworm_mito_genome.fasta

# Scaffold
cp ../17_final_scaffolds/CAJRHG02.fasta.gz .
gunzip CAJRHG02.fasta.gz
ragtag.py scaffold -t 24 -u -o ragtag_output CAJRHG02.fasta mealworm_nuclear_genome.fasta

# Add back in the mitochondrion
cat ragtag_output/ragtag.scaffold.fasta mealworm_mito_genome.fasta > mealworm_cleaned.fasta

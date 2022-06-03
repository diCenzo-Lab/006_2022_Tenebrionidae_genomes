# Create maximum likelihood phylogeny from orthofinder single copy orthogroups
mkdir orthofinder_output/Results_Apr16/Single_Copy_Orthologue_Sequences2/ # Make directory
perl scripts/rename_sequences.pl # Rename protein sequences with a consistent naming convention
mkdir orthofinder_output/Results_Apr16/mafft/ # Make directory
mkdir orthofinder_output/Results_Apr16/trimal/ # Make directory
perl scripts/align_trim.pl # Align and trim each set of orthologs
mkdir orthofinder_output/Results_Apr16/trimal_modified/ # Make directory
perl scripts/modify_trimal.pl # Modify the trimmed alignment file for concatenation
mkdir orthofinder_output/Results_Apr16/maximum_likelihood_phylogeny/ # Make directory
echo Asbolus_verrucosus > orthofinder_output/Results_Apr16/species_list.txt # Make a file of species names
echo Tenebrio_molitor_new >> orthofinder_output/Results_Apr16/species_list.txt # Make a file of species names
echo Tenebrio_molitor_old >> orthofinder_output/Results_Apr16/species_list.txt # Make a file of species names
echo Tribolium_castaneum >> orthofinder_output/Results_Apr16/species_list.txt # Make a file of species names
echo Tribolium_madens >> orthofinder_output/Results_Apr16/species_list.txt # Make a file of species names
echo Zophobas_morio >> orthofinder_output/Results_Apr16/species_list.txt # Make a file of species names
perl scripts/combine_alignments.pl > orthofinder_output/Results_Apr16/maximum_likelihood_phylogeny/final_alignment.fasta # Combine the trimmed alignments as one super alignment
cd orthofinder_output/Results_Apr16/maximum_likelihood_phylogeny/ # Change directory
raxmlHPC-HYBRID-AVX2 -T 28 -s final_alignment.fasta -N 5 -n test_phylogeny -f a -p 12345 -x 12345 -m PROTGAMMAAUTO # Determine the best model for the phylogeny
raxmlHPC-HYBRID-AVX2 -T 28 -s final_alignment.fasta -N 100 -n final_phylogeny -f a -p 12345 -x 12345 -m PROTGAMMAJTTF # Create a maximum likelihood phylogeny
cd ../../ # Change directory

# Create maximum likelihood phylogeny from BUSCO single copy complete genes
busco --in Genome_sequences/Asbolus_verrucosus.fna --out Asbolus_verrucosus.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 24 -f # Run BUSCO
busco --in Genome_sequences/Tenebrio_molitor.fna --out Tenebrio_molitor.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 24 -f # Run BUSCO
busco --in Genome_sequences/Tenebrio_molitor_v2.fna --out Tenebrio_molitor_v2.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 24 -f # Run BUSCO
busco --in Genome_sequences/Tenebrio_molitor_v3.fna --out Tenebrio_molitor_v3.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 24 -f # Run BUSCO
busco --in Genome_sequences/Tribolium_castaneum.fna --out Tribolium_castaneum.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 24 -f # Run BUSCO
busco --in Genome_sequences/Tribolium_confusum.fna --out Tribolium_confusum.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 24 -f # Run BUSCO
busco --in Genome_sequences/Tribolium_freemani.fna --out Tribolium_freemani.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 24 -f # Run BUSCO
busco --in Genome_sequences/Tribolium_madens.fna --out Tribolium_madens.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 24 -f # Run BUSCO
busco --in Genome_sequences/Zophobas_atratus.fna --out Zophobas_atratus.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 24 -f # Run BUSCO
busco --in Genome_sequences/Zophobas_morio.fna --out Zophobas_morio.eukaryota.busco --lineage_dataset eukaryota_odb10 --mode genome --cpu 24 -f # Run BUSCO
mv *.eukaryota.busco BUSCO_output/ # Move BUSCO output
mkdir BUSCO_output/gene_lists/ # Make directory
ls -1 BUSCO_output/Asbolus_verrucosus.eukaryota.busco/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences/ > BUSCO_output/gene_lists/Asbolus_verrucosus.txt # Get single-copy gene list
ls -1 BUSCO_output/Tenebrio_molitor.eukaryota.busco/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences/ > BUSCO_output/gene_lists/Tenebrio_molitor.txt # Get single-copy gene list
ls -1 BUSCO_output/Tenebrio_molitor_v2.eukaryota.busco/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences/ > BUSCO_output/gene_lists/Tenebrio_molitor_v2.txt # Get single-copy gene list
ls -1 BUSCO_output/Tenebrio_molitor_v3.eukaryota.busco/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences/ > BUSCO_output/gene_lists/Tenebrio_molitor_v3.txt # Get single-copy gene list
ls -1 BUSCO_output/Tribolium_castaneum.eukaryota.busco/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences/ > BUSCO_output/gene_lists/Tribolium_castaneum.txt # Get single-copy gene list
ls -1 BUSCO_output/Tribolium_confusum.eukaryota.busco/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences/ > BUSCO_output/gene_lists/Tribolium_confusum.txt # Get single-copy gene list
ls -1 BUSCO_output/Tribolium_freemani.eukaryota.busco/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences/ > BUSCO_output/gene_lists/Tribolium_freemani.txt # Get single-copy gene list
ls -1 BUSCO_output/Tribolium_madens.eukaryota.busco/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences/ > BUSCO_output/gene_lists/Tribolium_madens.txt # Get single-copy gene list
ls -1 BUSCO_output/Zophobas_atratus.eukaryota.busco/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences/ > BUSCO_output/gene_lists/Zophobas_atratus.txt # Get single-copy gene list
ls -1 BUSCO_output/Zophobas_morio.eukaryota.busco/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences/ > BUSCO_output/gene_lists/Zophobas_morio.txt # Get single-copy gene list
comm -12 <(comm -12 <(sort BUSCO_output/gene_lists/Zophobas_atratus.txt) <(sort BUSCO_output/gene_lists/Zophobas_morio.txt)) <(comm -12 <(comm -12 <(comm -12 <(sort BUSCO_output/gene_lists/Asbolus_verrucosus.txt) <(sort BUSCO_output/gene_lists/Tenebrio_molitor.txt)) <(comm -12 <(sort BUSCO_output/gene_lists/Tenebrio_molitor_v2.txt) <(sort BUSCO_output/gene_lists/Tenebrio_molitor_v3.txt))) <(comm -12 <(comm -12 <(sort BUSCO_output/gene_lists/Tribolium_castaneum.txt) <(sort BUSCO_output/gene_lists/Tribolium_confusum.txt)) <(comm -12 <(sort BUSCO_output/gene_lists/Tribolium_freemani.txt) <(sort BUSCO_output/gene_lists/Tribolium_madens.txt)))) > BUSCO_output/gene_lists/common_genes.txt # Find single-copy genes present in all insects
mkdir BUSCO_output/renamed_genes/ # Make directory
mkdir BUSCO_output/renamed_genes/Asbolus_verrucosus/ # Make directory
mkdir BUSCO_output/renamed_genes/Tenebrio_molitor/ # Make directory
mkdir BUSCO_output/renamed_genes/Tenebrio_molitor_v2/ # Make directory
mkdir BUSCO_output/renamed_genes/Tenebrio_molitor_v3/ # Make directory
mkdir BUSCO_output/renamed_genes/Tribolium_castaneum/ # Make directory
mkdir BUSCO_output/renamed_genes/Tribolium_confusum/ # Make directory
mkdir BUSCO_output/renamed_genes/Tribolium_freemani/ # Make directory
mkdir BUSCO_output/renamed_genes/Tribolium_madens/ # Make directory
mkdir BUSCO_output/renamed_genes/Zophobas_atratus/ # Make directory
mkdir BUSCO_output/renamed_genes/Zophobas_morio/ # Make directory
perl scripts/rename_sequences2.pl Asbolus_verrucosus # Rename gene sequences to have consistent naming
perl scripts/rename_sequences2.pl Tenebrio_molitor # Rename gene sequences to have consistent naming
perl scripts/rename_sequences2.pl Tenebrio_molitor_v2 # Rename gene sequences to have consistent naming
perl scripts/rename_sequences2.pl Tenebrio_molitor_v3 # Rename gene sequences to have consistent naming
perl scripts/rename_sequences2.pl Tribolium_castaneum # Rename gene sequences to have consistent naming
perl scripts/rename_sequences2.pl Tribolium_confusum # Rename gene sequences to have consistent naming
perl scripts/rename_sequences2.pl Tribolium_freemani # Rename gene sequences to have consistent naming
perl scripts/rename_sequences2.pl Tribolium_madens # Rename gene sequences to have consistent naming
perl scripts/rename_sequences2.pl Zophobas_atratus # Rename gene sequences to have consistent naming
perl scripts/rename_sequences2.pl Zophobas_morio # Rename gene sequences to have consistent naming
mkdir BUSCO_output/combined_genes/ # Make directory
perl scripts/concatenate_orthologs.pl BUSCO_output/gene_lists/common_genes.txt # Concatenate orthologs
mkdir BUSCO_output/mafft/ # Make directory
mkdir BUSCO_output/trimal/ # Make directory
perl scripts/align_trim2.pl # Align and trim each set of orthologs
mkdir BUSCO_output/trimal_modified/ # Make directory
perl scripts/modify_trimal2.pl # Modify the trimmed alignment file for concatenation
mkdir BUSCO_output/maximum_likelihood_phylogeny/ # Make directory
grep '>' BUSCO_output/mafft/996662at2759_mafft.fasta | sed 's/>//' > BUSCO_output/species_list.txt # Make a file of species names
perl scripts/combine_alignments2.pl > BUSCO_output/maximum_likelihood_phylogeny/final_alignment.fasta # Combine the trimmed alignments as one super alignment
cd BUSCO_output/maximum_likelihood_phylogeny/ # Change directory
raxmlHPC-HYBRID-AVX2 -T 28 -s final_alignment.fasta -N 5 -n test_phylogeny -f a -p 12345 -x 12345 -m PROTGAMMAAUTO # Determine the best model for the phylogeny
raxmlHPC-HYBRID-AVX2 -T 28 -s final_alignment.fasta -N 100 -n final_phylogeny -f a -p 12345 -x 12345 -m PROTGAMMAJTTF # Create a maximum likelihood phylogeny
cd ../../ # Change directory

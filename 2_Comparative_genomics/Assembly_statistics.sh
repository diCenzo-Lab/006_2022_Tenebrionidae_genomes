# Check assembly completeness
mkdir BUSCO_output
busco --offline --in Genome_sequences/Asbolus_verrucosus.fna --out Asbolus_verrucosus.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 28 -f
busco --offline --in Genome_sequences/Tenebrio_molitor.fna --out Tenebrio_molitor.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 32 -f
busco --offline --in Genome_sequences/Tenebrio_molitor_v2.fna --out Tenebrio_molitor_v2.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 32 -f
busco --offline --in Genome_sequences/Tenebrio_molitor_v3.fna --out Tenebrio_molitor_v3.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 32 -f
busco --offline --in Genome_sequences/Tribolium_castaneum.fna --out Tribolium_castaneum.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 32 -f
busco --offline --in Genome_sequences/Tribolium_confusum.fna --out Tribolium_confusum.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 32 -f
busco --offline --in Genome_sequences/Tribolium_freemani.fna --out Tribolium_freemani.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 32 -f
busco --offline --in Genome_sequences/Tribolium_madens.fna --out Tribolium_madens.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 32 -f
busco --offline --in Genome_sequences/Zophobas_atratus.fna --out Zophobas_atratus.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 32 -f
busco --offline --in Genome_sequences/Zophobas_morio.fna --out Zophobas_morio.endopterygota.busco --lineage_dataset endopterygota_odb10 --mode genome --cpu 32 -f
mv *.endopterygota.busco BUSCO_output/

# Get other statistics
stats.sh -Xmn1G Genome_sequences/Asbolus_verrucosus.fna
stats.sh -Xmn1G Genome_sequences/Tenebrio_molitor.fna
stats.sh -Xmn1G Genome_sequences/Tenebrio_molitor_v2.fna
stats.sh -Xmn1G Genome_sequences/Tenebrio_molitor_v3.fna
stats.sh -Xmn1G Genome_sequences/Tribolium_castaneum.fna
stats.sh -Xmn1G Genome_sequences/Tribolium_confusum.fna
stats.sh -Xmn1G Genome_sequences/Tribolium_freemani.fna
stats.sh -Xmn1G Genome_sequences/Tribolium_madens.fna
stats.sh -Xmn1G Genome_sequences/Zophobas_atratus.fna
stats.sh -Xmn1G Genome_sequences/Zophobas_morio.fna

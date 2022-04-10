## Extract long contigs and sort

# Tenebrio molitor
pullseq -i Genome_sequences/Tenebrio_molitor.fna -m 1000000 > Genome_sequences/Tenebrio_molitor.large.fna
seqkit sort --by-length --reverse Genome_sequences/Tenebrio_molitor.large.fna > Genome_sequences/Tenebrio_molitor.large.sorted.fna

# Tenebrio molitor v3
pullseq -i Genome_sequences/Tenebrio_molitor_v3.fna -m 1000000 > Genome_sequences/Tenebrio_molitor_v3.large.fna
seqkit sort --by-length --reverse Genome_sequences/Tenebrio_molitor_v3.large.fna > Genome_sequences/Tenebrio_molitor_v3.large.sorted.fna

# Tribolium castaneum
pullseq -i Genome_sequences/Tribolium_castaneum.fna -m 1000000 > Genome_sequences/Tribolium_castaneum.large.fna
seqkit sort --by-length --reverse Genome_sequences/Tribolium_castaneum.large.fna > Genome_sequences/Tribolium_castaneum.large.sorted.fna

#Tribolium confusum
pullseq -i Genome_sequences/Tribolium_confusum.fna -m 1000000 > Genome_sequences/Tribolium_confusum.large.fna
seqkit sort --by-length --reverse Genome_sequences/Tribolium_confusum.large.fna > Genome_sequences/Tribolium_confusum.large.sorted.fna

#Tribolium freemani
pullseq -i Genome_sequences/Tribolium_freemani.fna -m 1000000 > Genome_sequences/Tribolium_freemani.large.fna
seqkit sort --by-length --reverse Genome_sequences/Tribolium_freemani.large.fna > Genome_sequences/Tribolium_freemani.large.sorted.fna

# Tribolium madens
pullseq -i Genome_sequences/Tribolium_madens.fna -m 1000000 > Genome_sequences/Tribolium_madens.large.fna
seqkit sort --by-length --reverse Genome_sequences/Tribolium_madens.large.fna > Genome_sequences/Tribolium_madens.large.sorted.fna

# Zophobas atratus
pullseq -i Genome_sequences/Zophobas_atratus.fna -m 1000000 > Genome_sequences/Zophobas_atratus.large.fna
seqkit sort --by-length --reverse Genome_sequences/Zophobas_atratus.large.fna > Genome_sequences/Zophobas_atratus.large.sorted.fna

# Zophobas morio
pullseq -i Genome_sequences/Zophobas_morio.fna -m 1000000 > Genome_sequences/Zophobas_morio.large.fna
seqkit sort --by-length --reverse Genome_sequences/Zophobas_morio.large.fna > Genome_sequences/Zophobas_morio.large.sorted.fna


## Create dot dotplots

# Uploaded files to the D-Genies webserver for dotplot construction.

## Identify SNPs

# Tenebrio molitor new versus Tenebrio molitor v3
/datadisk1/Bioinformatics_programs/mummer/bin/nucmer --minmatch 100 --mincluster 1000 --diagfactor 10 --banded --diagdiff 5 -t 24 -p dotplots/TmolNew_Tmol3 Genome_sequences/Tenebrio_molitor_v3.fna Genome_sequences/Tenebrio_molitor.fna
/datadisk1/Bioinformatics_programs/mummer/bin/mummerplot --fat --filter --small --postscript -p dotplots/TmolNew_Tmol3 dotplots/TmolNew_Tmol3.delta
ps2pdf dotplots/TmolNew_Tmol3.ps dotplots/TmolNew_Tmol3.pdf
/datadisk1/Bioinformatics_programs/mummer/bin/show-snps -Clr dotplots/TmolNew_Tmol3.delta > dotplots/TmolNew_Tmol3.snps

# Zophobas morio new versus Zophobas atratus
/datadisk1/Bioinformatics_programs/mummer/bin/nucmer --minmatch 100 --mincluster 1000 --diagfactor 10 --banded --diagdiff 5 -t 24 -p dotplots/Zmor_Zatr Genome_sequences/Zophobas_atratus.fna Genome_sequences/Zophobas_morio.fna
/datadisk1/Bioinformatics_programs/mummer/bin/mummerplot --fat --filter --small --postscript -p dotplots/Zmor_Zatr dotplots/Zmor_Zatr.delta
ps2pdf dotplots/Zmor_Zatr.ps dotplots/Zmor_Zatr.pdf
/datadisk1/Bioinformatics_programs/mummer/bin/show-snps -Clr dotplots/Zmor_Zatr.delta > dotplots/Zmor_Zatr.snps

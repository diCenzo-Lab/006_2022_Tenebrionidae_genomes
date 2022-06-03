#!usr/bin/perl
use 5.010;

while(<>) {
  chomp;
  $_ =~ s/\.faa//;
  system("cat BUSCO_output/renamed_genes/*/$_.faa > BUSCO_output/combined_genes/$_.fa");
}

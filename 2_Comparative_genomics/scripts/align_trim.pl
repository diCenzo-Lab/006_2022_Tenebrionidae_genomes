#!usr/bin/perl
use File::Find;
use 5.010;
use Cwd;

$pwd = cwd();
$parent = "$pwd/orthofinder_output/Results_Apr16/Single_Copy_Orthologue_Sequences2";
$output1 = "$pwd/orthofinder_output/Results_Apr16/mafft";
$output2 = "$pwd/orthofinder_output/Results_Apr16/trimal";

find( \&search_all_folder, $parent );

sub search_all_folder {
	chomp $_;
	return if $_ eq '.' or $_ eq '..';
	&read_files ($_) if (-f);
	say($_);
}

sub read_files {
	($filename) = @_;
	$filename2 = substr($filename,0,-3);
	$filename3 = $filename2;
	$filename2 .= '_mafft.fasta';
	$filename3 .= '_trimal.fasta';
	system("mafft --localpair --thread 6 --quiet $parent/$filename > $output1/$filename2");
	system("trimal -in $output1/$filename2 -out $output2/$filename3 -fasta -automated1");
}

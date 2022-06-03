#!usr/bin/perl
use File::Find;
use 5.010;
use Cwd;

$insect = @ARGV[0];
$pwd = cwd();
$parent = "$pwd/BUSCO_output/$insect.eukaryota.busco/run_eukaryota_odb10/busco_sequences/single_copy_busco_sequences";
$output = "$pwd/BUSCO_output/renamed_genes/$insect";

find( \&search_all_folder, $parent );

sub search_all_folder {
	chomp $_;
	return if $_ eq '.' or $_ eq '..';
	&read_files ($_) if (-f);
}

sub read_files {
	($filename) = @_;
	open($in, '<', "$parent/$filename");
	open($out, '>', "$output/$filename");
	$count = 0;
	while(<$in>) {
		if(/>/) {
			say $out (">$insect");
		}
		else {
			print $out ("$_");
		}
	}
	close($in);
	close($out);
}

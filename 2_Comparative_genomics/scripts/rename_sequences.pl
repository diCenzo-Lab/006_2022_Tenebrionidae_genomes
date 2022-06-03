#!usr/bin/perl
use File::Find;
use 5.010;
use Cwd;

$pwd = cwd();
$parent = "$pwd/orthofinder_output/Results_Apr16/Single_Copy_Orthologue_Sequences";
$output = "$pwd/orthofinder_output/Results_Apr16/Single_Copy_Orthologue_Sequences2";

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
			$count++;
			if($count == 1) {
				say $out (">Asbolus_verrucosus");
			}
			elsif($count == 2) {
				say $out (">Tenebrio_molitor_new");
			}
			elsif($count == 3) {
				say $out (">Tenebrio_molitor_old");
			}
			elsif($count == 4) {
				say $out (">Tribolium_castaneum");
			}
			elsif($count == 5) {
				say $out (">Tribolium_madens");
			}
			elsif($count == 6) {
				say $out (">Zophobas_morio");
			}
		}
		else {
			print $out ("$_");
		}
	}
	close($in);
	close($out);
}

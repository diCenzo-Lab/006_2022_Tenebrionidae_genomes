#!usr/bin/perl
use 5.010;

while(<>) {
	if(/>/) {
		chomp;
		print("\n$_\t");
	}
	else {
		chomp;
		print("$_");
	}
}

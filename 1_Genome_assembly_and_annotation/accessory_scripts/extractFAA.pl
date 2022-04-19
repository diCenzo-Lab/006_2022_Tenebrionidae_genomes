#!usr/bin/perl
use 5.010;

while(<>) {
  chomp;
  if(/protein_id/) {
    $_ =~ s/\"//g;
    $_ =~ s/\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \/protein_id\=//;
    say(">$_");
  }
  elsif(/translation/) {
    $_ =~ s/\"//g;
    $_ =~ s/\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \/translation\=//;
    say("$_");
    $test = 1;
  }
  elsif($test == 1 && /\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ /) {
    $_ =~ s/\"//g;
    $_ =~ s/\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ //;
    say("$_");
  }
  else {
    $test = 0;
  }
}

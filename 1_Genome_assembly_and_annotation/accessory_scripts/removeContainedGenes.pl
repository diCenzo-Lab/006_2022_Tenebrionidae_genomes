#!usr/bin/perl
use 5.010;

$scaffold_old = 0;
$stard_old = 0;
$end_old = 0;
$direction_old = 0;

while(<>) {
  if(/ID=gene/) {
    @line = split("\t", $_);
    $scaffold = @line[0];
    $start = @line[3];
    $end = @line[4];
    $direction = @line[6];
    if($scaffold eq $scaffold_old && $direction eq $direction_old && $start < $end_old) {
      print($_);
    }
    $scaffold_old = $scaffold;
    $stard_old = $start;
    $end_old = $end;
    $direction_old = $direction;
  }
}

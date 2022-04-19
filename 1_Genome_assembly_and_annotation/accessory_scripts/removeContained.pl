#!usr/bin/perl
use 5.010;

while(<>) {
  @line = split("\t", $_);
  if(/###/) {
    @cds_start = [0];
    @cds_end = [0];
    @exon_start = [0];
    @exon_end = [0];
    print($_);
  }
  elsif(@line[2] eq 'CDS') {
    $test = 0;
    $n = -1;
    for $i (@cds_start) {
      $n++;
      if(@line[3] >= $i) {
        if(@line[4] <= @cds_end[$n]) {
          $test = 1;
        }
      }
    }
    if($test == 0) {
      print($_);
    }
    push(@cds_start, @line[3]);
    push(@cds_end, @line[4]);
  }
  elsif(@line[2] eq 'exon') {
    $test = 0;
    $n = -1;
    for $i (@exon_start) {
      $n++;
      if(@line[3] >= $i) {
        if(@line[4] <= @exon_end[$n]) {
          $test = 1;
        }
      }
    }
    if($test == 0) {
      print($_);
    }
    push(@exon_start, @line[3]);
    push(@exon_end, @line[4]);
  }
  else {
    print($_);
  }
}

$count = 0;
while(<>) {
  if(/>/) {
    $count = $count + 1;
    if(/mitochondrion/) {
      print(">contig_$count [mitochondrion]\n");
    }
    else {
      print(">contig_$count\n");
    }
  }
  else {
    print("$_");
  }
}

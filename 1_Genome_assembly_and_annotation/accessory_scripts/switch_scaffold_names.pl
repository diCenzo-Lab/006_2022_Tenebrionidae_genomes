$count = 0;
while(<>) {
  if(/>/) {
    $count = $count + 1;
    if(/mitochondrion/) {
      print(">scaffold_$count [mitochondrion]\n");
    }
    else {
      print(">scaffold_$count\n");
    }
  }
  else {
    print("$_");
  }
}

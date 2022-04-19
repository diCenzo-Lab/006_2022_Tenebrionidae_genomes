$test = 0;
while(<>) {
  if(/>/) {
    if(/mitochondrion/) {
      $test = 1;
    }
    else {
      $test = 0;
      print($_);
    }
  }
  else {
    if($test == 0) {
      print($_);
    }
  }
}

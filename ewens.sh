for i in {0..10}; do bob=$(( $i*10 + 1 )); perl ewens.pl 12 $bob; cat ewens.out | awk '{ SUM += $1} END { print SUM  }'; done
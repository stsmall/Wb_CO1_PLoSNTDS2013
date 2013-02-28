#!/usr/bin/perl -w
use strict;
use lib "/export/apps/local/lib/perl/ /opt/rocks/lib/perl5/5.8.8/";

use Math::BigFloat;

## main code
my $theta = $ARGV[0];
my $popsize = $ARGV[1];
my %stirling1_hash=(); # hash for caclulating stirling numbers of the first kind
open OUT, ">ewens.out";
&ewens;
close OUT;

##"ewens.pl 5 10"

##subroutines
sub kay{
	my $k=1;
	for my $qn (1..$popsize-1){
		$k+=$theta/($theta+$qn);
	}
	$k=int($k+0.5);
	return $k;
}

sub ratio_factorial{
	my $fact = Math::BigFloat -> new(1);
	my $a=shift();
	my $b=shift();
	if ($a==0){$fact=1; return $fact;}
	for my $i ($a-$b+1..$a){
		$fact*=$i;
	}
	return $fact;
}

sub stirling1 {
	my $n = int(shift());
	my $m = int(shift());
	$stirling1_hash{$n,$m} = ($n < $m) ? 0: ($n == $m) ? 1: ($m == 1) ? Math::BigFloat -> new(-($n-1)*&stirling1($n-1,$m)): Math::BigFloat -> new(&stirling1($n-1,$m-1)-($n-1)*&stirling1($n-1,$m)) unless exists($stirling1_hash{$n,$m});
  return $stirling1_hash{$n,$m};
}

sub ewens { # sets up the ewens sampling distribution, draws alleles from it according to $popsize and $theta.i
	my $N=$popsize;
	my $k=&kay;
	my $allele=0; my $freq=0;
	my @ewens=(); $ewens[0]=0;

	for my $j (1..$N){
		my $a = &stirling1($N-$j,$k-1);
		my $b = &stirling1($N,$k);
		my $c = (&ratio_factorial($N,$j))/$j;
		my $z = abs ($a * $c / $b);
		$ewens[$j]=$z+$ewens[$j-1];
		print OUT "$z\n";
	}
}


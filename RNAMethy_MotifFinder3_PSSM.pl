#!/usr/bin/perl -w
use strict;
use warnings;

if(@ARGV != 2) {
    print STDERR "Usage: RNAMethy_MotifFinder3_PSSM.pl inf(fasta) outf\n";
    exit(0);
}

my ($inf, $outf)= @ARGV;

my %hash;
my $seed;
my @seed;
my $min_score;
open(IN, $inf) or die "Cannot open $inf\n";
while(<IN>){
	if($_=~/^>Seed/i){
		$seed=<IN>;
		chomp($seed);
		@seed=split(//, $seed);
		for(my $i=10; $i<=10+$#seed; $i++){
			$hash{$i}{$seed[$i-10]}++;
			$hash{$i}{"total"}++;
		}
	}elsif($_=~/^>\w+ (\d+)/){
		$min_score=$1;
		my $seq=<IN>;
		chomp($seq);
		my ($i, $j) = &get_align_score($seq, $seed);
		my @seq=split(//,$seq);
		if($i==0){
			for(my $k=0; $k<=$#seq; $k++){
				my $index=10+$j+$k;
				$hash{$index}{$seq[$k]}++;
				$hash{$index}{"total"}++;
			}
		}elsif($j==0){
			for(my $k=0; $k<=$#seq; $k++){
				my $index=10-$i+$k;
				$hash{$index}{$seq[$k]}++;
				$hash{$index}{"total"}++;
			}
		}
	}
}
close IN;

my %print;
my %print_count;

my @rna=("A", "C", "G", "U");
my @index= sort{$a<=>$b} keys %hash;
#foreach my $index (sort{$a<=>$b} keys %hash){
for(my $i=0; $i<=$#index; $i++){
	$print{"index"}.="\t$index[$i]";
	$print_count{"index"}.="\t$index[$i]";
	my $total_ratio=0;
	foreach my $rna (@rna){
		if(!exists $hash{$index[$i]}{$rna}){
			$hash{$index[$i]}{$rna}=0;
		}
		$print_count{$rna}.="\t$hash{$index[$i]}{$rna}";
		my $ratio=sprintf("%.2f", $hash{$index[$i]}{$rna}/$hash{$index[$i]}{"total"});
		$total_ratio+=$ratio;
		if($total_ratio>1){
			$total_ratio=$total_ratio-$ratio;
			$ratio=1-$total_ratio;
			$total_ratio+=$ratio;
		}
		$print{$rna}.="\t$ratio";
	}
#	$print{"total"}.="\t".$hash{$index[$i]}{"total"};
}

open(OUT, ">$outf") or die "Cannot open $outf\n";
print OUT "Position".$print{"index"}."\n";
print OUT "A".$print{"A"}."\n";
print OUT "C".$print{"C"}."\n";
print OUT "G".$print{"G"}."\n";
print OUT "U".$print{"U"}."\n";
#print OUT "Total".$print{"total"}."\n";
close OUT;

my $outf2=$outf.".count";
open(OUT, ">$outf2") or die "Cannot open $outf\n";
print OUT "Position".$print_count{"index"}."\n";
print OUT "A".$print_count{"A"}."\n";
print OUT "C".$print_count{"C"}."\n";
print OUT "G".$print_count{"G"}."\n";
print OUT "U".$print_count{"U"}."\n";
#print OUT "Total".$print{"total"}."\n";
close OUT;

sub get_align_score{
	my ($k1, $k2)=@_;
#	my $score=100;
	if(length($k1)<length($k2)){
		($k1, $k2)= ($k2, $k1);
	}
	my @k1=split(//, $k1);
	my @k2=split(//, $k2);
	my @matrix;
	for(my $i=0; $i<=$#k1; $i++){
		for(my $j=0; $j<=$#k2; $j++){
			if($k1[$i] eq $k2[$j]){
				$matrix[$i][$j]=0;
			}else{
				$matrix[$i][$j]=1;
			}
		}
	}

	for(my $j=0; $j<=$#k2; $j++){
		my $tmp_score;
		for(my $k=0; $k+$j<=$#k2; $k++){
			$tmp_score+=$matrix[$k][$j+$k];
		}
		$tmp_score=$tmp_score+$j;
		if($tmp_score==$min_score){
			return (0, $j);
		}
	}

	for(my $j=1; $j<=$#k1; $j++){
		my $tmp_score;
		for(my $k=0; $k+$j<=$#k1; $k++){
			if($k>$#k2){
				last;
			}
			$tmp_score+=$matrix[$j+$k][$k];
		}
		if($#k2>$#k1-$j){
			$tmp_score=$tmp_score + $#k2 - ($#k1-$j);
		}
#		$tmp_score=$tmp_score+$j;
		if($tmp_score==$min_score){
			return ($j, 0);
		}
	}
}






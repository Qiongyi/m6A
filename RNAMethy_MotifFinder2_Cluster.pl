#!/usr/bin/perl -w
use strict;
use warnings;

if(@ARGV != 5) {
    print STDERR "Usage: RNAMethy_MotifFinder2_Cluster.pl k_cutoff(eg. 1) inf1 inf2 inf3 outf\n";
    exit(0);
}

my ($k_cutoff, $inf1, $inf2, $inf3, $outf)= @ARGV;

my @kmer; # record kmer
my %kmer; # record kmer

my %score_column; # add scores for each kmer column
my %score; # record all scores for each kmer pair;
open(IN, $inf1) or die "Cannot open $inf1\n";
while(<IN>){
	if($_=~/^Kmer/){
		next;
	}elsif($_=~/\w+/){
		my @info=split(/\t/, $_);
		push(@kmer, $info[0]);
		$kmer{$info[0]}=1;
	}
}
close IN;

open(IN, $inf2) or die "Cannot open $inf2\n";
while(<IN>){
	if($_=~/^Kmer/){
		next;
	}elsif($_=~/\w+/){
		my @info=split(/\t/, $_);
		push(@kmer, $info[0]);
		$kmer{$info[0]}=1;
	}
}
close IN;

open(IN, $inf3) or die "Cannot open $inf3\n";
while(<IN>){
	if($_=~/^Kmer/){
		next;
	}elsif($_=~/\w+/){
		my @info=split(/\t/, $_);
		push(@kmer, $info[0]);
		$kmer{$info[0]}=1;
	}
}
close IN;

open(OUT, ">$outf") or die "Cannot open $outf\n";
my $seed_count=1;
while(@kmer>=5){
	my $out_fasta="Kmer_sequences_seed$seed_count.fa";
	open(OUT1, ">$out_fasta") or die "Cannot open $out_fasta\n";
	my $tmp=scalar(@kmer);
	print STDERR "$tmp k-mer left for seed $seed_count\n";
	
### start analysis from here
for(my $i=0; $i<=$#kmer; $i++){
	for(my $j=0; $j<=$#kmer; $j++){
		my $score=0;
		if($kmer[$i] ne $kmer[$j]){
			$score = &get_align_score($kmer[$i], $kmer[$j]);
#			print OUT "$kmer[$i]\t$kmer[$j]\t$score\n";
		}
		$score_column{$kmer[$i]} += $score;
		$score{$kmer[$i]}{$kmer[$j]}=$score;
	}
}


# find the seed

my $seed_kmer;
foreach my $kmer (sort {$score_column{$a}<=>$score_column{$b}} keys %score_column){
	print OUT "$kmer\t$score_column{$kmer}\tSeed $seed_count\n";
	print OUT1 ">Seed$seed_count\n$kmer\n";
	$seed_count++;
	$seed_kmer=$kmer;
	delete $kmer{$kmer};
	last;
}

# extract all kmers whose distance is within k from the seed
my $count_tmp=0;
foreach my $kmer (keys %{$score{$seed_kmer}}){
	if($score{$seed_kmer}{$kmer}<=$k_cutoff){
		if($kmer ne $seed_kmer){
			$count_tmp++;
			print OUT "$kmer\t$score{$seed_kmer}{$kmer}\t \n";
			print OUT1 ">Kmer_$count_tmp $score{$seed_kmer}{$kmer}\n$kmer\n";
			delete $kmer{$kmer};
		}
	}
}
print OUT "\n";
close OUT1;

@kmer = keys(%kmer);
%score_column=();
%score=();
#print STDERR "$tmp k-mer left\n"; sleep 1;

}
close OUT;


sub get_align_score{
	my ($k1, $k2)=@_;
	my $score=100;
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
		if($tmp_score<$score){
			$score=$tmp_score;
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
		if($tmp_score<$score){
			$score=$tmp_score;
		}
	}

	return $score;
}






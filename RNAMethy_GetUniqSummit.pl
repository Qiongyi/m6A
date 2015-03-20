#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 2) {
    print STDERR "Usage: RNAMethy_GetUniqSummit.pl Naive_SummitPosition.xls Naive_UniqSummitPosition.xls\n";
    exit(0);
}
my ($inf, $outf)= @ARGV;

my $count_all;
my $count_uniq;
my %hash;
open(IN, $inf) or die "Cannot open $inf\n";
open(OUT, ">$outf") or die "Cannot open $outf\n";
while(<IN>){
	if($_=~/^ID/){
		print OUT $_;
		next;
	}
	my @info=split(/\t/,$_);
	my $tag="$info[1]\t$info[3]\t$info[7]\t$info[8]\t$info[9]\t$info[10]\t$info[11]\t$info[12]";
	if(exists $hash{$tag}){
		$count_all++;
		next;
	}else{
		print OUT $_;
		$hash{$tag}=1;
		$count_uniq++;
		$count_all++;
	}
}
close IN;
close OUT;

print STDERR "Number of all summit position in $inf: $count_all\n";
print STDERR "Number of unique summit position in $inf: $count_uniq\n";

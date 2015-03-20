#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 4) {
    print STDERR "Usage: RNAMethy_CombineUniqSummit.pl Naive_UniqSummitPosition.xls Context_UniqSummitPosition.xls FC_UniqSummitPosition.xls Naive_Context_FC_SummitPosition.xls\n";
    exit(0);
}
my ($inf1, $inf2, $inf3, $outf)= @ARGV;

my %hash;
open(OUT, ">$outf") or die "Cannot open $outf\n";
open(IN, $inf1) or die "Cannot open $inf1\n";
while(<IN>){
	if($_=~/^ID/){
		print OUT $_;
		next;
	}else{
		my @info=split(/\t/,$_);
		$hash{$info[1]}{$info[3]}=$_;
	}
}
close IN;

open(IN, $inf2) or die "Cannot open $inf2\n";
while(<IN>){
	if($_=~/^ID/){
#		print OUT $_;
		next;
	}else{
		my @info=split(/\t/,$_);
		$hash{$info[1]}{$info[3]}=$_;
	}
}
close IN;

open(IN, $inf3) or die "Cannot open $inf3\n";
while(<IN>){
	if($_=~/^ID/){
#		print OUT $_;
		next;
	}else{
		my @info=split(/\t/,$_);
		$hash{$info[1]}{$info[3]}=$_;
	}
}
close IN;

foreach my $chr (sort keys %hash){
	foreach my $pos (sort{$a<=>$b} keys %{$hash{$chr}}){
		print OUT "$hash{$chr}{$pos}";
	}
}
close OUT;


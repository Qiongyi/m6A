#!/usr/bin/perl -w
use strict;
use warnings;
use Statistics::R;

if(@ARGV != 2) {
    print STDERR "Usage: RNAMethy_Summit_TTest.pl Naive_Context_FC_UniqSummitPosition.xls Naive_Context_FC.ttest.xls\n";
    exit(0);
}
my ($inf, $outf)= @ARGV;

my ($c3, $c8, $c13, $c15, $c18, $c20)=(2501832, 5752214, 4153178, 4023330, 4491010, 4345696);
my $normalzation_factor=2000000;  ## normalized to 1M fragment (2M reads)
my $total_num=0;
# Create a communication bridge with R and start R
my $R = Statistics::R->new();
my %hash;
my %com1;
my %com2;
my %com3;
open(IN, $inf) or die "Cannot open $inf\n";
while(<IN>){
	chomp;
	if($_=~/^Peak_ID/){
		next;
	}else{
		my @info=split(/\t/,$_);
		my $c3_norm=$info[9]/$c3*$normalzation_factor;
		my $c8_norm=$info[10]/$c8*$normalzation_factor;
		my $c13_norm=$info[11]/$c13*$normalzation_factor;
		my $c15_norm=$info[12]/$c15*$normalzation_factor;
		my $c18_norm=$info[13]/$c18*$normalzation_factor;
		my $c20_norm=$info[14]/$c20*$normalzation_factor;		
=no need to set interger for t-test
		my $c3_norm=int($info[9]/$c3*$normalzation_factor+0.5);
		my $c8_norm=int($info[10]/$c8*$normalzation_factor+0.5);
		my $c13_norm=int($info[11]/$c3*$normalzation_factor+0.5);
		my $c15_norm=int($info[12]/$c3*$normalzation_factor+0.5);
		my $c18_norm=int($info[13]/$c3*$normalzation_factor+0.5);
		my $c20_norm=int($info[14]/$c3*$normalzation_factor+0.5);
=cut
		my $list1="$c3_norm,$c8_norm";
		my $list2="$c13_norm,$c15_norm";
		my $list3="$c18_norm,$c20_norm";

		$hash{$info[0]}="$c3_norm\t$c8_norm\t$c13_norm\t$c15_norm\t$c18_norm\t$c20_norm";

		my $p_value;
  		### Run R commands
  		# Naive vs. Context
  		if($c3_norm == $c8_norm && $c8_norm==$c13_norm && $c13_norm == $c15_norm){
  			$p_value=1;
  		}elsif($c3_norm == $c8_norm && $c13_norm == $c15_norm){
  			$c15_norm=$c15_norm+0.01;
  			$list2="$c13_norm,$c15_norm";
  			$R->run(qq`x <- t.test(c($list1), c($list2))`);
  			$p_value= $R -> get('x$p.value');
  		}else{
  			$R->run(qq`x <- t.test(c($list1), c($list2))`);
  			$p_value= $R -> get('x$p.value');
  		}
  		$hash{$info[0]}.="\t$p_value";
  		$com1{$info[0]}=$p_value;

  		# Naive vs. FC
  		if($c3_norm == $c8_norm && $c8_norm==$c18_norm && $c18_norm == $c20_norm){
  			$p_value=1;
  		}elsif($c3_norm == $c8_norm && $c18_norm == $c20_norm){
  			$c20_norm=$c20_norm+0.01;
  			$list3="$c18_norm,$c20_norm";
  			$R->run(qq`x <- t.test(c($list1), c($list3))`);
  			$p_value= $R -> get('x$p.value');
  		}else{
  			$R->run(qq`x <- t.test(c($list1), c($list3))`);
  			$p_value= $R -> get('x$p.value');
  		}
  		$hash{$info[0]}.="\t$p_value";
  		$com2{$info[0]}=$p_value;

  		# Context vs. FC
  		if($c13_norm == $c15_norm && $c15_norm==$c18_norm && $c18_norm == $c20_norm){
  			$p_value=1;
  		}elsif($c13_norm == $c15_norm && $c18_norm == $c20_norm){
  			$c20_norm=$c20_norm+0.01;
  			$list3="$c18_norm,$c20_norm";
  			$R->run(qq`x <- t.test(c($list2), c($list3))`);
  			$p_value= $R -> get('x$p.value');
  		}else{
  			$R->run(qq`x <- t.test(c($list2), c($list3))`);
  			$p_value= $R -> get('x$p.value');
  		}
  		$hash{$info[0]}.="\t$p_value";
  		$com3{$info[0]}=$p_value;
  		$total_num++;
  	}
}
close IN;

# False discovery rate (FDR) control is a statistical method used in multiple hypothesis testing to correct for multiple comparisons.
my $fdr;
my $order=0;
foreach my $peak (sort {$com1{$a}<=>$com1{$b}} keys %com1){
	$order++;
	$fdr = $com1{$peak}*$total_num/$order;
	if($fdr>1){
		$fdr=1;
	}
	$hash{$peak}.="\t$fdr";
}

$order=0;
foreach my $peak (sort {$com2{$a}<=>$com2{$b}} keys %com2){
	$order++;
	$fdr = $com2{$peak}*$total_num/$order;
	if($fdr>1){
		$fdr=1;
	}
	$hash{$peak}.="\t$fdr";
}

$order=0;
foreach my $peak (sort {$com3{$a}<=>$com3{$b}} keys %com3){
	$order++;
	$fdr = $com3{$peak}*$total_num/$order;
	if($fdr>1){
		$fdr=1;
	}
	$hash{$peak}.="\t$fdr";
}

open(IN, $inf) or die "Cannot open $inf\n";
open(OUT, ">$outf") or die "Cannot open $outf\n";
while(<IN>){
	chomp;
	if($_=~/^Peak_ID/){
		print OUT "$_\tdoc_norm_C3\tdoc_norm_C8\tdoc_norm_C13\tdoc_norm_C15\tdoc_norm_C18\tdoc_norm_C20\tP-value(Naive.vs.Context)\tP-value(Naive.vs.FC)\tP-value(Context.vs.FC)\tFDR(Naive.vs.Context)\tFDR(Naive.vs.FC)\tFDR(Context.vs.FC)\n";
	}else{
		my @info=split(/\t/,$_);
		if(!defined($hash{$info[0]})){
			print STDERR "Not defined for $info[0] hash value. Please check!\n";
		}
		print OUT "$_\t$hash{$info[0]}\n";
	}
}
close IN;
close OUT;


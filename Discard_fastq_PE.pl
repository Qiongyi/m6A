#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 6) {
    print STDERR "Usage: Discard_fastq_PE.pl sam inf1 inf2 outf1 outf2 outf_stats\n";
    exit(0);
}

my ($sam, $inf1, $inf2, $outf1, $outf2, $outf3)= @ARGV;

my $sampleID;
if($inf1=~/([^\/]+)_read1\.fastq\.gz/){
	$sampleID=$1;
}

my %hash;
open(IN, $sam) or die "Cannot open $!\n"; 
while(<IN>){
	my @info=split(/\t/,$_);
	$hash{$info[0]}=undef;
}
close IN;

open(OUT1, ">$outf1") or die $!;
open(OUT2, ">$outf2") or die $!;
if(-e $outf3){
	open(OUT3, ">>$outf3") or die $!;
}else{
	open(OUT3, ">$outf3") or die $!;
	print OUT3 "Sample ID\tTotal read pairs\tRead pairs after filter rRNA and PhiX\t%\n";
}

if ($inf1 =~ /\.gz$/) {
	open(IN1, "gunzip -c $inf1 |") || die "can't open pipe to $inf1";
}else {
	open(IN1, $inf1) || die "can't open $inf1";
}

if ($inf2 =~ /\.gz$/) {
	open(IN2, "gunzip -c $inf2 |") || die "can't open pipe to $inf2";
}else {
	open(IN2, $inf2) || die "can't open $inf2";
}

my $total=0;
my $filter=0;
while(my $fq1_line1=<IN1>){
	my $fq1_line2=<IN1>;
	my $fq1_line3=<IN1>;
	my $fq1_line4=<IN1>;

	my $fq2_line1=<IN2>;
	my $fq2_line2=<IN2>;
	my $fq2_line3=<IN2>;
	my $fq2_line4=<IN2>;

	$total++;

	if($fq1_line1=~/^\@([^\s]+)/){
		if(!exists($hash{$1})){
			my $name=$1;
			if($fq2_line1=~/^\@([^\s]+)/){
				if($name ne $1){
					print STDERR "read1 and read2 are not paired for $inf1 and $inf2!\n$fq1_line1$fq2_line1"; exit;
				}
			}
			if(length($fq1_line2) != length($fq1_line4)){
				print STDERR "The lengths of sequence and quality are not matched for $inf1!\n$fq1_line1$fq1_line2$fq1_line3$fq1_line4"; exit;
			}
			if(length($fq2_line2) != length($fq2_line4)){
				print STDERR "The lengths of sequence and quality are not matched for $inf2!\n$fq2_line1$fq2_line2$fq2_line3$fq2_line4"; exit;
			}
			if(length($fq1_line2)>=12 && length($fq2_line2)>=12){			
				$filter++;
				print OUT1 $fq1_line1.$fq1_line2.$fq1_line3.$fq1_line4;
				print OUT2 $fq2_line1.$fq2_line2.$fq2_line3.$fq2_line4;
			}
		}
	}else{
		print STDERR "Wrong line for $inf1:\n$fq1_line1"; exit;
	}
}
close IN1;
close IN2;
close OUT1;
close OUT2;

print OUT3 "$sampleID\t$total\t$filter\t";
printf OUT3 ("%.2f%%\n", $filter/$total*100);
close OUT3

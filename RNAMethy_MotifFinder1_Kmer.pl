#!/usr/bin/perl -w
use strict;
use warnings;
use Statistics::R;

if(@ARGV != 6) {
    print STDERR "Usage: RNAMethy_MotifFinder1_Kmer.pl len_kmer len_around_summit(eg. 25 for a 51bp window centred the summit) transcripts_seq Naive_Context_FC.ttest.GenomicLocation.xls Naive_Context_FC_control.GenomicLocation.xls m6A_motifs.xls\n";
    exit(0);
}

my ($len_kmer, $len, $seq, $peaks, $control, $outf)= @ARGV;

my %hash; # store all transcritp sequences
open(IN, $seq) or die "Cannot open $seq\n";
while(<IN>){
	if($_=~/^>(.+)/){
		my $id=$1;
		my $seq=<IN>;
		chomp($seq);
		$hash{$id}=$seq;
	}
}
close IN;

### read control
my %control;
my $control_total;
open(IN, $control) or die "Cannot open $control\n";
while(<IN>){
	chomp;
	my @info=split(/\t/,$_);
	if($info[0]=~/Peak_ID/){
		next;
	}else{
		my $id= "$info[1] $info[2]";
		if(exists $hash{$id}){
			my $centre= $info[6]+50-1;
			if($info[4] eq "-"){
				$centre= length($hash{$id})-$info[6]+1;
			}
			my $core_start=$centre-$len>0?($centre-$len):0;
			my $extract_len=$len*2+1;
			if($centre-$len<0){
				$extract_len=$extract_len+$centre-$len;
			}
			my $core_seq=uc(substr($hash{$id}, $core_start, $extract_len));
			$core_seq=~tr/T/U/;
			my $core_seq_len=length($core_seq);
			for(my $i=0; $i<$core_seq_len-$len_kmer+1; $i++){
				my $kmer=substr($core_seq, $i, $len_kmer);
				$control{$kmer}++;
				$control_total++;
			}				
		}
	}
}
close IN;

### read m6A peaks
my %case;
my $case_total;
open(IN, $peaks) or die "Cannot open $peaks\n";
while(<IN>){
		chomp;
		my @info=split(/\t/,$_);
		my $id= "$info[1] $info[2]";
		if(exists $hash{$id}){
			my $centre= $info[6]-1;
			if($info[4] eq "-"){
				$centre= length($hash{$id})-$info[6]+1;
			}			
			my $core_start=$centre-$len>0?($centre-$len):0;
			my $extract_len=$len*2+1;
			if($centre-$len<0){
				$extract_len=$extract_len+$centre-$len;
			}
			my $core_seq=uc(substr($hash{$id}, $core_start, $extract_len));
			$core_seq =~tr/T/U/;
			my $core_seq_len=length($core_seq);

			for(my $i=0; $i<$core_seq_len-$len_kmer+1; $i++){
				my $kmer=substr($core_seq, $i, $len_kmer);
				$case{$kmer}++;
				$case_total++;
			}
		}
}
close IN;

my $case_kmer=scalar(keys %case);
my $control_kmer=scalar(keys %control);

print STDERR "Total number of case kmer ($len_kmer): $case_kmer\n";
print STDERR "Total number of control kmer ($len_kmer): $control_kmer\n";
print STDERR "Total number of case kmer counts ($len_kmer): $case_total\n";
print STDERR "Total number of control kmer counts ($len_kmer): $control_total\n";

# Create a communication bridge with R and start R
my $R = Statistics::R->new();

open(OUT, ">$outf") or die "Cannot open $outf\n";
print OUT "Kmer\tP-value(Fisher's Exact Test)\tBonferonni-corrected P-value\tcase_prevalence\tcontrol_prevalence\tNumbers_for_FisherExactTest\n";
foreach my $kmer (keys %case){
	if(!exists $control{$kmer}){
		$control{$kmer}=0;
	}
	my $fisher_test=($case_total-$case{$kmer}).",".$case{$kmer}.",".($control_total-$control{$kmer}).",".$control{$kmer};
	$R->run(
qq`m<-matrix(c($fisher_test),2,2,byrow=T)`,
qq`x <- fisher.test(m)`,
);
	my $p_value= $R -> get('x$p.value');
	my $q_value= $p_value*$case_kmer>1?1:$p_value*$case_kmer;
	print OUT "$kmer\t$p_value\t$q_value\t".($case{$kmer}/$case_total)."\t".($control{$kmer}/$control_total)."\t"."$fisher_test\n";
}
close OUT;





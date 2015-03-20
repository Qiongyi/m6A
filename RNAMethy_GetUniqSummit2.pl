#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 4) {
    print STDERR "Usage: RNAMethy_GetUniqSummit2.pl cutoff Naive_Context_FC_SummitPosition.xls Naive_Context_FC_SummitPosition.group.xls Naive_Context_FC_UniqSummitPosition.xls\n";
    exit(0);
}
my ($cutoff, $inf, $group, $outf)= @ARGV;

my %hash;
my $chr="chr_0";
my $pos="-10000";
my @cluster;
my $count_peak=0;

my $total_count=0;
my $total_cluster=1;

my $deseq = "for_DESeq.txt";
my $stats = "peak_distance_stats.xls";

open(IN, $inf) or die "Cannot open $inf\n";
open(DE, ">$deseq") or die "Cannot open $deseq\n";
open(STAT, ">$stats") or die "Cannot open $stats\n";
open(OUT, ">$outf") or die "Cannot open $outf\n";
open(OUT2, ">$group") or die "Cannot open $group\n";
while(<IN>){
	if($_=~/^ID/){
		my @info=split(/\t/,$_);
		my $tmp=join"\t", @info[1..$#info];
		print OUT "Peak_ID\tgene_id\ttranscript_id\t$tmp";
#		print OUT $_;
		print OUT2 $_;
		next;
	}
	$total_count++;
	my @info=split(/\t/,$_);
	if($info[1] eq $chr && $info[3]-$pos<=$cutoff){
		print STAT ($info[3]-$pos)."\n";
		push(@cluster, $_);
		$chr=$info[1];
		$pos=$info[3];
	}else{
		if(@cluster>=1){
			my $tmp=join"",@cluster;
			print OUT2 "\#Group $total_cluster\n$tmp";
			&remove_redun(\@cluster);
			$total_cluster++;
		}
		@cluster=();
		push(@cluster, $_);
		$chr=$info[1];
		$pos=$info[3];
	}
}
close IN;

&remove_redun(\@cluster);
my $tmp=join"",@cluster;
print OUT2 "\#Group $total_cluster\n$tmp";
#$total_cluster++;

close OUT;
close OUT2;
close DE;
close STAT;

print STDERR "Number of all summit position in $inf: $total_count\n";
print STDERR "Number of unique summit position (after grouping) in $inf: $total_cluster\n";

sub remove_redun{
	my $cluster=shift;
	my @cluster=@{$cluster};
	if(@cluster == 1){
		$count_peak++;
		my @tmp=split(/\t/,$cluster[0]);
		my $tmp=join"\t", @tmp[1..$#tmp];
		if($tmp[0]=~/gene_id "([^"]+)"; transcript_id "([^"]+)";/){			
			print OUT "Peak_$count_peak\t$1\t$2\t$tmp";
		}
#		print OUT $cluster[0];
		chomp($cluster[0]);
		my @ele=split(/\t/,$cluster[0]);
		print DE "$ele[1]"."_"."$ele[3]\t$ele[7]\t$ele[8]\t$ele[9]\t$ele[10]\t$ele[11]\t$ele[12]\n";
	}elsif(@cluster>1){
		my $tag=0;
		chomp($cluster[0]);
		my @info=split(/\t/,$cluster[0]);
		my $total= $info[7]+$info[8]+$info[9]+$info[10]+$info[11]+$info[12];
		my $max_line=$cluster[0];
		for(my $i=1; $i<=$#cluster; $i++){
			chomp($cluster[$i]);
			@info=split(/\t/,$cluster[$i]);
				if(($info[7]+$info[8]+$info[9]+$info[10]+$info[11]+$info[12])>=$total){
					$max_line=$cluster[$i];
					$total= $info[7]+$info[8]+$info[9]+$info[10]+$info[11]+$info[12];
				}else{
					$count_peak++;
					my @tmp=split(/\t/,$max_line);
					my $tmp=join"\t", @tmp[1..$#tmp];
					if($tmp[0]=~/gene_id "([^"]+)"; transcript_id "([^"]+)";/){
						print OUT "Peak_$count_peak\t$1\t$2\t$tmp\n";
						$tag=1;
						my @ele=split(/\t/,$max_line);
						print DE "$ele[1]"."_"."$ele[3]\t$ele[7]\t$ele[8]\t$ele[9]\t$ele[10]\t$ele[11]\t$ele[12]\n";
						last;
					}
				}
		}
		if($tag==0){
			$count_peak++;
			my @tmp=split(/\t/,$max_line);
			my $tmp=join"\t", @tmp[1..$#tmp];
			if($tmp[0]=~/gene_id "([^"]+)"; transcript_id "([^"]+)";/){			
				print OUT "Peak_$count_peak\t$1\t$2\t$tmp\n";
				my @ele=split(/\t/,$max_line);
				print DE "$ele[1]"."_"."$ele[3]\t$ele[7]\t$ele[8]\t$ele[9]\t$ele[10]\t$ele[11]\t$ele[12]\n";
			}
		}
	}
}

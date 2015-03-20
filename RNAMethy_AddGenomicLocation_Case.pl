#!/usr/bin/perl -w
use strict;
use warnings;
use Statistics::R;

if(@ARGV != 3) {
	print STDERR "Usage: RNAMethy_AddGenomicLocation_Case.pl GTF Naive_UniqSummitPosition.xls Naive.final.GenomicLocation.xls\n";
    exit(0);
}
my ($gtf, $inf, $outf)=@ARGV;
my $cutoff=200;  # 200bp after TSS is treated as "TSS"; 200bp leftside and 200bp rightside from the stop codon is treated as "Stop";

my %exon;
my %cds;

open(IN, $gtf) or die "cannot open $gtf\n";
while(<IN>){
	my @info=split(/\t/,$_);
	if($info[8]=~/(gene_id \"[^"]+\"; transcript_id \"[^"]+\";)/){
		if($info[2] eq "exon"){
			$exon{$1}.=$_;
		}elsif($info[2] eq "CDS"){
			$cds{$1}.=$_;
		}
	}
}
close IN;

my %cds_start;
my %cds_end;
my %exon_start;
my %exon_end;

foreach my $gene (keys %exon){
	if(exists $cds{$gene}){
		my @t_pos;
		my @line=split(/\n/, $exon{$gene});
		foreach my $line (@line){
			my @info=split(/\t/,$line);
			for(my $i=$info[3]; $i<=$info[4]; $i++){
				push(@t_pos, $i);
			}
		}
#		my @exon1=split(/\t/, $line[0]);
#		my @exon2=split(/\t/, $line[$#line]);
		my ($exon_start, $exon_end)=($t_pos[0], $t_pos[$#t_pos]);

		@line=split(/\n/, $cds{$gene});
		my @cds1=split(/\t/, $line[0]);
		my @cds2=split(/\t/, $line[$#line]);
		my ($cds_start, $cds_end)=($cds1[3], $cds2[4]);

		for(my $i=0; $i<=$#t_pos; $i++){
			if($t_pos[$i]==$exon_start){
				$exon_start{$gene}=$i+1;
			}elsif($t_pos[$i]==$exon_end){
				$exon_end{$gene}=$i+1;
			}
			if($t_pos[$i]==$cds_start){
				$cds_start{$gene}=$i+1;
			}elsif($t_pos[$i]==$cds_end){
				$cds_end{$gene}=$i+1;
			}
		}
#		@line=split(/\n/, $stop{$gene});

#		if($cds1[6] eq "-"){
#			($cds_start, $cds_end)=($cds2[4], $cds1[3]);
#			@t_pos=reverse(@t_pos);
#		}


	}
}

my $count=0;
open(OUT, ">$outf") or die "cannot open $outf\n";
print OUT "Peak_ID\tgene_id\ttranscript_id\tChr\tstrand\tgenomic_pos\tsummit_pos_in_gene\tsummit_pos_in_rep1\tsummit_pos_in_rep2\tdoc_C3\tdoc_C8\tdoc_C13\tdoc_C15\tdoc_C18\tdoc_C20\tCategory(TSS/5'UTR/CDS/Stop/3'UTR)\n";
open(IN, $inf) or die "cannot open $inf\n";
while(<IN>){
	chomp;
	if($_=~/^gene_id/){
		my @info=split(/\t/, $_);
		if($info[0]=~/gene_id "([^"]+)"; transcript_id "([^"]+)";/){
			my $rest=join"\t", @info[1..$#info];
			my ($g_id, $t_id) = ($1, $2);
			$count++;
			if(exists $cds{$info[0]} && $info[2] eq "+"){
				if(!exists $cds_start{$info[0]} ||  !exists $cds_end{$info[0]}){
					print STDERR "Wrong gene:\n$_\n$cds_start{$info[0]}\n$cds_end{$info[0]}\n";
					next;
				}

				if($cds_start{$info[0]} == $exon_start{$info[0]} || $cds_end{$info[0]} == $exon_end{$info[0]}){
					print OUT "Peak_$count\t$g_id\t$t_id\t$rest\tlack of UTR info\n";
					next;
				}

				if($info[4]<=$cutoff && $info[4]<$cds_start{$info[0]}){
					print OUT "Peak_$count\t$g_id\t$t_id\t$rest\tTSS\n";
				}elsif($info[4]<$cds_start{$info[0]}){
					print OUT "Peak_$count\t$g_id\t$t_id\t$rest\t5'UTR\n";			
				}elsif($info[4]<=$cds_end{$info[0]}+3-$cutoff){ # in mm10.gtf cds_end doesn't include the stop codon
					print OUT "Peak_$count\t$g_id\t$t_id\t$rest\tCDS\n";
				}elsif($info[4]<=$cds_end{$info[0]}+3+$cutoff){
					print OUT "Peak_$count\t$g_id\t$t_id\t$rest\tStop\n";				
				}elsif($info[4]<=$exon_end{$info[0]}){
					print OUT "Peak_$count\t$g_id\t$t_id\t$rest\t3'UTR\n";
				}
			}elsif(exists $cds{$info[0]} && $info[2] eq "-"){
				if(!exists $cds_start{$info[0]} ||  !exists $cds_end{$info[0]}){
					print STDERR "$_\n$cds_start{$info[0]}\n$cds_end{$info[0]}\n";
					next;
				}
				if($cds_start{$info[0]} == $exon_start{$info[0]} || $cds_end{$info[0]} == $exon_end{$info[0]}){
					print OUT "Peak_$count\t$g_id\t$t_id\t$rest\tlack of UTR info\n";
					next;
				}
				if($info[4]<$cds_start{$info[0]}-3-$cutoff){
					print OUT "Peak_$count\t$g_id\t$t_id\t$rest\t3'UTR\n";			
				}elsif($info[4]<$cds_start{$info[0]}-3+$cutoff){
					print OUT "Peak_$count\t$g_id\t$t_id\t$rest\tStop\n";
				}elsif($info[4]<=$cds_end{$info[0]}){
					print OUT "Peak_$count\t$g_id\t$t_id\t$rest\tCDS\n";
				}elsif($info[4]>$exon_end{$info[0]}-200){
					print OUT "Peak_$count\t$g_id\t$t_id\t$rest\tTSS\n";
				}elsif($info[4]<=$exon_end{$info[0]}-200){
					print OUT "Peak_$count\t$g_id\t$t_id\t$rest\t5'UTR\n";
				}
			}else{
				print OUT "Peak_$count\t$g_id\t$t_id\t$rest\tNon-coding RNA\n";
			}
		}
	}
}
close IN;
close OUT;

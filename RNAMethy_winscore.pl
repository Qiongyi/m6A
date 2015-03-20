#!/usr/bin/perl -w
use strict;
use warnings;

use List::Util qw(sum);

if(@ARGV != 6) {
    print STDERR "Usage: RNAMethy_winscore.pl GTF gene_list IP input outf outf2(all)\n";
    exit(0);
}

# Analyses were limited to mouse genes with at least 40 reads per 1,000 nucleotides
my $cutoff=40; 
my $normalzation_factor=2000000;  ## normalized to 1M fragment (2M reads)
my ($gtf, $gene_list, $ip, $in, $outf, $outf2)= @ARGV;

my $sample_tag;
if($outf=~/(C\d+)\.winscore\.xls/){
	$sample_tag=$1."_cov";
}
system("mkdir -p $sample_tag");

my $normalization_factor_ip;
my $normalization_factor_in;

=all reads after filter
if($sample_tag =~ /^C3_cov/){
	$normalization_factor_ip= 17301942;
	$normalization_factor_in= 31783694;
}elsif($sample_tag =~ /^C8_cov/){
	$normalization_factor_ip= 19366306;
	$normalization_factor_in= 31936926;
}elsif($sample_tag =~ /^C13_cov/){
	$normalization_factor_ip= 18010072;
	$normalization_factor_in= 28329970;
}elsif($sample_tag =~ /^C15_cov/){
	$normalization_factor_ip= 20226724;
	$normalization_factor_in= 35066174;
}elsif($sample_tag =~ /^C18_cov/){
	$normalization_factor_ip= 17515350;
	$normalization_factor_in= 32109920;
}elsif($sample_tag =~ /^C20_cov/){
	$normalization_factor_ip= 17480324;
	$normalization_factor_in= 34450782;
}else{
	print STDERR "Wrong format for sample_tag: $sample_tag\n"; exit;
}
=cut

=mapped reads
if($sample_tag =~ /^C3_cov/){
	$normalization_factor_ip= 13281132;
	$normalization_factor_in= 23881203;
}elsif($sample_tag =~ /^C8_cov/){
	$normalization_factor_ip= 15367095;
	$normalization_factor_in= 25294237;
}elsif($sample_tag =~ /^C13_cov/){
	$normalization_factor_ip= 15087433;
	$normalization_factor_in= 22206713;
}elsif($sample_tag =~ /^C15_cov/){
	$normalization_factor_ip= 15742579;
	$normalization_factor_in= 26243342;
}elsif($sample_tag =~ /^C18_cov/){
	$normalization_factor_ip= 14390094;
	$normalization_factor_in= 23270195;
}elsif($sample_tag =~ /^C20_cov/){
	$normalization_factor_ip= 13889077;
	$normalization_factor_in= 25988853;
}else{
	print STDERR "Wrong format for sample_tag: $sample_tag\n"; exit;
}
=cut

if($sample_tag =~ /^C3_cov/){
	$normalization_factor_ip= 2501832;
	$normalization_factor_in= 5692246;
}elsif($sample_tag =~ /^C8_cov/){
	$normalization_factor_ip= 5752214;
	$normalization_factor_in= 3735602;
}elsif($sample_tag =~ /^C13_cov/){
	$normalization_factor_ip= 4153178;
	$normalization_factor_in= 5885242;
}elsif($sample_tag =~ /^C15_cov/){
	$normalization_factor_ip= 4023330;
	$normalization_factor_in= 4753470;
}elsif($sample_tag =~ /^C18_cov/){
	$normalization_factor_ip= 4491010;
	$normalization_factor_in= 3364516;
}elsif($sample_tag =~ /^C20_cov/){
	$normalization_factor_ip= 4345696;
	$normalization_factor_in= 6369686;
}else{
	print STDERR "Wrong format for sample_tag: $sample_tag\n"; exit;
}

my %transcript;

open(IN, $gtf) or die "Cannot open $gtf\n"; 
while(<IN>){
	my @info=split(/\t/, $_);
	if($info[2] eq "exon" && $info[0]=~/^chr[0-9XYMCxymc]+$/){
		if($info[8]=~/(gene_id "[^\"]+"\; transcript_id "[^\"]+";)/){
			$transcript{$1}.="$info[0]\t$info[3]\t$info[4]\t$info[6]\n";
		}else{
			print STDERR "Wrong for this line in GTF:\n$_\n"; exit;
		}
	}
}
close IN;

open(OUT, ">$outf") or die "Cannot open $outf\n"; 
open(OUT2, ">$outf2") or die "Cannot open $outf2\n"; 
print OUT "ID\tChr\tStrand\tWin_start-Win_end(gene)\tWin_start-Win_end(genome)\tMeanWin_IP\tMeanWin_In\tMedian_IP\tMedian_In\tWinscore\tWinscore_fdr\n";
print OUT2 "ID\tChr\tStrand\tWin_start-Win_end(gene)\tWin_start-Win_end(genome)\tMeanWin_IP\tMeanWin_In\tMedian_IP\tMedian_In\tWinscore\tWinscore_fdr\n";

open(IN, $gene_list) or die "Cannot open $gene_list\n"; 
while(<IN>){
	if($_=~/^ID/){
		next;
	}elsif($_=~/^gene_id/){
		chomp;
		my @info=split(/\t/, $_);
		### 12 and 13 for Naive IP; 14 and 15 for Context IP; 16 and 17 for FC IP
		if(($info[6]>=$cutoff && $info[7]>=$cutoff) || ($info[8]>=$cutoff && $info[9]>=$cutoff) || ($info[10]>=$cutoff && $info[11]>=$cutoff) || ($info[12]>=$cutoff && $info[13]>=$cutoff) || ($info[14]>=$cutoff && $info[15]>=$cutoff) || ($info[16]>=$cutoff && $info[17]>=$cutoff)){
#			print STDERR "$info[0]\n";
			my @gene=split(/\n/,$transcript{$info[0]});
			my @gene_pos;
			my @mean;
			my %hash_ip; ### record the read depths for each postion for IP
			my %hash_in; ### record the read depths for each postion for In
			foreach my $line (@gene){
				my @ele=split(/\t/,$line);
				for(my $i=$ele[1]; $i<=$ele[2]; $i++){
					push(@gene_pos, $i);
					$hash_ip{$i}=0; # default is 0
					$hash_in{$i}=0; # default is 0
				}
			}

			### get the sam file for this gene (IP)
			my $tmp_file1=$outf.".ip.sam";
			my $tmp_file2=$outf.".in.sam";
			my $start_sam=$info[3]-300>0?($info[3]-300):0;
			system("samtools view $ip $info[1]:$start_sam-$info[4] > $tmp_file1");
			system("samtools view $in $info[1]:$start_sam-$info[4] > $tmp_file2");
			### IP
			if(-s $tmp_file1){
				open(SAM, $tmp_file1) or die "Cannot open $!\n";
				while(<SAM>){
					my @sam=split(/\t/,$_);
					if($sam[8]=~/^(\d+)$/){
						for(my $i=$sam[3]; $i<=$sam[3]+$1-1; $i++){
							$hash_ip{$i}++;
						}
					}
				}
				close SAM;
			}
			### control
			if(-s $tmp_file2){
				open(SAM, $tmp_file2) or die "Cannot open $!\n";
				while(<SAM>){
					my @sam=split(/\t/,$_);
					if($sam[8]=~/^(\d+)$/){
						for(my $i=$sam[3]; $i<=$sam[3]+$1-1; $i++){
							$hash_in{$i}++;
						}
					}
				}
				close SAM;
			}
			system("rm -fr $tmp_file1");
			system("rm -fr $tmp_file2");
			my @gene_depth_ip;
			my @gene_depth_in;

			my $t_id;
			if($info[0]=~/transcript_id \"([^"]+)\"/){
				$t_id=$1;
			}
			my $input_cov=$t_id.".input.xls";
			my $ip_cov=$t_id.".ip.xls";
			open(OUT3, ">$sample_tag/$input_cov") or die "Cannot open $sample_tag/$input_cov\n";
			open(OUT4, ">$sample_tag/$ip_cov") or die "Cannot open $sample_tag/$ip_cov\n";

			for(my $i=0; $i<=$#gene_pos; $i++){
				push(@gene_depth_ip, $hash_ip{$gene_pos[$i]});
				push(@gene_depth_in, $hash_in{$gene_pos[$i]});
				my $pos=$i+1;
				print OUT3 "$pos\t$hash_in{$gene_pos[$i]}\t".(int($hash_in{$gene_pos[$i]}/$normalization_factor_in*$normalzation_factor+0.5))."\n";
				print OUT4 "$pos\t$hash_ip{$gene_pos[$i]}\t".(int($hash_ip{$gene_pos[$i]}/$normalization_factor_ip*$normalzation_factor+0.5))."\n";
			}
			close OUT3;
			close OUT4;

			### Median value for the gene
			my $Median_ip= &Median(\@gene_depth_ip);
			my $Median_in= &Median(\@gene_depth_in);
#			print STDERR "MedianGeneIP: $Median_ip; MedianGeneControl: $Median_in\n";
			### Mean value for each window 100bp with 50bp overlap
			for(my $i=0; $i<=$#gene_pos; $i+=50){
				my $j=$i+99;
				if($#gene_pos-$i+1<51){
					last;  ### no need for run fragment <=50bp because it has been included in the last 100bp window
				}
				if($#gene_pos-$i+1<100){
					#last window  if $j>=$#gene_pos
					$j=$#gene_pos;
				}
				my @gene_depth_win_ip=@gene_depth_ip[$i..$j];
				my @gene_depth_win_in=@gene_depth_in[$i..$j];
#				print STDERR @gene_depth_win_ip."\n";
#				print STDERR @gene_depth_win_in."\n";
				my $MeanWin_ip = &Mean(\@gene_depth_win_ip);
				my $MeanWin_in = &Mean(\@gene_depth_win_in);
#				print STDERR "MeanWinIP: $MeanWin_ip; MeanWinControl: $MeanWin_in\n";
				if($MeanWin_ip<1){
					$MeanWin_ip=1;
				}
				if($MeanWin_in<1){
					$MeanWin_in=1;
				}
				if($Median_ip<1){
					$Median_ip=1;
				}
				if($Median_in<1){
					$Median_in=1;
				}						
				my $winscore=log($MeanWin_ip*$Median_in/($Median_ip*$MeanWin_in))/log(2);
				my $winscore_fdr=log($MeanWin_in*$Median_ip/($Median_in*$MeanWin_ip))/log(2);
				if($winscore>=2 || $winscore_fdr>=2){
					print OUT "$info[0]\t$info[1]\t$info[2]\t".($i+1)."-".($j+1)."\t$gene_pos[$i]-$gene_pos[$j]\t$MeanWin_ip\t$MeanWin_in\t$Median_ip\t$Median_in\t$winscore\t$winscore_fdr\n";
				}
				print OUT2 "$info[0]\t$info[1]\t$info[2]\t".($i+1)."-".($j+1)."\t$gene_pos[$i]-$gene_pos[$j]\t$MeanWin_ip\t$MeanWin_in\t$Median_ip\t$Median_in\t$winscore\t$winscore_fdr\n";
#				print STDERR "Winscore: $winscore\n\n"; sleep 1;
			}
		}
	}
}
close IN;
close OUT;
close OUT2;

sub Mean{
	my $list=shift;
	my @list=@{$list};
	return sum(@list)/@list;
}

sub Median{
	my $list=shift;
	my @list= sort {$a<=>$b} @{$list};
	my $count=@list;
	if($count ==0){
		return undef;
	}
	if(($count%2)==1){
		return $list[int(($count-1)/2)];
	}else{
		return ($list[int(($count-1)/2)] + $list[int($count/2)])/2;
	}
}



#!/usr/bin/perl -w
use strict;
use warnings;
use Statistics::R;

if(@ARGV != 5) {
	print STDERR "Usage: RNAMethy_GenomicLocation.pl norm_file GTF Naive_SummitPosition.xls Naive_SummitPosition_GenomicLocation.xls Naive_SummitPosition_GenomicLocation_stat.xls\n";
    exit(0);
}
my ($norm, $gtf, $inf, $outf, $outf2)=@ARGV;
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

my %hash;

open(OUT, ">$outf") or die "cannot open $outf\n";
open(IN, $inf) or die "cannot open $inf\n";
while(<IN>){
	chomp;
	if($_=~/^ID/){
		print OUT "$_\tCategory(TSS/5'UTR/CDS/Stop/3'UTR)\n";
	}else{
		my @info=split(/\t/, $_);
		
		if(exists $cds{$info[0]} && $info[2] eq "+"){
			if(!exists $cds_start{$info[0]} ||  !exists $cds_end{$info[0]}){
				print "$_\n$cds_start{$info[0]}\n$cds_end{$info[0]}\n";
				next;
			}

			if($cds_start{$info[0]} == $exon_start{$info[0]} || $cds_end{$info[0]} == $exon_end{$info[0]}){
				print OUT "$_\tlack of UTR info\n";
				if(!exists $hash{$info[1].$info[2].$info[3]}){
					$hash{$info[1].$info[2].$info[3]}="lack of UTR info";
				}
				next;
			}

			if($info[4]<=$cutoff && $info[4]<$cds_start{$info[0]}){
				print OUT "$_\tTSS\n";
				if(!exists $hash{$info[1].$info[2].$info[3]}){
					$hash{$info[1].$info[2].$info[3]}="TSS";
				}elsif($hash{$info[1].$info[2].$info[3]} !~ /TSS/){
					$hash{$info[1].$info[2].$info[3]}.="/TSS";
				}
			}elsif($info[4]<$cds_start{$info[0]}){
				print OUT "$_\t5'UTR\n";
				if(!exists $hash{$info[1].$info[2].$info[3]}){
					$hash{$info[1].$info[2].$info[3]}="5'UTR";
				}elsif($hash{$info[1].$info[2].$info[3]} !~ /5'UTR/){
					$hash{$info[1].$info[2].$info[3]}.="/5'UTR";
				}				
			}elsif($info[4]<=$cds_end{$info[0]}+3-$cutoff){ # in mm10.gtf cds_end doesn't include the stop codon
				print OUT "$_\tCDS\n";
				if(!exists $hash{$info[1].$info[2].$info[3]}){
					$hash{$info[1].$info[2].$info[3]}="CDS";
				}elsif($hash{$info[1].$info[2].$info[3]} !~ /CDS/){
					$hash{$info[1].$info[2].$info[3]}.="/CDS";
				}	
			}elsif($info[4]<=$cds_end{$info[0]}+3+$cutoff){
				print OUT "$_\tStop\n";
				if(!exists $hash{$info[1].$info[2].$info[3]}){
					$hash{$info[1].$info[2].$info[3]}="Stop";
				}elsif($hash{$info[1].$info[2].$info[3]} !~ /Stop/){
					$hash{$info[1].$info[2].$info[3]}.="/Stop";
				}					
			}elsif($info[4]<=$exon_end{$info[0]}){
				print OUT "$_\t3'UTR\n";
				if(!exists $hash{$info[1].$info[2].$info[3]}){
					$hash{$info[1].$info[2].$info[3]}="3'UTR";
				}elsif($hash{$info[1].$info[2].$info[3]} !~ /3'UTR/){
					$hash{$info[1].$info[2].$info[3]}.="/3'UTR";
				}
			}
		}elsif(exists $cds{$info[0]} && $info[2] eq "-"){
			if(!exists $cds_start{$info[0]} ||  !exists $cds_end{$info[0]}){
				print "$_\n$cds_start{$info[0]}\n$cds_end{$info[0]}\n";
				next;
			}

			if($cds_start{$info[0]} == $exon_start{$info[0]} || $cds_end{$info[0]} == $exon_end{$info[0]}){
				print OUT "$_\tlack of UTR info\n";
				if(!exists $hash{$info[1].$info[2].$info[3]}){
					$hash{$info[1].$info[2].$info[3]}="lack of UTR info";
				}
				next;
			}

			if($info[4]<$cds_start{$info[0]}-3-$cutoff){
				print OUT "$_\t3'UTR\n";
				if(!exists $hash{$info[1].$info[2].$info[3]}){
					$hash{$info[1].$info[2].$info[3]}="3'UTR";
				}elsif($hash{$info[1].$info[2].$info[3]} !~ /3'UTR/){
					$hash{$info[1].$info[2].$info[3]}.="/3'UTR";
				}				
			}elsif($info[4]<$cds_start{$info[0]}-3+$cutoff){
				print OUT "$_\tStop\n";
				if(!exists $hash{$info[1].$info[2].$info[3]}){
					$hash{$info[1].$info[2].$info[3]}="Stop";
				}elsif($hash{$info[1].$info[2].$info[3]} !~ /Stop/){
					$hash{$info[1].$info[2].$info[3]}.="/Stop";
				}
			}elsif($info[4]<=$cds_end{$info[0]}){
				print OUT "$_\tCDS\n";
				if(!exists $hash{$info[1].$info[2].$info[3]}){
					$hash{$info[1].$info[2].$info[3]}="CDS";
				}elsif($hash{$info[1].$info[2].$info[3]} !~ /CDS/){
					$hash{$info[1].$info[2].$info[3]}.="/CDS";
				}
			}elsif($info[4]>$exon_end{$info[0]}-200){
				print OUT "$_\tTSS\n";
				if(!exists $hash{$info[1].$info[2].$info[3]}){
					$hash{$info[1].$info[2].$info[3]}="TSS";
				}elsif($hash{$info[1].$info[2].$info[3]} !~ /TSS/){
					$hash{$info[1].$info[2].$info[3]}.="/TSS";
				}
			}elsif($info[4]<=$exon_end{$info[0]}-200){
				print OUT "$_\t5'UTR\n";
				if(!exists $hash{$info[1].$info[2].$info[3]}){
					$hash{$info[1].$info[2].$info[3]}="5'UTR";
				}elsif($hash{$info[1].$info[2].$info[3]} !~ /5'UTR/){
					$hash{$info[1].$info[2].$info[3]}.="/5'UTR";
				}
			}
		}else{
			print OUT "$_\tNon-coding RNA\n";
			if(!exists $hash{$info[1].$info[2].$info[3]}){
					$hash{$info[1].$info[2].$info[3]}="Non-coding RNA";
			}elsif($hash{$info[1].$info[2].$info[3]} !~ /Non-coding RNA/){
					$hash{$info[1].$info[2].$info[3]}.="/Non-coding RNA";
			}
		}
	}
}
close IN;
close OUT;

my ($tss_count, $stop_count, $utr5_count, $utr3_count, $cds_count, $noncoding_count) = (0, 0, 0, 0, 0, 0);
# TSS > Stop > 5'UTR > 3'UTR > CDS if a peak locates in multiple categories
foreach my $key (keys %hash){
	if($hash{$key}=~/TSS/){
		$tss_count++;
	}
	elsif($hash{$key}=~/Stop/){
		$stop_count++;
	}
	elsif($hash{$key}=~/5'UTR/){
		$utr5_count++;
	}
	elsif($hash{$key}=~/3'UTR/){
		$utr3_count++;
	}
	elsif($hash{$key}=~/CDS/){
		$cds_count++;
	}
	elsif($hash{$key}=~/Non-coding RNA/){
		$noncoding_count++;
	}		
}
my $total= $tss_count+$stop_count+$utr5_count+$utr3_count+$cds_count;

my %norm;
my $num_total=0;
open(IN, $norm) or die "cannot open $norm\n";
while(<IN>){
	if($_=~/^Category/i){
		next;
	}else{
		chomp;
		my @info=split(/\t/, $_);
		$norm{$info[0]}=$info[1];
		$num_total+=$info[1];
	}
}
close IN;

my $norm_tss= $norm{"TSS"}/$num_total;
my $norm_utr5= $norm{"5'UTR"}/$num_total;
my $norm_cds= $norm{"CDS"}/$num_total;
my $norm_stop= $norm{"Stop"}/$num_total;
my $norm_utr3= $norm{"3'UTR"}/$num_total;

my $tss_by_chance=int($total*$norm_tss+0.5);
my $test_tss=($total-$tss_by_chance).",".$tss_by_chance.",".($total-$tss_count).",".$tss_count;

my $utr3_by_chance=int($total*$norm_utr3+0.5);
my $test_utr3=($total-$utr3_by_chance).",".$utr3_by_chance.",".($total-$utr3_count).",".$utr3_count;

my $utr5_by_chance=int($total*$norm_utr5+0.5);
my $test_utr5=($total-$utr5_by_chance).",".$utr5_by_chance.",".($total-$utr5_count).",".$utr5_count;

my $cds_by_chance=int($total*$norm_cds+0.5);
my $test_cds=($total-$cds_by_chance).",".$cds_by_chance.",".($total-$cds_count).",".$cds_count;

my $stop_by_chance=int($total*$norm_stop+0.5);
my $test_stop=($total-$stop_by_chance).",".$stop_by_chance.",".($total-$stop_count).",".$stop_count;

# Create a communication bridge with R and start R
my $R = Statistics::R->new();

$R->run(
qq`m<-matrix(c($test_tss),2,2,byrow=T)`,
qq`x <- fisher.test(m)`,
);
my $p_value_tss= $R -> get('x$p.value');

$R->run(
qq`m<-matrix(c($test_utr3),2,2,byrow=T)`,
qq`x <- fisher.test(m)`,
);
my $p_value_utr3= $R -> get('x$p.value');

$R->run(
qq`m<-matrix(c($test_utr5),2,2,byrow=T)`,
qq`x <- fisher.test(m)`,
);
my $p_value_utr5= $R -> get('x$p.value');

$R->run(
qq`m<-matrix(c($test_cds),2,2,byrow=T)`,
qq`x <- fisher.test(m)`,
);
my $p_value_cds= $R -> get('x$p.value');

$R->run(
qq`m<-matrix(c($test_stop),2,2,byrow=T)`,
qq`x <- fisher.test(m)`,
);
my $p_value_stop= $R -> get('x$p.value');

open(OUT, ">$outf2") or die "cannot open $outf2\n";
print OUT "Category\tnumber\t%\tenrichment\tP-value (Fisher's Exact Test)\tNumbers_for_test\n";
printf OUT ("%s\t%d\t%.2f%%\t%.2f\t","TSS", $tss_count, $tss_count/$total*100, $tss_count/$total/$norm_tss);
print OUT "$p_value_tss\t$test_tss\n";
printf OUT ("%s\t%d\t%.2f%%\t%.2f\t","5'UTR", $utr5_count, $utr5_count/$total*100,  $utr5_count/$total/$norm_utr5);
print OUT "$p_value_utr5\t$test_utr5\n";
printf OUT ("%s\t%d\t%.2f%%\t%.2f\t","CDS", $cds_count, $cds_count/$total*100, $cds_count/$total/$norm_cds);
print OUT "$p_value_cds\t$test_cds\n";
printf OUT ("%s\t%d\t%.2f%%\t%.2f\t","Stop", $stop_count, $stop_count/$total*100, $stop_count/$total/$norm_stop);
print OUT "$p_value_stop\t$test_stop\n";
printf OUT ("%s\t%d\t%.2f%%\t%.2f\t","3'UTR", $utr3_count, $utr3_count/$total*100, $utr3_count/$total/$norm_utr3);
print OUT "$p_value_utr3\t$test_utr3\n";
close OUT;

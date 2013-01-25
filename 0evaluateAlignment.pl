use warnings;
use strict;

my $dir = '/media/DataDrive/sutureGenomics/alignmentTest/';
my $timing = $dir . 'alignerTiming.out';
my @contacts = (["Carlia_N","Carlia_S"]);
my @progs = qw(bowtie bowtie2 bwa novo smalt soap stampy);
#my @progs = qw(bowtie2);
my $libfile = '/media/DataDrive/sutureGenomics/library';

my $qual = '/media/DataDrive/sutureGenomics/qualStats/qualityControl.out';
my $masterout = $dir . "alignEvaluation.out";

my %numReads;
open(IN, "<$qual");
while(<IN>) {
	chomp(my $line = $_);
	$line =~ s/"//g;
	my @d = split("\t",$line);
	$numReads{$d[1]}{'reads'} = $d[4];
	$numReads{$d[1]}{'pairs'} = $d[6];
	}
close(IN);	

#defines the library
my %lib;
open(IN, "<$libfile");
while (<IN>) {
        chomp(my $line = $_);
        my @d = split(/\t/, $line);
        my $lib = $d[2] . '/' . $d[0] . '/';
        $lib{$d[2]}{$d[0]} = {'dir' => $lib, 'index' => $d[1]}; 
        } 
close(IN);
	
#my %snp;
#similarity across snps & AFs
#foreach my $contact (@contacts) {
#	for (my $i = 0; $i < scalar(@progs); $i++) {
#	    my $prog1 = $progs[$i];
#	    for (my $j = $i + 1; $j < scalar(@progs); $j++) {
#		my $prog2 = $progs[$j];
#			if ($prog2 ne $prog1) {
#				print "comparing $prog1 to $prog2 now!\n";
#				my @sort = sort {$a cmp $b} ($prog1,$prog2);
#				my $compare = $sort[0] . "_" . $sort[1];
#		
#				my $name = $contact->[0] . '_' . $contact->[1];
#		
#				my $vcf1 = $dir . $name . "." . $prog1 . ".vcf";
#				my $vcf2 = $dir . $name . "." . $prog2 . ".vcf";
#		
#				my ($both,$uniq,$afcorr,$sameGeno,$diffGeno) = compare($vcf1,$vcf2);				
#				my $total = $both * 2 + $uniq;
#				my $totalGeno = $sameGeno + $diffGeno;
				
#				$snp{$name}{$compare}{'both'} = ($both * 2)/$total;
#				$snp{$name}{$compare}{'uniq'} = $uniq/$total;
#				$snp{$name}{$compare}{'af'} = $afcorr;	
#				$snp{$name}{$compare}{'sameGeno'} = $sameGeno/$totalGeno;
#				$snp{$name}{$compare}{'diffGeno'} = $diffGeno/$totalGeno;
#				}
#			}
#		}
#	}	
	
#my %values;	
#foreach my $name (keys %snp) {
#	foreach my $compare (keys %{$snp{$name}}) {
#		foreach my $key (keys %{$snp{$name}{$compare}}) {
#			$values{$key}++;
#			}
#		}
#	}

#my $out = $dir . "alignmentComparison_SNP.out";
#open(OUT, ">$out");
#my @head = sort {$a cmp $b} keys %values;
#print OUT "name\tcomparison\t"; print OUT join("\t",@head); print OUT "\n";
#foreach my $name (keys %snp) {
#	foreach my $compare (keys %{$snp{$name}}) {
#	    print OUT $name, "\t", $compare, "\t";
#		foreach my $key (@head) {
#			print OUT $snp{$name}{$compare}{$key}, "\t";
#			}
#		print OUT "\n";
#		}
#	}	
#close (OUT);	
			
my %master;
foreach my $prog (@progs) {
    my @bam = <$dir*$prog\.bam>;
	foreach my $bam (@bam) {
	    
	    my $lib;
	    if ($bam =~ m/(SS\d+)/) {
		$lib = $1;
	    }

		print "working on bam analysis of $lib now!\n";	
		my $call1 = system("samtools index $bam") unless(-f $bam . ".bai");
		my $pairs = paired($bam);
		my $reads = allReads($bam);
		my $time = timeAlign($timing,$lib,$prog);		
		$master{$lib}{$prog}{'pairs'} = $pairs/$numReads{$lib}{'pairs'};
		$master{$lib}{$prog}{'reads'} = $reads/$numReads{$lib}{'reads'};
		$master{$lib}{$prog}{'time'} = $time;
		}
	}
	
my %values1;	
foreach my $name (keys %master) {
	foreach my $program (keys %{$master{$name}}) {
		foreach my $key (keys %{$master{$name}{$program}}) {
			$values1{$key}++;
			}
		}
	}

my $out1 = $dir . "alignmentComparison_count.out";
open(OUT, ">$out1");
my @head1 = sort {$a cmp $b} keys %values1;
print OUT "name\tprogram\t"; print OUT join("\t",@head1); print OUT "\n";
foreach my $name (keys %master) {
	foreach my $prog (keys %{$master{$name}}) {
	    print OUT $name, "\t", $prog, "\t";
		foreach my $key (@head1) {
			print OUT $master{$name}{$prog}{$key}, "\t";
			}
	    print OUT "\n";
		}
	}	
close (OUT);		
		
##########################
# behold the subroutines #
##########################
	
#num of paired reads mapped
sub paired {
	my ($bam) = @_;
	my @call = `cat $bam | samtools view -f 0x0002 - | wc`;
	my $reads = $1 if $call[0] =~ m/^(\d+)/;
	$reads = 'NA' unless($reads);
	return($reads);
	}

#num of reads total mapped
sub allReads {
	my ($bam) = @_;
	my @call = `samtools idxstats $bam`;
	my $totReads;
	foreach my $line (@call) {
		my @d = split(/\t/,$line);
		$totReads += $d[2];
		}
	return($totReads);
	}
	
#time it takes to run
sub timeAlign {
	my ($file,$lib,$prog) = @_;
	
	my %time;
	open(IN, "<$file");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		$time{$d[0]}{$d[1]} = $d[2];
		}
	close(IN);	
	
	return($time{$prog}{$lib});
	}
	
sub compare {
	my ($vcf1,$vcf2) = @_;
	
	my %vcf; my %geno;
	open(IN, "<$vcf1");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		if ($d[0] =~ m/contig/) {
			$vcf{$d[0]}{$d[1]}{'1'}++ while $line =~ m/1\//g;
			$vcf{$d[0]}{$d[1]}{'1'}++ while $line =~ m/\/1/g;
			
			for (my $i = 9; $i < scalar(@d); $i++) {
				$geno{$d[0]}{$d[1]}{'1'}{$i} = $1 if $d[$i] =~ m/(\d+\/\d+)/;
				}
			}
		}	
	close(IN);
	
	open(IN, "<$vcf2");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		if ($d[0] =~ m/contig/) {
			$vcf{$d[0]}{$d[1]}{'2'}++ while $line =~ m/1\//g;
			$vcf{$d[0]}{$d[1]}{'2'}++ while $line =~ m/\/1/g;
			
			for (my $i = 9; $i < scalar(@d); $i++) {
				$geno{$d[0]}{$d[1]}{'2'}{$i} = $1 if $d[$i] =~ m/(\d+\/\d+)/;
				}
			}
		}	
	close(IN);
	
	my ($sameGeno,$diffGeno);
	foreach my $c (keys %geno) {
		foreach my $pos (keys %{$geno{$c}}) {
			if ($geno{$c}{$pos}{'1'} && $geno{$c}{$pos}{'2'}) {
				foreach my $indiv (keys %{$geno{$c}{$pos}{'1'}}) {
					if ($geno{$c}{$pos}{'2'}{$indiv}) {
						if ($geno{$c}{$pos}{'2'}{$indiv} eq $geno{$c}{$pos}{'1'}{$indiv}) {
							$sameGeno++;
							}
						else {
							$diffGeno++;
							}
						}
					}
				}
			}
		}	
	
	my $both; my $uniq; my $x = []; my $tracker = 0;
	foreach my $c (keys %vcf) {
		foreach my $pos (keys %{$vcf{$c}}) {
			if ($vcf{$c}{$pos}{'1'} && $vcf{$c}{$pos}{'2'}) {
				$both++;
				$x -> [$tracker][1] = $vcf{$c}{$pos}{'1'};
				$x -> [$tracker][2] = $vcf{$c}{$pos}{'2'};
				$tracker++;
				}
			else {
				$uniq++;
				}
			}
		}		
	my $afcorr = correlation($x);
	return($both,$uniq,$afcorr,$sameGeno,$diffGeno);	
	}	
	

###################################################################################
# Pearson correlation code from                                                   #  
# http://davetang.org/muse/2010/11/29/calculating-pearson-correlation-using-perl/ #	
###################################################################################	

sub mean {
	my ($x)=@_;
	my $num = scalar(@{$x}) - 1;
	my $sum_x = '0';
	my $sum_y = '0';
	for (my $i = 0; $i < scalar(@{$x}); ++$i){
		$sum_x += $x->[$i][1];
		$sum_y += $x->[$i][2];
   		}
   	my $mu_x = $sum_x / $num;
   	my $mu_y = $sum_y / $num;
   	return($mu_x,$mu_y);
	}

sub ss {
	my ($x,$mean_x,$mean_y,$one,$two)=@_;
	my $sum = '0';
	for (my $i=0;$i<scalar(@{$x});++$i){
		$sum += ($x->[$i][$one]-$mean_x)*($x->[$i][$two]-$mean_y);
		}
	return $sum;
	}
 
sub correlation {
	my ($x) = @_;
	my ($mean_x,$mean_y) = mean($x);
	my $ssxx=ss($x,$mean_x,$mean_y,1,1);
	my $ssyy=ss($x,$mean_x,$mean_y,2,2);
	my $ssxy=ss($x,$mean_x,$mean_y,1,2);
	my $xcorrel;
	if ($ssxx == 0 || $ssyy == 0 || $ssxy == 0) {
		$xcorrel = 'NA';
		}
	else {	
		my $correl=correl($ssxx,$ssyy,$ssxy);
		$xcorrel=sprintf("%.4f",$correl);
		}
	return($xcorrel);
	}
 
sub correl {
	my($ssxx,$ssyy,$ssxy)=@_;
	my $sign=$ssxy/abs($ssxy);
	my $correl=$sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy));
	return $correl;
	}		
		

		

use warnings;
use strict;
use List::Util qw[min max];

#this isn't really reciprocal blast homology! there is no reciprocal blast to it. sadness...but i don't think that is necessary

my $main_dir = '/Users/singhal/Desktop/genomics/seqfiles/';
my $reads_dir = '/Users/singhal/Desktop/genomics/seqfiles/';
my $carlia_n = $main_dir . 'Carlia_N_trinity.fa.final.annotated';
my $carlia_s = $main_dir . 'Carlia_S_trinity.fa.final.annotated';
my $lampro_n = $main_dir . 'Lampro_N_trinity.fa.final.annotated';
my $lampro_c = $main_dir . 'Lampro_C_trinity.fa.final.annotated';
my $lampro_s = $main_dir . 'Lampro_S_trinity.fa.final.annotated';
my $sapro_c = $main_dir . 'Sapro_C_trinity.fa.final.annotated';
my $sapro_s = $main_dir . 'Sapro_S_trinity.fa.final.annotated';
my $out = $main_dir . 'homologyComparison_redo.out';
my $np = 8;
my $evalue = 1e-20;
my $cov = '5'; #trust the other contig if it is at least at 5x coverage

my %compare = (
'CarliaN_CarliaS' => {'seq1' => $carlia_n, 'seq2' => $carlia_s },
'LamproN_LamproS' => {'seq1' => $lampro_n, 'seq2' => $lampro_s },
'LamproC_LamproS' => {'seq1' => $lampro_c, 'seq2' => $lampro_s },
'SaproC_SaproS' => {'seq1' => $sapro_c, 'seq2' => $sapro_s },
'CarliaN_LamproN' => {'seq1' => $carlia_n, 'seq2' => $lampro_n },
'CarliaN_SaproC' => {'seq1' => $carlia_n, 'seq2' => $sapro_c },
'SaproC_LamproN' => {'seq1' => $sapro_c, 'seq2' => $lampro_n }
);

open(FINAL, ">$out");
foreach my $compare (keys %compare) {
	my $seq1 = $compare{$compare}{'seq1'};
	my $seq2 = $compare{$compare}{'seq2'};
	#my $s1 = parseSeq($seq1);
	#my $s2 = parseSeq($seq2);
	#compare($s1,$s2,$compare);
	
	my $s1b = parseSeq2($seq1);
	my $s2b = parseSeq2($seq2);
	recipBlastNuc($s1b,$s2b,$seq2,$compare);
	recipBlastProt($s1b,$s2b,$seq2,$compare);

#	my ($bam,$pileup) = mapping($s1,$seq1,$compare);
#	snpHomology($seq1,$pileup,$compare);
	}
close(FINAL);

sub snpHomology {
	my ($seq1,$pileup,$compare) = @_;
	
	my %l;
	open(SEQ, "<$seq1");
	while(<SEQ>) {
		chomp(my $line = $_);
		chomp(my $seq = <SEQ>);
		if ($line =~ m/>(\S+).*gs(\d+).*ge(\d+)/) {
			my $c = $1; my $gs = $2; my $ge = $3; my $cds = $3 - $2 + 1;
			$l{$c} = {'length' => length($seq), 'gs' => $gs, 'ge' => $ge, 'cds' => $cds};
			}
		}
	close(SEQ);	
	
	my %d; my $c = 'NA'; my %cds;
	open(IN,"<$pileup");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		if ($d[0] eq $c) {
			if (length($d[5]) >= $cov) {
				if ($d[1] >= $l{$d[0]}{'gs'} && $d[1] <= $l{$d[0]}{'ge'}) {
					$d{'align'}++; $cds{'align'}++;
					my @matches = $d[4] =~ m/([atgc])/ig;
					if (scalar(@matches)/length($d[5]) >= 0.5) {
						$d{'mis'}++; 
						$cds{'mis'}++;
						}
					my @indels = $d[4] =~ m/[+|-](\d+)/ig;	
					if (scalar(@indels)/length($d[5]) >= 0.5) {	
						$d{'gap'} += $indels[0];
						$cds{'gap'} += $indels[0];
						}
					}
				else {
					$d{'align'}++; 
					my @matches = $d[4] =~ m/([atgc])/ig;
					if (scalar(@matches)/length($d[5]) >= 0.5) {
						$d{'mis'}++; 
						}
					my @indels = $d[4] =~ m/[+|-](\d+)/ig;	
					if (scalar(@indels)/length($d[5]) >= 0.5) {	
						$d{'gap'} += $indels[0];
						}
					}
				}
			}
		else {	
			#accounting
			if (%d && $l{$c}{'cds'}) {
				$d{'gap'} = 0 unless($d{'gap'});
				$cds{'gap'} = 0 unless($cds{'gap'});
				$d{'mis'} = 0 unless($d{'mis'});
				$cds{'mis'} = 0 unless($cds{'mis'});
				$cds{'align'} = 0 unless($cds{'align'});
				$d{'peralign'} = $d{'align'}/$l{$c}{'length'};
				$cds{'peralign'} = $cds{'align'}/$l{$c}{'cds'};
				print FINAL "snp\t>$c\t$compare\tall\t$d{'align'}\t$d{'peralign'}\t$d{'mis'}\t$d{'gap'}\n";	
				print FINAL "snp\t>$c\t$compare\tcds\t$cds{'align'}\t$cds{'peralign'}\t$cds{'mis'}\t$cds{'gap'}\n";	
				}
			#going to be sloppy and ignore the first bp of every contig	
			%d = (); %cds = ();
			$c = $d[0];
			}
		}
	}
	
sub mapping {
	my ($s1,$seq1,$compare) = @_;

	my $contact = $1 if $seq1 =~ m/$main_dir([a-z]+_[a-z])/i;
	my $reads1gz = $reads_dir . $contact . '/' . $contact . '_1.fastq.gz';
	my $reads2gz = $reads_dir . $contact . '/' . $contact . '_2.fastq.gz';
	my $readsugz = $reads_dir . $contact . '/' . $contact . '_u.fastq.gz';
	
	my %seq1 = %$s1;
	my $seq = 'annotatedRef.fa';
	open(OUT, ">$seq");
	foreach my $c (keys %seq1) {
		print OUT "$seq1{$c}{'contig'}\n", $seq1{$c}{'seq'}, "\n";
		}	
	close(OUT);	
	
	
	my $reads1 = unzip($reads1gz);
	my $reads2 = unzip($reads2gz);
	my $readsu = unzip($readsugz);
	my ($bam,$pileup) = samtools($reads1,$reads2,$readsu,$compare,$seq);
	rezip($reads1); rezip($reads2); rezip($readsu);
	return($bam,$pileup);
	}

sub rezip {
	my ($file) = @_;
	my $call = system("gzip -1 $file");
	}
	
sub unzip {
	my ($filegz) = @_;
	my $call = system("gunzip $filegz");
	my $out = $1 if $filegz =~ m/(.*)\.gz/;
	return($out);
	}
	
sub samtools {
	my ($reads1,$reads2,$readsu,$compare,$seq1) = @_;

	my $out = $main_dir .  $compare . ".mpileup";
	my $bam = $main_dir . $compare . ".sorted";
	my $call1 = system("bowtie2-build $seq1 $seq1.bowtie2");
	my $call2 = system("bowtie2 -x $seq1.bowtie2 -1 $reads1 -2 $reads2 -S bowtie1.sam -5 5 -3 5 --sensitive -k 10 -X 300 -p $np");
	my $call3 = system("bowtie2 -x $seq1.bowtie2 $readsu -S bowtie2.sam -5 5 -3 5 --sensitive -k 10 -p $np");
	my $call4 = system("samtools view -bS bowtie1.sam > bowtie1.bam");
	my $call5 = system("samtools view -bS bowtie2.sam > bowtie2.bam");
	my $call6 = system("samtools merge bowtie.bam bowtie1.bam bowtie2.bam");
	my $call7 = system("samtools sort bowtie.bam $bam");
	my $call8 = system("samtools faidx $seq1");
	$bam = $bam . '.bam';
	my $call9 = system("samtools mpileup -f $seq1 $bam > $out");
	my $call10 = system("rm $seq1.bowtie2* bowtie1.sam bowtie2.sam bowtie1.bam bowtie2.bam bowtie.bam");

	return($bam,$out);
	}

sub compare {
	my ($s1,$s2,$compare) = @_;
	
	my %seq1 = %$s1; my %seq2 = %$s2;

	
	foreach my $c (keys %seq1) {
		if ($seq2{$c}) {
			my %entire; my %cds; my %l1; my %l2;
			
			#align entire thing
			open(TAR, ">target.fa");
			open(QUER, ">query.fa");
			print TAR ">$c" , "_1\n$seq1{$c}{'seq'}\n";
			print QUER ">$c" , "_2\n$seq2{$c}{'seq'}\n";
			close(TAR); close(QUER);		
			my @call = `blastn -query query.fa -subject target.fa -outfmt 6`;
			if (@call) {
				foreach my $line (@call) {
					my @d = split(/\t/,$line);
					for (my $i = $d[6]; $i <= $d[7]; $i++) {
						$l1{$i}++;
						}
					$entire{'mis'} += $d[4];
					$entire{'gap'} += $d[5];
					}
				$entire{'align'} = scalar(keys %l1);
				$entire{'peralign'} = $entire{'align'}/min(length($seq1{$c}{'seq'}),length($seq2{$c}{'seq'}));
				$entire{'peralign'} = 1 if $entire{'peralign'} > 1;
				print FINAL "annotated\t$seq1{$c}{'contig'}\t$compare\tall\t$entire{'align'}\t$entire{'peralign'}\t$entire{'mis'}\t$entire{'gap'}\n";	
				}
			
			#align just cds
			my $info1 = $seq1{$c}{'info'};
			my $gs1 = $1 if $info1 =~ m/gs(\d+)/;
			my $ge1 = $1 if $info1 =~ m/ge(\d+)/;
			my $length = $ge1 - $gs1 + 1;
			my $seq1 = substr $seq1{$c}{'seq'}, $gs1 - 1, $length;			
			open(TAR, ">target.fa");
			open(QUER, ">query.fa");
			print TAR ">$c" , "_1\n$seq1\n";
			print QUER ">$c" , "_2\n$seq2{$c}{'seq'}\n";
			close(TAR); close(QUER);
			my @call1 = `blastn -query query.fa -subject target.fa -outfmt 6`;
			if (@call1) {
				foreach my $line (@call1) {
					my @d = split(/\t/,$line);
					for (my $i = $d[6]; $i <= $d[7]; $i++) {
						$l2{$i}++;
						}
					$cds{'mis'} += $d[4];
					$cds{'gap'} += $d[5];
					}
				$cds{'align'} = scalar(keys %l2);	
				$cds{'peralign'} = $cds{'align'}/$length;
				$cds{'peralign'} = 1 if $cds{'peralign'} > 1;
				print FINAL "annotated\t$seq1{$c}{'contig'}\t$compare\tcds\t$cds{'align'}\t$cds{'peralign'}\t$cds{'mis'}\t$cds{'gap'}\n";	
				}
			}
		}	
		
	my $call = system("rm query.fa target.fa");
	}

sub parseSeq {
	my ($seq) = @_;
	
	my %seq;
	
	open(IN, "<$seq");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>/) {
			my @d = split(/\t/,$line);
			chomp(my $seq = <IN>);
			if (scalar(@d) > 1) {
				if ($seq{$d[2]}) {
					$seq{$d[2]} = {'seq' => $seq, 'info' => $d[1], 'contig' => $d[0]} if length($seq) > length($seq{$d[2]}{'seq'});
					}
				else {
					$seq{$d[2]} = {'seq' => $seq, 'info' => $d[1], 'contig' => $d[0]};
					}
				}	
			}
		}	
	close(IN);	
	return(\%seq);	
	}
	
sub parseSeq2 {
	my ($seq) = @_;
	
	my %seq;
	
	open(IN, "<$seq");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $c = $1;
			my @d = split(/\t/,$line);
			chomp(my $seq = <IN>);
			if (scalar(@d) > 1) {
				$seq{$c} = {'seq' => $seq, 'info' => $d[1]};
				}	
			else {
				$seq{$c} = {'seq' => $seq};
				}
			}
		}	
	close(IN);	
	return(\%seq);	
	}	
	
sub	newdb {
	my ($s) = @_;
	
	open(IN, "<$s");
	my $newseq = 'newseq.fa';
	open(OUT, ">$newseq");
	while(<IN>) {
		my $line = $_;
		if ($line =~ m/(>\S+)/) {
			print OUT $1, "\n";
			}
		else {
			print OUT $line;
			}
		}	
	close(IN);
	close(OUT);
	return($newseq);
	}
	
sub recipBlastProt {
	my ($s1,$s2,$seq2,$compare) = @_;
	
	my %l;
	my %seq1 = %$s1; my %seq2 = %$s2;
	my $blast = "recipBlastProt.out";
	open(OUT, ">annotatedRef.fa");
	foreach my $c (keys %seq1) {
		if ($seq1{$c}{'info'}) {
			my $info1 = $seq1{$c}{'info'};
			my $gs1 = $1 if $info1 =~ m/gs(\d+)/;
			my $ge1 = $1 if $info1 =~ m/ge(\d+)/;
			my $length = $ge1 - $gs1 + 1;
			my $seq1 = substr $seq1{$c}{'seq'}, $gs1 - 1, $length;
			$seq1{$c}{'cds'} = $seq1;
			print OUT ">$c\n", $seq1, "\n";
			}
		}
	my $newseq = newdb($seq2);	
	my $call1 = system("formatdb -i $newseq -p F");
	my $call2 = system("blastall -p blastn -d $newseq -i annotatedRef.fa -a $np -e $evalue -m 8 -o $blast -b 1");
	close(OUT);
	my $call3 = system("rm $newseq*");		
	my $call = system("rm annotatedRef.fa");	
	
	my %matches;
	open(IN, "<$blast");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		$matches{$d[0]} = $d[1];
		}
	close(IN);
	
	foreach my $c1 (keys %matches) {
		open(TAR, ">target.fa");
		open(QUER, ">query.fa");
		print TAR ">$c1\n$seq1{$c1}{'cds'}\n";
		print QUER ">$matches{$c1}\n$seq2{$matches{$c1}}{'seq'}\n";
		close(TAR); close(QUER);		
		my @call = `blastn -query query.fa -subject target.fa -outfmt 6`;
		my %entire; my %l1;
		if (@call) {
			foreach my $line (@call) {
				my @d = split(/\t/,$line);
				for (my $i = $d[6]; $i <= $d[7]; $i++) {
					$l1{$i}++;
					}
				$entire{'mis'} += $d[4];
				$entire{'gap'} += $d[5];
				}
			$entire{'align'} = scalar(keys %l1);
			$entire{'peralign'} = $entire{'align'}/min(length($seq1{$c1}{'cds'}),length($seq2{$matches{$c1}}{'seq'}));
			$entire{'peralign'} = 1 if $entire{'peralign'} > 1;
			print FINAL "blast\t$c1\t>", "$compare\tcds\t$entire{'align'}\t$entire{'peralign'}\t$entire{'mis'}\t$entire{'gap'}\n";	
			}
		}	
	
	unlink($blast);	
	}

sub recipBlastNuc {
	my ($s1,$s2,$seq2,$compare) = @_;
	my %seq1 = %$s1; my %seq2 = %$s2;	
	my $blast = "recipBlastNuc.out";	
	open(OUT, ">annotatedRef.fa");
	my %l;
	foreach my $c (keys %seq1) {
		if ($seq1{$c}{'info'}) {
			print OUT ">$c\n", $seq1{$c}{'seq'}, "\n";
			}
		}	
	close(OUT);	
	my $newseq = newdb($seq2);	
	my $call1 = system("formatdb -i $newseq -p F");
	my $call2 = system("blastall -p blastn -d $newseq -i annotatedRef.fa -a $np -e $evalue -m 8 -o $blast -b 1");	
	my $call3 = system("rm $newseq*");
	my $call = system("rm annotatedRef.fa");	
		
	my %matches;
	open(IN, "<$blast");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		$matches{$d[0]} = $d[1];
		}
	close(IN);
		
	foreach my $c1 (keys %matches) {
		open(TAR, ">target.fa");
		open(QUER, ">query.fa");
		print TAR ">$c1\n$seq1{$c1}{'seq'}\n";
		print QUER ">$matches{$c1}\n$seq2{$matches{$c1}}{'seq'}\n";
		close(TAR); close(QUER);		
		my %entire; my %l1;
		my @call = `blastn -query query.fa -subject target.fa -outfmt 6`;
		if (@call) {
			foreach my $line (@call) {
				my @d = split(/\t/,$line);
				for (my $i = $d[6]; $i <= $d[7]; $i++) {
					$l1{$i}++;
					}
				$entire{'mis'} += $d[4];
				$entire{'gap'} += $d[5];
				}
			$entire{'align'} = scalar(keys %l1);
			$entire{'peralign'} = $entire{'align'}/min(length($seq1{$c1}{'seq'}),length($seq2{$matches{$c1}}{'seq'}));
			$entire{'peralign'} = 1 if $entire{'peralign'} > 1;
			print FINAL "blast\t$c1\t>", "$compare\tall\t$entire{'align'}\t$entire{'peralign'}\t$entire{'mis'}\t$entire{'gap'}\n";	
			}
		}
	
	unlink($blast);	
	}
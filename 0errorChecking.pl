##################################################################################
# a script to look at error rates by comparing two reads                         #
# external dependencies: blat (v34), trimmomatic                                 #
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 23 Jan 2012           #
##################################################################################

use warnings;
use strict;

#home directory with all my files; make sure to end with a backslash
my $dir = '/media/DataDrive/sutureGenomics/';
#this file has all the info about my library, and creates a directory of all my files
my $lib = '/media/DataDrive/sutureGenomics/library';
#this file has a list of the adaptor sequences for my library, including the reverse complements.
my $adaptorFile = '/media/DataDrive/sutureGenomics/illuminaAdapters.fa';
#path to trimmomatic executable
my $trimmomatic = '/home/singhal/programs/trimmomatic/trimmomatic-0.17.jar';
#path to blat executable
my $blat = 'blat';
#how many reads to test to look for error; can take a long time, so don't recommend doing more than 1e5
my $sample = 100000;

#defines the library
open(IN, "<$lib");
my %lib;
while (<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/, $line);
	my $lib = $d[2] . '/' . $d[0] . '/';
	$lib{$d[0]} = {'dir' => $lib, 'index' => $d[1]};	
	}
close(IN);

#defines the adaptors
open(AD, "<$adaptorFile");
my %ad;
while(<AD>) {
	chomp(my $line = $_);
	if ($line =~ m/>(\S+)/) {
		my $id = $1;
		chomp (my $seq = <AD>);
		$ad{$id} = $seq;
		}
	}
close(AD);	

###########################
# run the subroutines     #
###########################	

foreach my $lib (sort {$a cmp $b} keys %lib) {
	my $subdir = $dir . $lib{$lib}{'dir'};
	if (-d $subdir) {
		my $out = $subdir . "errorChecking.out";
		unless(-f $out) {
			print "Starting to process files from $lib! So exciting.\n";
			my $file1gz = $subdir . $lib . '_1.fastq.gz';
			my $file2gz = $subdir . $lib . '_2.fastq.gz';
			my $file1   = $subdir . $lib . '_1.fastq';
			my $file2   = $subdir . $lib . '_2.fastq';

			my $start1 = time;
			unzip($file1gz,$file2gz);
			my $time1 = int((time - $start1)/60);
			print "Unzipped files in $lib in $time1 minutes! Now we are cooking.\n";

			my $start2 = time;
			my ($out1p,$out2p,$out1u,$out2u) = adaptorRemovalandTrim($subdir,$file1,$file2,\%ad);
			my $time2 = int((time - $start2)/60);
			print "Found adaptor and trimmed in $lib in $time2 minutes! Slowly but surely...\n";			
			my $start3 = time;
			my ($qual,$error,$pos,$rpos,$totAlign) = errorFind($out1p,$out2p);
			my $time3 = int((time - $start3)/60);
			print "Found errors in $time3 minutes! Slowly but surely...\n";
			
			my $start4 = time;
			eraseAndZip($out1p,$out2p,$out1u,$out2u,$file1,$file2);
			my $time4 = int((time - $start4)/60);
			print "Zipped and cleaned up in $lib in $time4 minutes! Done!\n";
	
			my %qual = %{$qual}; my %error = %{$error}; my %pos = %{$pos}; my %rpos = %{$rpos};
	
			open(OUT, ">$out");
			foreach my $qual (keys %qual) {
				print OUT "QUAL\t", $qual, "\t", $qual{$qual}, "\n";
				}
			foreach my $e1 (keys %error) {
				foreach my $e2 (keys %{$error{$e1}}) {
					print OUT "ERR\t$e1", "_", $e2, "\t", $error{$e1}{$e2}, "\n";
					}
				}
			foreach my $pos (keys %pos) {
				print OUT "POS\t$pos\t$pos{$pos}\n";
				}
			foreach my $rpos (keys %rpos) {
				print OUT "RPOS\t$rpos\t$rpos{$rpos}\n";
				}
			print OUT "totalign\tNA\t$totAlign\n";	
			close(OUT)
			}
		}
	}
	
###########################
# behold the subroutines  #
###########################

sub errorFind {
	my ($in1p,$in2p) = @_;
	
	my $totAlign; my (%qual,%error,%pos,%rpos);
	
	my %yes;
	my @call = `wc $in1p`;
	my $num = $1 if $call[0] =~ m/^\s+(\d+)/;
	$num = $num/4;
	for (my $i = 0; $i < $sample; $i++) {
		$yes{int(rand()*$num)}++;
		}
	
	my $tracker = 1;
	open(IN1, "<$in1p");
	open(IN2, "<$in2p");
	while(<IN1>) {
		chomp(my $id1 = $_);
		chomp(my $seq1 = <IN1>);
		chomp(my $qualid1 = <IN1>);
		chomp(my $qual1 = <IN1>);
		
		chomp(my $id2 = <IN2>);
		chomp(my $seq2 = <IN2>);
		chomp(my $qualid2 = <IN2>);
		chomp(my $qual2 = <IN2>);
		
		if ($yes{$tracker}) {
			open(SEQ1, ">seq1.fa");
			open(SEQ2, ">seq2.fa");
			print SEQ1 ">seq1" . "\n" . $seq1 . "\n";
			print SEQ2 ">seq2" . "\n" . $seq2 . "\n";
			close(SEQ1); close(SEQ2);
			my @q1 = split(//,$qual1); my @q2 = split(//,$qual2); 
				
			my $call = system("$blat seq1.fa seq2.fa blat.out -noHead -dots=100 -out=axt");

			if (-s "blat.out" > 0) {
				open(IN, "<blat.out");
				chomp(my $d1 = <IN>); chomp(my $d2 = <IN>); chomp(my $d3 = <IN>);
				$totAlign += length($d2);
				my @d = split(/\s+/,$d1);
				my @d1 = split(//,$d2); my @d2 = split(//,$d3);
				for (my $i = 0; $i < scalar(@d1); $i++) {
					if ($d1[$i] ne $d2[$i]) {
						#this is an error
						if ($d[7] eq '-') {
							my $pos1 = $d[2] + $i; my $relPos1 = sprintf("%.2f",$pos1/length($seq1));
							my $pos2 = length($seq2) - ($d[5] + $i) + 1; my $relPos2 = sprintf("%.2f",$pos2/length($seq2));
							my $qual1 = ord($q1[$pos1 - 1])-33; my $qual2 = ord($q2[$pos2 - 1])-33;
							if ($qual1 < $qual2) {
								#error is in seq1
								$qual{$qual1}++;
								$pos{$pos1}++;
								$rpos{$relPos1}++;
								$error{$d1[$i]}{$d2[$i]}++;
								}
							else {
								#error is in seq2
								$qual{$qual2}++;
								$pos{$pos2}++;
								$rpos{$relPos2}++;
								$error{$d2[$i]}{$d1[$i]}++;
								} 
							}
						}
					}
				close(IN);	
				}	
			}
		$tracker++;	
		}	
	close(IN1); close(IN2);	
	unlink("blat.out"); unlink("seq1.fa"); unlink("seq2.fa");
	return(\%qual,\%error,\%pos,\%rpos,$totAlign);
	}

sub eraseAndZip {
	my ($out1p,$out2p,$out1u,$out2u,$file1,$file2) = @_;
	
	#need to erase all the trimmed files
	unlink($out1p); unlink($out2p); unlink($out1u); unlink($out2u);
 
	#need to compress my original files
	my $call1 = system("gzip -1 $file1");
	my $call2 = system("gzip -1 $file2");
	}

sub unzip {
	my ($file1, $file2) = @_;
	my $call1 = system("gunzip $file1");
	my $call2 = system("gunzip $file2");
	}


sub adaptorRemovalandTrim {
	my ($subdir, $file1, $file2,$ad) = @_;
	my %ad = %{$ad};

	my ($out1p,$out2p,$out1u,$out2u);	
	
	if ($file1 =~ m/([A-Z|0-9]+)_[1|2]\./i) {	
		#need to define the file names
		my $lib = $1;
		$out1p = $subdir . $lib . '_1p_trimmed.fastq';
		$out2p = $subdir . $lib . '_2p_trimmed.fastq';
		$out1u = $subdir . $lib . '_1u_trimmed.fastq';
		$out2u = $subdir . $lib . '_2u_trimmed.fastq';

		#need to print a file with the Illumina adaptors in it	
		my $ad = $dir . 'adaptors.fa';
		open(OUT, ">$ad");
		my @ad = ('adapt' . $lib{$lib}{'index'}, 'universal', 'adapt' . $lib{$lib}{'index'} . '_rc');
		foreach (@ad) {
			print OUT ">", $_, "\n", $ad{$_}, "\n";
			}
		close(OUT);

		my $call = system("java -classpath $trimmomatic org.usadellab.trimmomatic.TrimmomaticPE -phred33 $file1 $file2 $out1p $out1u $out2p $out2u ILLUMINACLIP:$ad:2:40:12 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:36");
		}
	unlink($ad);
	return($out1p,$out2p,$out1u,$out2u);
	}


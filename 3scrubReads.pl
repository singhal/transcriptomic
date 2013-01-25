##################################################################################
# a script to clean up reads (adaptor/duplicate/contamination removal), trimming #
# external dependencies: blat (v34), flash, trimmomatic                          #
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 7 Dec 2011            #
##################################################################################

use warnings;
use strict;

#home directory with all my files; make sure to end with a backslash
my $dir = '/media/DataDrive/sutureGenomics/';
my $tmpdir = '/media/DataDrive/tmp/';
#this file has all the info about my library, and creates a directory of all my files
my $lib = '/home/singhal/sutureGenomics/library';
#this file has a list of the adaptor sequences for my library, including the reverse complements.
my $adaptorFile = '/home/singhal/sutureGenomics/illuminaAdapters.fa';
#path to trimmomatic executable
my $trimmomatic = '/home/singhal/programs/trimmomatic/trimmomatic-0.17.jar';
#path to blat executable
my $blat = 'blat';
#indexed contaminant genome; here this includes ecoli and rRNA
my $contam = '/home/singhal/sutureGenomics/contaminants.fa';
my $readLength = 100;
my $nper = 0.6; #will get rid of reads for which more than $nper of bases are NNs
my $aper = 0.5; #will get rid of reads with any runs of bases longer than $aper*$readLength
#path to flash executable
my $flash = '/home/singhal/programs/flash/flash';

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
		unless (-f $subdir . $lib . '_1p_final.fastq.gz') {
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
			my $dup = $subdir . 'duplicates.out';
			duplicates($file1, $file2, $dup);
			my $time2 = int((time - $start2)/60);
			print "Found duplicates in $lib in $time2 minutes! Now this is something...\n";

			my $start3 = time;
			my $low = $subdir . 'lowComplexity.out';
			removeLowComplexity($file1,$low); removeLowComplexity($file2,$low);
			my $time3 = int((time - $start3)/60);
			print "Found low complexity reads in $lib in $time3 minutes! Whew, almost there!\n";

			my $start4 = time;
			my ($out1p,$out2p,$out1u,$out2u) = adaptorRemovalandTrim($subdir,$file1,$file2,\%ad);
			my $time4 = int((time - $start4)/60);
			print "Found adaptor and trimmed in $lib in $time4 minutes! Slowly but surely...\n";

			my $start5 = time;
			my $contaminants = $subdir . 'contam.out';
			removeContamination($out1p,$out2p,$out1u,$out2u,$contaminants,$blat);
			my $time5 = int((time - $start5)/60);
			print "Removed contamination in $lib in $time5 minutes! It's going...\n";

			my $start6 = time;
			my $unpaired = $subdir . $lib . '_u_trimmed.fastq';
			my $flashout = runFlash($out1p,$out2p,$flash);
			processFlash($flashout,$out1p,$out2p,$out1u,$out2u,$unpaired);
			my $time6 = int((time - $start6)/60);
			print "Combined reads in $lib in $time6 minutes. Next to last step!\n";
		
			my $start7 = time;
			my $final = makeFinalFiles($subdir,$out1p,$out2p,$unpaired,$dup,$low,$contaminants);
			my $time7 = int((time - $start7)/60);
			print "Made final libraries in $lib in $time7 minutes! Getting there!\n";

			my $start8 = time;
			eraseAndZip($out1p,$out2p,$unpaired,$final,$file1,$file2);
			my $time8 = int((time - $start8)/60);
			print "Zipped and cleaned up in $lib in $time8 minutes! Done!\n";
			}
		}
	}


###########################
# behold the subroutines  #
###########################

sub processFlash {
	my ($files, $file1, $file2, $file3, $file4,$unpaired) = @_;
	my @flashout = @{$files};	

	#first need to rename my not combined files
	my $call1 = system("mv $flashout[1] $file1");
	my $call2 = system("mv $flashout[2] $file2");

	#now need to combine all my single reads
	
	my $call3 = system("cat $file3 $file4 $flashout[0] > $unpaired");
	unlink($file3); unlink($file4); unlink($flashout[0]);
	}

sub runFlash {
	my ($file1,$file2,$flash) = @_;
	my $call = system("$flash $file1 $file2 -x 0.02");
	my @newfiles = ('extendedFrags.fastq', 'notCombined_1.fastq', 'notCombined_2.fastq');
	return(\@newfiles);	
	}

sub makeFinalFiles {
	my ($subdir,$out1p,$out2p,$unpaired,$dup,$low,$contam) = @_;

	my %junk;
	my @junk = ($dup,$low,$contam);
	foreach my $j (@junk) {
		open(IN, "<$j");
		while(<IN>) {
			$junk{$1}++ if $_ =~ m/(\S+)/
			}
		close(IN);
		}
	
	my @final;
	my @paired = ($out1p, $out2p);
	my %id;
	foreach my $read (@paired) {
		open(IN, "<$read");
		while(<IN>) {
			my $line = $_;
			if ($line =~ m/^\@(HS\S+)/) {
				$id{$1}++;
				}
			}
		}	
		
	my @reads = ($out1p, $out2p, $unpaired);
	foreach my $read (@reads) {
		my $id = $1 if $read =~ m/([A-Z|0-9|_]+)_trimmed/i;
		my $paired = $1 if $read =~ m/([1|2])p_/;
		my %fr = ('1' => 'F', '2' => 'R');	
		my $out = $subdir . $id . "_final.fastq";
		push(@final,$out);
		open(IN, "<$read"); open(OUT, ">$out");
		while(<IN>) {
			my $line = $_;
			if ($line =~ m/^\@(HS\S+)/) {
				my $readName = $1;
				my $seq = <IN>;		
				my $qualID = <IN>;
				my $qual = <IN>;		
				unless($junk{$readName}) {
					if ($paired) {
						if ($id{$readName} == 2) {
							print OUT '@' . $readName . '_' . $fr{$paired} . "\n" . $seq . $qualID . $qual;
							}
						}
					else {
						print OUT '@' . $readName . "\n" . $seq . $qualID . $qual;
						}				
					}
				}
			}
		close(IN); close(OUT);
		}

	return(\@final);
	}

sub eraseAndZip {
	my ($out1p,$out2p,$unpaired,$final,$file1,$file2) = @_;
	
	#need to erase all the trimmed files
	unlink($out1p); unlink($out2p); unlink($unpaired);
 
	#need to compress my original files
	my $call1 = system("gzip -1 $file1");
	my $call2 = system("gzip -1 $file2");

	#need to compress my final files
	foreach my $file (@{$final}) {
		my $call = system("gzip -1 $file");
		}
	}

sub unzip {
	my ($file1, $file2) = @_;
	my $call1 = system("gunzip $file1");
	my $call2 = system("gunzip $file2");
	}

sub duplicates {
	my ($file1,$file2,$dup) = @_;	
	my $sorted1 = sortFile($file1);
	my $sorted2 = sortFile($file2);	
	my $dub_ref = getDuplicates($sorted1);
	removeDuplicates($dub_ref, $sorted2, $dup);
	unlink($sorted1); unlink($sorted2);
	}

sub sortFile {
	my ($file) = @_;
	my $out = $file . "2";
	open(OUT, ">$out");
	my $sorted = $file . ".sorted";
	open(IN, "<$file");
	while(<IN>){
		chomp(my $line = $_);
		if ($line =~ m/^@/) {
			my $id = $1 if $line =~ m/^(\S+)/;
			chomp(my $seq = <IN>);
			chomp(my $qualID = <IN>);
			chomp(my $qual = <IN>);
			print OUT $id, "\t", $seq, "\t", $qualID, "\t", $qual, "\n";
			}
		}
	close(IN); close(OUT);
	my $call = system("sort -k 2,2 -T $tmpdir $out > $sorted");
	unlink($out);
	return($sorted);
	}

sub getDuplicates {
	my ($sorted) = @_;
	open(SORT, "<$sorted");
	my ($old_id, $old_seq);
	my $tracker = 1;
	my %dup;
	while(<SORT>){
		chomp(my $line = $_);
		my @d = split(/\t/, $line);
		if ($old_id) {
			unless ($old_seq eq $d[1]) {
				delete $dup{$tracker} if scalar(@{$dup{$tracker}}) < 2;
				$tracker++;			
				}
			}
		push(@{$dup{$tracker}}, $d[0]);
		$old_id = $d[0]; $old_seq = $d[1];
		}
	close(SORT);
	return(\%dup);
	}	

sub removeDuplicates {
	my ($dupref, $sorted2, $dup) = @_;
	my %dup = %$dupref;
	my ($old_id, $old_seq);	
	my %dup_rev;

	foreach my $tracker (keys %dup) {
		foreach my $id (@{$dup{$tracker}}) {
			$id = $1 if $id =~ m/^\@(\S+)/;
			$dup_rev{$id} = $tracker;
			}
		}

	open(SORT, "<$sorted2");
	open(DUP, ">$dup");
	while(<SORT>){
		chomp(my $line = $_);
		my @d = split(/\t/, $line);
		if ($old_id) {
			$old_id = $1 if $old_id =~ m/^\@(\S+)/; $d[0] = $1 if $d[0] =~ m/^\@(\S+)/;
			if ($dup_rev{$old_id} && $dup_rev{$d[0]}) {	
				#these are duplicates in reverse direction	
				if ($old_seq eq $d[1]) {		
					#these are duplicates in the forward direction too				
					if ($dup_rev{$old_id} eq $dup_rev{$d[0]}) {
						print DUP $d[0], "\n";
						}
					}
				}
			}
		$old_id = $d[0]; $old_seq = $d[1];
		}
	close(SORT); close(DUP); undef(%dup); undef(%dup_rev);
	}

sub removeLowComplexity {
	my ($file,$low) = @_;
	open(IN, "<$file");
	open(OUT, ">>$low");
	while(<IN>) {
		chomp(my $line = $_);		
		if ($line =~ m/^@(\S+)/) {
			my $id = $1;
			chomp(my $seq = <IN>);
			my $n = int($nper*length($seq));
			my $a = int($aper*length($seq));
			my $ncounter = ($seq =~ m/N/g);
			if ($seq =~ m/[A]{$a}/i || $seq =~ m/[T]{$a}/i || $seq =~ m/[G]{$a}/i || $seq =~ m/[C]{$a}/i || $ncounter >= $n) {
				print OUT $id, "\n";
				}
			}
		}	
	close(IN); close(OUT);	
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

sub removeContamination {
	my ($out1p,$out2p,$out1u,$out2u,$contaminants,$blat) = @_;
	
	my @reads = ($out1p, $out2p, $out1u, $out2u);
	my $fasta = "reads.fasta";
	my $out = "blatTemp.out";
	foreach my $read (@reads) {
		open(IN, "<$read");
		open(OUT, ">$fasta");
		while(<IN>) {
			chomp(my $line = $_);
    		if ($line =~ m/^\@(HS\S+)/) {
				my $id = $1;
				my $seq = <IN>;
				print OUT ">", "$1", "\n", $seq;
    			}
			}
		close(IN); close(OUT);	

		my $call1 = system("$blat $contam $fasta $out -noHead");
		my $call2 = system("cut -f 10 $out >> $contaminants");
		unlink($fasta); unlink($out);
		}
	}

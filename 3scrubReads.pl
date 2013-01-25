##################################################################################
# a script to clean up reads (adaptor/duplicate/contamination removal), trimming #
# external dependencies: blat (v34), flash, trimmomatic                          #
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 7 Dec 2011            #
# updated by Tyler Linderoth, tylerp.linderoth@gmail.com, 10 Dec 2012             # 								                 
##################################################################################

use warnings;
use strict;
use Getopt::Std;

my $version = '0.1.3'; # as of Dec. 6 2012

die (qq/
#########################
# 3scrubReads v$version #
#########################
Usage: 3scrubReads.pl options

Options:
-d	CHAR	directory containing fastq files 
-o	CHAR	output directory 
-t	CHAR	temporary directory
-a	FILE	file containing adaptor sequences for my library including reverse compliments
-g	CHAR	path to trimmoatic executable
-b	CHAR	path to blat executable
-c	CHAR	contaminant file (ex. indexed contaminant genome including with ecoli and rRNA)
-f	CHAR	path to flash executable
-l	INT	read length [100]
-n	FLOAT	will get rid of reads for which more than INT*read_length bases are NNs [0.6]
-r	FLOAT	will get rid of reads with any runs of bases longer than FLOAT*read_length [0.5]
-h	Y|N	perform hard trimming if -h is set to Y. No hard trimming if -h is set to N [Y]
-w	Y|N	warn if entire sequences are removed by hard trimming if -w is set to Y [N]    
-1	INT	number of base pairs to hard trim from 5' end of read 1 [0]
-2	INT	number of base pairs to hard trim from 3' end of read 1 [0]
-3	INT	number of base pairs to hard trim from 5' end of read 2 [0]
-4	INT	number of base pairs to hard trim from 3' end of read 2 [0]

Option Examples:
-d '\/Users\/sarahhykin\/Desktop\/Anolis_FF_data\/Sample_JMSH001B_index7\/Sample_1B_catted\/'
-o '\/Users\/sarahhykin\/Desktop\/Anolis_FF_data\/Sample_JMSH001B_index7\/Sample_1B_catted\/scrubbed\/'
-t '\/Users\/sarahhykin\/Desktop\/Anolis_FF_data\/Sample_JMSH001B_index7\/Sample_1B_catted\/tmp\/'
-a '\/Users\/sarahhykin\/Desktop\/Anolis_FF_data\/Sample_JMSH001B_index7\/Sample_1B_catted\/adaptorseqs_index7.fa'
-g '\/Users\/sarahhykin\/bin\/Trimmomatic-0.20\/trimmomatic-0.20.jar'
-b '\/Users\/sarahhykin\/bin\/blatSrc\/blat\/blat'
-c '\/Users\/sarahhykin\/Desktop\/Anolis_FF_data\/bin\/contaminant.fa'
-f '\/Users\/sarahhykin\/bin\/FLASH\/flash'

Notes:
-Hard trimming only works on files with illumina Casava 1.8+ header format:
 example: \@HS2:241:C0W6TACXX:4:1101:1130:2215 1:N:0:ACAGTG
-Hard trimming is performed AFTER trimmomatic run
-Assumes a library naming convention of ONLY letters & numbers with a _run[1|2].txt ending
-Home directory with all my files; make sure to end with a backslash
-Assumes the FLASH version released on 14 May 2012
\n/) if !@ARGV;

my %opts = (
d => undef,
o => undef,
t => undef,
a => undef,
g => undef,
b => undef,
c => undef,
f => undef,
l => 100,
n => 0.6,
r => 0.5,
h => 'Y',
w => 'N',
1 => 0,
2 => 0,
3 => 0,
4 => 0
);

getopts('d:o:t:a:g:b:c:f:l:n:r:h:w:1:2:3:4:', \%opts);

my $dir = $opts{d};
my $outdir = $opts{o};
my $tmpdir = $opts{t};
my $adaptorFile = $opts{a};
my $trimmomatic = $opts{g};
my $blat = $opts{b};
my $flash = $opts{f};
my $contam = $opts{c};
my $readLength = $opts{l};
my $nper = $opts{n};
my $aper = $opts{r};

###########################
# run the subroutines     #
# edit at your risk       #
###########################

mkdir($outdir) unless (-d $outdir);
my @files = <$dir*run1*fastq>;
	
foreach my $file1 (@files) {
	my $file2 = $file1;
	$file2 =~ s/run1/run2/;
	my $lib = $1 if $file1 =~ m/([A-Z|0-9]+)_run1/i;

	my $start2 = time;	
	my $dup = $file1 . '.duplicates.out';
	duplicates($file1, $file2, $dup);
	my $time2 = int((time - $start2)/60);
	print "Found duplicates in $lib in $time2 minutes! Now this is something...\n";

	my $start3 = time;
	my $low = $file1 . '.lowComplexity.out';
	removeLowComplexity($file1,$low); removeLowComplexity($file2,$low);
	my $time3 = int((time - $start3)/60);
	print "Found low complexity reads in $lib in $time3 minutes! Whew, almost there!\n";

	my $start4 = time;
	my ($out1p,$out2p,$out1u,$out2u) = adaptorRemovalandTrim($outdir,$file1,$file2,$adaptorFile);
	my $time4 = int((time - $start4)/60);
	print "Found adaptor and trimmed in $lib in $time4 minutes! Slowly but surely...\n";

	if ($opts{h} =~ /Y/i) {
	my $start4b = time;
	($out1p, $out2p, $out1u, $out2u) = hardTrim($out1p, $out2p, $out1u, $out2u, $opts{1}, $opts{2}, $opts{3}, $opts{4});
	my $time4b = int((time-$start4b)/60);
	print "Hard trimmed in $lib in $time4b\n";
	}
	
	my $start5 = time;
	my $contaminants =  $file1 . '.contam.out';
	removeContamination($out1p,$out2p,$out1u,$out2u,$contaminants,$blat);
	my $time5 = int((time - $start5)/60);
	print "Removed contamination in $lib in $time5 minutes! It's going...\n";

	my $start6 = time;
	my $unpaired = $outdir . $lib . '_u_trimmed.fastq';
	my $flashout = runFlash($out1p,$out2p,$flash);
	processFlash($flashout,$out1p,$out2p,$out1u,$out2u,$unpaired);
	my $time6 = int((time - $start6)/60);
	print "Combined reads in $lib in $time6 minutes. Next to last step!\n";
		
	my $start7 = time;
	my $final = makeFinalFiles($outdir,$out1p,$out2p,$unpaired,$dup,$low,$contaminants);
	my $time7 = int((time - $start7)/60);
	print "Made final libraries in $lib in $time7 minutes! Getting there!\n";
	}


###########################
# behold the subroutines  #
###########################

sub processFlash {
	my ($files, $file1, $file2, $file3, $file4,$unpaired) = @_;
	my @flashout = @{$files};	

	#first need to rename my not combined files
	my $call1 = system("mv $flashout[1] $file1") if (-f $flashout[1]);
	my $call2 = system("mv $flashout[2] $file2") if (-f $flashout[2]);

	#now need to combine all my single reads
	if (-f $flashout[0]) {	
		my $call3 = system("cat $file3 $file4 $flashout[0] > $unpaired");
		}
	else {
		my $call3 = system("cat $file3 $file4 > $unpaired");
		}
	
	unlink($file3); unlink($file4); unlink($flashout[0]);
	}

sub runFlash {
	my ($file1,$file2,$flash) = @_;
	my $call = system("$flash $file1 $file2 -x 0.35 -r 159 -f 200");
	my @newfiles = ('out.extendedFrags.fastq', 'out.notCombined_1.fastq', 'out.notCombined_2.fastq');
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
			if ($line =~ m/^\@(H\S+:\S+)/) {
				$id{$1}++;
				}
			}
		}	
		
	my @reads = ($out1p, $out2p, $unpaired);
	foreach my $read (@reads) {
		my $id = $1 if $read =~ m/([A-Z|0-9|_]+)_trimmed/i;
		my $paired = $1 if $read =~ m/([1|2])p/;
		my %fr = ('1' => 'F', '2' => 'R');	
		my $out = $subdir . $id . "_final.fastq";
		push(@final,$out);
		open(IN, "<$read"); open(OUT, ">$out");
		while(<IN>) {
			my $line = $_;
			if ($line =~ m/^\@(H\S+:\S+)/) {
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
			$id = $1 if $id =~ m/^\@(H\S+)/;
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

	my ($out1p,$out2p,$out1u,$out2u);	
	
	if ($file1 =~ m/([A-Z|0-9]+)_run1/i) {	
		#need to define the file names
		my $lib = $1;
		$out1p = $subdir . $lib . '_1p_trimmed.fastq';
		$out2p = $subdir . $lib . '_2p_trimmed.fastq';
		$out1u = $subdir . $lib . '_1u_trimmed.fastq';
		$out2u = $subdir . $lib . '_2u_trimmed.fastq';

		my $call = system("java -classpath $trimmomatic org.usadellab.trimmomatic.TrimmomaticPE -phred33 $file1 $file2 $out1p $out1u $out2p $out2u ILLUMINACLIP:$ad:2:40:12 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:36");
		}
	return($out1p,$out2p,$out1u,$out2u);
	}
	
sub hardTrim {
	my @fastqfiles = @_[0 .. 3]; #array w/ fastq files
	my ($f5, $f3, $r5, $r3) = @_[4 .. 7]; # number of bases to trim	
	my @newfiles;
	foreach (@fastqfiles) {
	# make sure file exists and has size != 0
	die ("$_ does not exist!\n") if (!-e $_);
	die ("$_ has size zero!\n") if (-z $_);
	my ($filename, $oldfile) = ($_, $_);
	# open files for reading/writing
		open(FASTQ, '<', $_);	
		$filename =~ s/trimmed\.fastq$/trimmedH\.fastq/; # new fastq file name
		push @newfiles, $filename; # store new file names in this array	
		open(OUT, '>', $filename);
		# check whether file contains read 1 or read 2 and set trimming parameters
		my $head;
		do {$head = <FASTQ> } until ($head =~ /^\@\S+:\d+:/);
		my $direction = $1 if ($head =~ /\s+([1|2]):[A-Z]:\d+:[A-Z]+\s*$/i);
		my $cut5 = ($direction == 1) ? $f5 : $r5;
		my $cut3 = ($direction == 1) ? $f3 : $r3;
		# process input fastq file	
		seek (FASTQ, 0, 0); # rewind fastq file before processing it
		while (<FASTQ>) {
			# ensure that you start at the header		
			if ($_ =~ /^\@\S+:\d+:/) {			
				$head = $_;
				# trim and print sequence
				my $seq;			
				do { $seq = <FASTQ> } until ($seq =~ /^[A-Z]+/i);
				# trim 5' end
				if ($seq =~ /^[A-Z]{$cut5}([A-Z]+)/i) {
					$seq = $1; 
				} else {
					# toss out read warning
					if ($opts{w} =~ /Y/i) {
						chomp $head;
						print STDERR "$head removed by hard trimming!\n";
					}
					next;
				}
				# trim 3' end 
				if ($seq =~ /^([A-Z]+)[A-Z]{$cut3}\s*$/i) {
					$seq = $1;
				} else {
					if ($opts{w} =~ /Y/i) {
						chomp $head;
						print STDERR "$head removed by hard trimming!\n";
					}		
					next;
				}   
				print OUT "$head$seq\n"; # print sequence
				# print optional sequence identifier
				my $id2;
				do { $id2 = <FASTQ> } until ($id2 =~ /^\+/);
				print OUT $id2;
				# trim and print quality scores
				my $qual;
				do { $qual = <FASTQ> } until ($qual =~ /^\S+/);
				$qual = $1 if $qual =~ /^\S{$cut5}(\S+)/; # trim 5' end
				$qual = $1 if $qual =~ /^(\S+)\S{$cut3}\s*$/;; # trim 3' end 	
				print OUT "$qual\n";		
			} else {
				next;
			} 
		}
		close FASTQ;
		close OUT; 
		unlink $oldfile; # delete previous version of fastq files
	}
	return @newfiles;
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
    		if ($line =~ m/^\@(H\S+)/) {
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

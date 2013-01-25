##################################################################################
# a script to prepare sequence files to use for assembly		                 #
# external dependencies: none													 #
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 13 Dec 2011           #
##################################################################################

use warnings;
use strict;

#home directory with all my files; make sure to end with a backslash
my $dir = '/media/DataDrive/sutureGenomics/';
#this file has all the info about my library, and creates a directory of all my files
my $lib = '/home/singhal/sutureGenomics/library';
#how many libraries do you have per lineage? will not make these files unless all the libraries are there.
my $numLib = 5; 

#defines the library
open(IN, "<$lib");
my %lib;
while (<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/, $line);
	my $lib = $d[2] . '/' . $d[0] . '/';
	$lib{$d[2]}{$d[0]} = $lib; 
	}
close(IN);

###########################
# run the subroutines     #
###########################

foreach my $lib (keys %lib) {
	my (@p1, @p2, @un);
	foreach my $name (keys %{$lib{$lib}}) {
		my $subdir = $dir . $lib{$lib}{$name};
		if (-d $subdir) {
			push(@p1, $subdir . $name . '_1p_final.fastq.gz');
			push(@p2, $subdir . $name . '_2p_final.fastq.gz');
			push(@un, $subdir . $name . '_u_final.fastq.gz');
			}
		}	
	if (scalar(@p1) == $numLib) {
		print "Making files in $lib!\n";
		makeCombo(\@p1,\@p2,\@un,$lib);
		interlace(\@p1,\@p2,\@un,$lib);
		makeUnpaired(\@un,$lib);
		makeAll($lib);
		}
	}

###########################
# behold the subroutines! #
###########################

sub makeAll {
	my ($lib) = @_;

	my $seq1 = $dir . $lib . '/' . $lib . '_1.fastq.gz';
	my $seq2 = $dir . $lib . '/' . $lib . '_2.fastq.gz';
	my $seq3 = $dir . $lib . '/' . $lib . '_u.fastq.gz';
	my $all = $dir . $lib . '/' . $lib . '.fastq.gz';
	
	my $call1 = system("gzip -dc -1 $seq1 $seq2 $seq3 | gzip -c -1 > $all") unless (-f $all);	
	}

sub makeUnpaired {
	my ($un,$lib) = @_;
	
	my $seq3 = $dir . $lib . '/' . $lib . '_u.fastq';
	my $seq3_short = $dir . $lib . '/' . $lib . '_u_short.fastq';

	unless(-f $seq3) {
		open(OUT, ">$seq3");
		open(OUTSHORT, ">$seq3_short");
		foreach my $seqgz (@{$un}) {
			my $call = system("gunzip $seqgz");
			my $seqfile = $1 if $seqgz =~ m/(\S+)\.gz/;			
			open(IN, "<$seqfile");
			while(<IN>) {
				chomp(my $id = $_);
				chomp(my $seq = <IN>);
				chomp(my $qualid = <IN>);
				chomp(my $qual = <IN>); 
				print OUT $id . "\n" . $seq . "\n" . $qualid . "\n" . $qual . "\n";
				if (length($seq) > 126) {
					$seq = substr $seq, 0, 126;
					$qual = substr $qual, 0, 126;
					}
				print OUTSHORT $id . "\n" . $seq . "\n" . $qualid . "\n" . $qual . "\n";
				}			
			close(IN);
			my $call2 = system("gzip -1 $seqfile");
			}
		close(OUT); close(OUTSHORT);
		}

	my $call1 = system("gzip -1 $seq3");
	my $call2 = system("gzip -1 $seq3_short");
	}

#subroutine modified from shuffleSequences_fastq.pl provided by Zerbino with Velvet
sub interlace {
	my ($p1,$p2,$un,$lib) = @_;

	my @p1 = @{$p1};
	my @p2 = @{$p2};
	my $interGZ = $dir . $lib . '/' . $lib . '_p.fastq.gz';	
	my $inter   = $dir . $lib . '/' . $lib . '_p.fastq';	
	
	unless(-f $interGZ) {
		open(OUT, ">$inter");
		for (my $i = 0; $i < scalar(@p1); $i++) {
			my $fAgz = $p1[$i];
			my $fBgz = $p2[$i];

			my $fA = $1 if $fAgz =~ m/(\S+)\.gz/;
			my $fB = $1 if $fBgz =~ m/(\S+)\.gz/;

			my $call1 = system("gunzip $fAgz"); 
			my $call2 = system("gunzip $fBgz");

			open(FA, "<$fA");
			open(FB, "<$fB");

			while(<FA>) {
				print OUT $_;
				$_ = <FA>;
				print OUT $_; 
				$_ = <FA>;
				print OUT $_; 
				$_ = <FA>;
				print OUT $_; 
	
				$_ = <FB>;
				print OUT $_; 
				$_ = <FB>;
				print OUT $_;
				$_ = <FB>;
				print OUT $_;
				$_ = <FB>;
				print OUT $_;
				}
			close(FA); close(FB);
			my $call3 = system("gzip -1 $fA");
			my $call4 = system("gzip -1 $fB");
			}
		close(OUT);
		my $call5 = system("gzip -1 $inter");
		}
	}

sub makeCombo {
	my ($p1,$p2,$un,$lib) = @_;

	my $files1 = join(" ", @{$p1}); 
	my $files2 = join(" ", @{$p2}); 
	my $files3 = join(" ", @{$un});

	my $seq1 = $dir . $lib . '/' . $lib . '_1.fastq.gz';
	my $seq2 = $dir . $lib . '/' . $lib . '_2.fastq.gz';
	
	my $call1 = system("gzip -dc -1 $files1 | gzip -c -1 > $seq1") unless (-f $seq1);	
	my $call2 = system("gzip -dc -1 $files2 | gzip -c -1 > $seq2") unless (-f $seq2);	
	}

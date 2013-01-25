use warnings;
use strict;

#home directory with all my files; make sure to end with a backslash
my $dir = '/media/DataDrive/sutureGenomics/';
#this file has all the info about my library, and creates a directory of all my files
my $lib = $dir . 'library';

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

my $out = $dir . 'qualityControl.out';
open(OUT, ">$out");
print OUT "lib\torigReads\torigBases\tfinalReads\tfinalBases\tpairedReads\tpairedBases\tcontam\tdupl\tlow\n";
foreach my $lib (sort {$a cmp $b} keys %lib) {
	my $subdir = $dir . $lib{$lib}{'dir'};
	if (-d $subdir) {
	    print "Doing $lib now!\n";
    	my $file1gz = $subdir . $lib . '_1.fastq.gz';        
        my $file1pgz = $subdir . $lib . '_1p_final.fastq.gz';  
        my $file2pgz = $subdir . $lib . '_2p_final.fastq.gz';  
        my $fileugz = $subdir . $lib . '_u_final.fastq.gz';
        my $contam = $subdir . 'contam.out';
        my $dupl = $subdir . 'duplicates.out';
        my $low = $subdir . 'lowComplexity.out';
        
        #print out: orig num reads, orig num bases, final num reads, final num bases, final num pairs, contam reads, dupl reads, low complex reads
        
        my $file1 = unzip($file1gz);
        my $file1p = unzip($file1pgz);
        my $file2p = unzip($file2pgz);
        my $fileu = unzip($fileugz);
        
        my $origReads = wc($file1) / 2;
     	my $origBases = $origReads * 100;
     	
     	my ($paired1Reads, $paired1Bases) = final($file1p); 
     	my ($paired2Reads, $paired2Bases) = final($file2p); 
     	my ($unpairedReads, $unpairedBases) = final($fileu);
     	my $finalReads = $paired1Reads + $paired2Reads + $unpairedReads;
     	my $finalBases = $paired1Bases + $paired2Bases + $unpairedBases;
     	my $pairedReads = $paired1Reads + $paired2Reads;
     	my $pairedBases = $paired1Bases + $paired2Bases;
     	
     	my $contamWc = wc($contam) * 2;
     	my $duplWc = wc($dupl) * 2;
     	my $lowWc = wc($low) * 2;
     	
     	print OUT $lib, "\t", $origReads, "\t", $origBases, "\t", $finalReads, "\t", $finalBases, "\t", $pairedReads, "\t", $pairedBases, "\t", $contamWc, "\t", $duplWc, "\t", $lowWc, "\n";

	my @files=($file1,$file1p,$file2p,$fileu);    
	    zip(\@files);
     	}
     }	
close(OUT); 

sub zip {
    my ($files) = @_;
    my @f = @{$files};
    foreach my $file (@f) {
	my $call = system("gzip -1 $file");
    }
}

sub final {
	my ($file) = @_;
	my ($reads, $bases);
	open(IN, "<$file");
	while(<IN>) {
		chomp(my $seqID = $_);
        chomp(my $seq = <IN>);
    	chomp(my $qualid = <IN>);
        chomp(my $qual = <IN>);
		$reads++;
		$bases += length($seq);
		}
	close(IN);
	return($reads,$bases);
	}
		
sub wc {
	my ($file) = @_;
	my @call = `wc $file`;
	my $count = $1 if $call[0] =~ m/^\s+(\d+)/;
	return($count);
	}		
		
sub unzip {
	my ($filegz) = @_;
	my $call = system("gunzip $filegz");
	my $out = $1 if $filegz =~ m/(.*)\.gz/;
	return($out);
	}		
     

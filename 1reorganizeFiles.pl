############################################################################
# a script to take the files as they come off the sequencer and reorganize #
# external dependencies: none                                              # 
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 5 Dec 2011      #
############################################################################

use warnings;
use strict;

#home directory with all my files
my $dir = '/media/DataDrive/sutureGenomics/';

#this file has all the info about my library
my $lib = '/home/singhal/sutureGenomics/library';

#create hash with library info
open(IN, "<$lib");
while (<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/, $line);

	#the name as it comes off the machine
	my $mach_name = $1 . '_index' . $d[1] if $d[3] =~ m/(\S{4})/;

	#the directory for the lineage
	my $lin_dir = $dir . $d[2] . '/';
	mkdir($lin_dir) unless(-d $lin_dir);	

	#the directory for the individual
	my $old_dir = $dir . $mach_name;
	my $new_dir = $lin_dir . $d[0] . '/';
	#move the files from the machine-named directory to this one
	my $call1 = system("mv $old_dir $new_dir") if (-d $old_dir);

	if (-d $new_dir) {
		my $fwd_seq = $new_dir . $d[0] . '_1.fastq.gz';
		my $rev_seq = $new_dir . $d[0] . '_2.fastq.gz';

		#combines files and removes individual files
		unless (-f $fwd_seq) {
			my $call2 = system("gzip -1 -dc $new_dir*R1* | gzip -1 -c > $fwd_seq");
			my $call3 = system("rm $new_dir*R1*");
			}
		unless (-f $rev_seq) {
			my $call4 = system("gzip -1 -dc $new_dir*R2* | gzip -1 -c > $rev_seq");
			my $call5 = system("rm $new_dir*R2*");
			}		
		}
	}
close(IN);

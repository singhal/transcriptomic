############################################################################
# a script to take raw reads and evaluate their quality			           #
# external dependencies: FastQC                                            # 
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 5 Dec 2011      #
############################################################################

use warnings;
use strict;

#home directory with all my files
my $dir = '/media/DataDrive/sutureGenomics/';
#how many threads can you dedicate to this adventure?
my $thread = '2'; 
my $fastqc_dir = '/home/singhal/programs/FastQC/';

my @dir = <$dir*>;
foreach my $dir1 (@dir) {
	$dir1 = $dir1 . '/';
	my @subdir = <$dir1*>;
	foreach my $dir2 (@subdir) {
		$dir2 = $dir2 . '/';
		my @files = <$dir2*>;
		foreach my $file (@files) {
			if ($file =~ m/_[1|2].fastq.gz/) {
				my $qual = $dir2 . 'eval';
				unless (-d $qual) {
					my $call = system($fastqc_dir . "fastqc -t $thread $file");
					}
				}
			}
		}
	}


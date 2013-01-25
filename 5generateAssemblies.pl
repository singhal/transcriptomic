##############################################################################################################################
# a script to prepare shell scripts to use for assembly with multiple de novo assemblers on TACC                             #
# external dependencies: to run the script, none but ... to use the scripts: abyss, velvet, trinity, oases, soapdenovo-trans #
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 13 Dec 2011                                                       #
##############################################################################################################################

use warnings;
use strict;

#home directory with all my files; make sure to end with a backslash
my $scratch = '/scratch/01302/ssinghal/';
my $home = '/work/01302/ssinghal/';
my $psc_scratch = '/brashear/ssinghal/';
my $email = 'sonal.singhal1@gmail.com';
my $insert = 30; #how long is the expected distance between the two reads

#parameters for soap
my @soapkmer = qw(21 31 41 51 61 71 81 91); #max 96
my @soapcov = qw(2 10);
my $soap_dir = '/usr/users/3/ssinghal/soap';

#parameters for abyss
my @abysskmer = qw(21 31 41 51 61 71 81 91); #max 96
my @abysscov = qw(2 6 10);

#parameters for velvet
my @velvetkmer = qw(31 41 51 61 71 81 91); #max 96
my $velvetdir = '/home/01302/ssinghal/bin/'; #directory where velvet executables are

my @lib = qw(Carlia_N Carlia_S Lampro_N Lampro_C Lampro_S Sapro_C Sapro_S);

###########################
# run the subroutines     #
###########################

foreach my $lib (@lib) {
	my $file_1 = $home . $lib . '_1.fastq';
	my $file_2 = $home . $lib . '_2.fastq';
	my $file_u = $home . $lib . '_u.fastq';
	my $file_1s = $scratch . $lib . '_1.fastq';
	my $file_2s = $scratch . $lib . '_2.fastq';
	my $file_us = $scratch . $lib . '_u.fastq';
	if (-f $file_1) {
		makeABYSS($file_1,$file_2,$file_u,$lib);
		}
	if (-f $file_us) {
		makeSOAP($file_1s,$file_2s,$file_us,$lib);
		}
	$file_u = $scratch . $lib . '_u.fastq';
	my $file = $scratch . $lib . '.fastq';
	if (-f $file) {		
		makeVelvet($file,$file_u,$lib);
		}
	my $file2 = $psc_scratch . $lib . '.fastq';
	if (-f $file2) {
		makeTrinity($file2,$lib);
		}
	}

###########################
# behold the subroutines! #
###########################

sub makeTrinity {
	my ($file,$lib) = @_;

	my $resultsDir =  $psc_scratch . 'trinityResults';
	my $runfileDir =  $psc_scratch . 'trinityScripts';
	mkdir($resultsDir) unless (-d $resultsDir);
	mkdir($runfileDir) unless (-d $runfileDir);
	my $libResults = $resultsDir . '/' . $lib;
	mkdir($libResults) unless (-d $libResults);

	my $job = 'trinity' . "_" . $lib;
	open(OUT, ">$runfileDir" . "/" . $job . ".sh");
	
	print OUT "#!/bin/bash\n";
	print OUT "#PBS -l ncpus=32\n";
	print OUT "#PBS -l walltime=96:00:00\n";
	print OUT "#PBS -j oe\n";
	print OUT "#PBS -q batch\n";
	print OUT "#PBS -m e\n";
	print OUT "#PBS -M $email\n";
	print OUT "#PBS -N $job\n";

	print OUT "set -x\n";
	print OUT "source /usr/share/modules/init/bash\n";

	print OUT "module load trinity/r2011-11-26\n";
	print OUT "module load java\n";

	print OUT "ulimit -u unlimited\n";
	
	print OUT "ja\n";

	print OUT "cd $libResults\n";
	print OUT "export OMP_NUM_THREADS=10\n";

	print OUT "Trinity.pl --seqType fq --single $file --output \$SCRATCH_RAMDISK --CPU 10 --bflyHeapSpace 10G > $libResults" . "/trinity.log\n";
	print OUT "cp -r \$SCRATCH_RAMDISK/Trinity.fasta $libResults\n";
	print OUT "cp -r \$SCRATCH_RAMDISK/\* $libResults\n";
	print OUT "ja -chlst\n";
	close(OUT);
	}


sub makeSOAP {
	my ($file_1,$file_2,$file_u,$lib) = @_;

	my $resultsDir =  $scratch . 'soapResults';
	my $runfileDir =  $scratch . 'soapScripts';
	mkdir($resultsDir) unless (-d $resultsDir);
	mkdir($runfileDir) unless (-d $runfileDir);
	foreach my $k (@soapkmer){
		foreach my $cov (@soapcov) {
			my $job = 'soap' . "_" . $lib . "_k" . $k . "_c" . $cov;
			open(CFG, ">$runfileDir" . "/" . $job . ".config");
			open(OUT, ">$runfileDir" . "/" . $job . ".sh");

			print OUT "#!/bin/bash\n"; 
			print OUT "#\$ -N $job\n"; 
			print OUT "#\$ -j y\n"; 
			print OUT "#\$ -P hpc\n"; 
			print OUT "#\$ -o $resultsDir" , "/", $job, ".out\n"; 
			print OUT "#\$ -pe 1way 16\n"; 
			print OUT "#\$ -q largemem\n"; 
			print OUT "#\$ -l h_rt=08:00:00\n"; 
			print OUT "#\$ -M $email\n"; 
			print OUT "#\$ -m e\n";
			print OUT "#\$ -cwd\n"; 
			print OUT "#\$ -V\n";
			print OUT $velvetdir . "SOAPdenovo-Trans-127mer all -s $runfileDir" . "/" . $job . ".config -K $k -o $resultsDir" . "/" . "$job -p 16 -d $cov -D $cov\n";		
			
			print CFG "max_rd_len=200\n";	
			print CFG "\[LIB\]\n";
			print CFG "avg_ins=$insert\n";
			print CFG "asm_flags=3\n";
			print CFG "rank=1\n";
			print CFG "q1=$file_1\n";
			print CFG "q2=$file_2\n";
			print CFG "q=$file_u\n";

			close(CFG); close(OUT);
			}
		}
    }

sub makeVelvet {
	my ($file,$file_u,$lib) = @_;
	
	my $runfileDir = $scratch . 'velvetScripts/';
	mkdir($runfileDir) unless (-d $runfileDir);
	my $allresults = $scratch . 'velvetResults/';
	mkdir($allresults) unless (-d $allresults);

	foreach my $k (@velvetkmer) {
		my $jobH = 'velveth' . "_" . $lib . "_k" . $k; 
		my $jobG = 'velvetg' . "_" . $lib . "_k" . $k; 
		my $jobO = 'oases' . "_" . $lib . "_k" . $k; 
		my $merge = 'mergeVelvet' . "_" . $lib; 

		my $resultsDir = $allresults . 'velvet_' . $lib . "_k" . $k; 
		my $finalResults = $allresults . 'velvet_' . $lib;
		mkdir($resultsDir) unless (-d $resultsDir);
		mkdir($finalResults) unless (-d $finalResults);
		
		open(OUTH, ">$runfileDir" . $jobH . ".sh");
		open(OUTG, ">$runfileDir" . $jobG . ".sh");
		open(OUTO, ">$runfileDir" . $jobO . ".sh");
		open(MERGE, ">$runfileDir" . $merge . ".sh");

		print OUTH "#!/bin/bash\n"; 
		print OUTH "#\$ -N $jobH\n"; 
		print OUTH "#\$ -j y\n"; 
		print OUTH "#\$ -P hpc\n"; 
		print OUTH "#\$ -o $resultsDir" , "/", $jobH, ".out\n"; 
		print OUTH "#\$ -pe 1way 8\n"; 
		print OUTH "#\$ -q largemem\n"; 
		print OUTH "#\$ -l h_rt=08:00:00\n"; 
		print OUTH "#\$ -M $email\n"; 
		print OUTH "#\$ -m e\n";
		print OUTH "#\$ -cwd\n"; 
		print OUTH "#\$ -V\n";
		print OUTH "$velvetdir" . 'velveth ' . $resultsDir . " $k " . '-shortPaired -fastq ' . $file . ' -short2 ' . $file_u . "\n";

		print OUTG "#!/bin/bash\n"; 
		print OUTG "#\$ -N $jobG\n"; 
		print OUTG "#\$ -j y\n"; 
		print OUTG "#\$ -P hpc\n"; 
		print OUTG "#\$ -o $resultsDir" , "/", $jobG, ".out\n"; 
		print OUTG "#\$ -pe 1way 8\n"; 
		print OUTG "#\$ -q largemem \n"; 
		print OUTG "#\$ -l h_rt=08:00:00\n"; 
		print OUTG "#\$ -M $email\n"; 
		print OUTG "#\$ -m e\n";
		print OUTG "#\$ -cwd\n"; 
		print OUTG "#\$ -V\n";
		print OUTG "$velvetdir" . 'velvetg ' . $resultsDir . ' -read_trkg yes -ins_length ' . $insert .  "\n";

		print OUTO "#!/bin/bash\n"; 
		print OUTO "#\$ -N $jobO\n"; 
		print OUTO "#\$ -j y\n"; 
		print OUTO "#\$ -P hpc\n"; 
		print OUTO "#\$ -o $resultsDir" , "/", $jobO, ".out\n"; 
		print OUTO "#\$ -pe 1way 8\n"; 
		print OUTO "#\$ -q largemem \n"; 
		print OUTO "#\$ -l h_rt=08:00:00\n"; 
		print OUTO "#\$ -M $email\n"; 
		print OUTO "#\$ -m e\n";
		print OUTO "#\$ -cwd\n"; 
		print OUTO "#\$ -V\n";
		print OUTO "$velvetdir" . 'oases ' . $resultsDir . ' -ins_length ' . $insert .  "\n";

		print MERGE "#!/bin/bash\n"; 
		print MERGE "#\$ -N $merge\n"; 
		print MERGE "#\$ -j y\n"; 
		print MERGE "#\$ -P hpc\n"; 
		print MERGE "#\$ -o $finalResults" , "/", $merge, ".out\n"; 
		print MERGE "#\$ -pe 1way 8\n"; 
		print MERGE "#\$ -q largemem \n"; 
		print MERGE "#\$ -l h_rt=24:00:00\n"; 
		print MERGE "#\$ -M $email\n"; 
		print MERGE "#\$ -m e\n";
		print MERGE "#\$ -cwd\n"; 
		print MERGE "#\$ -V\n";
		print MERGE "$velvetdir" . "velveth $finalResults 27 -long " . $finalResults . '_*/transcripts.fa' .  "\n";
		print MERGE "$velvetdir" . "velvetg $finalResults -read_trkg yes -conserveLong yes" . "\n";
		print MERGE "$velvetdir" . "oases $finalResults -merge" . "\n";

		close(OUTG); close(OUTH); close(OUTO); close(MERGE);
		}
	}
	
sub makeABYSS {
	my ($file_1,$file_2,$file_u,$lib) = @_;
	my $resultsDir =  $home . 'abyssResults';
	my $runfileDir =  $home . 'abyssScripts';
	mkdir($resultsDir) unless (-d $resultsDir);
	mkdir($runfileDir) unless (-d $runfileDir);
	foreach my $k (@abysskmer) {			
	    my $job = 'abyss' . "_" . $lib . "_k" . $k;
	    my $dir = $resultsDir . '/' . $job;
	    mkdir($dir) unless (-d $dir);
	    open(OUT, ">/$runfileDir" . "/" . $job . ".sh");
	    print OUT "#!/bin/bash\n"; 
	    print OUT "#\$ -N $job\n"; 
	    print OUT "#\$ -j y\n"; 
	    print OUT "#\$ -o $dir" , "/", $job, ".out\n"; 
	    print OUT "#\$ -pe 16way 128\n"; 
	    print OUT "#\$ -q normal\n"; 
	    print OUT "#\$ -l h_rt=20:00:00\n"; 
	    print OUT "#\$ -M $email\n"; 
	    print OUT "#\$ -m e\n";
	    print OUT "#\$ -cwd\n"; 
	    print OUT "#\$ -V\n";
	    print OUT "cd $dir\n";
	    print OUT "abyss-pe k=$k mpirun=/opt/apps/pgi7_2/openmpi/1.3/bin/mpirun np=32 n=5 s=200 in=\'$file_1 $file_2\' se=$file_u name=$dir/$job\n";
		}
	}


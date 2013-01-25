use warnings;
use strict;

use Bio::DB::Fasta;

my @contacts = (["Lampro_N","Lampro_C","Lampro_S","Sapro_C","Sapro_S"]);
my $dir = '/media/DataDrive/sutureGenomics/';
my $finaldir = $dir . 'lineageSNPs/';
my $tmpdir = $dir . 'tmp/';
my $lib = $dir . 'library';
my $np = 2;

#defines the library
my %lib;
open(IN, "<$lib");
while (<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/, $line);
	my $lib = $d[2] . '/' . $d[0] . '/';
	$lib{$d[2]}{$d[0]} = {'dir' => $lib, 'index' => $d[1]};	
	} 
close(IN);

mkdir($finaldir) unless(-d $finaldir);
mkdir($tmpdir) unless(-d $tmpdir);
#loop-de-loop!
for (my $i = 0; $i < scalar(@contacts); $i++) {
	
	my $name = 'allLampro';
 	
	#make file with just annotated transcripts
	my $seq = makeSeq($contacts[$i][0]);
	
	#map reads back to transcripts
	my $call1 = system("bowtie2-build -q $seq $seq") unless(-f $seq . ".1.bt2");
	#contact1
	my $out1 = mapping($seq,\%lib,$contacts[$i][0],$dir);
	#contact2
	my $out2 = mapping($seq,\%lib,$contacts[$i][1],$dir);
	#contact3 
	my $out3 = mapping($seq,\%lib,$contacts[$i][2],$dir);
	#outgroups
	my $out4 = mapping($seq,\%lib,$contacts[$i][3],$dir);
	my $out5 = mapping($seq,\%lib,$contacts[$i][4],$dir);

	#get depth from contacts and mpileup for the last
	my $pile1 = runDepth($seq,$dir,$contacts[$i][0],$out1,$finaldir);
	my $pile2 = runDepth($seq,$dir,$contacts[$i][1],$out2,$finaldir);
	my $pile3 = runDepth($seq,$dir,$contacts[$i][2],$out3,$finaldir);
	my $pile4 = runPileup($seq,$dir,'saproOutgroup',$out4,$out5,$finaldir);

	#get sophisticated snps from ingroups
	my $vcf1 = sophSNP($seq,$dir,'jointLampro',$out1,$out2,$out3,$finaldir);
	}
rmdir($tmpdir);	
	
###########################
# run the subroutines     #
###########################		

sub mapping {
	my ($seq, $libHash, $contact, $dir) = @_;
	
	my %lib = %{$libHash};
	my @outfiles;
	
    foreach my $lib (sort {$a cmp $b} keys %{$lib{$contact}}) {		
		my $subdir = $dir . $lib{$contact}{$lib}{'dir'};
	    my $file1gz = $subdir . $lib . '_1p_final.fastq.gz';
	    my $file2gz = $subdir . $lib . '_2p_final.fastq.gz';
	    my $fileugz = $subdir . $lib . '_u_final.fastq.gz';

		my $out = $subdir . $lib . ".sorted.bam";
		unless(-f $out) {
		    my $file1 = unzip($file1gz);
		    my $file2 = unzip($file2gz);
		    my $fileu = unzip($fileugz);
			my $out = runMapping($seq,$file1,$file2,$fileu,$subdir,$lib,$dir) unless (-f $out);
			my @files = ($file1,$file2,$fileu);
			zip(\@files);
			}
		push(@outfiles,$out);	
		}
	
	return(\@outfiles);
	}

sub makeSeq {
	my ($contact) = @_;
	my $file = $dir . "seqfiles/" . $contact . "_trinity.fa.final.annotated";
	my $out = $dir . "seqfiles/" . $contact . "_annotated.fa";
	unless(-f $out) {
		open(IN, "<$file");
		open(OUT, ">$out");
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/>/) {
				chomp(my $seq = <IN>);
				if ($line =~ m/ENSA/) {
					print OUT $line, "\n", $seq, "\n";
					}
				}
			}
		close(OUT); close(IN);	
		}
	return($out);
	}

sub runPileup {	
	my ($seq,$dir,$contact,$outfiles1,$outfiles2,$finaldir) = @_;
	my @out1 = @{$outfiles1};
	my @out2 = @{$outfiles2};
	my @out = (@out1,@out2);
	my $call1 = system("samtools faidx $seq") unless(-f $seq . ".fai");
	my $out = $finaldir . $contact . ".mpileup.out";
	unless (-f $out) {
		my $files = join("\t", @out);
		my $call2 = system("samtools mpileup -d 20000 -A -I -B -f $seq $files > $out");	
		}
	return($out);
	}


sub runDepth {	
	my ($seq,$dir,$contact,$outfiles,$finaldir) = @_;
	my @out = @{$outfiles};
	my $call1 = system("samtools faidx $seq") unless(-f $seq . ".fai");
	my $out = $finaldir . $contact . ".depth.out";
	unless (-f $out) {
		my $files = join("\t", @out);
		my $call2 = system("samtools depth $files > $out");	
		}
	return($out);
	}

sub sophSNP {	
	my ($seq,$dir,$contact,$outfiles1,$outfiles2,$outfiles3,$finaldir) = @_;
	my @out1 = @{$outfiles1};
	my @out2 = @{$outfiles2};
	my @out3 = @{$outfiles3};
	my @out = (@out1,@out2,@out3);
	my $call1 = system("samtools faidx $seq") unless(-f $seq . ".fai");
	my $out = $finaldir . $contact . ".vcf.out";
	unless (-f $out) {
		my $files = join("\t", @out);
		my $call2 = system("samtools mpileup -d 20000 -A -E -I -uf $seq $files | bcftools view -bvcg - > var_out.raw.bcf");
		my $call3 = system("bcftools view var_out.raw.bcf | vcfutilsTL02.pl varFilter -D 30000 -w 0 -e 0.00000001 > $out");
		}
	return($out);
	}
	
sub runMapping {
	my ($seq,$file1,$file2,$fileu,$subdir,$lib,$dir) = @_;
	my $out = $subdir . $lib . ".sorted";
	print $out, "\n";
	my $call1 = system("bowtie2 -x $seq -1 $file1 -2 $file2 -S bowtie1_out.sam -5 5 -3 5 --sensitive -k 10 -X 300 -p $np");
	my $call2 = system("bowtie2 -x $seq $fileu -S bowtie2_out.sam -5 5 -3 5 --sensitive -k 10  -p $np");
	my $call3 = system("samtools view -bS bowtie1_out.sam > bowtie1_out.bam");
	my $call4 = system("samtools view -bS bowtie2_out.sam > bowtie2_out.bam");
	my $call5 = system("samtools merge bowtie_out.bam bowtie1_out.bam bowtie2_out.bam");
	my $call6 = system("samtools sort bowtie_out.bam $out");
	my $call7 = system("rm bowtie1_out.sam bowtie2_out.sam bowtie_out.bam bowtie1_out.bam bowtie2_out.bam");	
	return($out);
	}

sub zip {
    my ($file) = @_;
    foreach my $_ (@$file) {
		my $call = system("gzip -1 $_");
    	}
	}

sub unzip {
    my ($file) = @_;
    my $unzip = $1 if $file =~ m/(\S+)\.gz/;
    my $call = system("gunzip $file");
    return($unzip);
	}	

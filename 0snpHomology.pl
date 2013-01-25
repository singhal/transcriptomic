use warnings;
use strict;

my %compare = (
'CarliaN_CarliaS' => {'seq1' => $carlia_n, 'seq2' => $carlia_s },
'LamproN_LamproS' => {'seq1' => $lampro_n, 'seq2' => $lampro_s },
'LamproC_LamproS' => {'seq1' => $lampro_c, 'seq2' => $lampro_s },
'SaproC_SaproS' => {'seq1' => $sapro_c, 'seq2' => $sapro_s },
'CarliaN_LamproN' => {'seq1' => $carlia_n, 'seq2' => $lampro_n },
'CarliaN_SaproC' => {'seq1' => $carlia_n, 'seq2' => $sapro_c },
'SaproC_LamproN' => {'seq1' => $sapro_c, 'seq2' => $lampro_n }
);

my $reads_dir '/media/DataDrive/sutureGenomics/';
my $main_dir = '/Users/singhal/Desktop/genomics/';
my $carlia_n = $main_dir . 'Carlia_N_trinity.fa.final.annotated';
my $carlia_s = $main_dir . 'Carlia_S_trinity.fa.final.annotated';
my $lampro_n = $main_dir . 'Lampro_N_trinity.fa.final.annotated';
my $lampro_c = $main_dir . 'Lampro_C_trinity.fa.final.annotated';
my $lampro_s = $main_dir . 'Lampro_S_trinity.fa.final.annotated';
my $sapro_c = $main_dir . 'Sapro_C_trinity.fa.final.annotated';
my $sapro_s = $main_dir . 'Sapro_S_trinity.fa.final.annotated';
my $out = $main_dir . 'snpHomology.out';
my $np = 2;

open(OUT, ">$out");
foreach my $compare (keys %compare) {
	my $seq1 = $compare{$compare}{'seq1'};
	my $seq2 = $compare{$compare}{'seq2'};
	
	my $contact = $1 if $seq1 =~ m/$main_dir([a-z]+_[a-z])/i;
	my $reads1gz = $reads_dir . $contact . '/' . $contact . '_1.fastq.gz';
	my $reads2gz = $reads_dir . $contact . '/' . $contact . '_2.fastq.gz';
	my $readsugz = $reads_dir . $contact . '/' . $contact . '_u.fastq.gz';
	
	my $reads1 = unzip($reads1gz);
	my $reads2 = unzip($reads2gz);
	my $readsu = unzip($readsugz);
	my $pileup = mapping($reads1,$reads2,$readsu,$compare,$seq2);
	rezip($reads1); rezip($reads2); rezip($readsu);
	compareSeq($seq2,$pileup);
	}
	
sub compareSeq {
	
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
	
sub mapping {
	my ($reads1,$reads2,$readsu,$compare,$seq2) = @_;

	my $out = $main_dir .  $compare . ".mpileup";
	my $call1 = system("bowtie2-build $seq2 $seq2.bowtie2");
	my $call2 = system("bowtie2 -x $seq2.bowtie2 -1 $reads1 -2 $reads2 -S bowtie1.sam -5 5 -3 5 --sensitive -k 10 -X 300 -p $np");
	my $call3 = system("bowtie2 -x $seq2.bowtie2 $readsu -S bowtie2.sam -5 5 -3 5 --sensitive -k 10 -p $np");
	my $call4 = system("samtools view -bS bowtie1.sam > bowtie1.bam");
	my $call5 = system("samtools view -bS bowtie2.sam > bowtie2.bam");
	my $call6 = system("samtools merge bowtie.bam bowtie1.bam bowtie2.bam");
	my $call7 = system("samtools sort bowtie.bam bowtie_sorted");
	my $call8 = system("samtools faidx $seq2");
	my $call9 = system("samtools mpileup -uf $seq2 bowtie_sorted.bam | bcftools view -bvcg - > var.raw.bcf");
	my $call10 = system("bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 -a 2 -w 0 > $out");
	my $call10 = system("rm $seq2.bowtie2* bowtie1.sam bowtie2.sam bowtie1.bam bowtie2.bam bowtie.bam bowtie_sorted* var.raw.bcf");

	return($out);
	}
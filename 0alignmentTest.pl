use warnings;
use strict;

my @contact = qw(Carlia_N Carlia_S);
my $rootdir = '/media/DataDrive/sutureGenomics/';
my $outdir = $rootdir . 'alignmentTest/';
mkdir($outdir) unless(-d $outdir);
my $np = 4;
my $time_out = $outdir . "alignerTiming.out"; 
open(OUT, ">>$time_out");
my $libfile = $rootdir . 'library';
my @aligners = qw(bowtie bowtie2 bwa soap smalt novo stampy);

#defines the library
my %lib;
open(IN, "<$libfile");
while (<IN>) {
        chomp(my $line = $_);
        my @d = split(/\t/, $line);
        my $lib = $d[2] . '/' . $d[0] . '/';
        $lib{$d[2]}{$d[0]} = {'dir' => $lib, 'index' => $d[1]}; 
        } 
close(IN);

foreach my $lineage (@contact) {
	foreach my $lib (keys %{$lib{$lineage}}) {    
	    if ($lib =~ m/SS56/) {
	    my $subdir  = $rootdir . $lib{$lineage}{$lib}{'dir'};
	    my $reads1gz = $subdir . $lib . '_1p_final.fastq.gz';
	    my $reads2gz = $subdir . $lib . '_2p_final.fastq.gz';
	    my $readsugz = $subdir . $lib . '_u_final.fastq.gz';

#		my $reads1 = unzipFiles($reads1gz);
#		my $reads2 = unzipFiles($reads2gz);
#		my $readsu = unzipFiles($readsugz);
	
#		my $seq = $rootdir . 'seqfiles/' . $contact[0] . '_annotated2.fa';
#		my $files = runAligners($seq,$np,$reads1,$reads2,$readsu,$lib);
	
#		rezipFiles($reads1);
#		rezipFiles($reads2);
#		rezipFiles($readsu);	    
	    }
	}
}

foreach my $align (@aligners) {
  my $lib;
   foreach my $lineage (@contact) {
	my @lib;
	foreach my $lib (sort {$a cmp $b} keys %{$lib{$lineage}}) {
	    my $file = $outdir . $lib . "." . $align . ".bam";
	    push(@lib,$file);
	}
	$lib .= join(" ",@lib); $lib .= ' ';
 } 
    my $out = $outdir . $contact[0] . '_' . $contact[1] . "." . $align . ".vcf";
    my $seq = $rootdir . 'seqfiles/' . $contact[0] . '_annotated.fa';
    my $call1 = system("samtools faidx $seq") unless(-f $seq . ".fai");
    unless(-f $out) {
	my $call2 = system("samtools mpileup -d 20000 -A -E -I -uf $seq $lib | bcftools view -bvcg - > var_out.raw.bcf");
	my $call3 = system("bcftools view var_out.raw.bcf | vcfutilsTL02.pl varFilter -D 30000 -w 0 -e 0.00000001 > $out");
	my $call4 = system("rm var_out.raw.bcf");
   }
}

sub rezipFiles {
	my ($file) = @_;	
	my $call = system("gzip -1 $file");
	}

sub runAligners {
	my ($assembly,$np,$forward,$reverse,$unpaired,$lib) = @_;
	
	#BWA
#	unless (-f $outdir . $lib . ".bwa.bam") {
#		my $call1 = system("bwa index $assembly");
#		my $start1 = time;
 #     	my $call2 = system("bwa aln -n 0.05 -t $np $assembly $forward > $forward" . ".sai");
#		my $call3 = system("bwa aln -n 0.05 -t $np $assembly $reverse > $reverse" . ".sai");
#		my $call4 = system("bwa aln -n 0.05 -t $np $assembly $unpaired > $unpaired" . ".sai");
#		my $call5 = system("bwa sampe -n 10 $assembly $forward.sai $reverse.sai $forward $reverse > bwa1.sam");
#		my $call6 = system("bwa samse -n 10 $assembly $unpaired.sai $unpaired > bwa2.sam");
#		my $end1 = int((time - $start1)/60);
#		print OUT "bwa\t$lib\t$end1\n";
#		my $call7 = system("samtools view -bS bwa1.sam > bwa1.bam");
#		my $call8 = system("samtools view -bS bwa2.sam > bwa2.bam");
#		my $tmp = $outdir . $lib . '.bam';
#		my $final = $outdir . $lib . '.bwa';
#		my $call9 = system("samtools merge $tmp bwa1.bam bwa2.bam");
#		my $call10 = system("samtools sort $tmp $final");
#		my $call11 = system("rm *sam bwa1.bam bwa2.bam $tmp $forward*sai $reverse*sai");
#		}
	
	
	#novoalign
#	unless (-f $outdir . $lib . ".novo.bam") {
#		my $call11 = system("novoindex $assembly.novo $assembly");
#		my $start2 = time;
#		my $call12 = system("novoalign -d $assembly.novo -F ILM1.8 -o SAM -f $forward $reverse > novo1.sam");
#		my $call13 = system("novoalign -d $assembly.novo -F ILM1.8 -o SAM -f $unpaired > novo2.sam");
#		my $end2 = int((time - $start2)/60);
#		print OUT "novo\t$lib\t$end2\n";
#		my $call14 = system("samtools view -bS novo1.sam > novo1.bam");
#		my $call15 = system("samtools view -bS novo2.sam > novo2.bam");
#		my $tmp = $outdir . $lib . '.novotmp.bam';
#		my $final = $outdir . $lib . '.novo';
#		my $call16 = system("samtools merge $tmp novo1.bam novo2.bam");
#		my $call10 = system("samtools sort $tmp $final");
#		my $call17 = system("rm *sam novo1.bam novo2.bam $tmp");		
#}


	#smalt
#	unless (-f $outdir . $lib . ".smalt.bam") {
#		my $call18 = system("smalt index $assembly.smalt $assembly");
#		my $start3 = time;
#		my $call19 = system("smalt map -f sam -l pe -n $np -o smalt1.sam $assembly.smalt $forward $reverse");
#		my $call20 = system("smalt map -f sam -n $np -o smalt2.sam $assembly.smalt $unpaired");
#		my $end3 = int((time - $start3)/60);
#		print OUT "smalt\t$lib\t$end3\n";
		
#		my $calla = system("samtools faidx $assembly");
#		my $call21 = system("samtools view -bS -t $assembly.fai smalt1.sam > smalt1.bam");
#		my $call22 = system("samtools view -bS -t $assembly.fai smalt2.sam > smalt2.bam");
#		my $tmp = $outdir . $lib . '.bam';
#		my $final = $outdir . $lib . '.smalt';
#		my $call23 = system("samtools merge $tmp smalt1.bam smalt2.bam");
#		my $call10 = system("samtools sort $tmp $final");
#		my $call24 = system("rm *sam smalt1.bam smalt2.bam $tmp");
#		}

	#bowtie1
#	unless (-f $outdir . $lib . ".bowtie.bam") {
#		my $call33 = system("bowtie-build $assembly $assembly.bowtie1");
#		my $start5 = time;
#		my $call34 = system("bowtie $assembly.bowtie1 -1 $forward -2 $reverse -5 5 -3 5 -X 300 -a --best -e 200 -S bowtie1.sam");
#		my $call35 = system("bowtie $assembly.bowtie1 $unpaired -5 5 -3 5 -a --best -e 200 -S bowtie2.sam");
#		my $end5 = int((time - $start5)/60);
#		print OUT "bowtie\t$lib\t$end5\n";
#		my $call36 = system("samtools view -bS bowtie1.sam > bowtie1.bam");
#		my $call37 = system("samtools view -bS bowtie2.sam > bowtie2.bam");
#		my $tmp = $outdir . $lib . '.bam';
#		my $final = $outdir . $lib . '.bowtie';
#		my $call38 = system("samtools merge $tmp bowtie1.bam bowtie2.bam");
#		my $call10 = system("samtools sort $tmp $final");
#		my $call39 = system("rm *sam bowtie1.bam bowtie2.bam $tmp");
#		}

	#bowtie2
#	unless (-f $outdir . $lib . ".bowtie2.bam") {
#		my $call40 = system("bowtie2-build $assembly $assembly.bowtie2");
#		my $start6 = time;
#		my $call41 = system("bowtie2 -x $assembly.bowtie2 -1 $forward -2 $reverse -S bow2tie1.sam -5 5 -3 5 --sensitive -k 10 -X 300 -p $np --no-mixed");
#		my $call42 = system("bowtie2 -x $assembly.bowtie2 $unpaired -S bow2tie2.sam -5 5 -3 5 --sensitive -k 10 -p $np");
#		my $end6 = int((time - $start6)/60);
#		print OUT "bowtie2\t$lib\t$end6\n";
#		my $call43 = system("samtools view -bS bow2tie1.sam > bow2tie1.bam");
#		my $call44 = system("samtools view -bS bow2tie2.sam > bow2tie2.bam");
#		my $tmp = $outdir . $lib . '.bam';
#		my $final = $outdir . $lib . '.bowtie2';
#		my $call45 = system("samtools merge $tmp bow2tie1.bam bow2tie2.bam");
  #    	my $call10 = system("samtools sort $tmp $final");
 #     	my $call46 = system("rm *sam bow2tie1.bam bow2tie2.bam $tmp");
#		}

	#stampy -- cannot run on that machine
	unless (-f $outdir . $lib . ".stampy.bam") {
		my $call25 = system("/home/singhal/programs/stampy-1.0.17/stampy.py --species=lizard --assembly=trinity -G $assembly $assembly");
		my $start4 = time;
		my $call26 = system("/home/singhal/programs/stampy-1.0.17/stampy.py -g $assembly -H $assembly");
		my $call27 = system("/home/singhal/programs/stampy-1.0.17/stampy.py -g $assembly -h $assembly -o stampy1.sam --substitutionrate=0.05 --insertsize=110 -M $forward $reverse");
		my $call28 = system("/home/singhal/programs/stampy-1.0.17/stampy.py -g $assembly -h $assembly -o stampy2.sam --substitutionrate=0.05 -M $unpaired");
		my $end4 = int((time - $start4)/60);
		print OUT "stampy\t$lib\t$end4\n";
		my $call29 = system("samtools view -bS stampy1.sam > stampy1.bam");
		my $call30 = system("samtools view -bS stampy2.sam > stampy2.bam");
		my $tmp = $outdir . $lib . '.stamptmp.bam';
		my $final = $outdir . $lib . '.stampy';
		my $call31 = system("samtools merge $tmp stampy1.bam stampy2.bam");
		my $call10 = system("samtools sort $tmp $final");
		my $call32 = system("rm *sam stampy1.bam stampy2.bam $tmp");
		}

	#soap2
#	unless (-f $outdir . $lib . ".soap.bam") {
#		my $call47 = system("2bwt-builder $assembly");
#		my $start7 = time;
#		my $call48 = system("soap -a $forward -b $reverse -D $assembly.index -o soap1 -2 out.sam -m 0 -x 500 -r 2 -v 4 -p $np");
#		my $call49 = system("soap -a $unpaired -D $assembly.index -o soap2 -r 2 -v 4 -p $np");
#		my $callb = system("soap2sam.pl soap1 > soap1.sam");
#		my $callc = system("soap2sam.pl soap2 > soap2.sam");
#		my $end7 = int((time - $start7)/60);
#		print OUT "soap2\t$lib\t$end7\n";
#		my $calla = system("samtools faidx $assembly");
#		my $call50 = system("samtools view -bS -t $assembly.fai soap1.sam > soap1.bam");
#		my $call51 = system("samtools view -bS -t $assembly.fai soap2.sam > soap2.bam");
#		my $tmp = $outdir . $lib . '.bam';
#		my $final = $outdir . $lib . '.soap';
#		my $call52 = system("samtools merge $tmp soap1.bam soap2.bam");
#		my $call10 = system("samtools sort $tmp $final");
#		my $call53 = system("rm *sam soap1.bam soap2.bam soap1 soap2 $tmp");
#		}
	}

sub unzipFiles {
	my ($file) = @_;
	my $call = system("gunzip $file");
	my $file2 = $1 if $file =~ m/(.*)\.gz/;
	return($file2);
	}

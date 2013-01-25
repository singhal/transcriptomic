use warnings;
use strict;

my %angsd = ('0'=>'A/A','1'=>'A/C','2'=>'A/G','3'=>'A/T','4'=>'C/C','5'=>'C/G','6'=>'C/T','7'=>'G/G','8'=>'G/T','9'=>'T/T');

my @contacts = (["Lampro_N","Lampro_C"],["Sapro_C","Sapro_S"],["Carlia_N","Carlia_S"],["Lampro_C","Lampro_S"]);
#my @contacts = (["Carlia_N","Carlia_S"]);
my $dir = '/media/DataDrive/sutureGenomics/';
my $finaldir = $dir . 'snpTesting/';
my $tmpdir = $dir . 'tmp/';
my $lib = $dir . 'library';
my $np = 2;

my $mincov = 6;
my $minreads = 3;
my $qual = 15;
my $pval = 0.05;

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
#runVarScan();
#runSamtools();
#runFreebayes();
#runSOAP();
#runBrute();
runANGSD();
rmdir($tmpdir);

########################
# subroutines level 1  #
########################

sub runBrute {
	#loop-de-loop!
    for (my $i = 0; $i < scalar(@contacts); $i++) {
	
		my $name = $contacts[$i][0] . "_" . $contacts[$i][1];
 	
		#make file with just annotated transcripts
		my $seq = makeSeq($contacts[$i][0]);

		my $call1 = system("bowtie2-build -q $seq $seq") unless(-f $seq . ".1.bt2");
		#contact1                                                                                                  
                my $out1 = mapping($seq,\%lib,$contacts[$i][0],$dir);
                #contact2                                                                                                  
                my $out2 = mapping($seq,\%lib,$contacts[$i][1],$dir);

		my ($pile1,$pile2) = runPileup($seq,$dir,$name,$out1,$out2,$finaldir);   
		#get depth from contacts and mpileup for the last
		runBruteSNP($name,$pile1,'1',$finaldir);
		runBruteSNP($name,$pile2,'2',$finaldir);
		}
	}

sub runVarScan {
	#loop-de-loop!
    for (my $i = 0; $i < scalar(@contacts); $i++) {
	
		my $name = $contacts[$i][0] . "_" . $contacts[$i][1];
 	
		#make file with just annotated transcripts
		my $seq = makeSeq($contacts[$i][0]);
	
		#map reads back to transcripts
		my $call1 = system("bowtie2-build -q $seq $seq") unless(-f $seq . ".1.bt2");
		#contact1
		my $out1 = mapping($seq,\%lib,$contacts[$i][0],$dir);
		#contact2
		my $out2 = mapping($seq,\%lib,$contacts[$i][1],$dir);
 
		#get depth from contacts and mpileup for the last
		my ($pile1,$pile2) = runPileup($seq,$dir,$name,$out1,$out2,$finaldir);
		runVar($name,$pile1,$pile2,$finaldir);
		}
	}
	
sub runFreebayes {
	#loop-de-loop!
    for (my $i = 0; $i < scalar(@contacts); $i++) {
	
		my $name = $contacts[$i][0] . "_" . $contacts[$i][1];
 	
		#make file with just annotated transcripts
		my $seq = makeSeq($contacts[$i][0]);
	
		#map reads back to transcripts
		my $call1 = system("bowtie2-build -q $seq $seq") unless(-f $seq . ".1.bt2");
		#contact1
		my $out1 = mapping($seq,\%lib,$contacts[$i][0],$dir);
		#contact2
		my $out2 = mapping($seq,\%lib,$contacts[$i][1],$dir);
 
		#get depth from contacts and mpileup for the last
		runFree($seq,$dir,$name,$out1,$out2,$finaldir);
		}
	}	
	
sub runSamtools {
	#loop-de-loop!
    for (my $i = 0; $i < scalar(@contacts); $i++) {
	
		my $name = $contacts[$i][0] . "_" . $contacts[$i][1];
		#make file with just annotated transcripts
		my $seq = makeSeq($contacts[$i][0]);
	
		#map reads back to transcripts
		my $call1 = system("bowtie2-build -q $seq $seq") unless(-f $seq . ".1.bt2");
		#contact1
		my $out1 = mapping($seq,\%lib,$contacts[$i][0],$dir);
		#contact2
		my $out2 = mapping($seq,\%lib,$contacts[$i][1],$dir);
 
		#get depth from contacts and mpileup for the last
		runSam($seq,$dir,$name,$out1,$out2,$finaldir);
		}
	}	

	
sub runANGSD {
	#loop-de-loop!
    for (my $i = 0; $i < scalar(@contacts); $i++) {
	
		my $name = $contacts[$i][0] . "_" . $contacts[$i][1];
		#make file with just annotated transcripts
		my $seq = $finaldir . $name . '.fa';
	
		#map reads back to transcripts
		my $call1 = system("bowtie2-build -q $seq $seq") unless(-f $seq . ".1.bt2");
		#contact1
		my $out1 = mapping($seq,\%lib,$contacts[$i][0],$dir,$name);
		#contact2
		my $out2 = mapping($seq,\%lib,$contacts[$i][1],$dir,$name);
		runANGSDsnp($out1,'1',$finaldir,$name,$seq);
		runANGSDsnp($out2,'2',$finaldir,$name,$seq);
		}
	}	

	
sub runSOAP {
	#loop-de-loop!
    for (my $i = 0; $i < scalar(@contacts); $i++) {
	
		my $name = $contacts[$i][0] . "_" . $contacts[$i][1];
		#make file with just annotated transcripts
		my $seq = makeSeq($contacts[$i][0]);
	
		#map reads back to transcripts
		my $call1 = system("2bwt-builder $seq") unless(-f $seq . ".index.amb");
		#contact1
		my $out1 = mappingSOAP($seq,\%lib,$contacts[$i][0],$dir);
		#contact2
		my $out2 = mappingSOAP($seq,\%lib,$contacts[$i][1],$dir);
 
		#get depth from contacts and mpileup for the last
		runSOAPsnp($seq,$dir,$name,$out1,$out2,$finaldir);
		}
	}		

###########################
# run the subroutines     #
###########################

sub runANGSDsnp {
	my ($bam,$id,$finaldir,$name,$seq) = @_;
	
	my @bamfiles = @{$bam};
	my $bamfiles = join("\t",@bamfiles);
	
	my $out = $finaldir . $name . "." . $id;
	
#	my $call1 = system("/home/singhal/programs/angsd0.03/angsd.g++ -anc $seq -realSFS 1 -doMaf 2 -doGeno 4 -outfiles $out mpileup -f $seq -d 20000 -A -I -g $bamfiles > NULL");
#	my $call2 = system("rm NULL");

	my $sfsout = $out . '.sfs';
#	my $call2 = system("/home/singhal/programs/angsd0.03/misc/optimSFS.gcc -nChr 10 -binput $sfsout");
	my $posterior = $out . '_posterior_probabilities.txt';
#	my $sfsml = $out . '.sfs.ml';
#	my $call3 = system("/home/singhal/programs/angsd0.03/misc/sfstools.g++ -nChr 10 -priorFile $sfsml -sfsFile $sfsout -dumpBinary 0 > $posterior");

	my $geno = $out . ".geno";
	my $actual = $finaldir . $name . ".angsd.vcf" . $id . ".out";
#	my $call = system("mv $geno $actual");

	my %var;
	my $pos = $sfsout . ".pos";
	open(POS, "<$pos");
	open(PP, "<$posterior");
	while(<POS>) {
	    chomp(my $posline = $_);
	    chomp(my $ppline = <PP>);
	    my $var = $1 if $ppline =~ m/^(\S+)/;
	    if ($var <= 0.05) {
		#this is a legit variant!
		my @d = split(/\t/,$posline);
		$var{$d[0]}{$d[1]}++;
	    }
	}
	close(PP); close(POS);

	open(GENO, "<$geno");
	open(OUT, ">$actual");
	my $junk = <GENO>;
	while(<GENO>){
	    chomp(my $line = $_);
	    my @d = split(/\t/,$line);
	    if ($var{$d[0]}{$d[1]}) {
		print OUT $d[0], "\t", $d[1], "\t";
		for (my $i = 2; $i < 7; $i++) {
		    print OUT $angsd{$d[$i]}, "\t";
		}
		print OUT "\n";
	    }
	}
	close(OUT); close(GENO);
}
	
sub runBruteSNP {
	my ($name,$pile,$index,$finaldir) = @_;
	
	my $out = $finaldir . $name . ".brute.vcf" . $index . ".out";
	open(OUT, ">$out");
	
	open(IN, "<$pile");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		
		my @geno; my $output = 0; my $winSNP;
		
		for (my $i = 4; $i < 17; $i = $i + 3) {
			my $snp = $d[$i];
			
			while ($snp =~ m/(\+\d+|\-\d+)/g) {
				my $match = $1;
				my $num = $1 if $match =~ m/(\d+)/;
				$match = '\\' . $match; 
				$snp =~ s/$match[atgcn]{$num}//i;
                }
                
        
            my $ref = 0; $ref++ while $snp =~ m/([\.|\,])/g;
            my %snp; my $alt = 0;
            while ($snp =~ m/([a|t|c|g])/ig) {
            	$snp{uc($1)}++;
            	$alt++;
            	}
            my $totcov = $alt + $ref;	
            my @snps = sort {$snp{$b} <=> $snp{$a}} keys %snp;	            	    
            my $geno = './.';	                	
		 
		 			
            if ($totcov >= $mincov) {
            	if ($alt > 0) {
            		if ($snp{$snps[0]} >= $minreads) {
            			if ($ref >= $minreads) {
            				$geno = '0/1';
            				$output = 1;
					$winSNP = $snps[0];
            				}
            			else {	
            				$geno = '1/1';
            				$output = 1;
					$winSNP = $snps[0];
            				}
            			}
            		else {	
            			$geno = '0/0';
            			}
            		}
            	else {
            		$geno = '0/0';
            		}
            	}	
       		push(@geno,$geno);
       		}
       		
       	if ($output) {
       		my $gt = join("\t",@geno);
       		print OUT $d[0], "\t", $d[1], "\t", $d[2], "\t", $winSNP, "\t", $gt, "\n";
       		}      		
       	}	
    close(OUT);	
	}

sub runSOAPsnp {
	my ($seq,$dir,$name,$outfiles1,$outfiles2,$finaldir) = @_;
	my @out1 = @{$outfiles1};
	my @out2 = @{$outfiles2};
	my @out = (@out1,@out2);
	
	foreach my $out (@out) {
		my $id = $1 if $out =~ m/([A-Z|0-9]+)\.sorted/;
		my $outfile = $finaldir . $name . "." . $id . ".soap.out";
		unless (-f $outfile) {	
			my $call = system("soapsnp -i $out -d $seq -o $outfile -r 0.0005 –e 0.001 -t -u -L 200 -q");
			}
		}
	}

sub runFree {
	my ($seq,$dir,$name,$outfiles1,$outfiles2,$finaldir) = @_;
	my @out1 = @{$outfiles1};
	my @out2 = @{$outfiles2};
	my $out1 = $finaldir . $name . ".freebayes.vcf1.out";
	my $out2 = $finaldir . $name . ".freebayes.vcf2.out";
	
	unless (-f $out1) {	
		my @new;
		foreach my $file (@out1) {
			my $id = $1 if $file =~ m/([A-Z|0-9]+)\.sorted/;
			my $new = $file . '.picard';
			my $call = system("java -jar /home/singhal/programs/picard-tools-1.57/AddOrReplaceReadGroups.jar I=$file O=$new SO=coordinate RGID=$id RGLB=$id RGPL=Illumina RGPU=$id RGSM=$id") unless (-f $new);
			push(@new,$new);
			}
		my $files = join("\t", @new);
		my $call = system("freebayes -v $out1 -f $seq -P $pval -i -4 --min-coverage $mincov $files");
		}

	unless (-f $out2) {	
		my @new;
		foreach my $file (@out2) {
			my $id = $1 if $file =~ m/([A-Z|0-9]+)\.soap/;
			my $new = $file . '.picard';
			my $call = system("java -jar /home/singhal/programs/picard-tools-1.57/AddOrReplaceReadGroups.jar I=$file O=$new SO=coordinate RGID=$id RGLB=$id RGPL=Illumina RGPU=$id RGSM=$id") unless (-f $new);
			push(@new,$new);
			}
		my $files = join("\t", @new);
		my $call = system("freebayes -v $out2 -f $seq -P $pval -i -4 --min-coverage $mincov $files");
		}
	}

sub runVar {
    my ($name,$pile1,$pile2,$finaldir) = @_;
    my $out1 = $finaldir . $name . '.varscan.vcf1.out';
    my $out2 = $finaldir . $name . '.varscan.vcf2.out';
    my $call1 = system("java -jar /home/singhal/programs/VarScan.jar mpileup2snp $pile1 --min-coverage $mincov --min-reads2 $minreads --min-avg-qual $qual --output-vcf 1 --p-value $pval > $out1") unless(-f $out1);
    my $call2 = system("java -jar /home/singhal/programs/VarScan.jar mpileup2snp $pile2 --min-coverage $mincov --min-reads2 $minreads --min-avg-qual $qual --output-vcf 1 --p-value $pval > $out2") unless(-f $out2);
	}
	
sub runSam {	
	my ($seq,$dir,$name,$outfiles1,$outfiles2,$finaldir) = @_;
	my @out1 = @{$outfiles1};
	my @out2 = @{$outfiles2};
	my $call1 = system("samtools faidx $seq") unless(-f $seq . ".fai");
	my $out1 = $finaldir . $name . ".samtools.vcf1.out";
	my $out2 = $finaldir . $name . ".samtools.vcf2.out";
	unless (-f $out1) {
		my $files = join("\t", @out1);
		my $call2 = system("samtools mpileup -d 20000 -A -E -I -uf $seq $files | bcftools view -bvcg - > var_out.raw.bcf");
		my $call3 = system("bcftools view var_out.raw.bcf | vcfutilsTL02.pl varFilter -D 30000 -w 0 -e 0.00000001 > $out1");
		my $call4 = system("rm var_out.raw.bcf");
		}
	unless (-f $out2) {
		my $files = join("\t", @out2);
		my $call2 = system("samtools mpileup -d 20000 -A -E -I -uf $seq $files | bcftools view -bvcg - > var_out.raw.bcf");
		my $call3 = system("bcftools view var_out.raw.bcf | vcfutilsTL02.pl varFilter -D 30000 -w 0 -e 0.00000001 > $out2");
		my $call4 = system("rm var_out.raw.bcf");
		}	
	}	
	
sub mappingSOAP {
	my ($seq, $libHash, $contact, $dir) = @_;
	
	my %lib = %{$libHash};
	my @outfiles;
	
	foreach my $lib (sort {$a cmp $b} keys %{$lib{$contact}}) {		
	    my $subdir = $dir . $lib{$contact}{$lib}{'dir'};
	    my $file1gz = $subdir . $lib . '_1p_final.fastq.gz';
	    my $file2gz = $subdir . $lib . '_2p_final.fastq.gz';
	    my $fileugz = $subdir . $lib . '_u_final.fastq.gz';

		my $out = $subdir . $lib . ".soap.sorted";
		unless(-f $out) {
		    my $file1 = unzip($file1gz);
		    my $file2 = unzip($file2gz);
		    my $fileu = unzip($fileugz);
			my $out = runMappingSOAP($seq,$file1,$file2,$fileu,$subdir,$lib,$dir) unless (-f $out);
			my @files = ($file1,$file2,$fileu);
			zip(\@files);
			}
		push(@outfiles,$out);	
		}
	
	return(\@outfiles);
	}
	
sub runMappingSOAP {
	my ($seq,$file1,$file2,$fileu,$subdir,$lib,$dir) = @_;
	my $out = $subdir . $lib . ".soap.sorted";
	print $out, "\n";
	$seq = $seq . ".index";
	my $call1 = system("soap -a $fileu -D $seq -o soap_u.out -r 2 -p $np -v 3");
	my $call2 = system("soap -a $file1 -b $file2 -D $seq -o soap_p1.out -2 soap_p2.out -m 0 -x 400 -r 2 -p $np -v 3");
	my $call3 = system("cat soap_u.out soap_p1.out soap_p2.out > soap.out");
	my $call4 = system("sort -k8,8 -k9,9n soap.out > $out");
	my $call5 = system("rm soap_u.out soap_p1.out soap_p2.out soap.out");	
	return($out);
	}	
	
sub mapping {
	my ($seq, $libHash, $contact, $dir,$name) = @_;
	
	my %lib = %{$libHash};
	my @outfiles;
	
	my $outdir = $finaldir . $name . '/';
	mkdir($outdir) unless(-d $outdir);

	foreach my $lib (sort {$a cmp $b} keys %{$lib{$contact}}) {		
	    my $subdir = $dir . $lib{$contact}{$lib}{'dir'};
	    my $file1gz = $subdir . $lib . '_1p_final.fastq.gz';
	    my $file2gz = $subdir . $lib . '_2p_final.fastq.gz';
	    my $fileugz = $subdir . $lib . '_u_final.fastq.gz';

		my $out = $outdir . $lib . ".sorted.bam";
		unless(-f $out) {
		    my $file1 = unzip($file1gz);
		    my $file2 = unzip($file2gz);
		    my $fileu = unzip($fileugz);
			my $out = runMapping($seq,$file1,$file2,$fileu,$subdir,$lib,$dir,$out) unless (-f $out);
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
	my $call1 = system("samtools faidx $seq") unless(-f $seq . ".fai");
	my $out1 = $finaldir . $contact . ".mpileup1.out";
	my $out2 = $finaldir . $contact . ".mpileup2.out";
	unless (-f $out1) {
		my $files = join("\t", @out1);
		my $call2 = system("samtools mpileup -d 20000 -A -I -B -f $seq $files > $out1");	
		}
	unless (-f $out2) {
		my $files = join("\t", @out2);
		my $call2 = system("samtools mpileup -d 20000 -A -I -B -f $seq $files > $out2");	
		}
	return($out1,$out2);
	}
	
sub runMapping {
	my ($seq,$file1,$file2,$fileu,$subdir,$lib,$dir,$out) = @_;

	$out =~ s/\.bam//;

	my $call1 = system("bowtie2 -x $seq -1 $file1 -2 $file2 -S bowtie1_out1.sam -5 5 -3 5 --sensitive -k 10 -X 300 -p $np");
	my $call2 = system("bowtie2 -x $seq $fileu -S bowtie2_out1.sam -5 5 -3 5 --sensitive -k 10  -p $np");
	my $call3 = system("samtools view -bS bowtie1_out1.sam > bowtie1_out1.bam");
	my $call4 = system("samtools view -bS bowtie2_out1.sam > bowtie2_out1.bam");
	my $call5 = system("samtools merge bowtie_out1.bam bowtie1_out1.bam bowtie2_out1.bam");
	my $call6 = system("samtools sort bowtie_out1.bam $out");
	my $call7 = system("rm bowtie1_out1.sam bowtie2_out1.sam bowtie_out1.bam bowtie1_out1.bam bowtie2_out1.bam");	
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

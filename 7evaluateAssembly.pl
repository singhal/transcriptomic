###########################################################################################
# a script to evaluate the quality of your assemblies by a number of metrics              #
# external dependencies: blat (v34), blastall/formatdb (>2.217), bowtie/bowtie-build (v1) #
# note that uniqueness measures % of contigs in evaluated assembly that are unique        #
# relative to the assembly under comparison												  #
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 3 January 2012                 #
###########################################################################################

use warnings;
use strict;

#home directory with all my files; make sure to end with a backslash
my $dir = '/media/DataDrive/sutureGenomics/';
#the drive with the assemblies
my $aDir = 'assemblies/';
#the libraries
my @lib = qw(Carlia_N Carlia_S Lampro_C Lampro_S Sapro_C Sapro_S Lampro_N);
#how many reads do you want to test?; used in readsUsed()
my $randomReads = 100000; 
#how many contigs do you want to test; used in orfFinder() & annotable() & contiguity() & compareAsemblies -- i.e., all blatting
my $contigs = 1000; 
#to see if the unique stuff is any good; used in compareAssemblies() and accuracy() -- i.e., all blasting
my $testSet = 100;
#database to blastx against
my $prot_db = '/media/DataDrive/sutureGenomics/genomes/Anolis_pep.fa';
#database to blat against
my $cdna_db = '/media/DataDrive/sutureGenomics/genomes/Anolis_cDNA.fa';
#number of processors you can use for blasting
my $proc = 4;
#evalue significance level when comparing two different assemblies
my $evalue1 = 1e-100;
#evalue significance level when comparing to Anolis (cDNA or prot)
my $evalue2 = 1e-20;
#location of cegma file; available from http://korflab.ucdavis.edu/Datasets/cegma/
my $cegma = '/media/DataDrive/sutureGenomics/genomes/cegma.fa';
my $finalOut = '/media/DataDrive/sutureGenomics/assemblyEvaluation.out';

###########################
# run the subroutines     #
###########################

my %master;
foreach my $lib (@lib) {
	my $subdir = $dir . $lib . '/';
	my $file1gz = $subdir . $lib . '_1.fastq.gz';
	my $file2gz = $subdir . $lib . '_2.fastq.gz';
	my $file1   = $subdir . $lib . '_1.fastq';
	my $file2   = $subdir . $lib . '_2.fastq';
	my $assemblyDir = $subdir . '/' . $aDir;
	
	if (-d $assemblyDir) {
		my @assemblies = <$assemblyDir*final>; 
		
		my $nLines = $randomReads * 4;
		my $call1 = system("gunzip $file1gz");
		my $call2 = system("gunzip $file2gz");
		my $reads1 = $subdir . "testReads1.fastq";
		my $reads2 = $subdir . "testReads2.fastq";
		my $call3 = system("head -n $nLines $file1 > $reads1");
		my $call4 = system("head -n $nLines $file2 > $reads2");
		my $call5 = system("gzip -1 $file1");
		my $call6 = system("gzip -1 $file2");
		
		foreach my $seq (@assemblies) {	
			print "Hey hey hey! I am evaluating $seq right now.\n";	
			my $id = $1 if $seq =~ m/([^\/]+)$/i;
			assemblyStats($seq,$id);
			readsUsed($seq,$id,$reads1,$reads2);
			my $random1 = bootstrap($seq,$contigs);	
			ORFfinder($id,$random1);
			foreach my $seq2 (@assemblies) {
				if ($seq2 ne $seq) {
					compareAssemblies($id,$seq2,$seq,$random1); 
					}
				}
			annotatable($id,$random1);
			my $random2 = bootstrap($seq,$testSet);
			accuracy($id,$random2);		
			my $random3 = bootstrap($cdna_db,$contigs);
			contiguity($seq,$id,$random3);
			cegma($seq,$id,$cegma);
			}
		unlink($reads1); unlink($reads2);	
		}	
	}	

###########################
# report the results      #
###########################	

open(OUT, ">$finalOut");
my %d;
foreach my $id (keys %master) {
	foreach (keys %{$master{$id}}) {
		$d{$_}++;
		}
	}
print OUT "contact\tassembly\t";	
foreach (sort {$a cmp $b} keys %d) {
	print OUT $_, "\t";
	}
print OUT "\n";	
foreach my $id (sort {$a cmp $b} keys %master) {
	my $contact = $1 if $id =~ m/^([a-z]+_[a-z]+)/i;
	my $assembly = $1 if $id =~ m/_([a-z]+)\.fa/i;
	print OUT $contact, "\t", $assembly, "\t";
	foreach my $d (sort {$a cmp $b} keys %d) {
		if ($master{$id}{$d}) {
			print OUT sprintf("%.3f",$master{$id}{$d}), "\t";
			}
		else {
			print OUT "NA\t";
			}
		}
	print OUT "\n";
	}
close(OUT);
	
###########################
# behold the subroutines  #
###########################

sub cegma {
    my ($seq,$id,$cegma) = @_;
    
    my $out = $seq . '.cegma.blast.out';
    my $call1 = system("formatdb -i $seq -p F");
    my $call2 = system("blastall -p tblastn -d $seq -i $cegma -e $evalue2 -m 8 -o $out -b 1");
    my @call3 = `cut -f 1 $out | uniq | wc`;
    
    my $count = $1 if $call3[0] =~ m/^\s+(\d+)/;
    $master{$id}{'cegma'} = $count;

    unlink($out); my $call4 = system("rm $seq" . ".n*");
}

sub contiguity {
	my ($seq, $id, $random) = @_;
	
	my $temp1 = 'query.fa';
	my $temp2 = 'blatout';
	
	my %seq;
	open(OUT, ">$temp1");
	for (my $i = 0; $i < scalar(@{$random}); $i++) { 
		print OUT ">cdna", $i, "\n", $random->[$i], "\n";
		$seq{'cdna' . $i}{'length'} = length($random->[$i]);
		}
	close(OUT);
	
	my %contiguity; my %complete;
	my $call = system("blat $seq $temp1 $temp2 -noHead -out=blast8");
	open(IN, $temp2);
	while(<IN>) {
		my @d = split(/\t/,$_);
		if ($d[10] < $evalue2) {
			my $min = (sort { $a <=> $b } ($d[6], $d[7]))[0];
			my $max = (sort { $a <=> $b } ($d[6], $d[7]))[1];
			unless ($contiguity{$d[0]}) {
				for (my $i = $min; $i <= $max; $i++) {
					$contiguity{$d[0]}{$i}++;
					}
				}
			for (my $i = $min; $i <= $max; $i++) {
				$complete{$d[0]}{$i}++;
				}	
			}	
		}
		
	my $complete; my $contiguity;	
	foreach my $c (%seq) {
		if ($complete{$c}) {
			$complete += scalar( keys %{$complete{$c}} ) / $seq{$c}{'length'};
			$contiguity += scalar( keys %{$contiguity{$c}} ) / $seq{$c}{'length'};
			}
		}
	$complete = $complete/scalar(keys %complete);
	$contiguity = $contiguity/scalar(keys %complete);
	
	$master{$id}{'complete'} = $complete;		
	$master{$id}{'contiguity'} = $contiguity;
	
	unlink($temp1); unlink($temp2);
	}
	
sub compareAssemblies {
	my ($id,$a1,$a2,$random) = @_;
	
	my %seq;
	open(IN, "<$a1");
	my $c;
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			$c = $1;
			}
		else {
			$seq{$c} .= $line;
			}
		}
	close(IN);	

	my %rand;
	my $temp1 = "randomSeq.fa";
	open(OUT, ">$temp1");
	for (my $i = 0; $i < scalar(@{$random}); $i++) { 
		print OUT ">contig", $i, "\n", $random->[$i], "\n";
		$rand{'contig' . $i} = $random->[$i];
		}
	close(OUT);

	my $db = $1 if $a1 =~ m/_([a-z|0-9]+).fa/i;
	my $out = $id . $db . ".out";
	my $blatCall = system("blat $a1 $temp1 $out -noHead -out=blast8");
	open(IN, "<$out");
	my %match;
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		#match already exists
		if ($match{$d[0]}) {
			#continuing the best part match
			if ($d[1] eq $match{$d[0]}{'match'}) {
				$match{$d[0]}{'error'} += ($d[4] + $d[5]);
				$match{$d[0]}{'al'} += $d[3];
				}	
			}
		#either a new match or a contig for which there isn't a good match	
		else {
			if ($d[10] < $evalue1) {
				$match{$d[0]} = { 'match' => $d[1] , 'error' => ($d[4] + $d[5]), 'al' => $d[3] }; 
				}
			}
		}
	close(IN);	
	
	my $accuracy; my $totLength; my $totAl;	my $numSeq = scalar( keys %{$seq{$a2}} );
	foreach my $c (keys %match) {
		$accuracy += $match{$c}{'error'};
		$totLength += $match{$c}{'al'};	
		$totAl += ($match{$c}{'al'} * 2) / ( length($seq{$match{$c}{'match'}}) + length($rand{$c}));
		delete($rand{$c});
		}
	$totAl = $totAl / ( scalar(keys %match) );
	$accuracy = 1 - ($accuracy / $totLength);
	my $unique = ($contigs - scalar(keys %match)) / $contigs;

	my $legitUniq = 'NA';
	unless (scalar(keys %rand) < 10) {
		my $testUnique = 'testUnique.fa';
		open(OUT, ">$testUnique");
		my @uniq = keys %rand;
		for (my $i = 0; $i < $testSet; $i++) {
			print OUT ">contig" . "$i" . "\n" . $rand{$uniq[int rand($#uniq)]} . "\n";
			}
		close(OUT);

		unless (-f $prot_db . ".pin") {
			my $call = system("formatdb -i $prot_db -t $prot_db");
			}
		my @call = `blastall -p blastx -d $prot_db -i $testUnique -m 8 -a $proc -b 1 -e $evalue2 | cut -f 1 | uniq`;
		$legitUniq = scalar(@call)/$testSet;
		unlink($testUnique); 
		}	

	$master{$id}{$db . '.totAl'} = $totAl;
	$master{$id}{$db . '.accuracy'} = $accuracy;
	$master{$id}{$db . '.unique'} = $unique;
	$master{$id}{$db . '.legitUniq'} = $legitUniq;
	
	unlink($out); unlink($temp1);
	}

sub bootstrap {
	my ($seq,$num) = @_;
	
	my @seq;
	my $c = -1;
	
	open(SEQ, "<$seq");
	while(<SEQ>) {
		chomp(my $line = $_);
		if ($line =~ m/^>(\S+)/) {
			$c++;
			}
		else { 
			$seq[$c] .= $line;
			}
		}
	close(SEQ);
		
		
	my @random;
	for (my $i = 0; $i < $num; $i++) {
		$random[$i] = $seq[int($#seq*rand())];
		}		
	return \@random;
	}		

sub ORFfinder {
	my ($id,$array_ref) = @_;
	$master{$id}{'orf'} = '0';
	foreach (@$array_ref) {
		my $seq = $_;
		my @longest;
		my $revcomp = reverse($seq);
		$revcomp =~ tr/ATGCatgc/GACTgact/; 
		my %seq = ('1' => $seq, '2' => substr($seq,1), '3' => substr($seq,2), '-1' => $revcomp, '-2' => substr($revcomp,1), '-3' => substr($revcomp,2));
		foreach(keys %seq) {
			my @orfs = translate($seq{$_}) =~ m/([A-Z]+)/g;
			@orfs = sort{length($a) <=> length($b)} @orfs;
			push(@longest, length($orfs[-1]));
			}
		@longest = sort{$a <=> $b} @longest;
		$master{$id}{'orf'} += $longest[-1]*3/length($seq); 
		}
	$master{$id}{'orf'} = sprintf("%.3f", $master{$id}{'orf'}/$contigs);
	}

sub translate {
	my ($string) = @_;
	$string = uc($string);
	my @codons = $string =~ m/(\S\S\S)/g;
	my %codons = (	'ATG'=>'M','ACG'=>'T','CTG'=>'L','CCG'=>'P','GTG'=>'V','GCG'=>'A','TTG'=>'L','TCG'=>'S',
					'ATA'=>'I','ACA'=>'T','CTA'=>'L','CCA'=>'P','GTA'=>'V','GCA'=>'A','TTA'=>'L','TCA'=>'S',
					'ATC'=>'I','ACC'=>'T','CTC'=>'L','CCC'=>'P','GTC'=>'V','GCC'=>'A','TTC'=>'F','TCC'=>'S',
					'ATT'=>'I','ACT'=>'T','CTT'=>'L','CCT'=>'P','GTT'=>'V','GCT'=>'A','TTT'=>'F','TCT'=>'S',
					'AGG'=>'R','AAG'=>'K','CGG'=>'R','CAG'=>'Q','GGG'=>'G','GAG'=>'E','TGG'=>'W','TAG'=>'*',
					'AGA'=>'R','AAA'=>'K','CGA'=>'R','CAA'=>'Q','GGA'=>'G','GAA'=>'E','TGA'=>'*','TAA'=>'*',
					'AGC'=>'S','AAC'=>'N','CGC'=>'R','CAC'=>'H','GGC'=>'G','GAC'=>'D','TGC'=>'C','TAC'=>'Y',
					'AGT'=>'S','AAT'=>'N','CGT'=>'R','CAT'=>'H','GGT'=>'G','GAT'=>'D','TGT'=>'C','TAT'=>'Y');
	my $translate;
	foreach(@codons) {
		if ($codons{$_}) {
			$translate = $translate . $codons{$_};
			}
		else {
			$translate = $translate . 'X';
			}
		}
	return($translate);
	}

sub readsUsed {
	my ($seq,$id,$reads1,$reads2) = @_;
	
	my $call1 = system("bowtie-build $seq $id -q");
	
	my $out = "bowtie.out";
	my $call2 = system("bowtie $id -1 $reads1 -2 $reads2 -k 1 -5 5 -3 5 --un un.fastq > $out");
	my $totDist;
	my $matches;
	open(OUT, "<$out"); 
	while(<OUT>) {
		$matches++;
		my $line1 = $_;
		my $line2 = <OUT>;
		my @d1 = split(/\s+/, $line1);
		my @d2 = split(/\s+/, $line2);
		$totDist += abs($d1[3] - $d2[3]); 	
		}
	close(OUT);
	$master{$id}{'numPairsMatched'} = $matches / $randomReads;
	$master{$id}{'avgDist'} = $totDist/ $matches;
	$master{$id}{'numReadsMatched'} = $matches*2;
	
	my $un1 = 'un_1.fastq';
	my $un2 = 'un_2.fastq';
	my $un_out1 = 'un1.out';
	my $un_out2 = 'un2.out';
	
	my $call3 = system("bowtie $id $un1 -k 3 --best -5 5 -3 5 > $un_out1");
	my $call4 = system("bowtie $id $un2 -k 3 --best -5 5 -3 5 > $un_out2");
	
	my %r1; my %d;
	open(IN, "<$un_out1");
	while(<IN>) {
		my @d = split(/\t/,$_);
		my $id1 = $1 if $d[0] =~ m/(\S+)_/;
		$r1{$id1}{$d[2]}++;
		$d{$d[0]}++;
		}
	close(IN);
	my %r2;
	open(IN, "<$un_out2");
	while(<IN>) {
		my @d = split(/\t/,$_);
		my $id1 = $1 if $d[0] =~ m/(\S+)_/;
		$r2{$id1}{$d[2]}++;
		$d{$d[0]}++;	
		}
	close(IN);
	
	my $chimera;
	foreach my $id (keys %r1) {
		#potential chimera
		my $match = 0;
		if ($r2{$id}) {
			foreach my $c (keys %{$r1{$id}}) {	
				$match++ if $r2{$id}{$c}
				}
			$chimera++ if $match == 0;	
			}
		}	
	
	$master{$id}{'numReadsMatched'} = $master{$id}{'numReadsMatched'} + scalar(keys %d);
	$master{$id}{'chimera'} = $chimera / ($master{$id}{'numReadsMatched'} / 2);
	$master{$id}{'numReadsMatched'} = $master{$id}{'numReadsMatched'} / ( $randomReads * 2 );
	
	my $call5 = system("rm $id"."*ebwt"); 
	unlink($un1,$un2,$un_out1,$un_out2,$out);
	}

sub assemblyStats {
	my ($seq,$id) = @_;
	
	my (%seq, $c);
	open(SEQ, "<$seq");
	while(<SEQ>) {
		chomp(my $line = $_);
		if ($line =~ m/^>(\S+)/) {
			$c = $1;
			}
		else { 
			$seq{$c} .= $line;
			}
		}
	close(SEQ);
	
	my $totLength;
	my @length;
	foreach my $c (sort {$a cmp $b} keys %seq) {
		push(@length, length($seq{$c}));
		$master{$id}{'totLength'} += length($seq{$c});
		}
		
	@length = sort {$a <=> $b} (@length);
	$master{$id}{'mean'} = int($master{$id}{'totLength'}/scalar(@length));
	$master{$id}{'median'} = $length[int($#length/2)];
	$master{$id}{'totContigs'} = scalar(@length);
	$master{$id}{'max'} = $length[$#length];
	my $track;
	for (my $i = 0; $i < scalar(@length); $i++) {
		$track += $length[$i];
		if ($length[$i] > 2000) {
			$master{$id}{'n2000'}++;
			$master{$id}{'n1000'}++;
			$master{$id}{'n500'}++;
			}
		elsif ($length[$i] > 1000) {
			$master{$id}{'n1000'}++;
			$master{$id}{'n500'}++;
			}
		elsif ($length[$i] > 500) {
			$master{$id}{'n500'}++;
			}
		if ($track/$master{$id}{'totLength'} > 0.5) {
			$master{$id}{'n50'} = $length[$i] unless $master{$id}{'n50'};			
			}	
		}				
	}
	
sub annotatable {
	my ($id, $random) = @_;
	
	my $temp1 = 'query.fa';
	my $temp2 = 'blatout';
	
	open(OUT, ">$temp1");
	for (my $i = 0; $i < scalar(@{$random}); $i++) { 
		print OUT ">contig", $i, "\n", $random->[$i], "\n";
		}
	close(OUT);
	
	my $call = system("blat $cdna_db $temp1 $temp2 -noHead -t=dnax -q=dnax -out=blast8");
	
	my %a; my $matches;
	open(IN, "<$temp2");
	while(<IN>){
		chomp(my $line = $_);
		my @d = split(/\t/, $line);
		unless($a{$d[0]}) {
			if ($d[10] < $evalue2) {
				#it is a match
				#how long is the match
				my $length = length($random->[$1]) if $d[0] =~ m/contig(\d+)/;
				$a{$d[0]}{'length'} = $d[3] / $length;
				#how accurate is the match
				$a{$d[0]}{'accuracy'} = $d[2];
				#how many matches do we have
				$matches++;
				}
			}
		}
	my $length; my $accuracy;
	foreach my $c (keys %a) {
		$length += $a{$c}{'length'};
		$accuracy += $a{$c}{'accuracy'};
		}
	$length = $length / $matches;
	$accuracy = $accuracy / $matches;
	
	$master{$id}{'anno'} = $matches/scalar(@{$random});
	$master{$id}{'lengthAnno'} = $length;
	$master{$id}{'accuAnno'} = $accuracy;
	
	unlink($temp1); unlink($temp2);
	}

sub accuracy {
	my ($id,$random) = @_;		
	
	my $temp1 = 'query.fa';
	open(OUT, ">$temp1");
	for (my $i = 0; $i < scalar(@{$random}); $i++) { 
		print OUT ">contig", $i, "\n", $random->[$i], "\n";
		}
	close(OUT);

	unless (-f $prot_db . ".pin") {
		my $call = system("formatdb -i $prot_db -t $prot_db");
		}
	
	my @call = `blastall -p blastx -d $prot_db -i $temp1 -m 8 -a $proc -b 1 -e $evalue2`;
	
	my %a; my $gap = 0; my $stop = 0;
	foreach my $line (@call) {
		my @d = split(/\t/,$line);
		unless($a{$d[0]}) {
			$a{$d[0]}++; my $subseq;
			if ($d[6] < $d[7]) {
				my $start = $d[6] - 1;
				my $length = $d[7] - $d[6] + 1;
				$subseq = substr $random->[$1], $start, $length if $d[0] =~ m/contig(\S+)/;
				}
			else {
				my $start = $d[7] - 1;
				my $length = $d[6] - $d[7] + 1;
				$subseq = substr $random->[$1], $start, $length if $d[0] =~ m/contig(\S+)/;
				$subseq = reverse($subseq); 
				$subseq =~ tr/atgcATGC/tacgTACG/;
				}			
			#check to see if there looks to be a gap
			$gap++ if length($subseq) % 3; 
			#check to see if there is a nonsense mutation
			my $stop1 = (translate($subseq) =~ tr/\*/\*/);
			$stop++ unless($stop1 == 0);
			}
		}
		
	$master{$id}{'nonsenseMu'} = $stop/scalar(keys %a); 
	$master{$id}{'gapMu'} = $gap/scalar(keys %a);
	unlink($temp1);
	}
		
		
		
		

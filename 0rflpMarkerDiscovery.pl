use warnings;
use strict;

use Bio::DB::Fasta;

my @contacts = (["Carlia_N","Carlia_S"]);
my $dir = '/media/DataDrive/sutureGenomics/';
my $finaldir = $dir . 'rflp/';
my $tmpdir = $dir . 'tmp/';
my $lib = $dir . 'library';
my $np = 2;
my $numInd = 5;
my $minInd = 4;
my $cov = 30;
my $af = 3;
my $dbP = $dir . "genomes/anolis_Prot.fa";
my $dbG = $dir . "genomes/anolis_DNA.fa";

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

my %renz = (
'AatII' => 'GACGTC',
'Acc651' => 'GGTACC',
'AccI' => 'GT[A|C][G|T]AC',
'AflII' => 'CTTAAG',
'AhdI' => 'GAC[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]GTC',
'AleI' => 'CAC[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]GTG',
'AluI' => 'AGCT',
'AlwI' => 'GGATC',
'AlwNI' => 'CAG[A|T|G|C][A|T|G|C][A|T|G|C]CTG',
'ApaI' => 'GGGCCC',
'ApaLI' => 'GTGCAC',
'ApoI' => '[A|G]AATT[C|T]',
'AscI' => 'GGCGCGCC',
'AseI' => 'ATTAAT',
'AsiSI' => 'GCGATCGC',
'AvaI' => 'C[C|T]CG[A|G]G',
'AvaII' => 'GG[A|T]CC',
'BaeGI' => 'G[G|T]GC[A|C]C',
'BamHI' => 'GGATCC',
'BanI' => 'GG[C|T][A|G]CC',
'BanII' => 'G[A|G]GC[C|T]C',
'BccI' => 'CCATC',
'BclI' => 'TGATCA',
'BfaI' => 'CTAG',
'BfuCI' => 'GATC',
'BglI' => 'GCC[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]GGC',
'BglII' => 'AGATCT',
'BlpI' => 'GCT[A|T|G|C]AGC',
'BsaAI' => '[C|T]ACGT[A|G]',
'BsaBI' => 'GAT[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]ATC',
'BsaHI' => 'G[A|G]CG[C|T]C',
'BsaI' => 'GGTCTC',
'BsaJI' => 'CC[A|T|G|C][A|T|G|C]GG',
'BsiEI' => 'CG[A|G][C|T]CG',
'BsiHKAI' => 'G[A|T]GC[A|T]C',
'BslI' => 'CC[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]GG',
'BsmAI' => 'GTCTC',
'BsmI' => 'GAATGC',
'BsoBI' => 'C[C|T]CG[A|G]G',
'BspDI' => 'ATCGAT',
'BspEI' => 'TCCGGA',
'BspHI' => 'TCATGA',
'BspQI' => 'GCTCTTC',
'BsrBI' => 'CCGCTC',
'BsrFI' => '[A|G]CCGG[C|T]',
'BsrGI' => 'TGTACA',
'BsrI' => 'ACTGG',
'BssHII' => 'GCGCGC',
'BstBI' => 'TTCGAA',
'BstEII' => 'GGT[A|T|G|C]ACC',
'BstNI' => 'CC[A|T]GG',
'BstUI' => 'CGCG',
'BstXI' => 'CCA[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]TGG',
'BstYI' => '[A|G]GATC[C|T]',
'BstZ17I' => 'GTATAC',
'Bsu36I' => 'CCT[A|T|G|C]AGG',
'BtgI' => 'CC[A|G][C|T]GG',
'BtsCI' => 'GGATG',
'BtsI' => 'GCAGTG',
'CfoI' => 'GCGC',
'ClaI' => 'ATCGAT',
'CspCI' => 'CAA[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]GTGG',
'CviQI' => 'GTAC',
'DdeI' => 'CT[A|T|G|C]AG',
'DpnI' => 'GATC',
'DpnII' => 'GATC',
'DraI' => 'TTTAAA',
'DraIII' => 'CAC[A|T|G|C][A|T|G|C][A|T|G|C]GTG',
'EagI' => 'CGGCCG',
'EarI' => 'CTCTTC',
'Eco53kI' => 'GAGCTC',
'EcoNI' => 'CCT[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]AGG',
'EcoO109I' => '[A|G]GG[A|T|G|C]CC[C|T]',
'EcoP15I' => 'CAGCAG',
'EcoRI' => 'GAATTC',
'EcoRII' => 'CC[A|T]GG',
'EcoRV' => 'GATATC',
'FokI' => 'GGATG',
'FspI' => 'TGCGCA',
'HaeII' => '[A|G]GCGC[C|T]',
'HaeIII' => 'GGCC',
'HhaI' => 'GCGC',
'HincII' => 'GT[C|T][A|G]AC',
'HindIII' => 'AAGCTT',
'HinfI' => 'GA[A|T|G|C]TC',
'HinP1I' => 'GCGC',
'HpaI' => 'GTTAAC',
'HpaII' => 'CCGG',
'HphI' => 'GGTGA[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]',
'Hpy166II' => 'GT[A|T|G|C][A|T|G|C]AC',
'Hpy188I' => 'TC[A|T|G|C]GA',
'Hpy188III' => 'TC[A|T|G|C][A|T|G|C]GA',
'HpyCH4IV' => 'ACGT',
'KpnI' => 'GGTACC',
'MboI' => 'GATC',
'MfeI' => 'CAATTG',
'MluI' => 'ACGCGT',
'MlyI' => 'GAGTC',
'MnlI' => 'CCTC',
'MseI' => 'TTAA',
'MslI' => 'CA[C|T][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|G]TG',
'MspA1I' => 'C[A|C]GC[G|T]G',
'MspI' => 'CCGG',
'MwoI' => 'GC[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]GC',
'NaeI' => 'GCCGGC',
'NarI' => 'GGCGCC',
'Nb.BbvCI' => 'CCTCAGC',
'Nb.BsmI' => 'GAATGC',
'Nb.BsrDI' => 'GCAATG',
'Nb.BtsI' => 'GCAGTG',
'NciI' => 'CC[C|G]GG',
'NcoI' => 'CCATGG',
'NgoMIV' => 'GCCGGC',
'NheI' => 'GCTAGC',
'NlaIII' => 'CATG',
'NotI' => 'GCGGCCGC',
'NruI' => 'TCGCGA',
'NsiI' => 'ATGCAT',
'Nt.AlwI' => 'GGATC',
'Nt.BbvCI' => 'CCTCAGC',
'Nt.BsmAI' => 'GTCTC',
'Nt.BspQI' => 'GCTCTTC',
'Nt.BstNBI' => 'GAGTC',
'PaeR71' => 'CTCGAG',
'PflFI' => 'GAC[A|T|G|C][A|T|G|C][A|T|G|C]GTC',
'PflMI' => 'CCA[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]TGG',
'PleI' => 'GAGTC',
'PmeI' => 'GTTTAAAC',
'PmlI' => 'CACGTG',
'PpuMI' => '[A|G]GG[A|T]CC[C|T]',
'PshAI' => 'GAC[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]GTC',
'PspGI' => 'CC[A|T]GG',
'PspOMI' => 'GGGCCC',
'PstI' => 'CTGCAG',
'PvuI' => 'CGATCG',
'PvuII' => 'CAGCTG',
'RsaI' => 'GTAC',
'RsrII' => 'CGG[A|T]CCG',
'SacI' => 'GAGCTC',
'SacII' => 'CCGCGG',
'SalI' => 'GTCGAC',
'Sau3AI' => 'GATC',
'Sau96I' => 'GG[A|T|G|C]CC',
'Sbf1' => 'CCTGCAGG',
'ScaI' => 'AGTACT',
'ScrFI' => 'CC[A|T|G|C]GG',
'SfiI' => 'GGCC[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]GGCC',
'SfoI' => 'GGCGCC',
'SgrAI' => 'C[A|G]CCGG[C|T]G',
'SmaI' => 'CCCGGG',
'SmlI' => 'CT[C|T][A|G]AG',
'SnaBI' => 'TACGTA',
'SpeI' => 'ACTAGT',
'SphI' => 'GCATGC',
'SspI' => 'AATATT',
'SstI' => 'CCGCGG',
'StuI' => 'AGGCCT',
'StyI' => 'CC[A|T][A|T]GG',
'SwaI' => 'ATTTAAAT',
'TaqI' => 'TCGA',
'TfiI' => 'GA[A|T]TC',
'TliI' => 'CTCGAG',
'Tsp509I' => 'AATT',
'TspRI' => '[A|T|G|C][A|T|G|C]CA[C|G]TG[A|T|G|C][A|T|G|C]',
'Tth111I' => 'GAC[A|T|G|C][A|T|G|C][A|T|G|C]GTC',
'XbaI' => 'TCTAGA',
'Xcm1' => 'CCA[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]TGG',
'Xho1' => 'CTCGAG',
'Xma1' => 'CCCGGG',
'Xmn1' => 'GAA[A|T|G|C][A|T|G|C][A|T|G|C][A|T|G|C]TTC',
	);

mkdir($finaldir) unless(-d $finaldir);
mkdir($tmpdir) unless(-d $tmpdir);
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

	#get mpileup from both
	my $pile1 = runPileup($seq,$dir,$contacts[$i][0],$out1,$finaldir);
	my $pile2 = runPileup($seq,$dir,$contacts[$i][1],$out2,$finaldir);
	
	#define fixed SNPs
	#get intron-exon boundaries to help design primers
	my $annoref = makeHash($dbP); 
	fixedSNPs($name,$tmpdir,$pile1,$pile2,$cov,$af,$seq,$finaldir,\%renz,$dbG,$dbP,$annoref);
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
	my ($seq,$dir,$contact,$outfiles,$finaldir) = @_;
	my @out = @{$outfiles};
	my $call1 = system("samtools faidx $seq") unless(-f $seq . ".fai");
	my $out = $finaldir . $contact . ".mpileup.out";
	unless (-f $out) {
		my $files = join("\t", @out);
		my $call2 = system("samtools mpileup -A -f $seq $files > $out");	
		}
	return($out);
	}
	
sub runMapping {
	my ($seq,$file1,$file2,$fileu,$subdir,$lib,$dir) = @_;
	my $out = $subdir . $lib . ".sorted";
	print $out, "\n";
	my $call1 = system("bowtie2 -x $seq -1 $file1 -2 $file2 -S bowtie1.sam -5 5 -3 5 --sensitive -k 10 -X 300 -p $np");
	my $call2 = system("bowtie2 -x $seq $fileu -S bowtie2.sam -5 5 -3 5 --sensitive -k 10  -p $np");
	my $call3 = system("samtools view -bS bowtie1.sam > bowtie1.bam");
	my $call4 = system("samtools view -bS bowtie2.sam > bowtie2.bam");
	my $call5 = system("samtools merge bowtie.bam bowtie1.bam bowtie2.bam");
	my $call6 = system("samtools sort bowtie.bam $out");
	my $call7 = system("rm bowtie1.sam bowtie2.sam bowtie.bam bowtie1.bam bowtie2.bam");	
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

sub fixedSNPs {
	my ($name,$tmpdir,$pileup1,$pileup2,$cov,$af,$seq,$finaldir,$renz,$dbG,$dbP,$annoref) = @_;
	
	splitPileup($pileup1,$tmpdir,1);	
	splitPileup($pileup2,$tmpdir,2);

	my $out = $finaldir . $name . '.RFLPsnps.out';
	open(OUT, ">$out"); 
	my $refseq = parseSeq($seq);
	my %seq = %$refseq;
	foreach my $c (keys %seq) {
		my $tmp1 = $tmpdir . $c . '_1.out';
		my $tmp2 = $tmpdir . $c . '_2.out';
		if (-f $tmp1 && -f $tmp2) {
			my $p1 = parsePileup($tmp1,$cov);
			my $p2 = parsePileup($tmp2,$cov);
			
			my $snp = comparePileup($p1,$p2,$af);
			$snp = annotateSNP($refseq,$snp);
	
			my %snp = %$snp;
			foreach my $c (sort {$a cmp $b} keys %snp) {
				foreach my $pos (sort {$a <=> $b} keys %{$snp{$c}}) {
					if ($snp{$c}{$pos}{'type'} eq "fixed") {
						my $anonHash = $snp{$c}{$pos};
						my $re = cuttable($c,$name,$anonHash,$pos,$renz,$seq{$c}{'seq'},$dbG,$dbP,$annoref);
						if (scalar(@{$re}) > 0) {
							print OUT $c, "\t", $pos, "\t", $snp{$c}{$pos}{'a1'}, "\t", $snp{$c}{$pos}{'a2'}, "\t", $snp{$c}{$pos}{'type'}, "\t", $snp{$c}{$pos}{'loc'}, "\t", $snp{$c}{$pos}{'coding'}, "\t";
							print OUT join("\t", @{$re});
							print OUT "\n";
							}
						}
					}
				}	
			}	
		}
	close(OUT);		
	}	
	
sub cuttable {
	my ($c,$name,$out,$pos,$reref,$seq,$databaseG,$databaseP,$annoref) = @_;
	my %renz = %$reref;
	
	my %snp = %{$out};

	#need to introduce mutation...
	my $start = substr $seq, 0, $pos - 1;
	my $end = substr $seq, $pos;
	#need to introduce my mutation1
	my $seq1 = $start . $snp{'a1'} . $end;
	#need to introduce my mutation2
	my $seq2 = $start . $snp{'a2'} . $end;
				
	my $seq1r = reverse($seq1);
	my $seq2r = reverse($seq2);
	$seq1r =~ tr/ATGCatgc/TACGtacg/;
	$seq2r =~ tr/ATGCatgc/TACGtacg/;
	
	my @renz;

	foreach my $re (keys %renz) {
		#need to make 4 pairwise comparisons
		my $match1 = $seq1 =~ m/($renz{$re})/g;	
		my $match1r = $seq1r =~ m/($renz{$re})/g;	
		my $match2 = $seq2 =~ m/($renz{$re})/g;	
		my $match2r = $seq2r =~ m/($renz{$re})/g;	
		if ($match1 ne $match2 || $match1r ne $match2r) {
			push(@renz,$re);
			}
		}

	if (scalar(@renz)) {
		my $seqout = $finaldir . $name . "." . $c . "_" . $pos . ".fa";		
		open(SEQ, ">$seqout");
	
		my $info = annotateProt($annoref,$databaseP,$databaseG, $seq1);
		print SEQ ">seq1\n", $seq1, "\n", ">seq2\n", $seq2, "\n";
		foreach(@$info) {
			print SEQ $_;
			}
		close(SEQ);	
		}
	return(\@renz);
	}

sub annotateProt {
	my ($anno,$databaseP,$databaseG,$seq1) = @_;
	
	my $query = "query.fa";
	open(OUT1, ">$query");
	print OUT1 ">seq1\n", $seq1, "\n"; 
	close(OUT1);
	
	my %seq;
	my @call = `blastall -p blastx -d $databaseP -i $query -m 8`;	
	if (scalar(@call) > 0) {	
 		my @d = split(/\t/,$call[0]);
		%seq = ('gene' => $d[1], 'eval' => $d[10], 'id' => $d[0]);
		}

	my %anno = %$anno;
		
	#need to use exonerate to define utr etc; call to external program
	my $dbP = Bio::DB::Fasta->new($databaseP);
	my $dbG = Bio::DB::Fasta->new($databaseG);

	my @final;	

	if ($seq{'id'}) {				
		#need to run exonerate			
		my $anolisID = $seq{'gene'};			
		my $DNA = $dbG->get_Seq_by_id($anno{$anolisID}{'contig'})->subseq($anno{$anolisID}{'start'} => $anno{$anolisID}{'end'});
		my $seq = $dbP->get_Seq_by_id($anolisID)->seq;	

		my $target2 = "target2.fa";
	    my $query2 = "query2.fa";
	    	
	    open(T2, ">$target2");
	    open(Q2, ">$query2");
		print T2 ">$anolisID\n$DNA\n";
		print Q2 ">protein\n$seq\n";
		my @call = `exonerate --model protein2genome $query2 $target2 --showvulgar no --showalignment no --showtargetgff yes`;
	   	my %cds;
		foreach(@call){
			if ($_ =~ m/orientation/) {
		    	$cds{'orient'} = $1 if $_ =~ m/([\+|\-])\s+\./;
				}
			elsif ($_ =~ m/exon\s+/) {
				my $bounds = $1 if $_ =~ m/exon\t([0-9]+\t[0-9]+)/;;
				my @bounds = split("\t", $bounds);
		    	$cds{'exon'}{$bounds[0]}=$bounds[1];
				}
			}
				
		
		my $line1 = ">" . $anolisID . "\n" . $seq . "\n";
		my $line2 = ">" . $anolisID . "\n" . $DNA . "\n";			
		push(@final,$line1,$line2);

		#need to allow for a failure to find cds 
		if (keys %cds){
			my $tracker = 1;
			if ($cds{'orient'} eq '+') {
				foreach(sort {$a <=> $b} keys %{$cds{'exon'}}){
				    my $start = $_;
		   			my $end = $cds{exon}{$start};
		   			my $length = $end - $start + 1;
		   			$start = $start - 1;
		   			my $sub = substr $DNA, $start, $length;
				    my $loc = $anno{$anolisID}{'contig'};
				    my $loc_s = $anno{$anolisID}{'start'} + $start;
				    my $loc_e = $loc_s + $length - 1;
				    $loc = $loc . '_' . $loc_s . '_' . $loc_e;

	   				my $line3 = ">" . $anolisID . "exon" . $tracker . "\t" . $loc . "\n" . $sub . "\n";
					push(@final,$line3);
    				$tracker++;
					}
				}
			else {
				foreach(sort {$b <=> $a} keys %{$cds{'exon'}}){
				    my $start = $_;
	    			my $end = $cds{exon}{$start};
	    			my $length = $end - $start + 1;
	    			$start = $start - 1;
	    			my $sub = substr $DNA, $start, $length;
	    			$sub = reverse($sub);
	    			$sub =~ tr/ATGCatgc/TACGtacg/;

				    my $loc = $anno{$anolisID}{'contig'};
				    my $loc_s =$anno{$anolisID}{'start'} + $start;
					my $loc_e = $loc_s + $length - 1;                       
					$loc = $loc. '_' .$loc_s . '_' . $loc_e;

					my $line3 = ">" . $anolisID . "exon" . $tracker . "\t" . $loc . "\n" . $sub . "\n";
					push(@final,$line3);
		   			$tracker++;
					}
				}
			}
		}		
	return(\@final);
	}
	
	
sub splitPileup {
	my ($pile, $tmpdir,$id) = @_;
	
	my $c = 'NA';
	open(IN, "<$pile");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/(contig\S+)/) {
			if ($1 eq $c) {
				print OUT $line, "\n";
				}
			else {
				#new file
				if ($c eq 'NA') {
					$c = $1;
					my $file = $tmpdir . $c . "_" . $id . '.out';
					open(OUT, ">$file");
					print OUT $line, "\n";
					}
				else {
					$c = $1;
					close(OUT);
					my $file = $tmpdir . $c . "_" . $id . '.out';
					open(OUT, ">$file");
					print OUT $line, "\n";
					}
				}
			}
		}
	close(IN);
	}	
	
sub annotateSNP {
	my ($seq,$snp) = @_;
	
	#define location of SNP
	my %snp = %$snp; 
	my %seq = %$seq;
	foreach my $c (keys %snp) {
		foreach my $pos (keys %{$snp{$c}}) {
			#is this in cds
			my $loc;
			if ($pos <= $seq{$c}{'ge'} && $pos <= $seq{$c}{'ge'}) {
				$loc = 'cds';
				}
			elsif ($pos > $seq{$c}{'ge'}) {
				if ($seq{$c}{'3u'}) {
					$loc = '3u';
					}
				else {
					$loc = 'undef';
					}
				}	
			elsif ($pos < $seq{$c}{'gs'}) {
				if ($seq{$c}{'5u'}) {
					$loc = '5u';
					}
				else {
					$loc = 'undef';
					}
				}		
			
			my $coding;
			if ($loc eq 'cds') {
				my $start = substr $seq{$c}{'seq'}, 0, $pos - 1;
				my $end = substr $seq{$c}{'seq'}, $pos;
				#need to introduce my mutation1
				my $orf1 = $start . $snp{$c}{$pos}{'a1'} . $end;
				#need to introduce my mutation2
				my $orf2 = $start . $snp{$c}{$pos}{'a2'} . $end;
				
				#this could lead to a SNP change -- let's see!
				my $gs = $seq{$c}{'gs'} - 1;
				my $ln = $seq{$c}{'ge'} - $seq{$c}{'gs'} + 1;
				$orf1 = substr $orf1, $gs, $ln;
				$orf2 = substr $orf2, $gs, $ln;
				my $aa1 = translate($orf1);
				my $aa2 = translate($orf2);
				if ($aa1 eq $aa2) {
					$coding = 'syn';
					}
				else {
					$coding = 'ns';
					}		
				}
			else {
				$coding = 'noncoding';
				}
			$snp{$c}{$pos}{'loc'} = $loc;
			$snp{$c}{$pos}{'coding'} = $coding;
			}
		}	
	return(\%snp);	
	}
		
sub parseSeq {
	my ($s) = @_;
	my %s;
	open(IN, "<$s");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+).*ENS.*/) {
			my $c = $1; my @d = split(/\t/,$line);
			chomp(my $seq = <IN>);
			$s{$c}{'seq'} = $seq;
			$s{$c}{'match'} = $d[2];
			
			if ($d[1] =~ m/5u(\d+)/) {
				$s{$c}{'5u'} = $1;
				}
			if ($d[1] =~ m/gs(\d+)/) {
				$s{$c}{'gs'} = $1;
				}
			if ($d[1] =~ m/ge(\d+)/) {
				$s{$c}{'ge'} = $1;
				}
			if ($d[1] =~ m/3u(\d+)/) {
				$s{$c}{'3u'} = $1;
				}	
			}
		}
	close(IN);	
	return(\%s);	
	}	

sub comparePileup {
	my ($p1,$p2,$af) = @_;
	my %p1 = %$p1; my %p2 = %$p2;
	
	my %snp;
	
	foreach my $c (sort {$a cmp $b} keys %p1) {
		foreach my $pos (sort {$a <=> $b} keys %{$p1{$c}}) {
			if ($p2{$c}{$pos}) {
				foreach my $ref (keys %{$p1{$c}{$pos}}) {		
				 	my $a1 = parseReads($ref,$p1{$c}{$pos}{$ref},$af);
				 	my $a2 = parseReads($ref,$p2{$c}{$pos}{$ref},$af);
				 	
				 	my @a1 = @{$a1}; my @a2 = @{$a2};
				 	my %a1; my %a2;
				 	for (my $i = 0; $i < scalar(@a1); $i++) {
				 		$a1{$a1[$i][0]}++; $a1{$a1[$i][1]}++;	
						}
					for (my $i = 0; $i < scalar(@a2); $i++) {
						$a2{$a2[$i][0]}++; $a2{$a2[$i][1]}++;
						}	
				 		 	
				 	if (scalar(keys %a1) == 1 && scalar(keys %a2) == 1) {
				 		#are they the same base?
				 		my @b1 = keys %a1; my @b2 = keys %a2;
				 		#if not, no SNP
				 		if ($b1[0] ne $b2[0]) {
				 			$snp{$c}{$pos} = {'type' => 'fixed', 'a1'=>$b1[0],'a2'=>$b2[0]};
				 			}
				 		}
				 	elsif (scalar(keys %a1) > 1 && scalar(keys %a2) == 1) {	
				 		my @b1 = sort {$a1{$b} <=> $a1{$a}} keys %a1;
				 		$snp{$c}{$pos} = {'type'=> 'poly1', 'a1'=>$b1[0],'a2'=>$b1[1]};
				 		}
				 	elsif (scalar(keys %a1) == 1 && scalar(keys %a2) > 1) {	
				 		my @b2 = sort {$a2{$b} <=> $a2{$a}} keys %a2;
				 		$snp{$c}{$pos} = {'type'=> 'poly2', 'a1'=>$b2[0],'a2'=>$b2[1]};
				 		}	
				 	else {
				 		my $shared;
				 		foreach (keys %a1) {
				 			$shared++ if $a2{$_};
				 			}
				 		if ($shared > 1) {
				 			my @b1 = sort {$a1{$b} <=> $a1{$a}} keys %a1;
				 			$snp{$c}{$pos} = {'type'=>'shared','a1'=>$b1[0],'a2'=>$b1[1]};
				 			}
				 		}
					}
				}
			}
		}
	return(\%snp);	
	}
	
sub parseReads {
	my ($ref,$r,$af) = @_;
	my @reads = @{$r}; my @a;
	for (my $i = 0; $i < scalar(@reads); $i++) {
		my %a;
		while ($reads[$i] =~ /([atgcATGC])/g) {
			$a{uc($1)}++;
			}
		foreach my $a (keys %a) {
			if ($a{$a} >= $af) {
				my @ref = $reads[$i] =~ m/([\.|\,])/g;
				if (scalar(@ref) >= $af) {
					#looks to be polymorphic
					$a[$i][0] = $ref; $a[$i][1] = $a;
					}
				else {
					#looks to be fixed
					$a[$i][0] = $a; $a[$i][1] = $a;
					}
				}
			else {
				$a[$i][0] = $ref; $a[$i][1] = $ref;
				}
			}
		unless ($a[$i][0]) {
			$a[$i][0] = $ref; $a[$i][1] = $ref;
			}
		}
	return(\@a);	
	}
	
sub parsePileup {
	my ($file,$cov) = @_;
	
	my %snp;
	open(IN, "<$file");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		my $c = $d[0]; my $pos = $d[1]; my $ref = $d[2];
		#checking to see what the coverage is at each individual
		my $pass = 0; my @reads;
		for (my $i = 0; $i < $numInd; $i++) {
			$pass++ if $d[3 + 3*$i] >= $cov;
			my $read = $d[4 + 3*$i];
			while ($read =~ /([\+|-](\d+))/g) {
				my $n = $2;
				$read =~ s/[\+|\-]\d+[atgcn]{$n}//i;
				}
			push(@reads, $read);
			}
		#have high enough coverage at the number of required individuals	
		if ($pass >= $minInd) {	
			$snp{$c}{$pos}{$ref} = \@reads;
			}
		}
	
	close(IN);

	return(\%snp);
	}
	
sub translate {
	my $string = shift;
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
	
sub makeHash {
	my ($database) = @_;
		
	my %anno;
	open(ANNO, "<$database");
	while(<ANNO>){
		chomp(my $line = $_);
		if ($line =~ m/^>(\S+)/) {
			my $id = $1;
			my @a = split(/\s/,$line);
			my @d = split(/:/,$a[2]);
			my $contig = $d[2];
			my $start = $d[3];
			my $end = $d[4];
			$anno{$id} = {'contig'=>$contig, 'start'=>$start, 'end'=>$end};
			}
		}	
	close(ANNO);

	return(\%anno)
	}	

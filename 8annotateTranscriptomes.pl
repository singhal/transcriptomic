###########################################################################################
# a script to annotate your assemblies via a "weak" recip blast method & to identify and  #
# break up chimeric contigs                                                               #
# external dependencies: blastall/formatdb (>2.2.17), exonerate, cd-hit-est, framedp      #
# this script requires framedp, which needs to have a formatted cfg file and whose        #
#	location needs to be identified via the bash/csh/etc profile			  #
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 29 January 2012                #
###########################################################################################

use warnings;
use strict;
use Bio::DB::Fasta;
use List::Util qw[min max];

#how many processors can you use for blasting?
my $np = '4';
#the assemblies you hope to annotate
my @assemblies = </media/DataDrive/sutureGenomics/*/assemblies/*trinity*final>;
#the database we will use to annotate results
my $dbP = '/media/DataDrive/sutureGenomics/protDB/annoProt/anolisProt.fa';
#not a match unless higher than this
my $evalue = '1e-20';
my $dir = '/media/DataDrive/sutureGenomics/';
#directory where framedp is 
my $framedir = '/home/singhal/framedp/';

########################
# run the subroutines! #
########################

#formats the protein database unless it already has been done
unless (-f $dbP . '.pin') {
	my $call = system("formatdb -i $dbP");
	}

foreach my $assembly (@assemblies) {
    my $seqout = $assembly . ".annotated";
    print "Doing assembly $assembly now!\n";
    unless (-f $seqout) {
		my $outfile1 = blastProteins($assembly); 
	        print "Done blasting 1 assembly $assembly now!\n";
		my $outfile2 = chimericTest($assembly,$outfile1);
		my $outfile3 = blastProteins($outfile2);

		print "Done blasting 2 assembly $assembly now!\n";
		my $recip = recipBlast($outfile2);
		print "Done recip blasting assembly $assembly now!\n";
	       	
		my ($seqref,$annoref) = makeHash($outfile2,$dbP); 
		my $seq = annotateProt($outfile3,$seqref,$annoref,$recip);
		print "Done annotating assembly $assembly now!\n";
		$seq = orfFinesse($seq);

########################
# report the results!  #
########################

		my $seqout = $assembly . ".annotated";
		open (SEQOUT, ">$seqout");
		my %seq = %$seq;
		foreach my $id (sort {$a cmp $b} keys %seq) {
		    print SEQOUT ">", $id, "\t";
		    if ($seq{$id}{'info'}) {
				print SEQOUT $seq{$id}{'info'}, "\t" if $seq{$id}{'info'};
				print SEQOUT $seq{$id}{'match'}[0]{'gene'}, "\t" if $seq{$id}{'match'};
				print SEQOUT $seq{$id}{'desc'}, "\t" if $seq{$id}{'desc'};
				print SEQOUT $seq{$id}{'match'}[0]{'eval'}, "\t" if $seq{$id}{'match'};
			    }
		    print SEQOUT "\n";
		    print SEQOUT $seq{$id}{'seq'}, "\n";
		}
		my $call = system("rm formatdb.log error.log");
    }
}

###########################
# behold the subroutines  #
###########################

sub recipBlast {
	my ($seq) = @_;
	
	my $recip = $seq . '.recipBlast.out';
	my $call1 = system("formatdb -i $seq -p F");
	my $call2 = system("blastall -p tblastn -d $seq -i $dbP -a $np -e $evalue -m 8 -o $recip -b 10");
	
	my %r; my $score;
	open(IN, "<$recip");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		if ($r{$d[0]}) {
			if ($d[11] == $score) {
				$r{$d[0]}{$d[1]}++;
				}
			}
		else {
			$r{$d[0]}{$d[1]}++;
			$score = $d[11];
			}
		}
	close(IN);
	unlink($recip);
	my $call3 = system("rm $seq.n*");
	
	return(\%r);
	}

sub makeHash {
	my ($assembly,$database) = @_;
	my (%seq, $id);
	
	open(IN, "<$assembly");
	while(<IN>){
		chomp(my $line = $_);
		if ($line =~ m/>(.*)/){
			$id = $1;
			}
		else {
			$seq{$id}{'seq'} .= $line;
			}
		}
	close(IN);
	
	my %anno;
	open(ANNO, "<$database");
	while(<ANNO>){
		chomp(my $line = $_);
		if ($line =~ m/^>(\S+)/) {
			my $id = $1;
			my ($gene,$info);
			if ($line =~ m/^>EN\S+\s+\S+_([A-Z]+)/) {
				$gene = $1;
				}
			else {
				$gene = 'NA';
				}
			if ($line =~ m/^>EN\S+\s+\S+\s+(\S+.*)$/) {
				$info = $1;
				}
			else {
				$info = 'NA';
				}
			$anno{$id} = {'gene' => $gene, 'info'=> $info};
			}
		}	
	close(ANNO);	

	return(\%seq,\%anno)
	}

sub blastProteins {
	my ($assembly) = @_;
		
	my $masterout = $assembly . '.blast.out';
	my $call = system("blastall -p blastx -d $dbP -i $assembly -a $np -e $evalue -m 8 -o $masterout -b 30");
	
	return($masterout);	
	}

sub parseBlast {
	my ($seq, $out) = @_;
	
	my %seq = %$seq;

	my $tracker;
	open(IN, "<$out");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		#if the line is for the same contig as before or not
		if ($seq{$d[0]}{'match'}) {
			#want to consider other high scoring matches
			unless ($d[1] eq $seq{$d[0]}{'match'}[$tracker]{'gene'}) {				
				#only consider addt'l matches if evalue is good and if there isn't a sharp decline in quality of match
				$d[10] = 1e-200 if $d[10] =~ m/^0.0$/;
				if ($d[10]/$seq{$d[0]}{'match'}[0]{'eval'} <= 10000) {
				    if (scalar(@{$seq{$d[0]}{'match'}}) < 10){
						push @{$seq{$d[0]}{'match'}}, {'gene' => $d[1], 'eval' => $d[10]};
						$tracker++;
				    	}
					}	
				}			
			}
		#this is a new contig	
		else {
			#will be dividing by this later, so it cannot equal 0!
			$d[10] = 1e-200 if $d[10] =~ m/^0.0$/; 
			push @{$seq{$d[0]}{'match'}}, {'gene' => $d[1], 'eval' => $d[10]};	
			$tracker = 0;
			}
		}											
	close(IN);

	return(\%seq);
	}
		
sub annotateProt {
	my ($out, $seq, $anno,$recip) = @_;
	
	$seq = parseBlast($seq,$out);	

	my %seq = %$seq; my %anno = %$anno; my %recip = %$recip; 
		
	#need to use exonerate to define utr etc; call to external program
	my $db = Bio::DB::Fasta->new($dbP);
		
	foreach my $id (keys %seq) {
		if ($seq{$id}{'match'}) {
			if ($recip{ $seq{$id}{'match'}[0]{'gene'} }{ $id }) {
				#exonerate first
				my $query = "query.fa";
				my $target = "target.fa";
				open(QUERY, ">$query");
				open(TARGET, ">$target");
				print QUERY ">$id\n$seq{$id}{'seq'}\n";
				my $seq = $db->get_Seq_by_id(${$seq{$id}{'match'}}[0]->{'gene'})->seq;
				my $protlength = $db->get_Seq_by_id(${$seq{$id}{'match'}}[0]->{'gene'})->length;
				my $contiglength = length($seq{$id}{'seq'});			
				print TARGET ">gene\n$seq\n";
				close(QUERY); close(TARGET);					
				my @call3 = `exonerate $target $query -m protein2genome --showalignment no --showcigar 0`;
				unlink($query); unlink($target);
				my $info;
				my @d = split(/\s+/,$call3[2]);
				if ($d[8] eq '+' || $d[8] eq '-') {
					if ($d[8] eq '+') {
						#this is in 5->3
						#identify gene start and stop
						my $start = $d[6] + 1;
						my $end = $d[7];
						$info = 'gs' . $start . '_ge' . $end;
						if ($d[3]/$protlength > 0.9 || $protlength - $d[3] < 11) {
							#yes, i am going to call it a  3' utr
							my $utr3 = $end+1;
							$info = $info . '_3u' . $utr3;
							}
						if ($d[2]/$protlength < 0.1 || $d[2] < 11) {
							#yes, i am going to call it a 5' utr
							my $utr5 = $start - 1;
							$info = '5u' . $utr5 . "_" . $info;
							}
						}		
					else {
						#this is in 3->5;
						$seq{$id}{'seq'} = reverse($seq{$id}{'seq'});
						$seq{$id}{'seq'} =~ tr/ATGC/TACG/;	
						
						my $gs = $contiglength - $d[6] + 1;
						my $ge = $contiglength - $d[7];
						$info = 'gs' . $gs . '_ge' . $ge;
					
						if ($d[3]/$protlength > 0.9 || $protlength - $d[3] < 11) {
							#yes, i am going to call it a  3' utr
							my $utr3 = $ge+1;
							$info = $info . '_3u' . $utr3;
							}
						if ($d[2]/$protlength < 0.1 || $d[2] < 11) {
							#yes, i am going to call it a 5' utr
							my $utr5 = $gs - 1;
							$info = '5u' . $utr5 . '_' . $info;
							}				
						}
					}
				#a blast match but no exonerate match?	
				else {
					print "Huh, this is odd. BLAST hit but not exonerate hit for $id?\n";
					}

				$seq{$id}{'info'} = $info;	
				#now define the gene name	
				foreach my $hashref (@{$seq{$id}{'match'}}) {
					my $gene =  $hashref->{'gene'};
					if( $anno{$gene} ) {
						unless ( $anno{$gene}{'info'} eq 'NA' ) {								
							$seq{$id}{'desc'} = $anno{$gene}{'info'};
							last;
							}
						}	
					}	
				}
			}				
		}
	return(\%seq)	
	}	
	
sub chimericTest {
	my ($seq,$out) = @_;
	
	my %seq; my $id;
	open(IN, "<$seq");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			$id = $1;
			}
		else {
			$seq{$id} .= $line;
			}
		}
	close(IN);

	#parses the blast output
	open(OUT, "<$out");
	my %match;
	my %chimera;
	my $t;
	while(<OUT>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
	
		#it is a gene that i already have	
		if ($match{$d[0]}) {	
			#it is a match to that gene i already have
			if ($match{$d[0]}[$t]{'gene'} eq $d[1]) {
				push @{$match{$d[0]}[$t]{'start'}}, $d[6];
				#ends
				push @{$match{$d[0]}[$t]{'end'}}, $d[7];
				}
			#it is a new match to that gene
			else {
				$t++;
				push @{$match{$d[0]}[$t]{'start'}}, $d[6];
				#ends
				push @{$match{$d[0]}[$t]{'end'}}, $d[7];
				$match{$d[0]}[$t]{'gene'} = $d[1];
				}
			}
		#it is a new gene
		else {
			$t = 0;
			#starts
			push @{$match{$d[0]}[$t]{'start'}}, $d[6];
			#ends
			push @{$match{$d[0]}[$t]{'end'}}, $d[7];
			$match{$d[0]}[$t]{'gene'} = $d[1];
			}
		}
	close(OUT);	
	
	#figures out which contigs appear to be chimeric	
	foreach my $c (sort {$a cmp $b} keys %match) {
		my %l;
		my @m = @{$match{$c}};
	
		my $start = min(@{$m[0]{'start'}},@{$m[0]{'end'}});
		my $end = max(@{$m[0]{'start'}},@{$m[0]{'end'}});
		for (my $j = $start; $j <= $end; $j++) {
			$l{$c}{$j}++;
			}
		
		for (my $i = 1; $i < scalar(@m); $i++) {	
			my $start = min(@{$m[$i]{'start'}},@{$m[$i]{'end'}});
			my $end = max(@{$m[$i]{'start'}},@{$m[$i]{'end'}});
			my $chimera = 1;
			for (my $j = $start; $j <= $end; $j++) {
				$chimera = 0 if $l{$c}{$j};	
				}
			if ($chimera) {	
				$chimera{$c}{$m[0]{'gene'}}++;
				$chimera{$c}{$m[$i]{'gene'}}++;
				for (my $j = $start; $j <= $end; $j++) {
					$l{$c}{$j}++;
					}
				}
			}
		}	
	
	#determines if these chimeras are real
	my $db = Bio::DB::Fasta->new($dbP);
	my %newseq; my $tracker = 1;
	foreach my $c (keys %chimera) {	
		my $query = "query.fa";
		my $target = "target.fa";
		open(QUERY, ">$query");
		open(TARGET, ">$target");
		print QUERY ">$c\n$seq{$c}\n";
		foreach my $c2 (keys %{$chimera{$c}}) {
			my $seq = $db->get_Seq_by_id($c2)->seq;
			print TARGET ">$c2\n$seq\n";
			}
		close(QUERY); close(TARGET);					
	
		my @call = `exonerate $target $query -m protein2genome --showalignment no --showcigar 0`;
		my @s; my @e;
		if (@call) {
			my %m;
			for (my $i = 2; $i < scalar(@call) - 1; $i++) {
				my @d = split(/\s+/,$call[$i]);
				unless ($m{$d[1]}) {
					if ($d[8] eq '+') {
						push @s, $d[6];
						push @e, $d[7];
						}
					else {
						push @s, $d[7];
						push @e, $d[6];
						}
					$m{$d[1]}++;
					}
				}		
			}	
		#a blast match but no exonerate match?	
		else {
			print "Huh, this is odd. BLAST hit but no exonerate hit for $id?\n";	
			}
	
		#ensure that the two annotated parts do not overlap
		my $overlap = 0;
		for (my $i = 0; $i < scalar(@s); $i++) {
			for (my $j = 0; $j < scalar(@e); $j++) {
				unless ($i == $j) {
					if (($s[$i] < $s[$j] &&   $s[$j] < $e[$i]) || ($s[$i] < $e[$j] &&  $e[$j] < $e[$i] )) {
						$overlap++;
						}
					}
				}
			}			
		if ($overlap) {
			my $max = 0;
			my $s; my $e;
			for (my $i = 0; $i < scalar(@s); $i++) {
				if ($e[$i] - $s[$i] > $max) {
					$s = $s[$i];
					$e = $e[$i];
					}
				@s = ($s); @e = ($e);
				}			
			}

		#separate the chimeric contig out to separate contigs
		for (my $i = 0; $i < scalar(@s); $i++) {
			my $ln = $e[$i] - $s[$i] + 1;
			$newseq{$tracker} = substr $seq{$c}, $s[$i], $ln; 
			$tracker++;
			}
		delete $seq{$c} if $seq{$c};	
		unlink($target); unlink($query);
		}
	
	my $outfile = $seq . ".noChimera";
	open(OUT, ">$outfile");
	foreach my $c (keys %newseq) {
		print OUT ">contig", $c, "\n", $newseq{$c}, "\n";
		}
	foreach my $c (keys %seq) {
		print OUT ">contig", $tracker, "\n", $seq{$c}, "\n";
		$tracker++;
		}
	close(OUT);	
	
	my $call = system("cd-hit-est -i $outfile -o temp.fa -c 1.00");
	open(IN, "<temp.fa");
	open(OUT2, ">temp2.fa"); 
	my $tracker2 = 1;
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>/) {
			print OUT2 ">contig", $tracker2, "\n";
			$tracker2++;
			}
		else {	
			print OUT2 $line, "\n";
			}
		}
	close(IN); close(OUT); 
	my $call2 = system("rm temp.fa*");
	my $call3 = system("mv temp2.fa $outfile");
	
	return($outfile);
	}
	
sub orfFinesse {
	my ($seq) = @_;
	
	my %seq = %$seq;
	my $frameProt = $dir . 'interruptedORFs.fa';
	open(OUT, ">$frameProt");
	foreach my $c (keys %seq) {
		if ($seq{$c}{'info'}) {
			my $info = $seq{$c}{'info'};
			my $gs = $1 if $info =~ m/gs(\d+)/;
			my $ge = $1 if $info =~ m/ge(\d+)/;
			my $gs0 = $gs - 1;
			my $length = $ge - $gs + 1;
		
			my $s = $seq{$c}{'seq'};
		
			my $def_orf = substr $s, $gs0, $length;
			my $pot_orf = substr $s, $gs0;
		
			my $def_aa = translate($def_orf);
			my $pot_aa = translate($pot_orf);
		
			if ($pot_aa =~ m/\*/) {
				$pot_aa = $1 if $pot_aa =~ m/^([A-Z]+)\*/;
				}
			if ($pot_aa =~ m/^\*/) {
			    $pot_aa = 1;
			}
		
			if (length($pot_aa) >  length($def_aa)) {
			    if (length($pot_aa) - length($def_aa) < 20 || ((length($pot_aa) - length($def_aa))/length($def_aa)) < 0.2) {
				my $newlength = 3 * length($pot_aa);
				my $newend = $gs + $newlength - 1;
				$info =~ s/ge$ge/ge$newend/;
			
				if ($info =~ m/_3u(\d+)/) {
					my $u3 = $1;
					my $u3new = $newend + 1;
					$info =~ s/_3u$u3/_3u$u3new/;
					}
				$seq{$c}{'info'} = $info;		
				}
			}
			elsif (length($pot_aa) < length($def_aa)) {
				#short enough to chop
				if (length($def_aa) - length($pot_aa) < 10 || (length($def_aa) - length($pot_aa))/length($def_aa) < 0.1) {
					my $newlength = 3 * length($pot_aa);
					my $newend = $gs + $newlength - 1;
					$info =~ s/ge$ge/ge$newend/;
			
					if ($info =~ m/_3u(\d+)/) {
						my $u3 = $1;
						my $u3new = $newend + 1;
						$info =~ s/_3u$u3/_3u$u3new/;
						}
					$seq{$c}{'info'} = $info;	
					}
				else {	
				    print OUT ">", $c, "\n", $seq{$c}{'seq'}, "\n";
					}
				}	
			}	
		}
	close(OUT);
	
	#now need to run framedp
	my $outdir = $dir . '/framedpTmp/'; mkdir($outdir) unless(-d $outdir);
	my $call = system("$framedir" . "bin/FrameDP.pl --cfg $framedir" . "cfg/FrameDP.cfg --infile $frameProt --outdir $outdir");
	my @pepdb = <$outdir*pepdb.fa>;
	my @seqdb = <$outdir*seqdb.fa>;
	
	my %orf;
	open(IN, "<$pepdb[0]");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>/) {
			my @d = split(/:/,$line);
			if ($d[-1] eq '+') {
				my $c = $1 if $d[0] =~ m/>(\S+)/;
				my $ln = $d[2] - $d[1] + 1;
				if ($orf{$c}) {
					$orf{$c} = {'ln' => $ln, 'gs' => $d[1], 'ge' => $d[2] - 3} if $ln > $orf{$c}{'ln'};
					}
				else {	
					$orf{$c} = {'ln' => $ln, 'gs' => $d[1], 'ge' => $d[2] - 3};  
					}
				}
			}
		}	
	close(IN);
	
	my %frameseq;
	open(IN, "<$seqdb[0]");
	my $c;
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			$c = $1;
			}
		else {
			$frameseq{$c} .= $line;
			}
		}
	close(IN);	
		
	foreach my $c (keys %orf) {	
		$seq{$c}{'seq'} = $frameseq{$c};
		my $oldinfo = $seq{$c}{'info'};
		my $oldln = $2 - $1 + 1 if $oldinfo =~ m/gs(\d+).*ge(\d+)/;
		my $newln = $orf{$c}{'ge'} - $orf{$c}{'gs'} + 1;
		
    	if ($newln > 0.75 * $oldln) {
			$oldinfo =~ s/gs(\d+)/gs$orf{$c}{'gs'}/;
			$oldinfo =~ s/ge(\d+)/ge$orf{$c}{'ge'}/;
			if ($oldinfo =~ m/_3u(\d+)/) {
				my $u3 = $1;
				my $u3new = $orf{$c}{'ge'} + 1;
				$oldinfo =~ s/_3u$u3/_3u$u3new/;
				}
			if ($oldinfo =~ m/5u(\d+)/) {
				my $u5 = $1;
				my $u5new = $orf{$c}{'gs'} - 1;
				$oldinfo =~ s/5u$u5/5u$u5new/;
				}			
			$seq{$c}{'info'} = $oldinfo; 
			}			
		}
	
	my $callRm = system("rm -r $outdir");	

	my $db = Bio::DB::Fasta->new($dbP);
	#need to check ORF again
	foreach my $c (keys %seq) {
		if ($seq{$c}{'info'}) {
			my $info = $seq{$c}{'info'};
			my $gs = $1 if $info =~ m/gs(\d+)/;
			my $ge = $1 if $info =~ m/ge(\d+)/;
			my $gs0 = $gs - 1;
			my $length = $ge - $gs + 1;
		
			my $s = $seq{$c}{'seq'};
		
			my $def_orf = substr $s, $gs0, $length;
			my $pot_orf = substr $s, $gs0;
		
			my $def_aa = translate($def_orf);
			my $pot_aa = translate($pot_orf);
		
			if ($pot_aa =~ m/\*/) {
				$pot_aa = $1 if $pot_aa =~ m/^([A-Z]+)\*/;
				}
			if ($pot_aa =~ m/^\*/) {
			    $pot_aa = 1;
			}
		
			if (length($pot_aa) < length($def_aa)) {
				#short enough to chop
				if (length($def_aa) - length($pot_aa) < 10 || (length($def_aa) - length($pot_aa))/length($def_aa) < 0.1) {
					my $newlength = 3 * length($pot_aa);
					my $newend = $gs + $newlength - 1;
					$info =~ s/ge$ge/ge$newend/;
			
					if ($info =~ m/_3u(\d+)/) {
						my $u3 = $1;
						my $u3new = $newend + 1;
						$info =~ s/_3u$u3/_3u$u3new/;
						}
					$seq{$c}{'info'} = $info;	
					}
				else {	
					#re-exonerate
					my $query = "query.fa";
					my $target = "target.fa";
					print "re exonerating for contig $c\n";
					open(QUERY, ">$query");
					open(TARGET, ">$target");
					print QUERY ">$c\n$seq{$c}{'seq'}\n";
					my $seq = $db->get_Seq_by_id($seq{$c}{'match'}[0]{'gene'})->seq;
					my $protlength = $db->get_Seq_by_id($seq{$c}{'match'}[0]{'gene'})->length;
					my $contiglength = length($seq{$c}{'seq'});			
					print TARGET ">gene\n$seq\n";
					close(QUERY); close(TARGET);					
					my @call3 = `exonerate $target $query -m protein2genome --showalignment no --showcigar 0`;
					unlink($query); unlink($target);
					my $info;
					my @d = split(/\s+/,$call3[2]);
					if ($d[8] eq '+' || $d[8] eq '-') {
						if ($d[8] eq '+') {
							#this is in 5->3
							#identify gene start and stop
							my $start = $d[6] + 1;
							my $end = $d[7];
							$info = 'gs' . $start . '_ge' . $end;
							if ($d[3]/$protlength > 0.9 || $protlength - $d[3] < 11) {
								#yes, i am going to call it a  3' utr
								my $utr3 = $end+1;
								$info = $info . '_3u' . $utr3;
								}
							if ($d[2]/$protlength < 0.1 || $d[2] < 11) {
								#yes, i am going to call it a 5' utr
								my $utr5 = $start - 1;
								$info = '5u' . $utr5 . "_" . $info;
								}
							$seq{$c}{'info'} = $info;
							}		
						}
					#a blast match but no exonerate match?	
					else {
						print "Huh, this is odd. BLAST hit but not exonerate hit for $c?\n";
						delete($seq{$c}{'info'});
					}
				}
			}
		}
	}
	return(\%seq);	
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
#			print "ERROR: ILLEGAL PASS TO CODON TRANSLATION: $_ is not a codon!\n";
			$translate = $translate . 'X';
			}
		}
	return($translate);
	}	

use warnings;
use strict;
use Bio::DB::Fasta;
use List::Util qw[min max];

#figure out incidence of chimerism -- never did this!

my $out = '/media/DataDrive/sutureGenomics/chimericIncidence.out';  
#how many processors can you use for blasting?
my $np = '2';
#the assemblies you hope to annotate
my @assemblies = </media/DataDrive/sutureGenomics/*/assemblies/*trinity*final>;
#the database we will use to annotate results
my $dbP = '/media/DataDrive/sutureGenomics/protDB/annoProt/anolisProt.fa';
#not a match unless higher than this
my $evalue = '1e-20';

########################
# run the subroutines! #
########################

#formats the protein database unless it already has been done
unless (-f $dbP . '.pin') {
	my $call = system("formatdb -i $dbP");
	}

open(OUT, ">$out");
foreach my $assembly (@assemblies) {
    print "Doing assembly $assembly now!\n";
	my $outfile1 = blastProteins($assembly); 
	print "Done lasting 1 assembly $assembly now!\n";
	my $outfile2 = chimericTest($assembly,$outfile1);
	my $call = system("rm formatdb.log error.log");
	}
close(OUT);	

sub blastProteins {
	my ($assembly) = @_;
		
	my $masterout = $assembly . '.blast.out';
#	my $call = system("blastall -p blastx -d $dbP -i $assembly -a $np -e $evalue -m 8 -o $masterout -b 30");
	
	return($masterout);	
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
	open(IN, "<$out");
	my %match;
	my %chimera;
	my $t;
	while(<IN>) {
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
	close(IN);	
	
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
			print OUT $seq, "\t", $c, "\t", length($seq{$c}), "\t", $tracker, "\t", $ln, "\n";
			}
		delete $seq{$c} if $seq{$c};	
		unlink($target); unlink($query);
		}
       
	}	

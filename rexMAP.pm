package rexMAP;

use strict;
use Usu;
use Digest::MD5 qw(md5_hex);

sub prepref {
	
	# check and filter for toplevel
	# extract chr names and warn if non unique
	my $exp        = $_[0];
	my $database   = $exp->get_db;
	my $sqlhost    = $exp->get_sqlhost;
	my $sqluser    = $exp->get_sqluser;
	my $sqlport    = $exp->get_sqlport;
	my $out        = $exp->get_refs;
	my $srt        = $exp->annotab('seq_region');
	my ($i, $read, $written);
	my %names;
	
	my $dbh = DBI->connect("DBI:mysql:$database:$sqlhost:$sqlport","$sqluser","") or die "DB connection error: $DBI::errstr"; 
	my $sth = $dbh->prepare( qq( SELECT name FROM $srt ) );
	$sth->execute();
	$sth->bind_columns( \$i );
	while ( $sth->fetch() ) { $names{$i}++ }
	$dbh->disconnect();
	foreach my $n (keys %names) { if ($names{$n} > 1) { die "\nFATAL ERROR: Chromosome $n is not unique in annotation table $srt ($names{$n})\n" } }
	my $cleanout = $out . ".clean";
	my $dirtyout = $out . ".dirty";
	open(FIN,  "<$out") or die "cannot find $out";
	open(FOUT, ">$cleanout") or die;
	my $printswitch;
	while(<FIN>) {
		$read++;
		if (/^>(\S+)/) {
			if (defined $names{$1}) { $printswitch = 1; $names{$1}--; }
			else                    { $printswitch = 0; $names{$1} = -1 }
		}
		if ($printswitch) {
			print FOUT;
			$written++;
		}
	}
	close(FIN);
	close(FOUT);
	# writeback not found identifiers
	my $oops = 0;
	print STDERR "Annotation/Reference Sequence Crossmatch:\n";
	foreach my $n (sort keys %names) {
		if     ($names{$n} == 0) { print STDERR "\tOK PRESENT                 \t$n\n" }
		elsif ($names{$n} == -1) { print STDERR "\t!! NO_ANNOTATION (EXCLUDE?)\t$n\n"; $oops++ }
		else                     { print STDERR "\t?? NOT_IN_REFERENCE        \t$n\n"; $oops++ }
	}
	if ($oops == 0 || ($read == $written)) {
		system("rm -f $cleanout");
	}
	else {
		my $choice;
		do {
			print STDERR "\nDO YOU WANT TO EXCLUDE CHROMOSOMES WITHOUT ANNOTATION (Y)es (N)o (A)bort? ";
			$choice = <STDIN>;
		} while ($choice !~ /(y|n|a)/i);
		if ($choice =~ /y/i) {
			system("mv $out $dirtyout");
			system("gzip $dirtyout");
			system("mv $cleanout $out");
		}
		elsif ($choice =~ /n/i) {
			system("rm -f $cleanout");
			print STDERR "Completing Index...";
			my $db = DBI->connect("DBI:mysql:$database:$sqlhost:$sqlport","$sqluser","") or die "DB connection error: $DBI::errstr"; 
			foreach my $n (sort keys %names) {
				if ($names{$n} == -1) {
					my $sql = "INSERT IGNORE INTO $srt (name) VALUES (\"$n\")";
					$db->prepare( $sql )->execute();
					print STDERR "\.";
				}
			}
			$db->disconnect();
			print STDERR "DONE\n";			
		}
		else {
			die "ABORTED BY USER";
		}
	}
	# bowtie
	if (-s ($out . '.1.ebwt')) { warn "WARNING: BOWTIE index already present. Won't update!\n" }
	else { system("bowtie-build $out $out") }
	return $out;
}

sub genomemap {
	my $remap     = $_[1];
	my $expid   = $_[0]->get_exid;
	my $reads   = $_[0]->get_read;
	my $reflib  = $_[0]->get_refs;
	
	my $bwtk    = $_[0]->get_bwtk;
	my $bwtm    = $_[0]->get_bwtm;
	my $bwte    = $_[0]->get_bwte;
	my $bwts    = $_[0]->get_bwts;
	my $bwtl    = $_[0]->get_bwtl;
	my $bwtb    = $_[0]->get_bwtb;

	my $mapfile    = $_[0]->mapping;
	my $unmapfile  = $_[0]->unmapped_fa;
	my $repeatfile = $_[0]->repeat_fa;
	my $pnum       = $_[0]->get_threads;
	
	# map
	if (-s $mapfile && !($remap)) {
		warn "WARNING: Reads ($reads) already mapped ($mapfile). Won't update!\n";
	}
	else {
		my $exec = "bowtie -q -t -p $pnum ".($bwts ? "--solexa1.3-quals" : "--phred33-quals")." ".($bwtb ? "--best --strata " : "")." --chunkmbs 256 --un $unmapfile --max $repeatfile -l $bwtl -k $bwtk -m $bwtm -e $bwte  $reflib $reads";
		print STDERR " -> Running BOWTIE as follows: $exec\n";
		##############################################################
		my $processedlines;
		open(MAP, "$exec | cut -f1-4,7-8 | tee $mapfile |") or die;
		while(<MAP>) {
			if (/^(\S+)/) {
				++$processedlines;
				if ($processedlines % 10000 == 0) {
					print STDERR "\r" . $1 . "   ";
				}
			}
			else { die }
		}
		print STDERR "\rProcessed $processedlines mappings from bowtie                           \n";
		close(MAP);
	}
	return;
}

sub splicemap {
	my $exp       = $_[0];
	my $remap     = $_[1];
	my $bwtlib    = $exp->annot("splicigarseq");
	my $r         = $exp->get_read;
	
	my $bwtv      = $exp->get_bwtv;
	my $bwte      = $exp->get_bwte;
	my $bwtk      = $exp->get_bwtk;
	my $bwtm      = $exp->get_bwtm;
	my $bwtb      = $exp->get_bwtb;
	my $bwtl      = $exp->get_bwtl;
	my $bwts      = $exp->get_bwts;

	my $mapfile   = $exp->splicedbwtout;
	my $pnum      = $exp->get_threads;
	
	# map
	if (-s $mapfile && !($remap)) {
		warn "WARNING: Reads ($r) already mapped ($mapfile). Won't update!\n";
	}
	else {
		my $maximalRedundancy = $exp->getSpliceRedundancy;
		if ($bwtk <= $maximalRedundancy) { $bwtk += $maximalRedundancy }
		if ($bwtm <= $maximalRedundancy) { $bwtm += $maximalRedundancy }
		my $exec = "bowtie -q -t -p $pnum ".($bwts ? "--solexa1.3-quals" : "--phred33-quals")." ".($bwtb ? "--best --strata" : "")." --chunkmbs 256 -l $bwtl -e $bwte " . (($bwtv < 0) ? "" : "-v $bwtv") . " -k $bwtk -m $bwtm $bwtlib $r"; ### using max mismatch / no quals
		print STDERR " -> Running BOWTIE as follows: $exec\n";
		##############################################################
		my $processedlines;
		open(MAP, "$exec | cut -f1-4,7-8 | tee $mapfile |") or die;
		while(<MAP>) {
			if (/^(\S+)/) {
				++$processedlines;
				if ($processedlines % 10000 == 0) {
					print STDERR "\r" . $1 . "   ";
				}
			}
			else { die }
		}
		print STDERR "\rProcessed $processedlines mappings from bowtie                           \n";
		close(MAP);
	}
	return;
}

##########################################################################################################################################################################################
##### OLD ################################################################################################################################################################################
##########################################################################################################################################################################################

sub OLD_genomemap {
	my $remap     = $_[1];
	my $expid   = $_[0]->get_exid;
	my $reads   = $_[0]->get_read;
	my $reflib  = $_[0]->get_refs;

	my $bwtk    = $_[0]->get_bwtk;
	my $bwtm    = $_[0]->get_bwtm;
	my $bwte    = $_[0]->get_bwte;
	my $bwts    = $_[0]->get_bwts;
	my $bwtl    = $_[0]->get_bwtl;
	my $bwtb    = $_[0]->get_bwtb;

	my $mapfile    = $_[0]->mapping;
	my $unmapfile  = $_[0]->unmapped_fa;
	my $repeatfile = $_[0]->repeat_fa;
	my $pnum       = $_[0]->get_threads;

	# map
	if (-s $mapfile && !($remap)) {
		warn "WARNING: Reads ($reads) already mapped ($mapfile). Won't update!\n";
	}
	else {
		my $exec = "bowtie -q -t -p $pnum ".($bwts ? "--solexa1.3-quals" : "--phred33-quals")." ".($bwtb ? "--best --strata " : "")." --chunkmbs 256 --un $unmapfile --max $repeatfile -l $bwtl -k $bwtk -m $bwtm -e $bwte  $reflib $reads | cut -f1-4,7-8 > $mapfile";
		print STDERR " -> Running BOWTIE as follows: $exec\n";
		die if system($exec);
	}
	return;
}

sub OLD_splicemap {
	my $exp       = $_[0];
	my $remap     = $_[1];
	my $bwtlib    = $exp->annot("splicigarseq");
	my $r         = $exp->get_read;

	my $bwtv      = $exp->get_bwtv;
	my $bwte      = $exp->get_bwte;
	my $bwtk      = $exp->get_bwtk;
	my $bwtm      = $exp->get_bwtm;
	my $bwtb      = $exp->get_bwtb;
	my $bwtl      = $exp->get_bwtl;
	my $bwts      = $exp->get_bwts;

	my $mapfile   = $exp->splicedbwtout;
	my $pnum      = $exp->get_threads;

	# map
	if (-s $mapfile && !($remap)) {
		warn "WARNING: Reads ($r) already mapped ($mapfile). Won't update!\n";
	}
	else {
		my $maximalRedundancy = $exp->getSpliceRedundancy;
		if ($bwtk <= $maximalRedundancy) { $bwtk += $maximalRedundancy }
		if ($bwtm <= $maximalRedundancy) { $bwtm += $maximalRedundancy }
		my $exec = "bowtie -q -t -p $pnum ".($bwts ? "--solexa1.3-quals" : "--phred33-quals")." ".($bwtb ? "--best --strata" : "")." --chunkmbs 256 -l $bwtl -e $bwte " . (($bwtv < 0) ? "" : "-v $bwtv") . " -k $bwtk -m $bwtm $bwtlib $r | cut -f1-4,7-8 > $mapfile"; ### using max mismatch / no quals
		print STDERR " -> Running BOWTIE as follows: $exec\n";
		die if system($exec);
	}
	return;
}

1;

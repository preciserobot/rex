package rexARE;

use strict;
use rexExperiment;
use rexRemote;
use rexJUNCTION;
use Bio::SeqIO;
use Digest::MD5 qw(md5_hex);
use FileHandle;
#NG
sub build_reads_from_junctions_and_exons {
		my $exp     = $_[0];
		my $junc    = $_[1];
		my $readl       = $exp->get_readl;
		my $mockreaddb  = $exp->mrdb;
		my $sqldb    = $exp->get_db;
		my $sqlhost  = $exp->get_sqlhost;
		my $sqlport  = $exp->get_sqlport;
		my $sqluser  = $exp->get_sqluser;
		my $m1 = md5_hex(rand());
		my $m2 = md5_hex(rand());
		my $depth = 2;
		my $maxfork = $exp->get_threads;
		my @npids;
		my %fh; # the filehandles

		my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
		my $sth;
		print STDERR "\nPreparing Tables...";
		$dbh->prepare( qq{ DROP TABLE IF EXISTS $mockreaddb; }         )->execute();
		$dbh->prepare( qq{CREATE TABLE $mockreaddb (seq CHAR(255));} )->execute();
		print STDERR "DONE\n";
		$dbh->disconnect();

		# add sequences from exons 
		print STDERR "Building reads...";
		foreach my $db ($exp->annot("splicigarseq"), $exp->annot("exonseqs")) {
			my $seqnumber; # done sequence count
			print STDERR "\n            from $db";
			my $fin  = Bio::SeqIO->new(-file => $db ,-format => 'fasta');
			while ( my $seq = $fin->next_seq() ) {
				#FW
				for(my $i=0; $i < ($seq->length - $readl); $i++) {
					#FW
					my $id = substr($seq->seq, $i, $depth);
					fhprint(\%fh, $m1, $id, substr($seq->seq, $i, $readl));
				}
				#RV
				$seq = $seq->revcom;
				for(my $i=0; $i < ($seq->length - $readl); $i++) {
					#FW
					my $id = substr($seq->seq, $i, $depth);
					fhprint(\%fh, $m1, $id, substr($seq->seq, $i, $readl));
				}
				print STDERR "\r".++$seqnumber;
			}
		}

		# forked sort $maxfork
		my $current;
		my $merge = [];
		foreach my $f (sort keys %fh) {
			# childwaiter
			while (scalar @npids >= $maxfork) {
				print STDERR "\n\tWaiting for child process (" . scalar @npids . ")\n";
				foreach my $p (@npids) { print STDERR "\t$p" }
				my $finished = wait;
				for (my $i=0;$i < scalar @npids; $i++) {
					if ($npids[$i] == $finished) {
						splice @npids, $i, 1;
						print STDERR "\n\tChild PID $finished terminated";
						last;
					}
				}
			}
			my $partfile       = $f . $m1;
			my $sortedpartfile = "sorted_" . $f . $m1;
			push @{$merge}, $sortedpartfile;

			print STDERR "\nProcessing " . ++$current . " out of " . scalar(keys(%fh));
			my $kidpid;
			if (!defined($kidpid = fork())) { die "cannot fork: $!" }
			elsif ($kidpid == 0) {
				# child
				$fh{$f}->close;
				unless (-e $partfile) { die $partfile . " not found" }
				system("sort -u $partfile -o $sortedpartfile");
				system("rm -f $partfile");
				exit; # exit the child process
			}
			else {
				#parent
				push @npids, $kidpid;
				print STDERR "\t($kidpid)";
			}
		}
		# wait for all childs to finish
		while (@npids) {
			print STDERR "\n\tWaiting for child process (" . scalar @npids . ")\n";
			my $finished = wait;
			for (my $i=0;$i < scalar @npids; $i++) {
				if ($npids[$i] == $finished) {
					splice @npids, $i, 1;
					print STDERR "\n\tChild PID $finished terminated";
					last;
				}
			}
		}
		my $joint = join ' ', @{$merge};
		unless (system("cat $joint > $m2")) { system("rm -f $joint") };

		# KEEP STATIC COPY
		my $rawmock = $exp->rawmockreads;
		system("ln $m2 $rawmock");

		$exp->wait4database;
		$dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
		print STDERR "\nLoading...";
		if ($dbh->prepare( qq{LOAD DATA LOCAL INFILE \"$m2\" INTO TABLE $mockreaddb } )->execute()) { system("rm -f $m2") }
		else { die "Mockreaddb store failed" }
		print STDERR "DONE\n";
		$dbh->disconnect();
		return $exp;
}

# calc splice redundancy prototype
sub spliceSeqRedu {
	my $exp     = $_[0];
	my $junc    = $_[1];
	my $readl   = $exp->get_readl;
	my $ovrh    = $exp->get_ovrh;
	my $m1 = md5_hex(rand());
	my @npids;
	my %fh; # the filehandles

	# create the filehandles by generating the splice clusters
	my $splicigars = $exp->annot("splicigars");
	my $sqldb    = $exp->get_db;
	my $sqlhost  = $exp->get_sqlhost;
	my $sqlport  = $exp->get_sqlport;
	my $sqluser  = $exp->get_sqluser;
	my ( $sql_id, $sql_region, $sql_start, $sql_end );
	my (%lulu, %cluster);
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my $sth = $dbh->prepare( qq{ SELECT id, seq_region_id, MIN(start) - ($readl - $ovrh) + 1 , MAX(end) + ($readl - $ovrh) - 1 FROM $splicigars GROUP BY id } );
	$sth->execute();
	$sth->bind_columns( \$sql_id, \$sql_region, \$sql_start, \$sql_end);
	while( $sth->fetch() ) { push @{$lulu{$sql_region}{$sql_start}}, [$sql_end, $sql_id] }
	$dbh->disconnect();
	
	# build clusters
	my $clusternumber = 0;
	my $end; # strictly increasing / reset on new region
	foreach my $r (keys %lulu) {
		$end = 0; 
		foreach my $start (sort {$a <=> $b} keys %{$lulu{$r}}) {
			foreach my $s (@{$lulu{$r}{$start}}) {
				if ($start > $end) { ++$clusternumber } # nonoverlapping ->start new cluster
				$end = ($s->[0] > $end ) ? $s->[0] : $end ;
				push @{$cluster{$clusternumber}}, $s->[1];
			}
		}
	}
	
	# get all multisplices
	my %ms;
	$dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	$sth = $dbh->prepare( qq{ SELECT id, seq_region_id, start, end FROM $splicigars ORDER BY start } );
	$sth->execute();
	$sth->bind_columns( \$sql_id, \$sql_region, \$sql_start, \$sql_end);
	while( $sth->fetch() ) { push @{$ms{$sql_id}}, [$sql_start, $sql_end] }
	$dbh->disconnect();
	
	# find largest
	my $maxovp = 10;
	my $ccount = scalar(keys %cluster);
	foreach my $c (sort { scalar(@{$cluster{$b}}) <=> scalar(@{$cluster{$a}}) } keys %cluster) { # starting with biggest cluster should optimize speed
		next if (scalar(@{$cluster{$c}}) < $maxovp); # skip clusters that cannot have more overlaps because they are themselves too small...
		print STDERR "\r" . (" " x 120) . "\r" . --$ccount . "\t" . scalar(@{$cluster{$c}}) . " => $maxovp" ;
		# get hash with counts
		my $maps = _getClusterMappings(\%ms, \@{$cluster{$c}}, $readl-$ovrh);
		foreach my $cnt (sort { $maps->{$b} <=> $maps->{$a} } keys %{$maps}) {
			if ($maxovp < $maps->{$cnt}) {
				$maxovp = $maps->{$cnt};	
			}
			last;
		}
	}
	#update annotation table
	$exp->updateSpliceRedundancy($maxovp);	
	return $maxovp;
}

#NG#
sub WriteMockReads {
	my $exp  = $_[0];
	my $mockreaddb  = $exp->mrdb;
	my $MockReadSampleSize = 300000000; # static to be consistent in all genomes (max value to avoid memory overflow)
	my $sqldb    = $exp->get_db;
	my $sqlhost  = $exp->get_sqlhost;
	my $sqlport  = $exp->get_sqlport;
	my $sqluser  = $exp->get_sqluser;

	# open fresh connection as after indexing the server might have gone away due to timeouts
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	print STDERR "Creating FASTQ file...";
	# count reads
	my $sth = $dbh->prepare( qq{ SELECT COUNT(*) FROM $mockreaddb; } );
	$sth->execute();
	$sth->bind_columns( \$a );
	my $mrc = 0;
	while( $sth->fetch() ) { $mrc = $a }
	$dbh->disconnect();
	print STDERR "($mrc)";
	
	# CREATE FASTQ (write directly to file)
	open(MRDB, ">$mockreaddb") or die;

	my $localexecute = "mysql --quick --unbuffered -u $sqluser -h $sqlhost -P $sqlport $sqldb";
	my $sqlquery = qq{ SELECT seq FROM $mockreaddb; };
	my $sqlexec = "echo \"$sqlquery\" | $localexecute";
	my $num;
	my $fails;
	if (-s $exp->rawmockreads) {
		my $rawmock = $exp->rawmockreads;
		open(SQLQUERY, "$rawmock" ) or die;
		while(<SQLQUERY>) {
			if (/^((A|G|C|T|N)+)/) {
				if (++$num % 100000 == 0) { print STDERR "\rWriting FASTQ file for $mrc (max $MockReadSampleSize) sampled reads... $num" }
				if (rand(1) < ($MockReadSampleSize / $mrc)) {
					print MRDB "@" . $num . "\n";
					print MRDB $1 . "\n";
					print MRDB "+" . "\n";
					print MRDB substr($exp->get_quba, 0, length($1)) . "\n";
				}
			} else { if (++$fails > 1) { print $_ ; die "line not recognized" } }
		}
		close(SQLQUERY);
		unless ($fails) { system("rm -f $rawmock") } # remove rawmockreads once mockread.fastq is written
	}
	else {
		open(SQLQUERY, "$sqlexec |" ) or die;
		while(<SQLQUERY>) {
			if (/^((A|G|C|T|N)+)/) {
				if (++$num % 100000 == 0) { print STDERR "\rWriting FASTQ file for $mrc (max $MockReadSampleSize) sampled reads... $num" }
				if (rand(1) < ($MockReadSampleSize / $mrc)) {
					print MRDB "@" . $num . "\n";
					print MRDB $1 . "\n";
					print MRDB "+" . "\n";
					print MRDB substr($exp->get_quba, 0, length($1)) . "\n";
				}
			} else { if (++$fails > 1) { print $_ ; die "line not recognized" } }
		}
		close(SQLQUERY);
	}
	close(MRDB);
	print STDERR "\rWriting FASTQ file... ($num out of $mrc reads) DONE\n";
	return $exp;
}

#NG#
sub MockMap {
	my $exp       = $_[0];
	my $expid     = $exp->get_exid;
	# BWT param
	my $bwte      = $exp->get_bwte;
	my $maxmism   = $exp->get_bwtv;
	# mapping libs
	my $genomelib = $exp->get_refs;
	my $splicelib = $exp->get_splc ? $exp->annot("splicigarseq") : $exp->splicedbwt;
	# mock reads
	my $mockreads = $exp->mrdb;
	# outfiles & SQL tables
	my $genomeout = $exp->mockgenomebwtout;
	my $spliceout = $exp->mocksplicebwtout;
	my $srt       = $exp->annot("seq_region");
	my $introntab = $exp->annot("introns");
	# SQL
	my $sqldb     = $exp->get_db;
	my $sqlhost   = $exp->get_sqlhost;
	my $sqlport   = $exp->get_sqlport;
	my $sqluser   = $exp->get_sqluser;
	# eachside
	my $pnum      = $exp->get_threads;
	my $bwts      = $exp->get_bwts;
#get chromo list
	print STDERR "Creating Chromosome index...\n";
	my ($a, $b);
	my %region;
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my $sth = $dbh->prepare( qq{ SELECT seq_region_id, name FROM $srt ; } );
	$sth->execute();
	$sth->bind_columns( \$a, \$b );
	while( $sth->fetch() ) { $region{$b} = $a }
	$dbh->disconnect();
	
	my @rotor    = ("/ ","- ","\\ ","| ");
	my $rotorpos = 0;
	
	my $processedlines = 9999;
	my $bwtexec = "bowtie -q -t ". ($bwts ? "--solexa1.3-quals" : "--phred33-quals") ." --chunkmbs 256 -p $pnum -k 1 -m 1 -e $bwte $genomelib $mockreads";
	print STDERR "Mapping $mockreads to $genomelib ($bwtexec)...\n";
	open(MAPOUT, ">$genomeout") or die;
	open(MAP, "$bwtexec | cut -f1-4 |") or die;
	print STDERR "\n";
	while(<MAP>) {
		if (/^(\d+)\s+([\+|-])\s+(\S+)\s+(\d+)/) {
			++$processedlines;
			if ($processedlines % 10000 == 0) {
				print STDERR "\r" . $1 . "   ";
			}
			if ($processedlines % 1000 == 0) {
				print STDERR "\b" x length($rotor[$rotorpos]);
				$rotorpos = ++$rotorpos % scalar(@rotor);
				print STDERR $rotor[$rotorpos];
			}
			print MAPOUT $1 . "\t";
			print MAPOUT (($2 eq '-') ? -1 : (($2 eq '+') ? 1 : die)) . "\t";
			print MAPOUT (($region{$3}) ? $region{$3} : die "Unknown Region") . "\t";
			print MAPOUT ($4 + 1) . "\n"; # corrects for offset
		}
		else { die }
	}
	print STDERR "\n";
	close(MAP);
	close(MAPOUT);
	# calculate maximal redundancy
	my $maximalRedundancy = $exp->getSpliceRedundancy;
	$processedlines = 9999;
	# $bwtexec = "bowtie -q -t ". ($bwts ? "--solexa1.3-quals" : "--phred33-quals") ." -p $pnum -k $maximalRedundancy -m $maximalRedundancy -e $bwte -v $maxmism $splicelib $mockreads";
	$bwtexec = "bowtie -q -t ". ($bwts ? "--solexa1.3-quals" : "--phred33-quals") ." --chunkmbs 256 -p $pnum -k $maximalRedundancy -m $maximalRedundancy -e $bwte " . (($maxmism < 0) ? "" : "-v $maxmism") . " $splicelib $mockreads";
	print STDERR "Mapping $mockreads to $splicelib ($bwtexec)...\n"; # more permissive as splices can be redundant
	open(MAPOUT, ">$spliceout") or die;
	open(MAP, "$bwtexec | cut -f1-4 |") or die;
	while(<MAP>) {
		if (/^(\d+)\s+([\+|-])\s+(\S+)\s+(\d+)/) {
			++$processedlines;
			if ($processedlines % 10000 == 0) {
				print STDERR "\r" . $1 . "   ";
			}
			if ($processedlines % 1000 == 0) {
				print STDERR "\b" x length($rotor[$rotorpos]);
				$rotorpos = ++$rotorpos % scalar(@rotor);
				print STDERR $rotor[$rotorpos];
			}
			print MAPOUT $1 . "\t";
			print MAPOUT (($2 eq '-') ? -1 : (($2 eq '+') ? 1 : die)) . "\t";
			print MAPOUT $3 . "\t";
			print MAPOUT ($4 + 1) . "\n"; # corrects for offset
		}
		else { die }
	}
	print STDERR "\n";
	close(MAP);
	close(MAPOUT);
	return;
}

#NG## now runs remotely automatically
sub AmEsMs { 
	my $exp         = $_[0];
	my $rfexon      = $exp->rfexon;
	my $genomeout   = $exp->mockgenomebwtout($_[1]); # argument is identifier override
	my $spliceout   = $exp->mocksplicebwtout($_[1]); # argument is identifier override
	my $splicigars  = $exp->annot("splicigars");
	my $readlength  = $exp->get_readl;
	my $rexon       = 'amesout_' . md5_hex(rand());
	my $rcexon      = $exp->rcexon;
	my $rbexon      = $exp->rbexon; # brute without correction for multiople anchors
	my $sqldb    = $exp->get_db;
	my $sqlhost  = $exp->get_sqlhost;
	my $sqlport  = $exp->get_sqlport;
	my $sqluser  = $exp->get_sqluser;
	
	
	unless (-s $splicigars) {
		my $sqlhost = $exp->get_sqlhost;
		my $sqluser = $exp->get_sqluser;
		my $sqlport = $exp->get_sqlport;
		my $sqldb   = $exp->get_db;
		my $localextract = "mysql --skip-column-names -u $sqluser -h $sqlhost $sqldb";
		my $cigarextract = qq{ SELECT * FROM $splicigars };
		die "ERROR: splicigar extract failed from table $splicigars" if system("echo \"$cigarextract\" | $localextract > $splicigars");
	}
	
	my $remote  = rexRemote->new();
	my $amesbin = $remote->get_ames_bin;
	my $amesrun = "$amesbin $rfexon $genomeout $spliceout $splicigars $readlength $rexon";
	#die "\nDEBUG AMES-MS USING:\n\t$amesrun\n(then remove this line)\n";
	if ((-s $rfexon) && (-s $genomeout) && (-s $spliceout) && (-s $splicigars)) {
		# run ames
		print STDERR "$amesrun\n";
		# send
		$remote->sendFile($rfexon);
		$remote->sendFile($genomeout);
		$remote->sendFile($spliceout);
		$remote->sendFile($splicigars);
		# compile
		$remote->sendCompile("ames");
		# run 
		$remote->execute($amesrun);
		# fetch
		$remote->fetchFile($rexon);
		# cleanup
		$remote->cleanup($rfexon, $genomeout, $spliceout, $splicigars, $rexon);
		$remote->cleanup();
		
		# split file
		system("cut -f1-7,9 $rexon > $rbexon");
		system("cut -f1-7,8 $rexon > $rcexon");
		#system("rm -f $rexon");
		
		# load new exon file
		my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
		my $sqlstack = [qq(DROP TABLE IF EXISTS $rbexon;),
						qq(CREATE TABLE $rbexon (gene_id INT, exon_id INT, seq_region INT, seq_region_start INT, seq_region_end INT, seq_region_strand INT, tag INT, uniscore INT);),
						qq(LOAD DATA LOCAL INFILE \"$rbexon\" INTO TABLE $rbexon),
						qq(DROP TABLE IF EXISTS $rcexon;),
						qq(CREATE TABLE $rcexon (gene_id INT, exon_id INT, seq_region INT, seq_region_start INT, seq_region_end INT, seq_region_strand INT, tag INT, uniscore INT);),
						qq(LOAD DATA LOCAL INFILE \"$rcexon\" INTO TABLE $rcexon),
				];
		foreach my $sql (@{$sqlstack}) { print STDERR "\."; $dbh->prepare( $sql )->execute() }
		# compress mock stuff
		print STDERR "\.";
	}
	else {
		warn "\nERROR: cannot find required files for AmEs call ($amesrun)\n";
		print STDERR "$rfexon\t" . ((-s $rfexon) ? "OK" : "!!") . "\n";
		print STDERR "$genomeout\t" . ((-s $genomeout) ? "OK" : "!!") . "\n";
		print STDERR "$spliceout\t" . ((-s $spliceout) ? "OK" : "!!") . "\n";
		print STDERR "$splicigars\t" . ((-s $splicigars) ? "OK" : "!!") . "\n";
		die;
	}
	print STDERR "\n";
	return;
}

#NG#
sub PseudoAmEs {
	my $exp = $_[0];
	my $rfexon = $exp->rfexon;
	my $rdexon = $exp->rdexon;
	my $rzexon = $exp->rzexon;
	my ($li,$er);
	open(FIN,  "<$rfexon") or die;
	open(DOUT, ">$rdexon") or die;
	open(ZOUT, ">$rzexon") or die;
	while (<FIN>) {
		++$li;
		if (/(\S+\s+\S+\s+\S+\s+)(\S+)\s+(\S+)(\s+\S+\s+\S+)/) {
			print DOUT $1 . $2 . "\t" .$3 . $4 . "\t" . ($3 - $2 + 1) . "\n";
			print ZOUT $1 . $2 . "\t" .$3 . $4 . "\t" . 0 . "\n";
		}
		else { warn "RFexon file ($rfexon) is corrupt at line $li\n"; die "TOO MANY ERRORS" if (++$er > 10) }
	}
	close(FIN);
	close(DOUT);
	close(ZOUT);
	
	# store to sql
	my $sqldb    = $exp->get_db;
	my $sqlhost  = $exp->get_sqlhost;
	my $sqlport  = $exp->get_sqlport;
	my $sqluser  = $exp->get_sqluser;
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	
	foreach my $l ($rdexon, $rzexon) {
		$exp->wait4database;
		my $sqlstack = [ qq(DROP TABLE IF EXISTS $l;),
						qq(CREATE TABLE $l (gene_id INT, exon_id INT, seq_region INT, seq_region_start INT, seq_region_end INT, seq_region_strand INT, tag INT, uniscore INT);),
						qq(LOAD DATA LOCAL INFILE \"$l\" INTO TABLE $l)
		];
		foreach my $sql (@{$sqlstack}) { print STDERR "\."; $dbh->prepare( $sql )->execute() }
		print STDERR "\.";
	}
	print STDERR "\n";
	$dbh->disconnect();
	return;
}

##############
### HELPER ###
##############
#NG#
sub seqdbfetch {
	my $fetch = "xdget -n -a $_[2]->[0] -b $_[2]->[0] $_[0] $_[1] |";
	my ($name,$seq);
	open(SEQ, $fetch) or die;
	while(<SEQ>) {
		my $line = $_; chomp($line);
		if ($line =~ /^>/) { $name .= $line }
		else               { $seq .= $line }
	}
	close(SEQ);
	die unless ($seq && $name);
	my $seqobj = Bio::Seq->new(-name => $name, -seq => $seq);
	return $seqobj;
}

#NG#
sub fhprint {
	my $f = $_[0];
	my $b = $_[1];
	my $i = $_[2];
	my $s = $_[3];
	unless ($f->{$i}) {
		# open new filehandle
		$f->{$i} = FileHandle->new;
		$f->{$i}->open(">$i" . $b);
	}
	$f->{$i}->print("$s\n");
	return $f;
}

#NG#
sub _revcom {
	my $dna = shift;
	my $revcomp = reverse($dna);
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	return $revcomp;
}

#NG#
sub _chooseformat {
	my $limit = 5000;
	# check seq number
	my $chrcount = 0;
	open(CHK, "<$_[0]") or die;
	while (<CHK>) {
		if ($chrcount >= $limit) { last }
		if (/^>/) {print STDERR "\rNOTE: Counting reference sequences... ".$chrcount++ }
	}
	close(CHK);
	my $format = ($chrcount < $limit) ? 'largefasta' : 'fasta';
	print STDERR "\rNOTE: Using $format for sequence extraction\n";
	return $format;
}
#reducheck
sub _getClusterMappings {
	my %mapCount;
	my $msp       = $_[0];
	my $mID       = $_[1];
	my $length    = $_[2];
	# foreach cluster member
	for(my $i = 0; $i < scalar(@{$mID}); $i++) {
		my $ms = $msp->{$mID->[$i]};
		for(my $o=0; $o < $length; $o++) {
			$mapCount{__anchorChecksum([$o, $o + $length],$ms)}++;
		}
	}
	return \%mapCount;
}
#reducheck
sub __anchorChecksum {
	#my @anchors;
	my $o = $_[0]; # array ref to start and end offset
	my $m = $_[1]; # array ref of multisplice slices
	my $offsetsum = 0;
	my $cigar;
	for (my $a = 0; $a < scalar(@{$m}); $a++) {
		if (($o->[0] - $offsetsum) + $m->[$a]->[0] <= $m->[$a]->[1]) {
			my $i = 0;
			# set start
			$cigar .= ($o->[0] - $offsetsum) + $m->[$a]->[0];
			# find next anchors
			for (my $b = $a; $b < scalar(@{$m}); $b++) { # take up where $a left
				# start
				if ($i > 0) {
					$cigar .= $m->[$b]->[0];
				}
				# end
				if (($o->[1] - $offsetsum) + $m->[$b]->[0] <= $m->[$b]->[1]) {
					# add last and return
					$cigar .= ($o->[1] - $offsetsum) + $m->[$b]->[0];
					return md5_hex($cigar);
				} else {
					$cigar .= $m->[$b]->[1];
				}
				# update for next
				$offsetsum += $m->[$b]->[1] - $m->[$b]->[0] + 1;
				$i++;
			}
			# no anchors found for given offset->[1]
			die "out of range ($o->[0]|$o->[1])<";
		}
		$offsetsum += $m->[$a]->[1] - $m->[$a]->[0] + 1;
	}
	die "out of range >($o->[0]|$o->[1])";
}



##############
### LEGACY ###
##############

1;

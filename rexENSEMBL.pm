package rexENSEMBL;

use strict;
use Digest::MD5 qw(md5_hex);

# use ensannot for all ensembl related stuff

#NG#
sub selectCoordinateSystem {
	my $exp           = $_[0];
	my $local         = $exp->get_db;
	my $seqregion     = $exp->ensannot('seq_region');
	my $seqregionfull = $exp->ensannot('seq_region_full');
	my $coordsystem   = $exp->ensannot('coord_system');
	
	my $sqlhost = $exp->get_sqlhost;
	my $sqluser = $exp->get_sqluser;
	my $sqlport = $exp->get_sqlport;
	my $localexecute = "mysql -u $sqluser -h $sqlhost --skip-column-names $local";
	
	# refetch full table
	my $remote = $exp->get_edb;
	my $host   = 'ensembldb.ensembl.org';
	my $port   = 5306;
	my $user   = 'anonymous';
	my $tables = [ 'seq_region', 'coord_system' ];
	my %index  = ( seq_region => ['name'] );
	# pure mysqldump scheme
	print STDERR "\n\nFetching Ensembl gene data ($remote)\n";
	
	foreach my $t (@{$tables}) {
		my $structfile = $t . time() . '.sql';
		my $datafile   = $t . time() . '.data';
		# copy table struct
		print STDERR "\tDumping table structure\n";
		system("mysqldump -h $host -P $port -u $user -d $remote $t --single_transaction > $structfile");
		# fetch data
		my $fquery = "select * from $t;";
		print STDERR "\t$fquery ($host)\n";
		system("echo \"$fquery\" | mysql -h $host -u $user -P $port --skip-column-names $remote > $datafile");
		# drop table
		print STDERR "\tDrop table $t\n";
		system("echo \"DROP TABLE IF EXISTS $t\" | $localexecute");
		# create table
		print STDERR "\tCreating table $t\n";
		system("$localexecute < $structfile");
		# drop table
		my $t2 = $_[0]->ensannot($t);
		print STDERR "\tDrop table $t2\n";
		system("echo \"DROP TABLE IF EXISTS $t2\" | $localexecute");
		#rename table
		my $rquery = "RENAME TABLE $t TO $t2;";
		print STDERR "\t$rquery\n";
		system("echo \"$rquery\" | $localexecute");
		# insert into table
		my $iquery = 'LOAD DATA LOCAL INFILE \''. $datafile . '\' INTO TABLE '. $t2 . ";";
		print STDERR "\t$iquery\n";
		system("echo \"$iquery\" | $localexecute");
		# Additional Indexing
		if ($index{$t}) {
			print STDERR "\tIndexing table $t\n";
			foreach my $i (@{$index{$t}}) {
				my $iquery = "ALTER TABLE $t2 ADD INDEX idx_$i ($i);";
				print STDERR "\t$iquery\n";
				system("echo \"$iquery\" | $localexecute");
			}
			print STDERR "\n";
		}
		# remove tmpfile
		system("rm -f $datafile $structfile");
	}
	
	# Select coordinate system
	my @choices;
	my $sel = qq(SELECT DISTINCT version from $coordsystem);
	open(SEQ, "echo \"$sel\" | $localexecute |") or die;
	while (<SEQ>) {
		my $line = $_;
		chomp($line);
		if ($line !~ /(version|NULL)/i) { push @choices, $line }
	}
	close(SEQ);
	for (my $i = 0; $i < scalar(@choices); $i++) {
		print STDERR $i . ". " . $choices[$i] . "\n";
	}
	print STDERR "Please select version: ";
	my $c = <STDIN>;
	chomp($c);
	my $clone = qq(DROP TABLE IF EXISTS $seqregionfull; CREATE TABLE $seqregionfull LIKE $seqregion; INSERT INTO $seqregionfull SELECT * FROM $seqregion);
	my $clean = "DELETE FROM sr USING $seqregion as sr, $coordsystem as cs WHERE sr.coord_system_id = cs.coord_system_id AND cs.version != \'$choices[$c]\';";
	print STDERR "\n\tClone...";
	system("echo \"$clone\" | $localexecute");
	print STDERR "DONE\n\tClean...";
	system("echo \"$clean\" | $localexecute");
	print STDERR "DONE\n";
	return 1;
}

#NG#
sub fetchTables {
	my $exp = $_[0];
	my $remote = $exp->get_edb;
	my $local  = $exp->get_db;
	my $host   = 'ensembldb.ensembl.org';
	my $port   = 5306;
	
	my $sqlhost = $exp->get_sqlhost;
	my $sqluser = $exp->get_sqluser;

	my $localexecute = "mysql -u $sqluser -h $sqlhost --skip-column-names $local";
	
	my $user   = 'anonymous';
	my @tables = ('exon','exon_transcript','gene','gene_stable_id','transcript','transcript_stable_id');#, 'seq_region', 'coord_system');
	my %index  = (	exon       => ['seq_region_start', 'seq_region_end', 'seq_region_strand'],
					transcript => ['seq_region_start', 'seq_region_end', 'seq_region_strand'],
					gene       => ['seq_region_start', 'seq_region_end', 'seq_region_strand']); #, seq_region => ['name'] );


	# pure mysqldump scheme
	print STDERR "\n\nFetching Ensembl gene data ($remote)\n";
	foreach my $t (@tables) {
		my $structfile = $t . time() . '.sql';
		my $datafile   = $t . time() . '.data';
		# copy table struct
		print STDERR "\tDumping table structure\n";
		system("mysqldump -h $host -P $port -u $user -d $remote $t --single_transaction > $structfile");
		# fetch data
		my $fquery = "select * from $t;";
		print STDERR "\t$fquery ($host)\n";
		system("echo \"$fquery\" | mysql -h $host -u $user -P $port --skip-column-names $remote > $datafile");
		# drop table
		print STDERR "\tDrop table $t\n";
		system("echo \"DROP TABLE IF EXISTS $t\" | $localexecute");
		# create table
		print STDERR "\tCreating table $t\n";
		system("$localexecute < $structfile");
		# drop table
		my $t2 = $_[0]->ensannot($t);
		print STDERR "\tDrop table $t2\n";
		system("echo \"DROP TABLE IF EXISTS $t2\" | $localexecute");
		#rename table
		my $rquery = "RENAME TABLE $t TO $t2;";
		print STDERR "\t$rquery\n";
		system("echo \"$rquery\" | $localexecute");
		# insert into table
		my $iquery = 'LOAD DATA LOCAL INFILE \''. $datafile . '\' INTO TABLE '. $t2 . ";";
		print STDERR "\t$iquery\n";
		system("echo \"$iquery\" | $localexecute");
		# Additional Indexing
		if ($index{$t}) {
			print STDERR "\tIndexing table $t\n";
			foreach my $i (@{$index{$t}}) {
				my $iquery = "ALTER TABLE $t2 ADD INDEX idx_$i ($i);";
				print STDERR "\t$iquery\n";
				system("echo \"$iquery\" | $localexecute");
			}
			print STDERR "\n";
		}
		# remove tmpfile
		system("rm -f $datafile $structfile");
		#print STDERR "\n###############################################################################################################################";
		#print STDERR "\n######## PLEASE RESTRICT SEQ REGIONS TO COORD SYSTEM RANK 1 TO EXCLUDE OBSOLETE/LOWLEVEL COORDINATE SYSTEMS TO BE USED ########";
		#print STDERR "\n###############################################################################################################################\n";
	}
	return;
}

#NG#
sub buildJunctionFile {
	my $exp  = $_[0]; 
	my $sr = $exp->ensannot("seq_region");
	my $er = $exp->ensannot("exon_ranks");
	my $ex = $exp->ensannot("exon");
	my $et = $exp->ensannot("exon_transcript");
	my $jeff = $exp->ensembljunctions;
	my $sqlhost = $exp->get_sqlhost;
	my $sqluser = $exp->get_sqluser;
	my $sqlport = $exp->get_sqlport;
	my $sqldb   = $exp->get_db;
	my $localexecute = "mysql -u $sqluser -h $sqlhost $sqldb";
	my @s;

	# create seq_region_index
	my %seqregion;
	my $sel = qq(SELECT DISTINCT seq_region_id, name from $sr);
	open(SEQ, "echo \"$sel\" | $localexecute |") or die;
	while (<SEQ>) {
		if (/^(\d+)\s+(\S+)/) { $seqregion{$1} = $2 }
	}
	close(SEQ);
	
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "DB connection error: $DBI::errstr"; 
	my $jctid = 0;
	# drop/create exon rank table
	$s[0] = qq( DROP TABLE IF EXISTS $er; );
	$s[1] = qq( CREATE TABLE $er SELECT e.seq_region_id, e.seq_region_start, e.seq_region_end, e.seq_region_strand, et.rank, et.transcript_id FROM $et AS et, $ex AS e WHERE e.exon_id = et.exon_id; );
	#index the exon rank table
	$s[2] = qq( alter table $er add index idx_tid (transcript_id), add index idx_rank (rank), add index idx_strand (seq_region_strand), add unique index idx_unique (transcript_id,rank,seq_region_id,seq_region_start, seq_region_end); );
	my $t = qq( select distinct l.seq_region_id, l.seq_region_strand, l.seq_region_end+1, r.seq_region_start-1 from $er as l left join $er as r on l.transcript_id = r.transcript_id where (l.seq_region_strand = 1 and r.seq_region_strand = 1 and l.rank+1 = r.rank) OR (l.seq_region_strand = -1 and r.seq_region_strand = -1 and l.rank = r.rank+1););
	print STDERR "Extracting Ensembl Junctions.";
	for(my $i=0; $i < scalar(@s); $i++) {
		$dbh->prepare($s[$i])->execute();
		print STDERR "$i";
	}
	open(EJCT, ">$jeff") or die;
	my $sth = $dbh->prepare($t);
	$sth->execute();
	my ( $srid, $strand, $start, $end );
	$sth->bind_columns( \$srid, \$strand, \$start, \$end );
	my %notfound;
	while( $sth->fetch() ) {
		if (defined $seqregion{$srid}) {
			print EJCT ++$jctid . "\t";
			print EJCT $seqregion{$srid} . "\t";
			print EJCT $strand . "\t";
			print EJCT $start . "\t";
			print EJCT $end . "\t";
			print EJCT "Annotated\n";
		}
		else {
			$notfound{$srid}++; 
		}
	}
	foreach my $sss (sort {$a <=> $b} keys %notfound) {
		warn "WARNING: Did not find chromsome with seq_region_id $sss\n";
	}
	close(EJCT);	
	print STDERR "Storing";
	my $sqlcreate = [];
	push @{$sqlcreate}, qq{DROP TABLE IF EXISTS $jeff;};
	push @{$sqlcreate}, qq{CREATE TABLE $jeff (junction_id INT, region CHAR(32), strand INT, start INT, end INT, source CHAR(128));};
	push @{$sqlcreate}, qq{LOAD DATA LOCAL INFILE \"$jeff\" INTO TABLE $jeff };
	foreach my $query (@{$sqlcreate}) { print STDERR "\."; $dbh->prepare( $query )->execute() }
	print STDERR "DONE\n";
	$dbh->disconnect();
	print STDERR ".DONE\n";
	return;
}

#NG#
sub buildCoreExons {
	# Build core exons (present in all transcripts)
	# calculate transcript lengths
	
	my $tx       = $_[0]->ensannot("transcript");
	my $ex       = $_[0]->ensannot("exon");
	my $extx     = $_[0]->ensannot("exon_transcript");
	my $sr       = $_[0]->ensannot("seq_region");
	my $gs       = $_[0]->ensannot("gene_stable_id");
	
	my $efile    = $_[0]->ensemblexons;
	my $tmp_fh = *TMPFILE;
	
	my $sqldb    = $_[0]->get_db;
	my $sqlhost  = $_[0]->get_sqlhost;
	my $sqlport  = $_[0]->get_sqlport;
	my $sqluser  = $_[0]->get_sqluser;
	

	print STDERR "Defining core exons\n";
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my $annot_sql = qq(select distinct gs.stable_id, t.transcript_id, s.name, e.seq_region_start, e.seq_region_end, e.seq_region_strand from $ex as e, $extx as et, $tx as t, $sr as s, $gs as gs where e.exon_id = et.exon_id and et.transcript_id = t.transcript_id and e.seq_region_id = s.seq_region_id and t.gene_id = gs.gene_id and t.biotype = "protein_coding" order by gs.stable_id);
	my $sth = $dbh->prepare( $annot_sql );
	$sth->execute();
	my $gex = [];
	my ($try, $success) = (0,0);
	open($tmp_fh, ">$efile") or die;
	while( my $row = $sth->fetchrow_hashref ) {
		print STDERR "\r" . ++$try;
		if (!($gex->[0]) || $gex->[0]->{stable_id} eq $row->{stable_id} ) { push @{$gex}, $row }
		else {
			print STDERR "\t" . ++$success;
			_build_core($tmp_fh, $gex);
			$gex = [];
			push @{$gex}, $row;
		}
	}
	print STDERR "\r" . $try . "\t" . ++$success;
	_build_core($tmp_fh, $gex);
	close(TMPFILE);
	
	# load annotations into global annotation table (retro + core exons)
	print STDERR "\nStoring";
	my $sqlcreate = [];
	push @{$sqlcreate}, qq{DROP TABLE IF EXISTS $efile;};
	push @{$sqlcreate}, qq{CREATE TABLE $efile (gene_id INT, exon_id INT, seq_region INT, seq_region_start INT, seq_region_end INT, seq_region_strand INT, tag INT);};
	push @{$sqlcreate}, qq{LOAD DATA LOCAL INFILE \"$efile\" INTO TABLE $efile};
	foreach my $query (@{$sqlcreate}) { print STDERR "\."; $dbh->prepare( $query )->execute() }
	print STDERR "DONE\n";
	$dbh->disconnect();
	print STDERR "\n";
	return;
}

#NEW VERSION USING R?EXON
sub geco {
	my $exp      = $_[0];
	my $strand = ($_[1]) ? (($_[1] =~ /\+/) ? "fwd" : ( ($_[1] =~ /-/) ? "rev" : $_[1] )) : undef ;
	my $exon     = $exp->exonmodel;
	my $coverage = $exp->raparm_coverage($strand);
	my $geco_out = $exp->geco($strand);
	my $sqldb    = $exp->get_db;
	my $sqlhost  = $exp->get_sqlhost;
	my $sqlport  = $exp->get_sqlport;
	my $sqluser  = $exp->get_sqluser;
	
	if ($exp->get_strands) {
		warn "\nWARNING: GeCo will not run on strand separated coverage tracks -> SKIPPING\n";
		return;
	}
	
	unless (-s $exon) {
		die "\nERROR: Cannot filed exon file ($exon)\n";
	}
	unless (-s $coverage) {
		print STDERR " -> Fetching coverage\n";
		my $sql = qq(SELECT * FROM $coverage);
		system("mysql -u $sqluser -h $sqlhost $sqldb --execute=\"$sql\" > $coverage");
	}
	# execute geco 
	my $exec = "geco " . $exon ." ". $coverage ." ". $geco_out;
	print STDERR " -> Running GeCo: $exec\n";
	my $return = system("$exec");
	if ($return) { die "FATAL: GECO returned an error code ($return)\n" }
	
	# load data into sql
	my $sqlstack = [];
	push @{$sqlstack}, qq(DROP TABLE IF EXISTS $geco_out;);
	push @{$sqlstack}, qq(CREATE TABLE IF NOT EXISTS $geco_out (seq_region INT NOT NULL, gene_id INT NOT NULL, transcript_tag INT NOT NULL, length INT NOT NULL, avg_coverage FLOAT NOT NULL, max_coverage FLOAT NOT NULL, rpkm FLOAT NOT NULL, intercept FLOAT NOT NULL, slope FLOAT NOT NULL, coeff FLOAT NOT NULL););
	push @{$sqlstack}, qq(LOAD DATA LOCAL INFILE \"$geco_out\" IGNORE INTO TABLE $geco_out;);

	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	foreach my $sql (@{$sqlstack}) { print STDERR "\."; $dbh->prepare( $sql )->execute() }
	
#	print STDERR " -> Geco Score: " . join " ", @{geco_score($exp)} . "\n";
	
	return;	
}

#############
## HELPERS ##
#############
#NG#
sub _build_core {
	my $fh = $_[0];
	my $e  = $_[1];
	my $g  = $_[1]->[0]->{stable_id}; #geneid
	my $c  = $_[1]->[0]->{name};      #chromosome
	my $s  = $_[1]->[0]->{seq_region_strand}; #strand
	
	#union the exons and cound coverage for each base
	my %deltamap;
	for (my $i = 0; $i < scalar(@{$e}); $i++) {
		$deltamap{$e->[$i]->{seq_region_start}}++;
		$deltamap{$e->[$i]->{seq_region_end}+1}--;
	}

	my @newexons; #[start,end,score]
	my $pre_start = 0;
	my $pre_score = 0;
	my %scores; # count occurences
	foreach my $s (sort {$a <=> $b} keys %deltamap) {
		if ($deltamap{$s} == 0) { next } # can be skipped as there is no change in coverage
		if ($pre_start != 0 && $pre_score != 0) {
			#create slice
			push @newexons, [ $pre_start, $s-1, $pre_score];
			$scores{$pre_score}++;
		}
		$pre_start  = $s;
		$pre_score += $deltamap{$s};
	}
	die "ABORT TRAP: algoritm did not end properly" if ($pre_score != 0 || $pre_start == 0);

	my $maxscore = 0;
	my $midscore = 0; # the other cutoff
	my $esum = 0;
	# max is not always best (allow for max and max - 1 if max-1 has more members than max)
	foreach my $sc (sort {$b <=> $a} keys %scores) {
		if ($sc > $maxscore)       { $maxscore = $sc }
		if ($scores{$sc} >= $esum) { $midscore = $sc }
		$esum += $scores{$sc};
	}
	
	# determine what exon to tag (2: conservative, 1: old core method)
	my $exonid;
	foreach my $e (@newexons) {
		print $fh $g . "\t";
		print $fh $g . '.' . ++$exonid . "\t";
		print $fh $c . "\t";
		print $fh $e->[0] . "\t";
		print $fh $e->[1] . "\t";
		print $fh $s . "\t";
		if ($e->[2] == $maxscore)    { print $fh "2" }
		elsif ($e->[2] >= $midscore) { print $fh "1" }
		else                         { print $fh "0" }
		print $fh "\n";
	}
	return;
}

#LEGACY
sub _core_print {
	my $fh = $_[0];
	my $e  = $_[1];
	my $c  = [];
	my $tagged = 0;
	for (my $i = 0; $i < scalar(@{$e}); $i++) {
		$c->[$i]++; # sets to 1
		for (my $j = 0; $j < scalar(@{$e}); $j++) {
			next if ($i==$j);
			if (($e->[$i]->{seq_region_start} <= $e->[$j]->{seq_region_start})
				&& ($e->[$j]->{seq_region_end} <= $e->[$i]->{seq_region_end})) { $c->[$j]++ }
		}
	}
	# find max
	# max is not always best (allow for max and max - 1 if max-1 has more members than max)
	#OLD#foreach my $m (@{$c}) { $max = ($max < $m) ? $m : $max }
	my $max = -1;
	my $esum = 0; 
	my @classes;
	for(my $x = 0; $x < scalar @{$c}; $x++) { $classes[$c->[$x]]++ }
	for(my $c = scalar(@classes); 0 < $c; $c--) {
		next unless ($classes[$c]);
		if (($max == -1) || ($classes[$c] >= $esum)) {
			$max = $c;
			$esum += $classes[$c];
		}
		else {
			last;
		}
	}
	
	my %printed;
	for (my $a = 0; $a < scalar(@{$e}); $a++) {
		if ( ($c->[$a] >= $max) && !(defined $printed{$e->[$a]->{exon_id}}) ) {
			#print $fh $c->[$a] . ">$max>\t";
			print $fh $e->[$a]->{exon_id} . "\t";
			print $fh $e->[$a]->{seq_region_id} . "\t";
			print $fh $e->[$a]->{seq_region_start} . "\t";
			print $fh $e->[$a]->{seq_region_end} . "\t";
			print $fh $e->[$a]->{seq_region_strand} . "\t";
			print $fh $e->[$a]->{phase} . "\t";
			print $fh $e->[$a]->{end_phase} . "\t";
			print $fh $e->[$a]->{is_current} . "\t";
			print $fh $e->[$a]->{is_constitutive} . "\n";
			$printed{$e->[$a]->{exon_id}}++;
			$tagged++;
		}
	}
	die "Could not tag core exon" if ($tagged == 0);
	return;
}



1;
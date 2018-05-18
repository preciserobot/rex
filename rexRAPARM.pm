package rexRAPARM;

use strict;

use Usu;

#########################
##### REMOTE CONFIG #####
#########################
# => now in experiment module (static)
sub rebuildreadindex {
	my $exp     = $_[0];
	my $pre     = $exp->get_exid;
	my $sqldb   = $exp->get_db;
	my $readfile = $exp->get_read;
	my $readtab  = $exp->reads;
	my $sqlhost = $exp->get_sqlhost;
	my $sqlport = $exp->get_sqlport;
	my $sqluser = $exp->get_sqluser;
	# safe and slow procedure (replaces former sed/perl oneliners)
	print STDERR "\nWriting read index to disk...\n";
	my $seqwrt = 1;
	my $line = 0;
	my $length = 1;
	my $readID;
	my $prevtile = 1;
	if (-s $readtab) { print STDERR "NOTICE: Using Previous read identifier file ($readtab)\n" }
	else {
		open(RIN, "<$readfile") or die $readfile;
		open(ROT, ">$readtab") or die $readtab;
		while (<RIN>) {
			my $chp = $_;
			chomp($chp);
			if    (($chp =~ /^(A|T|G|C|N)+$/) && ($seqwrt == 2))   { die "2|$chp" unless ($seqwrt == 2); $seqwrt++; $length = length($chp) }
			elsif (($chp =~ /^[^\s]{$length}$/) && ($seqwrt == 4)) { die "4|$chp" unless ($seqwrt == 4); $seqwrt = 1; $length = 1 }
			elsif ($chp =~ /^@[^:]+:\d+:(\d+)/) {
				if ($prevtile > $1) { die "READ FILE IS NOT ORDERED BY TILES!" }
				else { $prevtile = $1 }
				die "1|$chp" unless ($seqwrt == 1);
				$chp =~ s/^@//;
				print ROT ++$readID . "\t$chp\n";
				$seqwrt++;
			}
			elsif ($chp =~ /^\+/)              { die "3|$chp" unless ($seqwrt == 3); $seqwrt++ }
			elsif ($chp =~ /^\s*$/)            { next } # empty line (will not count)
			else                               { die ">$chp<$line>$readfile<" } # something else
			if ((++$line % 40000) == 0) { print STDERR "\r ". ($line / 4) . " reads indexed" };
		}
		print STDERR "\r ". ($line / 4) ." reads indexed\n";
		close(ROT);
		close(RIN);
	}
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my $load_sql = [qq( DROP TABLE IF EXISTS \`$readtab\`; ),
					qq( CREATE TABLE IF NOT EXISTS \`$readtab\` (id INT, readid char(64) NOT NULL, PRIMARY KEY (id)); ),
					qq(LOAD DATA LOCAL INFILE \'$readtab\' IGNORE INTO TABLE \`$readtab\` (id, readid)),
					qq(ALTER TABLE \`$readtab\` ADD UNIQUE INDEX (readid);),
					];
	$exp->wait4database;
	foreach my $query (@{$load_sql}) {
		print STDERR "## SQL ## ". substr($query,0,120)." ...\n";
		$dbh->prepare( $query )->execute();
	}
	$dbh->disconnect();
	return;
}

# makes a live readindex (no sql load for build)
sub liveFetch {
	my $exp        = $_[0];
	my $storeindex = $_[1];
	my $pre        = $exp->get_exid;
	my $sqldb      = $exp->get_db;
	my $readfile   = $exp->get_read;
	my $readtab    = $exp->reads;
	my $sqlhost    = $exp->get_sqlhost;
	my $sqlport    = $exp->get_sqlport;
	my $sqluser    = $exp->get_sqluser;
	
	# safe and slow procedure (replaces former sed/perl oneliners)
	print STDERR "\nWriting read index to disk...\n";
	my $seqwrt = 1;
	my $line = 0;
	my $length = 1;
	my $readID;
	my $prevtile = 1;
	if (-s $readtab) { print STDERR "NOTICE: Using Previous read identifier file ($readtab)\n" }
	else {
		open(RIN, "<$readfile") or die $readfile;
		open(ROT, ">$readtab") or die $readtab;
		while (<RIN>) {
			my $chp = $_;
			chomp($chp);
			if    (($chp =~ /^(A|T|G|C|N)+$/) && ($seqwrt == 2))   { die "2|$chp" unless ($seqwrt == 2); $seqwrt++; $length = length($chp) }
			elsif (($chp =~ /^[^\s]{$length}$/) && ($seqwrt == 4)) { die "4|$chp" unless ($seqwrt == 4); $seqwrt = 1; $length = 1 }
			elsif ($chp =~ /^@[^:]+:\d+:(\d+)/) {
				if ($prevtile > $1) { die "READ FILE IS NOT ORDERED BY TILES!" }
				else { $prevtile = $1 }
				die "1|$chp" unless ($seqwrt == 1);
				$chp =~ s/^@//;
				print ROT ++$readID . "\t$chp\n";
				$seqwrt++;
			}
			elsif ($chp =~ /^\+/)              { die "3|$chp" unless ($seqwrt == 3); $seqwrt++ }
			elsif ($chp =~ /^\s*$/)            { next } # empty line (will not count)
			else                               { die ">$chp<$line>$readfile<" } # something else
			if ((++$line % 40000) == 0) { print STDERR "\r ". ($line / 4) . " reads indexed" };
		}
		print STDERR "\r ". ($line / 4) ." reads indexed\n";
		close(ROT);
		close(RIN);
	}
	if ($storeindex > 1) {
		my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
		my $load_sql = [qq( DROP TABLE IF EXISTS \`$readtab\`; ),
						qq( CREATE TABLE IF NOT EXISTS \`$readtab\` (id INT, readid char(64) NOT NULL, PRIMARY KEY (id)); ),
						qq(LOAD DATA LOCAL INFILE \'$readtab\' IGNORE INTO TABLE \`$readtab\` (id, readid))
						];
		if ($storeindex > 2) { push @{$load_sql}, qq(ALTER TABLE \`$readtab\` ADD UNIQUE INDEX (readid);) } ## index
		$exp->wait4database;
		foreach my $query (@{$load_sql}) {
			print STDERR "## SQL ## ". substr($query,0,120)." ...\n";
			$dbh->prepare( $query )->execute();
		}
		$dbh->disconnect();
	}

	# GET SEQ_REGION INDEX
	my $srt       = $exp->annotab('seq_region');
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my %seq_region;
	my ($sr,$chr);
	my $sth = $dbh->prepare( qq{ SELECT seq_region_id, name FROM $srt });
	$sth->execute();
	$sth->bind_columns( \$sr, \$chr );
	while( $sth->fetch ) { $seq_region{$chr} = $sr }
	$dbh->disconnect();

	# LIVE EXTRACT OF SPLICES
	my $gmap  = $exp->mapping;
	my $smap = $exp->splicedbwtout;
	my $genomemap  = $exp->rawmap;
	my $splicemap  = $exp->splicedmap;

	my $idx = *INDEX;
	my $readindex; #{tile}{readid}=readnumber
	my $currenttile = 1;
	open($idx, "<$readtab") or die "cannot open $readtab";
	open(GMAP, "<$gmap") or die "cannot open $gmap";
	open(FOUT, ">$genomemap") or die "cannot open $genomemap";
	#inital indexread
	$readindex = _readIndexByTile($idx, $readindex, $currenttile);
	while (<GMAP>) {
		#      READID......TILE..POS     STR     SRT     POS    RES    MISM
		if (/^([^:]+:\d+:)(\d+)(\S+)\s+([\+|-])\s+(\S+)\s+(\d+)\s+\S+(.*)/) {
			if ($currenttile < $2) {
				$currenttile = $2;
				$readindex = _readIndexByTile($idx, $readindex, $currenttile);
			}
			if ($seq_region{$5} && $readindex->{$2}->{$1.$2.$3}) {
				print FOUT $readindex->{$2}->{$1.$2.$3} . "\t";
				print FOUT (($4 eq '-') ? -1 : 1) . "\t";
				print FOUT $seq_region{$5} . "\t";
				print FOUT ($6+1) . "\t";
				my $mism = $7;
				my $count = 0;
				++$count while $mism =~ /:/g;
				print FOUT $count . "\n";
			}
			else {
				warn "  WARNING: seq_region_id or readid not found!";
				print STDERR "           seq_region_id: $5 => "      . $seq_region{$5} . "\n";
				print STDERR "           readid: ".$1.$2.$3 . " => " . $readindex->{$2}->{$1.$2.$3} . "\n";
				die;
			}
		}
		else { die "cannot parse line " . $_ }
	}
	close(FOUT);
	close(GMAP);
	close($idx);
	print STDERR "\n";
	$readindex = {}; #{tile}{readid}=readnumber
	$currenttile = 1;
	open($idx, "<$readtab") or die "cannot open $readtab";
	open(SMAP, "<$smap") or die "cannot open $smap";
	open(FOUT, ">$splicemap") or die "cannot open $splicemap";
	#inital indexread
	$readindex = _readIndexByTile($idx, $readindex, $currenttile);
	while (<SMAP>) {
		#      READID......TILE..POS     STR     SRT     POS    RES    MISM
		if (/^([^:]+:\d+:)(\d+)(\S+)\s+([\+|-])\s+(\S+)\s+(\d+)\s+\S+(.*)/) {
			if ($currenttile < $2) {
				$currenttile = $2;
				$readindex = _readIndexByTile($idx, $readindex, $currenttile);
			}
			if ($readindex->{$2}->{$1.$2.$3}) {
				print FOUT $readindex->{$2}->{$1.$2.$3} . "\t";
				print FOUT (($4 eq '-') ? -1 : 1) . "\t";
				print FOUT $5 . "\t";
				print FOUT ($6+1) . "\t";
				my $mism = $7;
				my $count = 0;
				++$count while $mism =~ /:/g;
				print FOUT $count . "\n";
			}
			else {
				die "  ERROR: readid for ".$1.$2.$3 . " not found!";
			}
		}
		else { die "cannot parse line " . $_ }
	}
	close(FOUT);
	close(SMAP);
	close($idx);
	print STDERR "\r" . (" " x 200) . "\r\n";
	_getRAPARMannotation($exp);
	
	#STORE MAPPINGS FINAL CLEANUP (IF STORED)
	if ($storeindex > 1) { # prevents indexing
		my $preventindexing = ($storeindex > 2) ? 0 : 1;
		rexDB::loadmappings($exp,$preventindexing);
		die if (system("rm -f $readtab"));
	}
	else {
		$exp->countReads; # uses readfile
	}
	return;
}

sub _readIndexByTile {
	my $ix = $_[0];
	my $ri = $_[1];
	my $ct = $_[2]; # current tile
	
	my $ra = 2; # ahead
	my $rb = 8; # back
	my $mt = 120;
	my $pt = 0; # previous tile
	my $lc = 0; # line counter
	my $line;
	my $rt = 0;
	# read line to get status
	unless($line = defined($ix) ? <$ix> : die ) { return $ri }
	if ($line =~ /^(\d+)\t([^:]+:\d+:)(\d+)(\S+)/) {
		$rt = $3;
		my $name = $2.$3.$4;
		my $id = $1;
		$ri->{$rt}->{$name} = $id;
	}
	else { die $line } 

	while($rt <= ($ct + $ra)) {
		#print STDERR (++$lc % 2) ? "\\" : "/";
		$line = <$ix>;
		unless ($line) { return $ri }
		if ($line =~ /^(\d+)\t([^:]+:\d+:)(\d+)(\S+)/) {
			$rt = $3;
			my $name = $2.$3.$4;
			my $id = $1;
			$ri->{$rt}->{$name} = $id;
			if ($rt > $pt) {
				
				my $dt = $rt % $mt;
				my $pre = $dt-1-$rb;
				print STDERR "\r" . ("." x $pre) . (":" x ( ($pre < 0) ? $rb + $pre : $rb)) . "|" . (":" x (($dt > $mt-$ra ) ? $mt - $dt : $ra)) . ("." x ($mt-$dt)) . " " .  $rt . "    \r";
				$pt = $rt;
			}
		}
		else { die "OOPS at:\n$line\nshould be " . '(\d+)\t([^:]+:\d+:)(\d+)(\S+)' }
	}
	# cleanup no longer needed tiles from readindex
	# everything smaller than ct-rb
	foreach my $tile (keys %{$ri}) {
		if ($tile < ($ct - $rb)) {
			#print STDERR "X".$tile."X";
			delete $ri->{$tile};
		}
	}
	#print STDERR "\n";
	return $ri;
}

sub _getRAPARMannotation {
	my $exp = $_[0];
	my $sqlhost    = $exp->get_sqlhost;
	my $sqlport    = $exp->get_sqlport;
	my $sqluser    = $exp->get_sqluser;
	my $sqldb      = $exp->get_db;
	my $exon       = $exp->exonmodel;
	my $splicigars = $exp->annot("splicigars");
	
	## ANNOTATION TABLES (EXONS,SPLICIGARS)
	my ( $id, $gid, $eid, $rg, $start, $end, $strand, $tag, $uniscore, $spliceid );
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	### cigars
	print STDERR " -> Extracting Exon Annotation\n";
	my $lico = 0;
	open(FOUT, ">$exon") or die;
	my $exonquery = ($exon =~ s/RZexon/RFexon/) ? qq{ SELECT gene_id, exon_id, seq_region, seq_region_start, seq_region_end, seq_region_strand, tag, 0 FROM $exon; } : qq{ SELECT gene_id, exon_id, seq_region, seq_region_start, seq_region_end, seq_region_strand, tag, uniscore FROM $exon; };
	my $sth = $dbh->prepare( $exonquery );
	#$exp->wait4database;
	$sth->execute();
	$sth->bind_columns( \$gid, \$eid, \$rg, \$start, \$end, \$strand, \$tag, \$uniscore );
	while( $sth->fetch ) { if (++$lico % 1000 == 0) { print STDERR "\r$lico" }; print FOUT (join "\t", ($gid, $eid, $rg, $start, $end, $strand, $tag, $uniscore)) . "\n" }
	print STDERR "\r$lico\n";
	close(FOUT);
	### model -> rc/rd/rz
	print STDERR " -> Extracting Splice Annotation\n";
	$lico = 0;
	open(FOUT, ">$splicigars") or die;
	$sth = $dbh->prepare( qq{ SELECT id, seq_region_id, strand, start, end, splice_id FROM $splicigars; } );
	#$exp->wait4database;
	$sth->execute();
	$sth->bind_columns( \$id, \$rg, \$strand, \$start, \$end, \$spliceid );
	while( $sth->fetch ) { if (++$lico % 1000 == 0) { print STDERR "\r$lico" }; print FOUT (join "\t", ($id, $rg, $strand, $start, $end, $spliceid)) . "\n" }
	print STDERR "\r$lico\n";
	close(FOUT);
	$dbh->disconnect();
	
	return;
}


#########################
#########################
#RAPARM-MULTISPLICE

sub quickFetchMappings {
	# don't need to be ordered anymore
	my $exp       = $_[0];
	my $swample   = $_[1];
	my $sqldb     = $exp->get_db;
	my $sqlhost   = $exp->get_sqlhost;
	my $sqlport   = $exp->get_sqlport;
	my $sqluser   = $exp->get_sqluser;
	my $srt       = $exp->annotab('seq_region');

	my $reads      = $exp->reads;
	my $genomemap  = $exp->rawmap;
	my $splicemap  = $exp->splicedmap;
	my $exon       = $exp->exonmodel;
	my $splicigars = $exp->annot("splicigars");
	
	my $sqlexec = "mysql --skip-column-names --quick -u $sqluser -h $sqlhost $sqldb";
		
	my %swampledIDs;
	if ($swample) {
		my $lico;
		#extract all mapped ids
		my $randomfile = Usu::randomfile("idstobeswampled");
		open(IDS, ">$randomfile") or die "Cannot open file $randomfile";
		foreach my $table ($genomemap,$splicemap) {
			my $query = qq{ SELECT r.id FROM $table as m, $reads as r WHERE r.readid = m.readid };
			print STDERR "\n### SQL ### $query\n";
			open(SQL, "echo \"$query\" | $sqlexec |") or die;
			while( <SQL> ) { 
				print IDS;
				if (++$lico % 100000 == 0) { print STDERR "\r$lico" };
			}
			close(SQL);
		}
		close(IDS);
		open(SWAMP, "swample $randomfile $swample |") or die "Cannot find SWAMPLE tool";
		while (<SWAMP>) {
			if (/(\d+)/) {
				$swampledIDs{$1} = 1;
			}
		}
		close(SWAMP);
		system("rm -f $randomfile");
		if (scalar(keys %swampledIDs) != $swample) {
			warn "WARNING: Swampled " . scalar(keys %swampledIDs) . " read but $swample were requested. Press any key to acknowledge\n";
			<STDIN>;
		}
	}
	
	
	# get seq_region_hashtable
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my %seq_region;
	my ($sr,$chr);
	my $sth = $dbh->prepare( qq{ SELECT seq_region_id, name FROM $srt });
	$sth->execute();
	$sth->bind_columns( \$sr, \$chr );
	while( $sth->fetch ) { $seq_region{$chr} = $sr }
	$dbh->disconnect();
	
	# extract 4 files (2X SPLICE, CIGARS, RCEXON)
	my ($lico,$expected) = (0,0);
	print STDERR " -> Extracting Genome Mappings\n";
	# count mappings to check
	my $qqq = qq(SELECT COUNT(*) FROM \\\`$genomemap\\\`);
	open(SQL, "echo \"$qqq\" | $sqlexec |") or die;
	while (<SQL>) {
		if (/^(\d+)/) {
			$expected = $1;
			$lico = $1;
			last;
		}
		else { die "Could not count expected mappings" }
	}
	close(SQL);
	# extract
	my $skipextract = 0;
	if (-s $genomemap) {
		open(LICO, "wc -l $genomemap |") or die;
		while (<LICO>) {
			if (/(\d+)/ && $1 == $expected && !($swample)) {
				print "\rNOTICE: Using previously extracted mapping file ($genomemap) as it appears to be complete ($expected lines)\n";
				++$skipextract;
			}
		}
		close(LICO);
	}
	unless ($skipextract) {
		my $query = qq( SELECT r.id, m.strand, m.chr, m.pos + 1, LENGTH(mismatch) - LENGTH(REPLACE(mismatch, ':', '')) FROM \\\`$reads\\\` as r, \\\`$genomemap\\\` as m WHERE m.readid = r.readid );
		$exp->wait4database;
		open(FOUT, ">$genomemap") or die;
		open(SQL, "echo \"$query\" | $sqlexec |") or die;
		while (<SQL>) {
			if (/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
				if (--$lico % 100000 == 0) { print STDERR "\r$lico " };
				if ($swample && !($swampledIDs{$1})) { next }
				if ($seq_region{$3}) {
					print FOUT $1 . "\t";
					print FOUT (($2 eq '-') ? -1 : 1) . "\t";
					print FOUT $seq_region{$3} . "\t";
					print FOUT $4 . "\t";
					print FOUT $5 . "\n";
				}
				else { warn "  WARNING: seq_region_id for $3 not found!" }
			}
			else { die }
		}
		print STDERR "\r". ($expected - $lico) ."\n";
		if (!($swample) && abs($lico/$expected) > 0.01) {
			warn "WARNING: SQL extract may be incomplete. Expected $expected lines and got " .($expected - $lico). ". Acknowledge by pressing <ENTER>\n";
			<STDIN>;
		}
		close(SQL);
		close(FOUT);
	}

	print STDERR " -> Extracting Splice Mappings\n";
	# count mappings to check
	($expected,$lico) = (0,0);
	$qqq = qq(SELECT COUNT(*) FROM \\\`$splicemap\\\`);
	open(SQL, "echo \"$qqq\" | $sqlexec |") or die;
	while (<SQL>) {
		if (/^(\d+)/) {
			$expected = $1;
			$lico = $1;
			last;
		}
		else { die "Could not count expected mappings" }
	}
	close(SQL);
	# extract
	$skipextract = 0;
	if (-s $splicemap) {
		open(LICO, "wc -l $splicemap |") or die;
		while (<LICO>) {
			if (!($swample) && /(\d+)/ && $1 == $expected) {
				print "\rNOTICE: Using previously extracted mapping file ($splicemap) as it appears to be complete ($expected lines)\n";
				++$skipextract;
			}
		}
		close(LICO);
	}
	unless ($skipextract) {
		my $query = qq{ SELECT r.id, m.strand, m.splice, m.pos + 1, LENGTH(mismatch) - LENGTH(REPLACE(mismatch, ':', ''))  FROM \\\`$reads\\\` as r, \\\`$splicemap\\\` as m WHERE m.readid = r.readid; } ;
		$exp->wait4database;
		open(FOUT, ">$splicemap") or die;
		open(SQL, "echo \"$query\" | $sqlexec |") or die;
		while (<SQL>) {
			if (/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
				if (--$lico % 100000 == 0) { print STDERR "\r$lico " };
				if ($swample && !($swampledIDs{$1})) { next }
				print FOUT $1 . "\t" . (($2 eq '-') ? -1 : 1) . "\t" . $3 . "\t" . $4 . "\t" . $5 . "\n";
			}
			else { die }
		}
		print STDERR "\r". ($expected - $lico) ."\n";
		if (!($swample) && abs($lico/$expected) > 0.01) {
			warn "WARNING: SQL extract may be incomplete. Expected $expected lines and got " .($expected - $lico). ". Acknowledge by pressing <ENTER>\n";
			<STDIN>;
		}
		close(SQL);
		close(FOUT);
	}
	_getRAPARMannotation($exp);
	$dbh->disconnect();
	$exp->updateReadCount;
	return;
}

sub runRaparm { 
	my $exp        = $_[0];
	my $local      = $_[1];
	my $override   = $_[2];
	my $swampled   = $_[3];
	my $details    = $_[4];
	my $um         = $_[5];
	my $fast       = $_[6];
	
	
	#input
	my $exon       = $exp->exonmodel;
	my $genomemap  = ($exp->get_paired) ? $exp->smackout($exp->rawmap)     : $exp->rawmap($override);
	my $splicemap  = ($exp->get_paired) ? $exp->smackout($exp->splicedmap) : $exp->splicedmap($override);
	my $splicigars = $exp->annot("splicigars");
	my $in_tag     = $exp->get_atag;
	my $readlength = $exp->get_readl;
	my $strsp      = $exp->get_strands; ##### store in protocol! #############################################################################
	my $extend     = $exp->get_retros;  ##### store in protocol! #############################################################################
	my $bestmatch  = $exp->get_bestmap; ##### store in protocol! #############################################################################
	my $strands    = 0; # replace by switch for strand specificity
	#output
	my $out_exp = $exp->raparm_rpkm;
	my $out_isl = $fast ? '-' : $exp->raparm_island;
	my $out_cov = $fast ? '-' : $exp->raparm_coverage;
	my $out_uni = $fast ? '-' : $exp->raparm_unicoverage;
	my $out_spf = $fast ? '-' : $exp->raparm_splicefreq;
	my $out_map = $details ? $exp->raparm_map : "-";
	
	# rmms (rename files)
	if ($um) {
		$out_exp =~ s/raparm/rmms/;
		$out_isl =~ s/raparm/rmms/;
		$out_cov =~ s/raparm/rmms/;
		$out_uni =~ s/raparm/rmms/;
		$out_spf =~ s/raparm/rmms/;
		$out_map =~ s/raparm/rmms/;
	}
	
	unless (-s $exon) {
		_getRAPARMannotation($exp);
	}
	
	if ($swampled) {
		# append number
		$out_exp .= "_sampled" . $swampled;
		$out_isl .= "_sampled" . $swampled;
		$out_cov .= "_sampled" . $swampled;
		$out_uni .= "_sampled" . $swampled;
		$out_spf .= "_sampled" . $swampled;
		if ($details) { $out_map .= "_sampled" . $swampled }
	}
	
	# extract splicecigar if not exists
	unless (-s $splicigars) {
		my $sqlhost = $exp->get_sqlhost;
		my $sqluser = $exp->get_sqluser;
		my $sqlport = $exp->get_sqlport;
		my $sqldb   = $exp->get_db;
		my $localextract = "mysql --skip-column-names -u $sqluser -h $sqlhost $sqldb";
		my $cigarextract = qq{ SELECT id, seq_region_id, strand, start, end, splice_id FROM $splicigars };
		die "ERROR: splicigar extract failed from table $splicigars" if system("echo \"$cigarextract\" | $localextract > $splicigars");
	}
	
	
	# remote
	my $remote    = rexRemote->new();
	my $raparmbin = $um ? $remote->get_rmms_bin : $remote->get_raparm_bin;
	# execute string
	my $io = "$exon $genomemap $splicemap $splicigars $readlength $in_tag $strsp $extend $bestmatch $out_map $out_exp $out_isl $out_cov $out_uni $out_spf";
	if ($local) {
		my $exec = $raparmbin . " " . $io;
		print STDERR "$exec\n";
		system("$exec");
	} 
	else {
		my $exec = $raparmbin . " " . $io;
		if ((-s $exon) && (-s $genomemap) && (-s $splicemap) && (-s $splicigars)) {
			# run string
			print STDERR "$exec\n";
			# compile
			$remote->sendCompile($um ? "rmms" : "raparm");
			# send
			$remote->sendFile($exon);
			$remote->sendFile($genomemap);
			$remote->sendFile($splicemap);
			$remote->sendFile($splicigars);
			# run 
			$remote->execute($exec);
			# fetch
			if ($details) { $remote->fetchFile($out_map) }
			$remote->fetchFile($out_exp);
			$remote->fetchFile($out_isl);
			$remote->fetchFile($out_cov);
			$remote->fetchFile($out_uni);
			$remote->fetchFile($out_spf);
			# cleanup
			$remote->cleanup();
		}
		else { die "cannot find required files for RAPARM call ($exec)" }
	}
	return;
}

sub storeRaparmOutput {
	my $exp      = $_[0];
	my $strand   = $_[1]; # or any other suffix
	my $charstrand = ($strand) ? (($strand =~ /\+/) ? "fwd" : ( ($strand =~ /-/) ? "rev" : $strand )) : undef ;
	
	my $sqldb    = $exp->get_db;
	my $sqlhost  = $exp->get_sqlhost;
	my $sqlport  = $exp->get_sqlport;
	my $sqluser  = $exp->get_sqluser;
	
	my $out_map = $exp->raparm_map($strand);
	my $out_exp = $exp->raparm_rpkm($strand);
	my $out_isl = $exp->raparm_island($strand);
	my $out_cov = $exp->raparm_coverage($strand);
	my $out_uni = $exp->raparm_unicoverage($strand);
	my $out_spl = $exp->raparm_splicefreq($strand);
	
	
	#sql tables
	my $sql_map = $exp->raparm_map($charstrand);
	my $sql_exp = $exp->raparm_rpkm($charstrand);
	my $sql_isl = $exp->raparm_island($charstrand);
	my $sql_cov = $exp->raparm_coverage($charstrand);
	my $sql_uni = $exp->raparm_unicoverage($charstrand);
	my $sql_spl = $exp->raparm_splicefreq($charstrand);
	
	my $sqlstack = [];
 	
	# expression
	push @{$sqlstack}, qq(DROP TABLE IF EXISTS $sql_exp;);
	push @{$sqlstack}, qq(CREATE TABLE IF NOT EXISTS $sql_exp
						(geneid INT NOT NULL, type CHAR(16), region INT NOT NULL, start INT NOT NULL, end INT NOT NULL, strand INT NOT NULL, 
						tagged_length INT NOT NULL, mappedreads INT NOT NULL, rpkm FLOAT NOT NULL););
	push @{$sqlstack}, qq(LOAD DATA LOCAL INFILE \"$out_exp\" IGNORE INTO TABLE $sql_exp;);
	push @{$sqlstack}, qq(ALTER TABLE $sql_exp ADD INDEX (geneid););

	# islands
	if (-s $out_isl) {
		push @{$sqlstack}, qq(DROP TABLE IF EXISTS $sql_isl;);
		push @{$sqlstack}, qq(CREATE TABLE IF NOT EXISTS $sql_isl
						(seq_region INT NOT NULL, source CHAR(32), type CHAR(32), start INT NOT NULL, end INT NOT NULL, score FLOAT NOT NULL, 
						strand CHAR(8) NOT NULL, phase CHAR(8) NOT NULL, attributes TEXT););
		push @{$sqlstack}, qq(LOAD DATA LOCAL INFILE \"$out_isl\" IGNORE INTO TABLE $sql_isl;);
	}
	
	# coverage
	if (-s $out_cov) {
		push @{$sqlstack}, qq(DROP TABLE IF EXISTS $sql_cov;);
		push @{$sqlstack}, qq(CREATE TABLE IF NOT EXISTS $sql_cov (region INT NOT NULL, strand INT, start INT NOT NULL, end INT NOT NULL, length INT NOT NULL, coverage FLOAT NOT NULL););
		push @{$sqlstack}, qq(LOAD DATA LOCAL INFILE \"$out_cov\" IGNORE INTO TABLE $sql_cov;);
	}
	
	# unique coverage
	if (-s $out_uni) {
		push @{$sqlstack}, qq(DROP TABLE IF EXISTS $sql_uni;);
		push @{$sqlstack}, qq(CREATE TABLE IF NOT EXISTS $sql_uni (region INT NOT NULL, strand INT, start INT NOT NULL, end INT NOT NULL, length INT NOT NULL, coverage FLOAT NOT NULL););
		push @{$sqlstack}, qq(LOAD DATA LOCAL INFILE \"$out_uni\" IGNORE INTO TABLE $sql_uni;);
	}
	
	#splicefreq
	if (-s $out_spl) {
		push @{$sqlstack}, qq(DROP TABLE IF EXISTS $sql_spl;);
		push @{$sqlstack}, qq(CREATE TABLE IF NOT EXISTS $sql_spl (splice INT NOT NULL, coverage FLOAT NOT NULL););
		push @{$sqlstack}, qq(LOAD DATA LOCAL INFILE \"$out_spl\" IGNORE INTO TABLE $sql_spl;);
	}

# DISABLED as too large data (and not really needed)
#	# mapping protocol
#	if (-s $out_map) {
#		push @{$sqlstack}, qq(DROP TABLE IF EXISTS $sql_map;);
#		push @{$sqlstack}, qq(CREATE TABLE IF NOT EXISTS $sql_map
#							(readid INT NOT NULL, flag INT NOT NULL, geneid INT NOT NULL, exonid INT NOT NULL, 
#							ovpfrac FLOAT NOT NULL, expfrac FLOAT NOT NULL, region INT NOT NULL, strand INT, cigarline TEXT););
#		push @{$sqlstack}, qq(LOAD DATA LOCAL INFILE \"$out_map\" IGNORE INTO TABLE $sql_map;);
#	}

	for(my $i=0; $i < scalar(@{$sqlstack}); $i++) { $sqlstack->[$i] =~ s/\n\s+/ /g }
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	foreach my $sql (@{$sqlstack}) { print STDERR "## SQL ## ". substr($sql,0,120). ((length($sql) > 120) ? " ...\n" : "\n"); $dbh->prepare( $sql )->execute() }
	return;
}

1;

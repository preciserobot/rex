package rexDB;

use strict;

use DBI;
use DBD::mysql;
use Digest::MD5 qw(md5_hex);
use Time::localtime;


#ANNOTATION OVERRIDE (still testing)
sub cloneLibrary {
	my $exp        = $_[0];
	my $templateid = $_[1];
	
	die "ERROR: Override identifier is mandatory!" unless $_[1];
	die "ERROR: Override and current AnnotationID are equal! -> Check control File\n" if ($_[1] eq $exp->annotationID);

	my $sqldb     = $exp->get_db;
	my $sqlhost   = $exp->get_sqlhost;
	my $sqlport   = $exp->get_sqlport;
	my $sqluser   = $exp->get_sqluser;
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	
	# clone SQL TABLE splicigars/seq_region/gene_stable_id
	my @clonetables = (
		[$exp->annot("splicigars", $templateid), $exp->annot("splicigars")],
		[$exp->annot("seq_region", $templateid), $exp->annot("seq_region"), ['seq_region_id','name']],
		[$exp->annot("gene_stable_id", $templateid), $exp->annot("gene_stable_id"), ['gene_id','stable_id']],
		);
	foreach my $t (@clonetables) {
		print STDERR " -> CLONING $t->[0] <=> $t->[1] ";
		$dbh->prepare( qq{ DROP TABLE IF EXISTS $t->[1] } )->execute(); print STDERR ".";
		$dbh->prepare( qq{ CREATE TABLE $t->[1] SELECT * FROM $t->[0] } )->execute(); print STDERR ".";
		if ($t->[2]) {
			foreach my $idx (@{$t->[2]}) {
				$dbh->prepare( qq{ ALTER TABLE $t->[1] ADD UNIQUE INDEX ($idx) } )->execute(); print STDERR ".";
			}
		}
		 print STDERR "\n";
	}
	# hardlink splice reference
	my $templatefile = $exp->annot("splicigarseq", $templateid);
	my $targetfile = $exp->annot("splicigarseq");
	opendir(DIR, ".");
	foreach my $file (grep(/\.ebwt$/,readdir(DIR))) {
		#print STDERR "FILE: $file $templatefile";
		if ($file =~ /$templatefile(\.(?:rev\.)?\d\.ebwt)/) {
			print STDERR " -> LINKING $file <=> $targetfile".$1."\n";
			unless (-s $targetfile.$1) { system("ln $file $targetfile" . $1) }
		}
	}
	closedir(DIR);
	
	# clone splice redundancy from override
	my $an = $exp->get_annotid;
	my $thisannotid = $exp->annotationID;
	$dbh->prepare( qq{ UPDATE $an as a, $an as b SET a.splice_seq_redundancy = b.splice_seq_redundancy WHERE a.annotationID = \"$thisannotid\" and b.annotationID = \"$templateid\" } )->execute();
	$dbh->disconnect();
	
	return;
}

#NG#
sub writeCombinedExonFile {
	# preapres core not core exons and retros
	my $exp      = $_[0];
	my $efile    = $exp->rfexon;
	my $sr       = $exp->annot("seq_region");
	my $gs       = $exp->annot("gene_stable_id");
	my $sqldb    = $exp->get_db;
	my $sqlhost  = $exp->get_sqlhost;
	my $sqlport  = $exp->get_sqlport;
	my $sqluser  = $exp->get_sqluser;
	my $cnt = 0; # counte the retrieved exons
	# OPEN and PREPARE DBI handle
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my ($sql,$sth);
	my $maxindex = 0;
	
	if (-s $exp->get_cexo) {
		# load into mysql, extract using seq_regions and gene_id
		my $cefile = $exp->get_cexo;
		my $ce = $exp->annot($cefile);
		my $loadcustom;
		push @{$loadcustom}, qq{DROP TABLE IF EXISTS $ce;};
		push @{$loadcustom}, qq{CREATE TABLE $ce (stable_id CHAR(64), exon_id INT AUTO_INCREMENT, exon_name CHAR(64), chr CHAR(32), start INT, end INT, strand INT, tag INT, PRIMARY KEY (exon_id));};
		push @{$loadcustom}, qq{LOAD DATA LOCAL INFILE \"$cefile\" INTO TABLE $ce (stable_id, exon_name, chr, start, end, strand, tag) };
		foreach my $query (@{$loadcustom}) { print STDERR "### SQL ### " . substr($query,0,120) . " ...\n"; $dbh->prepare( $query )->execute() }
		$sql = qq{SELECT DISTINCT gs.gene_id, c.exon_id, s.seq_region_id, c.start, c.end, c.strand, c.tag FROM $ce as c LEFT JOIN $sr as s ON c.chr = s.name LEFT JOIN $gs AS gs ON c.stable_id = gs.stable_id};
	}
	else { die "ERROR: Exon file not found" }

	print STDERR "### SQL ### " . substr($sql,0,120). " ...\n";
	open(EOUT, ">$efile") or die;
	$sth = $dbh->prepare( $sql );
	my( $aa, $bb, $cc, $dd, $ee, $ff, $gg );
	$sth->execute();
	$sth->bind_columns( \$aa, \$bb, \$cc, \$dd, \$ee, \$ff, \$gg );
	while( $sth->fetch() ) {
		if ($aa =~ /NULL/ || $cc =~ /NULL/) { die ">>> SQL QUERY ERROR: asymmetric table join <<<" }
		print STDERR "$cnt\r" if (++$cnt % 1000 == 0);
		if ($ff =~ /\d+/)  { }
		elsif ($ff eq '+') { $ff =  1 }
		elsif ($ff eq '-') { $ff = -1 }
		else               { die }
		print EOUT "$aa\t$bb\t$cc\t$dd\t$ee\t$ff\t$gg\n";
	}
	close(EOUT);
=pod	
	###### RETROS
	#my $gsid = $exp->annot("gene_stable_id");
	#my $g    = $exp->annot("gene");
	################
		################
			################
				#### LOAD RTAB!!!!!!
			################
		################
	################
	# calculate retro offset ID
	my $offset = 100000;
	while($maxindex > $offset) { die if ($offset > (2**24)); $offset += 100000 }
	print STDERR "Proprietary Retros have index offset of $offset\n";
	if ($exp->get_rtab ne "none" && $exp->get_rtab ne "") {
		# create reto unique ids
		print STDERR "Building Retro Index Table\n";
		my $requery = [];
		my $retroidtable = $exp->retroidtable;
		my $rtab = $exp->get_rtab;
		push @{$requery}, qq{DROP TABLE IF EXISTS $retroidtable;};
		push @{$requery}, qq{CREATE TABLE $retroidtable (internal INT AUTO_INCREMENT PRIMARY KEY, retro_id INT NULL, retro_name CHAR(32), seq_region_id INT, start INT, end INT, strand CHAR(2), parent CHAR(32), disab CHAR(5), INDEX (retro_id));};
#		push @{$requery}, qq{INSERT INTO $retroidtable (retro_name, seq_region_id, start, end, strand, parent, disab) select r.retro_id, srt.seq_region_id, r.start, r.end, r.strand, r.parent_ensg, r.disab from $rtab as r left join $gsid as gs on r.retro_ensg = gs.stable_id left join $g as g on gs.gene_id = g.gene_id, $sr as srt where (g.biotype IS NULL or g.biotype != "protein_coding") and srt.name = r.chrom ORDER BY r.retro_id;};
		push @{$requery}, qq{INSERT INTO $retroidtable (retro_name, seq_region_id, start, end, strand, parent, disab) select r.retro_id, srt.seq_region_id, r.start, r.end, r.strand, r.parent_ensg, r.disab from $rtab as r, $sr as srt where srt.name = r.chrom AND r.retro_ensg NOT LIKE "ENS%" ORDER BY r.retro_id;};
		push @{$requery}, qq{UPDATE $retroidtable set retro_id = internal+$offset};  ### calculate the retro ids
		
		foreach my $query (@{$requery}) { print STDERR "\."; $dbh->prepare( $query )->execute() }
		# extract retros for annotation
		$sql = qq{ select distinct retro_id, seq_region_id, start, end, strand from $retroidtable; };
		$sth = $dbh->prepare( $sql );
		my( $Q_id, $Q_sr, $Q_start, $Q_end, $Q_strand );
		$sth->execute();
		$sth->bind_columns( \$Q_id, \$Q_sr, \$Q_start, \$Q_end, \$Q_strand );
		while( $sth->fetch() ) {
			print STDERR "$cnt\r" if (++$cnt % 1000 == 0);
			if ($Q_strand =~ /\d+/)  { }
			elsif ($Q_strand eq '+') { $Q_strand =  1 }
			elsif ($Q_strand eq '-') { $Q_strand = -1 }
			else                     { die }
			print EOUT "$Q_id\t$Q_id\t$Q_sr\t$Q_start\t$Q_end\t$Q_strand\t1\n";
		}
	}
	else { print STDERR "NOTICE: File (" . $exp->get_rtab . ") not found -> skipping Retros\n" }
=cut	
	# load annotations into global annotation table (retro + core exons)
	print STDERR "Storing";
	my $sqlcreate = [];
	push @{$sqlcreate}, qq{DROP TABLE IF EXISTS $efile;};
	push @{$sqlcreate}, qq{CREATE TABLE $efile (gene_id INT, exon_id INT, seq_region INT, seq_region_start INT, seq_region_end INT, seq_region_strand INT, tag INT);};
	push @{$sqlcreate}, qq{LOAD DATA LOCAL INFILE \"$efile\" INTO TABLE $efile};
	foreach my $query (@{$sqlcreate}) { print STDERR "\."; $dbh->prepare( $query )->execute() }
	print STDERR "DONE\n";
	$dbh->disconnect();
	return;
}

#NG
sub buildIntegerIndex {
	my $e = $_[0];
	my $gs = $e->annot("gene_stable_id");
	my $sr = $e->annot("seq_region");
	my $sqldb   = $e->get_db;
	my $sqlhost = $e->get_sqlhost;
	my $sqluser = $e->get_sqluser;

	my %seq_region;
	my %gene;
	my $gene_tmp   = md5_hex(rand()) . "_gene";
	my $region_tmp = md5_hex(rand()) . "_region";
	
	open(EXIN, $e->get_cexo) or die "ERROR: File (". $e->get_cexo .") not found!";
	while(<EXIN>) {
		if (/^(\S+)\s+\S+\s+(\S+)/) {
			$gene{$1}++;
			$seq_region{$2}++;
		}
		else { warn "line not recognized" }
	}
	close(EXIN);

	open(REGION, ">$region_tmp") or die "Cannot write to file ($region_tmp)";
	foreach my $r (keys %seq_region) { print REGION "$r\n" }
	close(REGION);
	
	open(GENE, ">$gene_tmp") or die "Cannot write to file ($gene_tmp)";
	foreach my $g (keys %gene) { print GENE "$g\n" }
	close(GENE);
	
	my $sql = [
		qq(DROP TABLE IF EXISTS $gs;),
		qq(DROP TABLE IF EXISTS $sr;),
		qq(CREATE TABLE IF NOT EXISTS $gs (gene_id INT AUTO_INCREMENT, stable_id char(128) NOT NULL, PRIMARY KEY (gene_id));),
		qq(CREATE TABLE IF NOT EXISTS $sr (seq_region_id INT AUTO_INCREMENT, name char(40) NOT NULL, PRIMARY KEY (seq_region_id), UNIQUE INDEX idx_name (name));),
		qq(LOAD DATA LOCAL INFILE \'$gene_tmp\'   IGNORE INTO TABLE $gs (stable_id);),
		qq(LOAD DATA LOCAL INFILE \'$region_tmp\' IGNORE INTO TABLE $sr (name);),
		qq(ALTER TABLE $gs ADD UNIQUE INDEX (stable_id);),
		qq(ALTER TABLE $sr ADD UNIQUE INDEX (name);),
	];
	foreach my $i (@{$sql}) {
		print STDERR "## SQL ## ". substr($i,0,120)." ...\n";
		die if (system("echo \"$i\" | mysql -u $sqluser -h $sqlhost $sqldb"));
	}
	die if (system("rm -f $gene_tmp $region_tmp"));
	return;
}

#NG 
sub buildExonSeqs {
	# IMPROVEMENT OVER LEGACY Building also joint (overlapping?) exons (without linking junction);
	
	my $exp     = $_[0];
	my $exon    = $exp->rfexon;
	my $exonseq = $exp->annot("exonseqs");
	my $sr      = $exp->annot("seq_region");
	my $db      = $exp->get_refs;
	my $sqlhost = $exp->get_sqlhost;
	my $sqlport = $exp->get_sqlport;
	my $sqluser = $exp->get_sqluser;
	my $sqldb   = $exp->get_db;
	my $seqcounter;
	my %e;
	my $localexecute = "mysql --skip-column-names --quick -u $sqluser -h $sqlhost $sqldb";
	my %seqregion;
	open(SEQ, "echo \"SELECT DISTINCT seq_region_id, name from $sr\" | $localexecute |") or die;
	while (<SEQ>) {
		if (/^(\d+)\s+(\S+)/) { $seqregion{$1} = $2 }
	}
	close(SEQ);
	
	
	open(EIN, "echo \"SELECT gene_id, exon_id, seq_region, seq_region_start, seq_region_end, seq_region_strand from $exon\" | $localexecute |") or die;
	while(<EIN>) {
		if (/^(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)/) { # geneid, exonid, seq_region, start, end, strand
			if ($4 > $5) { warn "WARNING: Ignored exon with start > end ($2)\n" }
			else {
				push @{$e{$seqregion{$3}}{$1}}, [$2, $4, $5, $6];
			}
		}
		else {die}
	}
	close(EIN);
	
	# make a fused exon list
	print STDERR "\nDetecting Fusion Exons...\n";
	my %f; # $f{$seqregion{$3}}, [ids, start, end]
	my $fusionsum = 0;
	foreach my $r (keys %e) {
		print STDERR "\r  $r  \r";
		foreach my $g (keys %{$e{$r}}) {
			# print genes 
			push @{$f{$r}{$g}}, _resolve($e{$r}{$g});
			$fusionsum += scalar(@{$f{$r}{$g}});
		}
	}
	print STDERR "$fusionsum merged exons\nGenerating Exon database...\n";
	open(SPDB, ">$exonseq") or die "cannot write to $exonseq";
	my $stream  = Bio::SeqIO->new(-file => $db ,-format => rexDB::_chooseformat($db));
	while ( my $seq = $stream->next_seq() ) {
		if ($e{$seq->primary_id}) {
			#print STDERR "\n\t" . $seq->primary_id . " => " . scalar(keys %{$e{$seq->primary_id}}) . "\n";
			foreach my $r (keys %{$e{$seq->primary_id}}) {
				foreach my $x (@{$e{$seq->primary_id}{$r}}) {
					print STDERR "\r".(" " x 120)."\r\t" . $seq->primary_id . "\t" . ++$seqcounter . "       ";
					print SPDB ">".$x->[0]."\n";
					print SPDB ($x->[3] < 0) ? _revcom($seq->subseq($x->[1], $x->[2])) : $seq->subseq($x->[1], $x->[2]);
					print SPDB "\n";
				}
				#fusion exons
				foreach my $x (@{$f{$seq->primary_id}{$r}}) {
					print STDERR "\r".(" " x 120)."\r\t" . $seq->primary_id . "\t" . ++$seqcounter . "       ";
					print SPDB ">".$x->[0]."\n";
					print SPDB ($e{$seq->primary_id}{$r}->[0]->[3] < 0) ? _revcom($seq->subseq($x->[1], $x->[2])) : $seq->subseq($x->[1], $x->[2]);
					print SPDB "\n";
				}
			}
		}
	}
	close(SPDB);
	print STDERR "\n";
	return $exonseq;
}

#NG
sub loadJunctions {
	my $exp  = $_[0];
	my $jadd = $_[1];
	my $jn = $exp->annot("junctions");
	my $sr = $exp->annot("seq_region");
	my $sqlhost = $exp->get_sqlhost;
	my $sqluser = $exp->get_sqluser;
	my $sqlport = $exp->get_sqlport;
	my $sqldb   = $exp->get_db;
	my $localexecute = "mysql -u $sqluser -h $sqlhost $sqldb";
	my @s;
	# drop/create exon rank table
	$s[0] = qq(DROP TABLE IF EXISTS $jn;);
	# create intron table
	$s[1] = qq(CREATE TABLE $jn (id INT, seq_region_id INT NOT NULL, strand INT NOT NULL, start INT NOT NULL, end INT NOT NULL, type TEXT, PRIMARY KEY (id), UNIQUE INDEX idx_intron (seq_region_id, strand, start, end), INDEX idx_region (seq_region_id), INDEX idx_strand (strand), INDEX idx_start (start), INDEX idx_end (end) ););
	for(my $i=0; $i < scalar(@s); $i++) { print STDERR "## SQL ## ". substr($s[$i],0,64) . " ...\n"; die if (system("echo \"$s[$i]\" | $localexecute")) }

	# create seq_region_index
	my %seqregion;
	my $sel = qq(SELECT DISTINCT seq_region_id, name from $sr);
	open(SEQ, "echo \"$sel\" | $localexecute |") or die;
	while (<SEQ>) {
		if (/^(\d+)\s+(\S+)/) {
			$seqregion{$2} = $1;
		}
	}
	close(SEQ);
	
	# write temporary file and load
	my $jeff = "tmp_junctions_" . md5_hex(rand());
	open(FIN, "<$jadd") or die "cannot open $jadd";
	open(FOUT, ">$jeff") or die "cannot open $jeff";
	while(<FIN>) {
		if (/^(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)/) {
			unless ($seqregion{$2}) {
				die "Chromosome name $2 not recognized";
			}
			print FOUT $1 . "\t";
			print FOUT $seqregion{$2} . "\t";
			print FOUT $3 . "\t";
			print FOUT $4 . "\t";
			print FOUT $5 . "\t";
			print FOUT $6 . "\n";
		}
	}
	close(FIN);
	close(FOUT);
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "DB connection error: $DBI::errstr"; 
	$dbh->prepare( qq(LOAD DATA LOCAL INFILE \"$jeff\" INTO TABLE $jn) )->execute();
	$dbh->disconnect;
	die if system("rm -f $jeff");
	return;
}

#NG (for faster reanalysis without having the files, indexes mappings and reads)
sub indexmappings {
		my $exp = $_[0];
		my $sqldb   = $exp->get_db;
		my $sqlhost = $exp->get_sqlhost;
		my $sqluser = $exp->get_sqluser;
		my $sqlport = $exp->get_sqlport;
		
		# HEADER
		my ($sec,$min,$hour,$mday,$mon,$yr,$wday,$yday,$isdst)=@{localtime(time)};
		my $pre = "I";
		my $post = "Indexing mapping and read IDs";
		print STDERR $pre." > ";
		printf STDERR "%02d:%02d:%02d",$hour,$min,$sec;
		print STDERR " > " . $post ."\n";
		
		my $maptab    = $exp->rawmap; # the SQL table
		my $splicetab = $exp->splicedmap;
		my $readtab   = $exp->reads;
		my $sqlstack = [
			qq(ALTER IGNORE TABLE \`$maptab\` ADD INDEX (readid);),
			qq(ALTER IGNORE TABLE \`$splicetab\` ADD INDEX (readid);),
			qq(ALTER IGNORE TABLE \`$readtab\` ADD INDEX (readid);)
		];
		$exp->wait4database;
		foreach my $i (@{$sqlstack}) {
			print STDERR "## SQL ## ". substr($i,0,120).((length($i) > 120) ? " ...\n" : "\n");
			my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
			$dbh->prepare( $i )->execute() or die "Couldn't execute statement!";
			$dbh->disconnect();
		}
		return;
}

#NG
sub loadmappings {
	my $exp = $_[0];
	my $noindex = $_[1]; # SWITCH disables indexing (if tables will not be used with readid)
	my $sqldb   = $exp->get_db;
	my $sqlhost = $exp->get_sqlhost;
	my $sqluser = $exp->get_sqluser;
	my $sqlport = $exp->get_sqlport;
	my $localexecute = "mysql -u $sqluser -h $sqlhost $sqldb";
	
	my $m = $exp->mapping;
	my $s = $exp->splicedbwtout;
	
	######## GENOME
	if (-s $m) {
		my $maptab = $exp->rawmap; # the SQL table
		my $sqlstack = [
			qq(DROP TABLE if EXISTS \`$maptab\`;),
			qq(CREATE TABLE IF NOT EXISTS \`$maptab\` (id INT NOT NULL AUTO_INCREMENT, readid char(64) NOT NULL, strand char(1) NOT NULL, chr char(32) NOT NULL, pos INT NOT NULL, reserved INT, mismatch TEXT, PRIMARY KEY (id));),
			qq(LOAD DATA LOCAL INFILE \'$m\' IGNORE INTO TABLE \`$maptab\` (readid,strand,chr,pos,reserved,mismatch);),
			qq(ALTER TABLE \`$maptab\` ADD INDEX idx_readid (readid);)
		];
		$exp->wait4database;
		foreach my $i (@{$sqlstack}) {
			if ($noindex && $i =~ /ADD INDEX/) { next }
			print STDERR "## SQL ## ". substr($i,0,120)." ...\n";
			my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
			$dbh->prepare( $i )->execute() or die "Couldn't execute statement!";
			$dbh->disconnect();
		}
	}
	else { warn "WARNING: Genome Mapping not present!\n" }
	######## SPLICE
	if (-s $s) {
		my $splicetab = $exp->splicedmap;
		my $sqlstack = [
			qq(DROP TABLE if EXISTS \`$splicetab\`;),
			qq(CREATE TABLE IF NOT EXISTS \`$splicetab\` (id INT NOT NULL AUTO_INCREMENT, readid char(64) NOT NULL, strand char(1) NOT NULL, splice char(18) NOT NULL, pos INT NOT NULL, reserved INT, mismatch TEXT, PRIMARY KEY (id));),
			qq(LOAD DATA LOCAL INFILE \'$s\' IGNORE INTO TABLE \`$splicetab\` (readid,strand,splice,pos,reserved,mismatch);),
			qq(ALTER TABLE \`$splicetab\` ADD INDEX idx_readid (readid);)
		];
#		$exp->wait4database;
		foreach my $i (@{$sqlstack}) {
			if ($noindex && $i =~ /ADD INDEX/) { next }
			print STDERR "## SQL ## ". substr($i,0,120)." ...\n";
			my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
			$dbh->prepare( $i )->execute() or die "Couldn't execute statement!";
			$dbh->disconnect();
		}

	}
	else { warn "WARNING: Splice Mapping not present!\n" }
	$exp->updateMapCount;
	$exp->updateReadCount;
	return;
}

#NG#
sub exp_done {
	# checks if allmap is ready
	my $exp = $_[0];
	my $allmap = $exp->allmap;
	my $sqldb    = $exp->get_db;
	my $sqlhost  = $exp->get_sqlhost;
	my $sqlport  = $exp->get_sqlport;
	my $sqluser  = $exp->get_sqluser;
	#print STDERR " -> looking for $allmap\n";
	#unless (-s $allmap) { return 0 }
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	return 0 unless ($dbh->prepare( "SHOW TABLES LIKE \"$allmap\"" )->execute() > 0);
	my $sql = qq(SELECT COUNT(*) FROM $allmap;);
	#print STDERR "## SQL ## $sql \n";
	my $sth = $dbh->prepare( $sql );
	my $rc = $sth->execute();
	if ($rc) {
		my ( $aa );
		$sth->bind_columns( \$aa );
		while( $sth->fetch() ) { print STDERR "\nFound $aa mappings from previous run\n" }
		$dbh->disconnect();
		return $aa;
	}
	$dbh->disconnect();
	return $rc;
}

#NG#
sub loadseq {
	my $exp     = $_[0];
	my $pre     = $exp->get_exid;
	my $sqldb   = $exp->get_db;
	my $readfile = $exp->get_read;
	my $readtab  = $exp->reads;
	my $sqlhost = $exp->get_sqlhost;
	my $sqlport = $exp->get_sqlport;
	my $sqluser = $exp->get_sqluser;

	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 

	my $create_sql = [	qq( DROP TABLE IF EXISTS \`$readtab\`; ),
						qq( CREATE TABLE IF NOT EXISTS \`$readtab\` (id INT AUTO_INCREMENT, readid char(64) NOT NULL, PRIMARY KEY (id)); ) ];

	foreach my $query (@{$create_sql}) {
		print STDERR "## SQL ## ". substr($query,0,120)." ...\n";
		$dbh->prepare( $query )->execute();
	}
	$dbh->disconnect();
	
	# load data
	my $rifi = time() . "_" . int(rand(10000)). "_reads";
	
	# conversion to lineperline format (MySQL suitable)
	# safe and slow procedure (replaces former sed/perl oneliners)
	print STDERR "\nStoring sequences in central database";
	my $seqwrt = 1;
	my $line = 0;
	my $length = 1;
	open(RIN, "<$readfile") or die $readfile;
	open(ROT, ">$rifi") or die $rifi;
	while (<RIN>) {
		my $chp = $_;
		chomp($chp);
		if    (($chp =~ /^(A|T|G|C|N)+$/) && ($seqwrt == 2))   { die "2|$chp" unless ($seqwrt == 2); $seqwrt++; $length = length($chp) }
		elsif (($chp =~ /^[^\s]{$length}$/) && ($seqwrt == 4)) { die "4|$chp" unless ($seqwrt == 4); $seqwrt = 1; $length = 1 }
		elsif ($chp =~ /^@/)               { die "1|$chp" unless ($seqwrt == 1); $chp =~ s/^@//; print ROT "$chp\n"; $seqwrt++ }
		elsif ($chp =~ /^\+/)              { die "3|$chp" unless ($seqwrt == 3); $seqwrt++ }
		elsif ($chp =~ /^\s*$/)            { next } # empty line
		else                               { die ">$chp<$line>$readfile<" } # something else
		if ((++$line % 10000) == 0) { print STDERR "\rStoring sequences in central database\t". ($line / 4) };
	}
	print STDERR "\n";
	close(ROT);
	close(RIN);
	
	$dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	
	my $load_sql = [qq(LOAD DATA LOCAL INFILE \'$rifi\' IGNORE INTO TABLE \`$readtab\` (readid)),
					qq(ALTER TABLE \`$readtab\` ADD UNIQUE INDEX (readid);),
					];
					
	$exp->wait4database;
	foreach my $query (@{$load_sql}) {
		print STDERR "## SQL ## ". substr($query,0,120)." ...\n";
		$dbh->prepare( $query )->execute();
	}
	$dbh->disconnect();
	$exp->updateReadCount;
	die if (system("rm -f $rifi"));
	return $readtab;
}

###################################
## INTERNALS ######################
###################################

#NG
# resolves a set of potential fusion exons
sub _resolve {
	#order exons by start
	my %o;
	foreach my $a (@{$_[0]}) {
		unless ($o{$a->[1]} && $o{$a->[1]}->[0] > $a->[2]) {
			$o{$a->[1]} = [$a->[2], $a->[0]];
		}
	}
	# resolve
	my @fusion;
	my $end = -1;
	foreach my $a (sort {$a <=> $b} keys %o) {
		if ($end > $a) { next }
		$end = ($o{$a}->[0] > $end) ? $o{$a}->[0] : $end;
		while (my $newend = __another($end, \%o)) {
			$end = $newend;
		}
		unless ($end == $o{$a}->[0]) {
			push @fusion, [$o{$a}->[1]."+", $a, $end]
		}
	}
	return @fusion;	
}

#NG
sub __another {
	foreach my $e (sort {$a <=> $b} keys %{$_[1]}) {
		if (($e <= $_[0]+1) && ($_[0]+1 <= $_[1]->{$e}->[0])) {
			return $_[1]->{$e}->[0];
		}
	}
	return 0;
}

#NG#
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
	# max is not always best (allow for max and max - 1 if max-1 has more meber than max)
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

#NG#
sub _revcom {
	my $dna = shift;
	my $revcomp = reverse($dna);
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	return $revcomp;
}


#REX# chooseforamt is based on file size rather than on sequence count
sub _chooseformat {
	my $limit = 104857600; # 100 Mbytes
	my $format = 'fasta';
	# seq size
	if (-s $_[0] >= $limit) {
		$format = 'largefasta';
	}
	print STDERR "\rNOTE: Using $format for sequence extraction\n";
	return $format;
}

###################################
## LEGACY CODE ####################
###################################
sub LEGACY_chooseformat {
	my $limit = 100000;
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

1;

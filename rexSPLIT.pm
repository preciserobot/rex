package rexSPLIT;

use strict;

sub mapfilter {
	my $exp     = $_[0];
	my $local   = $_[1];
	my $pairs   = (-s $exp->matepairs) ? $exp->matepairs : _getPairs($exp);
	my $gmap    = $exp->rawmap;
	my $smap    = $exp->splicedmap;
	my $ciga    = $exp->annot("splicigars");
	my $introns = (-s $exp->annot("junctions")) ? $exp->annot("junctions") : _getIntrons($exp);
	my $readlen = $exp->get_readl;
	my $mindist = $exp->get_mindist;
	my $maxdist = $exp->get_maxdist;
	my $gmap_out = $exp->smackout($gmap);
	my $smap_out = $exp->smackout($smap);
	# remote
	my $remote   = rexRemote->new();
	my $smackbin = $remote->get_smack_bin;

	my $io = "$pairs $gmap $smap $ciga $introns $readlen $mindist $maxdist $gmap_out $smap_out";

	if ($local) {
		my $exec = "smack-ms " . $io;
		print STDERR " -> Running SMACK locally! ($exec)\n";
		my $return = system($exec);
		if ($return) { die "FATAL: SMACK! returned an error code ($return)\n" }
	}
	else {
		my $exec = $remote->get_smack_bin . " " . $io;
		if ((-s $pairs) && (-s $gmap) && (-s $smap) && (-s $ciga) && (-s $introns)) {
			# run smack
			print STDERR " -> Running SMACK remotely! ($exec)\n";
			# compile
			$remote->sendCompile("smack");
			# send
			$remote->sendFile($pairs);
			$remote->sendFile($gmap);
			$remote->sendFile($smap);
			$remote->sendFile($ciga);
			$remote->sendFile($introns);
			# run 
			$remote->execute($exec);
			# fetch
			$remote->fetchFile($gmap_out);
			$remote->fetchFile($smap_out);
			# cleanup
			$remote->cleanup();
		}
		else { die "cannot find required files for SMACK call ($exec)" }
	}
	return;
}

###############
## INTERNALS ##
###############

sub _getIntrons {
	my $e = $_[0];
	my $i = $e->annot("junctions");
	my $d = $e->get_db;
	my $h = $e->get_sqlhost;
	my $u = $e->get_sqluser;
	print STDERR " -> Getting introns...";
	my $getintrons = qq{SELECT seq_region_id, start, end from $i};
	die if (system("mysql --skip-column-names --quick -h $h -u $u $d --execute=\"$getintrons\" > $i"));
	print STDERR "\n";
	return $i;
}

sub _getPairs {
	my $e = $_[0];
	my $r = $e->reads;
	my $o = $e->matepairs;
	my $d = $e->get_db;
	my $h = $e->get_sqlhost;
	my $p = $e->get_sqlport;
	my $u = $e->get_sqluser;
	my ( $a, $b );
	my %pairs;
	my $lico;
	print STDERR " -> Getting mate pair IDs...\n";
	my $getreads = qq{SELECT id, SUBSTRING_INDEX(readid,'/', 1) FROM $r;};
	my $mysqlquery = "mysql --skip-column-names --quick -h $h -u $u $d --execute=\"$getreads\"";
	open(POUT, ">$o") or die;
	open(SQL, "$mysqlquery |") or die;
	while( <SQL> ) {
		print STDERR "\r    $lico" if (++$lico % 100000 == 0);
		if (/^(\d+)\s+(\S+)/) {
			if ($pairs{$2}) {
				print POUT $pairs{$2} . "\t" . $1 . "\n";
				delete $pairs{$2};
			}
			else {
				$pairs{$2} = $1;
			}
		}
	}
	# set unpaired reads to mate with read 0
	foreach my $r (keys %pairs) {
		print POUT $pairs{$r} . "\t" . 0 . "\n";
	}
	print STDERR "\r    $lico";
	close(SQL);
	close(POUT);
	# load into sql
	my @sqlstack = (
		qq{DROP TABLE IF EXISTS $o; },
		qq{CREATE TABLE $o (a INT, b INT); },
		qq{LOAD DATA LOCAL INFILE \"$o\" INTO TABLE $o; }
	);
	my $dbh = DBI->connect("DBI:mysql:$d:$h:$p","$u","") or die "Database connection not made: $DBI::errstr"; 
	foreach my $q (@sqlstack) {
		print STDERR "\.";
		$dbh->prepare( $q )->execute;
	}
	$dbh->disconnect();
	print STDERR "\n";
	return $o;
}


1;
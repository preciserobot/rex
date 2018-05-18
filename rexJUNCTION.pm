#  Package for managing splice reconstruction and read mapping to those
# (c) David Brawand. All right reserved. Absolutely no warranty!

package rexJUNCTION;

use strict;

use Usu;
use DBI;
use DBD::mysql;
use Digest::MD5 qw(md5_hex);

sub generateMultisplice {
	my $exp = $_[0];
	my $stringent = $exp->get_estr;
	my $ovrh    = $exp->get_ovrh;
	my $readl   = $exp->get_readl;
	my $sqlhost = $exp->get_sqlhost;
	my $sqluser = $exp->get_sqluser;
	my $sqlport = $exp->get_sqlport;
	my $sqldb   = $exp->get_db;
	my $localextract = "mysql --skip-column-names -u $sqluser -h $sqlhost $sqldb";
	# generate cigars
	my $splicigars = $exp->annot("splicigars");
	my $junctions  = $exp->annot("junctions");
	my $mintron    = $exp->get_mintron;
	my $tmpin  = 'splicigarinput_' . md5_hex(rand());
	my $tmpout = 'splicigaroutput_' . md5_hex(rand());
	
	my $tmpin2;

	if ($stringent) {
		my $efile = $exp->rfexon;
		$tmpin2  = 'splicigarinput2_' . md5_hex(rand());
		my $exonquery = qq{ SELECT seq_region, seq_region_start, seq_region_end, seq_region_strand FROM $efile } ;
		die if system("echo \"$exonquery\" | $localextract > $tmpin2");
	}
	else {
		$tmpin2 = '-';
	}
	my $sel = qq( SELECT DISTINCT id, seq_region_id, strand, start, end from $junctions ); ############## SELECT HERE
	die if system("echo \"$sel\" | $localextract > $tmpin");
	my $exec = "splicigar $tmpin $tmpin2 $ovrh $readl $mintron > $tmpout";
	print STDERR "Running SPLICIGAR: $exec\n";
	die if system($exec);
	# LOAD TO MYSQL
	my @s;
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "DB connection error: $DBI::errstr"; 
	$s[0] = qq(DROP TABLE IF EXISTS $splicigars;);

	# create intron table
	$s[1] = qq(CREATE TABLE $splicigars (id INT, centersplice int, seq_region_id INT NOT NULL, strand INT NOT NULL, start INT NOT NULL, end INT NOT NULL, splice_id INT, INDEX idx_region (seq_region_id)););
	$s[2] = qq(LOAD DATA LOCAL INFILE \"$tmpout\" INTO TABLE $splicigars);
	for(my $i=0; $i < scalar(@s); $i++) { print STDERR "## SQL ## ". substr($s[$i],0,128) . " ...\n"; $dbh->prepare( $s[$i] )->execute(); }
	system("mv $tmpout $splicigars");
	system("rm -f $tmpin");
	if ($stringent) { system("rm -f $tmpin2") }
	return $splicigars;
}

sub buildMultispliceSeqs {
	my $exp = $_[0];
	my $splicigars   = $exp->annot("splicigars");
	my $splicigarseq = $exp->annot("splicigarseq");
	my $splicigarseq_tmp = $splicigarseq . "\.unsorted";
	my $sr           = $exp->annot("seq_region");
	my $db           = $exp->get_refs;
	my $sqlhost = $exp->get_sqlhost;
	my $sqlport = $exp->get_sqlport;
	my $sqluser = $exp->get_sqluser;
	my $sqldb   = $exp->get_db;
	my $splicelength = $exp->splicelength;
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "DB connection error: $DBI::errstr"; 
	my $fetch = qq( SELECT s.id, r.name, s.strand, s.start, s.end from $splicigars as s, $sr as r where s.seq_region_id = r.seq_region_id);
	my $sth = $dbh->prepare( $fetch );
	my( $aa, $bb, $cc, $dd, $ee, $seqcounter );
	$sth->execute();
	$sth->bind_columns( \$aa, \$bb, \$cc, \$dd, \$ee );
	my %slice;
	while( $sth->fetch() ) {
		$slice{$bb}{$cc}{$aa}{$dd} = $ee; # name, strand, id, start, end
	}
	print STDERR "\nGenerating Splice database...\n";
	open(SPDB, ">$splicigarseq_tmp") or die "cannot write to $splicigarseq_tmp";
	my $stream  = Bio::SeqIO->new(-file => $db ,-format => rexDB::_chooseformat($db));
	while ( my $seq = $stream->next_seq() ) {
		if ($slice{$seq->primary_id}) {
			foreach my $strand (keys %{$slice{$seq->primary_id}}) {
				foreach my $id (keys %{$slice{$seq->primary_id}{$strand}}) {
					print STDERR "\r\t" . $seq->primary_id . " " . ++$seqcounter . "     ";
					print SPDB "$id\t";
					my $sequence;
					my $donors;
					foreach my $start ( sort { $a <=> $b } keys %{$slice{$seq->primary_id}{$strand}{$id}} ) {
						
						# redifines start if negative (and adds Ns)
						my $st;
						my $end = $slice{$seq->primary_id}{$strand}{$id}{$start};
						if ($start < 1) {
							#warn "WARNING: Annotation close to region start ($start < ". 1 .") -> completing with Ns";
							$sequence .= "N" x (abs($start) + 1);
							$st = 1;
						}
						else {
							$st = $start;
						}
						
						# checks if sequence has to be completed with Ns
						if ($end > $seq->length) {
							#warn "WARNING: Annotation close to region end (". $seq->length ." < $end) -> completing with Ns";
							$sequence .= $seq->subseq($st, $seq->length);
							$sequence .= "N" x ($end - $seq->length);
						}
						else {
							$sequence .= $seq->subseq($st, $end);
						}
						#$donors .= '...';
						#$donors .= ($strand < 0) ? _revcom($seq->subseq($st - 2, $st - 1)) : $seq->subseq($st - 2, $st - 1);
						#$donors .= '[]';
						#$donors .= ($strand < 0) ? _revcom($seq->subseq($end + 1, $end + 2)) : $seq->subseq($end + 1, $end + 2)
					}
					if (length($sequence) != $splicelength) {
						die "FATAL: Spliceseq length not consistent with readlength and overlap (".length($sequence).")";
					}
					print SPDB $sequence . "\n";
					#print STDERR "$id\t$donors\n";
				}
			}
		}
	}
	close(SPDB);
	$dbh->disconnect();
	# sort by id and write
	open(SPDB, ">$splicigarseq") or die "cannot write to $splicigarseq";
	open(SORTED, "sort -k1n $splicigarseq_tmp |") or die "cannot sort";
	while(<SORTED>) {
		if (/^(\d+)\s+(\S+)/) {
			print SPDB ">" . $1 . "\n" . $2 . "\n";
		}
		else {die}
	}
	close(SORTED);
	close(SPDB);
	print STDERR "\nBuilding Splice site library\n";
	die if system("bowtie-build -f $splicigarseq $splicigarseq > /dev/null");	
	return $splicigarseq;
}

sub spliceRedu {
	my $exp        = $_[0];
	my $splicigars = $exp->annot("splicigars");
	# SQL
	my $sqldb     = $exp->get_db;
	my $sqlhost   = $exp->get_sqlhost;
	my $sqlport   = $exp->get_sqlport;
	my $sqluser   = $exp->get_sqluser;
	# redundancy level
	my $redu = 1;
	my $rr;
	# calculate maximal redundancy
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my $sth = $dbh->prepare( qq{ SELECT count(splice_id) FROM $splicigars WHERE splice_id != 0 GROUP BY splice_id ORDER BY count(splice_id) DESC LIMIT 1; } );
	$sth->execute();
	$sth->bind_columns( \$rr );
	while( $sth->fetch() ) { $redu = $rr }
	$dbh->disconnect();
	return $redu;
}

#######################################
## HELPERS ############################
#######################################
sub _revcom {
	my $dna = shift;
	my $revcomp = reverse($dna);
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	return $revcomp;
}
#######################################
## LEGACY AND TEMPLATES ###############
#######################################

1;

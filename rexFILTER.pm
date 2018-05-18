package rexFILTER;

use strict;
use Usu;
use Digest::MD5 qw(md5_hex);

sub run {
	my $e = $_[0];
	#check read length
	readlchk($e);

#	# MAKE BACKUP
#	print STDERR " -> Making backup of sequences in file ".$e->get_read."\n";
#	my $r   = $e->get_read;
#	my $unfi = $r . '.unfiltered'; # unfiltered reads;
#	unless (-e $unfi.'.bz2') { system("cp $r $unfi && pbzip2 $unfi") }

	# DROP FILTER (BY TILE)
	if ($e->get_infodrop < $e->get_readl) {
		my $rj = dropFilter($e);
		print STDERR "\tRejected ".scalar(@{$rj})." tiles ";
		foreach my $i (@{$rj}) { print STDERR "($i)" }
		print STDERR "\n";
	}
	# INFOBASE FILTER
	if ($e->get_infobase) {
		printf STDERR "\tAccepted %.2f percent of reads\n", filterInfobase($e)*100;
	}
	#GARBAGE FILTER
	unless ($e->get_paired) { filterGarbage($e) }; # filters reads over 0.9 of (A|T|G|C) or over 0.1 N

	idconvert($e); # coverts read file to uniqe id prefixes
	
	return;
}

####################################################################################################
## FILTER FUNCTIONS ################################################################################
####################################################################################################
sub idconvert {
	# check if head and tail identifiers are identical
	my ($pre,$mid);
	open(HEAD, "head -n6 " . $_[0]->get_read . " |") or die "\nERROR: Read file (".$_[0]->get_read.") not found\n";
	while (<HEAD>) {
		if (/^@([^:]+:\d:)(\d{1,3})(:-?\d{1,6}:-?\d{1,6})/) {
		#if (/^@([^:]+:\d:)(\d{1,3})(:-?\d{1,6}:-?\d{1,6})$/) {
			$pre = $1;
			$mid = $2;
			last;
		}
	}
	close(HEAD);
	unless ($pre && $mid) { die "Cannot parse read identifier" }
	
	my $rewrite = 0;
	open(TAIL, "tail -n6 " . $_[0]->get_read . " |") or die "\nERROR: Read file (".$_[0]->get_read.") not found\n";
	while (<TAIL>) {
		if (/^@([^:]+:\d:)(\d{1,3})(:-?\d{1,6}:-?\d{1,6})/) {
		#if (/^@([^:]+:\d:)(\d{1,3})(:-?\d{1,6}:-?\d{1,6})$/) {
			if ($pre ne $1) {
				++$rewrite;
			}
			last;
		}
	}
	close(TAIL);
	
	# rewrite ids to have cnsecutive tile numbers and same baseID
	if ($rewrite) {
	    my $maxtile = 120;
	    my $seqcounter = 0;
	    my $nonext = 1;
	    my $secondfile = 0;
	    my $previoustile = 0;
	    my $previousid;
	    print STDERR "\r -> UNIFYING IDs    ";
	    my $infile  = $_[0]->get_read;
	    my $outfile = $_[0]->get_read . '.rewrite';
	    open(FOUT, ">$outfile") or die "cannot write to $outfile";
	    open(FIN,  "<$infile") or die "\nERROR: Read file ($infile) not found\n";
	    while (<FIN>) {
		if (/^([@|\+])([^:]+:\d:)(\d{1,3})(:-?\d{1,6}:-?\d{1,6})$/ && !($nonext)) {
		    $nonext = 1;
		    my $tile = $3;
		    if ($pre ne $2) { # not same as first
			if ($previousid ne $2) {
			    $secondfile++;
			}
			$tile = $secondfile * $maxtile + $3;
		    }
		    $previousid = $2;
		    # print
		    print FOUT $1 . $pre . $tile . $4 . "\n";


		    # check tile ordering
		    if ($tile < $previoustile) { die "Reads are not ordered by tile ($previoustile => $tile)!" }
		    else                       { $previoustile = $tile }


		    if (++$seqcounter % 1000 == 0)      { print STDERR "\r -> UNIFYING IDs    $tile $seqcounter" }
		}
		elsif (/^\+$/) {
		    print FOUT;
		    $nonext = 1;
		}
		else {
		    print FOUT;
		    $nonext = 0;
		}
	    }
	    close(FIN);
	    close(FOUT);
		system("mv -f $outfile $infile");
		print STDERR "\n";
	}
	return;
}


sub readlchk {
	my $checked;
	open(FIN, $_[0]->get_read) or die "\nERROR: Read file (".$_[0]->get_read.") not found\n";
	while (<FIN>) {
		if (/^((G|C|A|T|N)+)$/) {
			if (length($1) != $_[0]->get_readl) {
				die "\nFATAL ERROR: Read length error (".length($1)." instead of ".$_[0]->get_readl.")\n";
			}
			elsif (++$checked > 1000) {
				close(FIN);
				return;
			}
		}
	}
	close(FIN);
	return;
}
sub filterGarbage {
	my $r   = $_[0]->get_read;
	my $fi = $r . '.garbagefiltered'; # unfiltered reads;
	print STDERR "\r -> GARBAGEFILTER   ";
	system("filter_garbage -q < $r > $fi && mv -f $fi $r");
	return 1;
}

sub dropFilter {
	my $e = $_[0];
	my $infile = $e->get_read;
	open(RIN, "<$infile") or die;
	my $outfile = $infile . '.dropfilter';
	open(ROUT, ">$outfile") or die;
	my $seqwrt = 1;
	my $length = $e->get_readl;
	my $seq = [];
	my $t = 0; #the current tile
	my (@rejected,@tilestack,@tilestat);
	my $seqcounter = 0;
	my $line = 0;
	while (<RIN>) {
		++$line;
		my $chp = $_;
		chomp($chp);
		if    (($chp =~ /^(A|T|G|C|N){$length}$/) && ($seqwrt == 2))   { die "2|$chp" unless ($seqwrt == 2); $seq->[1] = $chp; $seqwrt++;   } #SEQ
		elsif (($chp =~ /^[^\s]{$length}$/) && ($seqwrt == 4))         { die "4|$chp" unless ($seqwrt == 4); $seq->[3] = $chp; $seqwrt=1;
			# add to tile stack
			push @tilestack, $seq;
			push @tilestat, c_dp($chp, $e->get_infoqual, ($e->get_bwts) ? 64 : 33);
			$seq = [];
			
		} #QUAL
		elsif ($chp =~ /^@/)                                           { die "1|$chp" unless ($seqwrt == 1); $seq->[0] = $chp; $seqwrt++;
			if ($chp =~ /^@([^:]+:[^:]+:)([^:]+):/) {
				if ($t != $2) {
					$t = $2;
					# print if accepted
					die unless (scalar(@tilestack) == scalar(@tilestat));
					if (scalar(@tilestat)) {
						my $m = _median(\@tilestat);
						if ($m >= $e->get_infodrop) {
							_fiPrint(*ROUT, \@tilestack) ;
						}
						else {
							push @rejected, $1.$2."($m)";
						}
					}
					# reset
					@tilestack = ();
					@tilestat = ();
				}
				if (++$seqcounter % 1000 == 0) {
					print STDERR "\r -> TILEFILTER      $t ".$seqcounter;
				}
			}
		} #ID
		elsif ($chp =~ /^\+/)                                          { die "3|$chp" unless ($seqwrt == 3); $seq->[2] = $chp; $seqwrt++;   } #ID+
		elsif ($chp =~ /^\s*$/)                                        { next } # empty line
		else                                                           { die ">$chp<$line>$infile<" } # something else
		
	}
	close(RIN);
	# do last
	die unless (scalar(@tilestack) == scalar(@tilestat));
	if (_median(\@tilestat) >= $e->get_infodrop) {
		_fiPrint(*ROUT, \@tilestack) ;
	}
	else {
		push @rejected, $1.$2.'('._median(\@tilestat).')';
	}
	close(ROUT);
	die if system("mv -f $outfile $infile");
	return \@rejected;
}

sub filterInfobase {
	my $e = $_[0];
	my $infile = $e->get_read;
	open(RIN, "<$infile") or die;
	my $outfile = $infile . '.infobasefilter';
	open(ROUT, ">$outfile") or die;
	my $seqwrt = 1;
	my $length = $e->get_readl;
	my $seq = [];
	my ($rejected, $seqcounter, $line) = (0, 0, 0);
	while (<RIN>) {
		++$line;
		my $chp = $_;
		chomp($chp);
		if    (($chp =~ /^(A|T|G|C|N){$length}$/) && ($seqwrt == 2))   { die "2|$chp" unless ($seqwrt == 2); $seq->[1] = $chp; $seqwrt++;   } #SEQ
		elsif (($chp =~ /^[^\s]{$length}$/) && ($seqwrt == 4))         { die "4|$chp" unless ($seqwrt == 4); $seq->[3] = $chp; $seqwrt=1;
			# add to tile stack
			if (c_dp($chp, $e->get_infoqual, ($e->get_bwts) ? 64 : 33) >= $e->get_infobase) {
				_fiPrint(*ROUT, [ $seq ]) ;
			}
			else {
				$rejected++;
			}
			$seq = [];
		} #QUAL
		elsif ($chp =~ /^@/)       { die "1|$chp" unless ($seqwrt == 1); $seq->[0] = $chp; $seqwrt++; print STDERR "\r -> INFOBASEFILTER  $seqcounter $rejected" if (++$seqcounter % 1000 == 0) } #ID
		elsif ($chp =~ /^\+/)      { die "3|$chp" unless ($seqwrt == 3); $seq->[2] = $chp; $seqwrt++;   } #ID+
		elsif ($chp =~ /^\s*$/)    { next } # empty line
		else                       { die ">$chp<$line>$infile<" } # something else
		
	}
	close(RIN);
	close(ROUT);
	die if system("mv -f $outfile $infile");
	return 1-($rejected/$seqcounter);
}

####################################################################################################
### C FUNCTIONS ####################################################################################
####################################################################################################
use Inline C => <<'END_C';

int c_dp(char* s, int ib, int phred) {
	int droppos = strlen(s);
	int i;
	for(i = 0; i < strlen(s); i++) {
		if ((int) s[i] - phred < ib) {
			if (droppos == strlen(s)) {
				droppos = i;
			}
		}
		else {
			droppos = strlen(s);
		}
	}
	return droppos;
}

END_C

####################################################################################################
### INTERNALS ######################################################################################
####################################################################################################

sub _fiPrint {
	my $fh = $_[0];
	foreach my $s (@{$_[1]}) {
		print $fh $s->[0] . "\n";
		print $fh $s->[1] . "\n";
		print $fh $s->[2] . "\n";
		print $fh $s->[3] . "\n";
	}
	return;
}

sub _median {
	return "" unless ($_[0] && ($_[0] =~ /ARRAY/));
	my $size = scalar(@{$_[0]});
	my @Array = sort {$a <=> $b} @{$_[0]};
	my $median;
	if (($size % 2) == 0) {
		#even -> average of middle 2
		$median = $Array[int($#Array / 2)] + (($Array[int($#Array / 2) + 1] - $Array[int($#Array / 2)]) / 2);
	}
	else {
		# odd -> center value
		$median = $Array[($#Array / 2)];
	}
	die unless defined $median;
	return $median;
}

1;

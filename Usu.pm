package Usu;

use strict;

use Getopt::Long;
use File::stat;
use Time::localtime;
use Digest::MD5 qw(md5 md5_hex md5_base64);

sub header {
    my $date_string = ctime(stat($_[0])->mtime);
	my $prog = $_[1];
	my $desc = $_[2];
	my $cprt = '(c) '. $_[3] .' (' . $date_string . ')';
	my $sp   = 3;
	my $maxl = _longest($prog, @{$desc}, $cprt);

	print STDERR "\n\n";
	print STDERR "\t+"  . '=' x ($sp+$maxl+$sp) . "+\n";
	print STDERR "\t\|" . ' ' x ($sp+$maxl+$sp) . "\|\n";
	print STDERR "\t\|" . (' ' x $sp) . $prog . (' ' x ($maxl - (length $prog) + $sp)) . "\|\n";
	print STDERR "\t\|" . ' ' x ($sp+$maxl+$sp) . "\|\n";
	print STDERR "\t+"  . '=' x ($sp+$maxl+$sp) . "+\n";
	print STDERR "\t\|" . ' ' x ($sp+$maxl+$sp) . "\|\n";
	for(my $i=0; $i < scalar @{$desc}; $i++) {
		print STDERR "\t\|" . (' ' x $sp) . $desc->[$i] . (' ' x ($maxl - (length $desc->[$i]) + $sp)) . "\|\n";
	}
	print STDERR "\t\|" . ' ' x ($sp+$maxl+$sp) . "\|\n";
	print STDERR "\t\|" . (' ' x $sp) . $cprt . (' ' x ($maxl - (length $cprt) + $sp)) . "\|\n";
	print STDERR "\t\|" . ' ' x ($sp+$maxl+$sp) . "\|\n";
	print STDERR "\t+"  . '-' x ($sp+$maxl+$sp) . "+\n";

	return 1;
}

sub legacy {
	my $d = $_[0] ? $_[0] : 0;
	warn "\n>>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<<\n";
	warn "\n>>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<<\n";
	warn "\n>>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<< >>>LEGACY<<<\n";
	sleep $d;
	return;
}

sub default {
	my $d = $_[2] ? $_[2] : 0;
	warn "\n>>>DEFAULT<<< >>>DEFAULT<<< >>>DEFAULT<<< >>>DEFAULT<<< >>>DEFAULT<<< >>>DEFAULT<<< >>>DEFAULT<<<\n";
	print STDERR "\tUsing $_[1] as default value for $_[0]";
	warn "\n>>>DEFAULT<<< >>>DEFAULT<<< >>>DEFAULT<<< >>>DEFAULT<<< >>>DEFAULT<<< >>>DEFAULT<<< >>>DEFAULT<<<\n";
	sleep $d;
	return;
}

sub confirm {
	my $confirm = substr(md5_hex(rand), 0, 3);
	print STDERR $_[0] if ($_[0]);
	print STDERR "\ninput triplet ($confirm) to continue: ";
	my $in = <STDIN>;
	chomp($in);
	return 1 if ($in eq $confirm);
	print STDERR "\n-> NOT CONFIRMED BY USER\n\n";
	exit(0);
}

sub alert {
	for (my $i=0; $i < 3; $i++) {
		print STDERR "\a";
		sleep 1;
	}
	print STDERR $_[0];
	print STDERR "\n";
	return;
}

sub compress {
	unless (-e $_[1]) { return } # return if file does not exist
	my $prog = $_[0];
	print STDERR "\rFile " . $_[1] . " is being compressed by " . $_[0] . "...";
	die if system("$_[0] $_[1]");
	print STDERR "\rFile " . $_[1] . " has been compressed by " . $_[0] . "      \n";
	return;
}

sub countdownRemove {
	unless (-e $_[1]) { return } # return if file does not exist
	my $duration = $_[0];
	my @delfiles;
	if (-s $_[1])            { push @delfiles, $_[1] }
	if (-s ($_[1] . '.gz'))  { push @delfiles, ($_[1] . '.gz') }
	if (-s ($_[1] . '.bz2')) { push @delfiles, ($_[1] . '.bz2') }
	return unless $duration =~ /^\d+$/;
	return unless @delfiles;
	for(my $i=0; $i < $duration; $i++) {
		print STDERR "\rFile(s) " . (join ', ', @delfiles) . " will be removed in " . ($duration - $i) . " seconds";
		sleep 1;
	}
	foreach my $f (@delfiles) { system("rm -f $_[1]") }
	print STDERR "\rFile(s) " . (join ', ', @delfiles) . " have been deleted                   \n";
	return;
}

#NG
sub TimeAndSwitch {
	my $pre = $_[0];
	my $post = $_[1];
	my $selector = $_[2];
	my ($sec,$min,$hour,$mday,$mon,$yr,$wday,$yday,$isdst)=@{localtime(time)};
	print STDERR $pre." > ";
	printf STDERR "%02d:%02d:%02d",$hour,$min,$sec;
	print STDERR " > " . $post ."\n";
	unless (defined $selector) { return 1 } # returns 1 if no jumping demanded

	if ($selector =~ /,/) { 
		my @j = split /,/, $selector;
		my $match = 0;
		foreach my $s (@j) {
			if ($pre == $s) {
				$match++;
			}
		}
		if ($match == 0) {
			print STDERR "\t-> SKIPPED\n";
			return 0;
		}
	}
	
	elsif ($selector =~ /(\d+)-(\d+)/) {
		# step range
		if ($1 > $pre) {
			print STDERR "\t-> SKIPPED\n";
			return 0;
		}
		elsif ($pre > $2) {
			print STDERR "\t-> SKIPPED\n";
			return 0;
		}
	}
	elsif ($selector =~ /(\d+)\+/) {
		# step and following
		if ($1 > $pre) {
			print STDERR "\t-> SKIPPED\n";
			return 0;
		}
	}
	elsif ($selector =~ /((?:-)?\d+)/) {
		# just the step
		if ($1 != $pre) {
			print STDERR "\t-> SKIPPED\n";
			return 0;
		}
	}
	return 1;
}

sub pseudorandom {
	return md5_hex(time() + rand(100000000));
}

sub randomfile {
	my $prefix = "_tmp_";
	my $name   = $_[0] ? $_[0] : time();
	my $random = "_" . substr(md5_hex(rand(100)), 0, 8);
	my $raf = $prefix . $name . $random;
	system("touch $raf");
	return $raf; 
}

sub tag {
	my $t = $_[0];
	my $d = $_[1] ? $_[1] : 0;
	my $sp = 1;
	my $maxl = length($t);
	print STDERR "\n";
	print STDERR "\t ##"  . '#' x ($sp+$maxl+$sp) . "##\n";
	print STDERR "\t###" . (' ' x $sp) . $t . (' ' x $sp) . "###\n";
	print STDERR "\t ##"  . '#' x ($sp+$maxl+$sp) . "##\n";
	return;
}

sub countdown {
	my $duration = $_[0];
	return unless $duration =~ /^\d+$/;
	for(my $i=0; $i < $duration; $i++) {
		print STDERR ($duration - $i);
		sleep 1;
		print STDERR "\b" x length($duration - $i);
	}
	print STDERR "FINISHED";
	return;
}

sub usage {
	my $prog = $_[0];
	my $desc = $_[1];
	print STDERR "\n";
	print STDERR "\n\t$prog";
	print STDERR "\n";
	for(my $i=0; $i < scalar @{$desc}; $i++) {
		print STDERR "\t" . $desc->[$i] . "\n";
	}
	print STDERR "\n";
	return 1;
}

sub start {
	my ($sec,$min,$hour,$mday,$mon,$yr)=localtime(time);
	printf STDERR "\n\tPID %6d started at %4d-%02d-%02d %02d:%02d:%02d\n", $_[0],$yr+1900,$mon+1,$mday,$hour,$min,$sec;
	return 1;
}

sub _longest {
	my $maxl = 0;
	my @lines = @_;
	foreach my $s (@lines) {
		$maxl = ((length $s) > $maxl ) ? (length $s) : $maxl
	}
	return $maxl;
}

1;
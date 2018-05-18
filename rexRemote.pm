package rexRemote;

use strict;
use warnings;
our $AUTOLOAD;
use Carp;
use Digest::MD5 qw(md5_hex);

my $user = "kguschan";

# Class data and methods, referring to ALL created objects by this class
# Subroutines are made available to the entire package scope even if they are conditionally isolated or in a seperate scope
{
	my $classid = substr(md5_hex(time()),0,4);
	my %_attribute_properties = (
		#_server          => [ "dee-serv02.vital-it.ch",									'read.write' ],
		_server          => [ "cig-serv01.vital-it.ch",									'read.write' ], # open SSH forward first: ssh -f -l dbrawand -L 1234:dee-serv02.vital-it.ch:22 dev.vital-it.ch -N ### OBSOLETE
		_serverport      => [ "22",														'read.write' ],
		_serveruser      => [ $user,													'read.write' ],
#		_scratch         => [ "/home/cig/kaessmann/" . $user ."/",						'read.write' ],
#		_scratch         => [ "/scratch/frt/daily/" . $user ."/",						'read.write' ],
		_scratch         => [ "/scratch/local/weekly/" . $user ."/",					'read.write' ],
#		_scratch         => [ "/scratch/ul/monthly/" . $user ."/",						'read.write' ],
		_compileroptions => [ '-Wall -I/mnt/common/DevTools/src/boost/boost_1_46_1',	'read.write'],
		_directory       => [ "run_".$classid."/",										'read.write' ],
		_raparm_src      => [ '~/rex/bin/src/raparm-ms/*.?pp',							'read.write' ],
		_raparm_bin      => [ "raparm_" . $classid,										'read.write' ],
		_rmms_src        => [ '~/rex/bin/src/rmms/*.?pp',								'read.write' ],
		_rmms_bin        => [ "rmms_" . $classid,										'read.write' ],
		_ames_src        => [ '~/rex/bin/src/ames-ms/*.?pp',							'read.write' ],
		_ames_bin        => [ "ames_" . $classid,										'read.write' ],
		_smack_src       => [ '~/rex/bin/src/smack-ms/*.?pp',							'read.write' ],
		_smack_bin       => [ "smack_" . $classid,										'read.write' ],
	);

	# attribute management
	sub _all_attributes {
		keys %_attribute_properties;
	}
	
	sub _permissions {
		my($self, $attribute, $permissions) = @_;
		$_attribute_properties{$attribute}[1] =~ /$permissions/;
	}
	
	sub _attribute_default {
		my($self, $attribute) = @_;
		$_attribute_properties{$attribute}[0];
	}
	
	# object management
	my $_count = 0;
	sub get_count {
		$_count;
	}
	sub _incr_count {
		++$_count;
	}
	sub _decr_count {
		--$_count;
	}
}

# remote execution
sub sendFile    {
	my ($e, $arg, $force) = @_;
	my $port      = $e->{_serverport};
	my $user      = $e->{_serveruser};
	my $remote    = $e->{_server};
	my $remotedir = $e->{_scratch} . $e->{_directory};
	my $exec = "scp -c arcfour -P $port $arg $user".'@'."$remote:$remotedir";
	system($exec);
	return;
}
sub sendCompile {
	my ($e, $arg) = @_;
	my $remote    = $e->{_server};
	my $user      = $e->{_serveruser};
	my $port      = $e->{_serverport};
	my $remotedir = $e->{_scratch} . $e->{_directory};
	my $opts      = $e->{_compileroptions};
	my $success = 0;
	foreach my $key ($e->_all_attributes()) {
		if ($key =~ /$arg/ && $key =~ /src/) {
			my $src = $e->{$key};
			my $notokay = system("scp -P $port $src $user".'@'."$remote:$remotedir");
			$success++ unless $notokay;
		}
	}
	if ($success) {
		# find binary flag
		foreach my $key ($e->_all_attributes()) {
			if ($key =~ /$arg/ && $key =~ /bin/) {
				my $bin = $e->{$key};
				print STDERR "Compiling...";
				my $creturn = system("ssh -p $port $user".'@'."$remote \"g++ ".($opts ? $opts : '')." ".$remotedir."main.cpp -o ".$remotedir.$bin."\"");
				print STDERR ($creturn ? "FAILED!\n" : "SUCCESS!\n");
				return;
			}
		}
	}
	else {
		warn "There was a file transfer problem";
	}
	return;
}
sub run {
	my ($e, $arg) = @_;
	my $remote    = $e->{_server};
	my $port      = $e->{_serverport};
	my $user      = $e->{_serveruser};
	my $remotedir = $e->{_scratch} . $e->{_directory};

	my $script = "runme.sh";
	open(SCRIPT, ">$script");
	print SCRIPT "cd " . $remotedir . "\n";
	print SCRIPT $arg . "\n";
	close(SCRIPT);
	
	system("scp -P $port $script $user".'@'."$remote:$remotedir/$script");
	my $return = system("ssh -n -s $user".'@'."$remote \"$script\"");
	if ($return) { die "FATAL: returned an error code ($return)\n" }
	return;
}

sub execute     {
	my ($e, $arg) = @_;
	my $remote    = $e->{_server};
	my $port      = $e->{_serverport};
	my $user      = $e->{_serveruser};
	my $remotedir = $e->{_scratch} . $e->{_directory};
	my @exec = split /\s+/, $arg;
	my $execstring;
	foreach my $x (@exec) {
		if ($x =~ /^\d+$/ || $x =~ /^-$/) {
			$execstring .= $x . " ";
		}
		else {
			$execstring .= $remotedir . $x . " ";
		}
	}
	print STDERR "\n>>>\n$execstring\n<<<\n";
	
	my $return = system("ssh -p $port $user".'@'."$remote \"$execstring\"");
	if ($return) { die "FATAL: returned an error code ($return)\n" }
	return;
}
sub fetchFile   {
	my ($e, $arg) = @_;
	my $remote    = $e->{_server};
	my $port      = $e->{_serverport};
	my $user      = $e->{_serveruser};
	my $remotedir = $e->{_scratch} . $e->{_directory};
	system("scp -P $port $user".'@'.$remote.":".$remotedir.$arg." .");
	return;
}
sub cleanup     {
	my ($e, @arg) = @_;
	my $remote    = $e->{_server};
	my $port      = $e->{_serverport};
	my $user      = $e->{_serveruser};
	my $remotedir = $e->{_scratch} . $e->{_directory};
	if (@arg) {
		foreach my $file (@arg) {
			system("ssh -p $port $user".'@'."$remote \"rm -f $remotedir".$file."\"");
		}
	}
	else {
		system("ssh -p $port $user".'@'."$remote \"rm -rf $remotedir\"");
	}
	return;
}

################################################################################################################################################
## CONSTRUCTORS ################################################################################################################################
################################################################################################################################################
# Intelligent constructor
sub new {
	my ($class, %arg) = @_;
	my $self = bless { }, $class;
	
	foreach my $attribute ($self->_all_attributes()) {
		my($argument) = ($attribute =~ /^_(.*)/);
		if (exists $arg{$argument}) {
			$self->{$attribute} = $arg{$argument};
		}
		elsif ($self->_permissions($attribute, 'required')) {
			croak("Missing mandatory argument ($argument)");
		}
		else {
			$self->{$attribute} = $self->_attribute_default($attribute);
		}
	}
	$class ->_incr_count();
	# create directory for job
	my $remotedir = $self->{_scratch} . $self->{_directory};
	my $remotesrv = $self->{_server};
	my $port      = $self->{_serverport};
	my $user      = $self->{_serveruser};
	system("ssh -p $port $user".'@'."$remotesrv \"mkdir -p $remotedir\"");
	return $self;
}
sub load {
	my ($class, $arg) = @_;
	my $self = bless { }, $class;
	foreach my $attribute ($self->_all_attributes()) {
		my($argument) = ($attribute);
		if (exists $self->{$argument})                      { }# print STDERR "($argument|".$self->{$argument}.")\n" } # argument defined
		elsif ($self->_permissions($attribute, 'required')) { croak("Missing mandatory argument ($argument)") }
		else                                                { $self->{$attribute} = $self->_attribute_default($attribute) }
	}
	$class ->_incr_count();
	return $self;
}
sub AUTOLOAD {
	my ($self, $newvalue)  = @_;
	my ($op, $attrib) = ($AUTOLOAD =~/(get|set)(\_\w+)$/);
	unless ($op && $attrib) { croak "Method name $AUTOLOAD is not in the recongnized form (get|set)_attribute\n" }
	unless (exists $self->{$attrib}) { croak "No such attribute '$attrib' exists in the class ", ref($self) }
	no strict 'refs';
	if ($op eq 'get') {
		unless($self->_permissions($attrib, 'read')) { croak("ERROR: Read permission denied ($attrib)") }	
		*{$AUTOLOAD} = sub {
			my($self) = @_;
			unless($self->_permissions($attrib, 'read')) { croak("ERROR: Read permission denied ($attrib)") }
			$self->{$attrib}
		};

	}
	elsif ($op eq 'set') {
		unless($self->_permissions($attrib, 'write')) { croak("ERROR: Write permission denied ($attrib)") }	
		$self->{$attrib} = $newvalue;
		*{$AUTOLOAD} = sub {
			my($self, $newvalue) = @_;
			unless($self->_permissions($attrib, 'write')) { croak("ERROR: Write permission denied ($attrib)") }	
			$self->{$attrib} = $newvalue;
		};
	}
	use strict 'refs';
	return $self->{$attrib};	
}
sub DESTROY {
	my ($self) = @_;
	$self->_decr_count();
}

1;

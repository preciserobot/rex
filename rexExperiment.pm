package rexExperiment;

# a synteny slice (by orthologies in ensembl)

use strict;
use warnings;
our $AUTOLOAD;
use Carp;
use Digest::MD5 qw(md5_hex);
use rexRemote;


my $defaultlength = 75;


# Class data and methods, referring to ALL created objects by this class
# Subroutines are made available to the entire package scope even if they are conditionally isolated or in a seperate scope
{
	my %_attribute_properties = (
		_exid => [ '',			 					'read.write.required'	],	#the experiment ID (must be unique!)
		_spec => [ '',			 					'read.write'			],	#imaginary
		_gend => [ '',			 					'read.write'			],	#asexual
		_tiss => [ '',			 					'read.write'			],	#some
		_type => [ '',			 					'read.write'			],	#polya
		_desc => [ '',			 					'read.write'			],	#test dataset for tophat/bowtie
		_date => [ '',			 					'read.write'			],	#25.11.08/11:00
## BOWTIE PARAMETERS
		_bwtk => [ '10',		 					'read.write'			],	# How may allowed mappings
		_bwte => [ '70',		 					'read.write'			],	# max mismatch qual score
		_bwtm => [ '11',		 					'read.write'			],	# How may allowed mappings
		_bwtv => [ '-1',		 					'read.write'			],	# Max mismatch in splice DB remapping (to score splice sites) (negative values turn off / default was 3)
		_bwtb => [ '1',			 					'read.write'			],	# best
		_bwtl => [ '28',		 					'read.write'			],	# seed length
		_bwts => [ '0',			 					'read.write'			],	# solexa quality 1.3+ phred64 (otherwise phred33) 
# infobase filtering
		_infobase => [ '0',							'read.write'			],	# number of informative bases requested
		_infoqual => [ '5',							'read.write'			],	# minimum quality of informative base
		_infodrop => [ '50',						'read.write'			],	# quality drop must be after base pos (evaluated per tile!)
# splice mapping
		_ovrh =>    [ '1',		 					'read.write'			],	# the overhang required in splice site remapping (readl-3 + readl-3)
		_mintron => [ '40',		 					'read.write'			],	# minimum intron size to be included in splicigar (default is 40)
## FILES AND DATABASES
		_db   => [ 'rex',	 						'read.write'			],	#SQL_database_name
		_splc => [ '',								'read.write'			],	#junction database
		_cexo => [ '',								'read.write'			],	#custom exons
		_estr => [ '1',								'read.write'			],	#exon junction strigency for multisplice
		_read => [ '',			 					'read.write'			],	#the read file (FASTQ)
		_refs => [ '',			 					'read.write'			],	#the reference genome
		_edb  => [ '',			 					'read.write'			],	#Ensembl_SQL_database_name
		_rtab => [ '',			 					'read.write'			],	#retro table name
# STANDARD FILENAMES                       
		_unmp => [ 'unmapped',	 					'read.write'			],	#unmapped prefix
		_unrp => [ 'repeat',	 					'read.write'			],	#repeated (maps too many times)
## ANALYSIS TAG
		_atag => [ '1',								'read.write'			],	# analysis tag to be used
		_emod => [ '1',								'read.write'			],	# exon duplication model (0: zero, 1:mockmap, 2:duplication)
## SERVER STUFF
		_sqlhost  => [ 'hkserv.local',				'read.write'			],	#SQL database host (usually default)
		_sqlport  => [ '3306',						'read.write'			],	#SQL database port (usually default)
		_sqluser  => [ 'dbrawand',					'read.write'			],	#SQL database user (usually default)
		_proto    => [ 'ANALYSIS',					'read.write'			], # the protocol table
		_annotid  => [ 'ANNOTATION',				'read.write'			], # the ANNOTATION ID TABLE
		_mapid    => [ 'MAPPING',					'read.write'			], # the MAPPING ID TABLE
## BASIC DATA (quality baseline/realdlength)
		_readl => [ $defaultlength,					'read.write'			], # qualbaseline
		_quba  => [ substr('aaaaaaaaaaaaa```````````````_______^^^^^^^]]]]]\\\\\[[[ZZZYYXXXXWWVVVUUUTTR',0,$defaultlength),'read.write'], # qualbaseline
## SPLIT AND PAIRED END
		_paired   => [ '0',							'read.write'			], # max distance between mappings to be preserved ()
		_mindist  => [ '0',							'read.write'			], # min distance between mappings to be preserved ()
		_maxdist  => [ '0',							'read.write'			], # max distance between mappings to be preserved ()
## RAPARM ANALYSIS
		_retros   => [ '0',							'read.write'			], # ADD RETROGENES TO ANNOTATION IF DETECTED
		_strands  => [ '0',							'read.write'			], # RESPECT STRANDS
		_bestmap  => [ '0',							'read.write'			], # USE ONLY MAPPING WITH LEAST MISMATCHES
		_threads  => [ '4',							'read.write'			], # number of threads to use
#SPECIALS SET ON LOAD
		_mappingid  => [ '',						'read.write'			], # the mapping id to avoid opening a filehandle each time
		_annotationid  => [ '',						'read.write'			], # the annotation id set by protocol restore
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
sub remote_server     { return "dee-serv02.vital-it.ch" }
sub remote_scratch    { return "/scratch/frt/weekly/dbrawand/"}
sub remote_raparm_src { return '~/bioperl/cprog/raparm-ng/*.?pp'}
sub remote_raparm_bin { return "raparm-ng-x86_64_linux"}
sub remote_ames_src   { return '~/bioperl/cprog/ames-ms/*.?pp'}
sub remote_ames_bin   { return "ames-ms-x86_64_linux"}

### DYNFILES
# file suffixes
my $il = 'islands';
my $rf = 'ref';
my $rd = 'read';
my $sp = 'splice';
my $uq = 'unique';
my $rp = 'raparm';
# file extensions
my $gff   = '.gff';
my $fasta = '.fa';
my $txt   = '.txt';
my $cns   = '.cns';
# Short species name
sub _short    { return "_" . _s_($_[0]) . "_" . ($_[1] ? $_[1] : annotationID($_[0]) ) }
sub _shortens { return "_" . _s_($_[0]) . "_" . $_[0]->{_edb} }
sub shortchk  { return _s_($_[0]) . "_" . annotationID($_[0]) }

sub _s_ {
	my $ret;
	if ($_[0]->{_edb} =~ /macaca/) {
		$ret = (($_[0]->{_edb} =~ /(.)\w+\_(.).(.)\w+\_core\_(\d+)/) ? $1 . $2 . $3 . $4 : die);
	}
	else {
		$ret = (($_[0]->{_edb} =~ /(.)\w+\_(..)\w+\_core\_(\d+)/) ? $1 . $2 . $3 : die ">> $_[0]->{_edb} <<");
	}
	return $ret;
}

sub annotationID {
	if ($_[0]->{_annotationid}) { return $_[0]->{_annotationid} }
	#my $a = (-s $_[0]->{_cexo}) ? $_[0]->{_cexo} : '';
	#my $b = (-s $_[0]->{_splc}) ? $_[0]->{_splc} : '';
	my $a = $_[0]->{_cexo};
	my $b = $_[0]->{_splc};
	my $c = $_[0]->{_estr};
	my $d = $_[0]->{_mintron};
	my $e = $_[0]->{_edb};
	my $f = $_[0]->{_refs};
	my $g = $_[0]->{_readl};
	my $h = (-s $_[0]->{_rtab}) ? $_[0]->{_rtab} : '';
	return substr(md5_hex($a . $b . $c . $d . $e . $f . $g . $h),0,8);
}
sub mappingID {
	if ($_[0]->{_mappingid} && !($_[1])) {
		# reuse 
	}
	else {
		my $a = annotationID($_[0]);
		my $b = $_[0]->{_edb};
		my $c = $_[0]->{_tiss};
		my $d = $_[0]->{_gend};
		my $e = $_[0]->{_read};
		my $f = headstring($_[0]->{_read});
		$_[0]->{_mappingid} = substr(md5_hex($a . $c . $d . $e . $f),0,8);
		# store to sql
		store_mappingID($_[0],$f);
	}
	return $_[0]->{_mappingid};
}

sub headstring{
	open(HEAD, "head $_[0] |");
	while (<HEAD>){
		if (/@([^:]+:\d):\d+:-?\d+:\d+/) {
			close(HEAD);
			return $1;
		}
	}
	die "\nERROR: Cannot parse readfile ($_[0]) -> set mapid manually\n";
}

# splice sequence length
sub splicelength { my ($self, $arg) = @_; return 2*($self->{_readl} - $self->{_ovrh}); }


sub annot          { my ($self, $arg, $arg2) = @_; $arg =~ s/\./_/g; return $arg . _short($self, $arg2) }
sub ensannot       { my ($self, $arg) = @_; $arg =~ s/\./_/g; return $arg . _shortens($self) }
sub annotab        { my ($self, $arg) = @_; $arg =~ s/\./_/g; return $arg . _short($self) }
sub exonmodel      { 
	my ($self, $arg) = @_;
	if    ($self->{_emod} == 0) { return $self->rzexon }
	elsif ($self->{_emod} == 1) { return $self->rcexon }
	elsif ($self->{_emod} == 2) { return $self->rbexon }
	elsif ($self->{_emod} == 3) { return $self->rdexon }
	warn("\nWARNING: Unknown Duplication model, using default (RC)\n");
	return $self->rcexon;
}
# Mapping files
sub rawmap         { my ($self, $arg) = @_; return ($arg ? $arg : mappingID($self)) . '_rawmap'      } # argument overrides id
sub splicedmap     { my ($self, $arg) = @_; return ($arg ? $arg : mappingID($self)) . '_splicedmap'  } # argument overrides id
# dyn routines
sub splcread       { my ($self, $arg) = @_; return $self->{_exid} . '_' . $sp . $txt }
sub uniq           { my ($self, $arg) = @_; return $self->{_exid} . '_uniq'  }
sub introntab      { my ($self, $arg) = @_; return $self->{_exid} . '_introns'  }
sub cln_introntab  { my ($self, $arg) = @_; return $self->{_exid} . '_introns_clean'  }
sub splicedbwtout  { my ($self, $arg) = @_; return mappingID($self) . "_spliced" . '.bwtout' }
sub mapping        { my ($self, $arg) = @_; return $self->mappingID() . "_full"    . '.bwtout' }  
sub cleanannotab   { my ($self, $arg) = @_; return substr($arg,0,1) . _short($self) }
sub custom         { my ($self, $arg) = @_; return $arg . _short($self) }
sub unmapped_fa    { my ($self, $arg) = @_; return mappingID($self) . "_" . $self->{_unmp} . $fasta } # fa unmapped reads
sub repeat_fa      { my ($self, $arg) = @_; return mappingID($self) . "_" . $self->{_unrp} . $fasta } # fa repeated reads
sub reads          { my ($self, $arg) = @_; return mappingID($self) . "_" . $rd }
sub matepairs      { my ($self, $arg) = @_; return mappingID($self) . "_" . 'matepairs' }

# ENSEMBL ANNOTATION FILES
sub ensembljunctions { my ($self, $arg) = @_; return $self->{_edb} . '_junctions' }
sub ensemblexons     { my ($self, $arg) = @_; return $self->{_edb} . '_exons' }

# RETRO
sub retroidtable { my ($self, $arg) = @_; return annotationID($self) . "_retroid" } 

# statics
sub rfexon { my ($self, $arg) = @_; return 'RFexon' . _short($self) } # without uniscore
sub rcexon { my ($self, $arg) = @_; return 'RCexon' . _short($self) } # mockmap evaluation
sub rbexon { my ($self, $arg) = @_; return 'RBexon' . _short($self) } # mockmap evaluation (without equilibrating splice anchors)
sub rzexon { my ($self, $arg) = @_; return 'RZexon' . _short($self) } # zero model
sub rdexon { my ($self, $arg) = @_; return 'RDexon' . _short($self) } # erange like duplicate model

# mock reads (unique)
sub mrdb             { my ($self, $arg) = @_; return 'mockreads_' . annotationID($self) }
sub rawmockreads     { my ($self, $arg) = @_; return 'rawmockreads_' . annotationID($self) }
sub mockgenomebwtout { my ($self, $arg) = @_; return 'mockreads_' . ($arg ? $arg : annotationID($self)) . '_genomemap' }
sub mocksplicebwtout { my ($self, $arg) = @_; return 'mockreads_' . ($arg ? $arg : annotationID($self)) . '_splicemap' }
# raparm_ng output files
sub raparm_map         { my ($self, $arg) = @_; return ($arg) ? $self->{_exid} . '_raparm_map_'      . $arg : $self->{_exid} . '_raparm_map'      ; } 
sub raparm_rpkm        { my ($self, $arg) = @_; return ($arg) ? $self->{_exid} . '_raparm_rpkm_'     . $arg : $self->{_exid} . '_raparm_rpkm'     ; }
sub raparm_island      { my ($self, $arg) = @_; return ($arg) ? $self->{_exid} . '_raparm_islands_'  . $arg : $self->{_exid} . '_raparm_islands'  ; }
sub raparm_coverage    { my ($self, $arg) = @_; return ($arg) ? $self->{_exid} . '_raparm_coverage_' . $arg : $self->{_exid} . '_raparm_coverage' ; }
sub raparm_unicoverage { my ($self, $arg) = @_; return ($arg) ? $self->{_exid} . '_raparm_unique_'   . $arg : $self->{_exid} . '_raparm_unique' ; }
sub raparm_splicefreq  { my ($self, $arg) = @_; return ($arg) ? $self->{_exid} . '_raparm_splicefreq_' . $arg : $self->{_exid} . '_raparm_splicefreq' ; }
# PE reads
sub smackout   { my ($self, $arg) = @_; return ($arg) ? ($arg . '_smacked') : $self->{_exid} . '_smacked'; }
#geco/allEns
sub allEns { my ($self, $arg) = @_; return $self->{_edb} . '_allensembl'  }
sub geco   { my ($self, $arg) = @_; return ($arg) ? $self->{_exid} . '_geco_' . $arg : $self->{_exid} . '_geco' ; }

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
	return $self;
}
sub load {
	my ($class, $arg) = @_;
	my $self = bless { }, $class;
	
	die "File not found" unless (-e $arg);
	open(FIN, $arg);
	while(<FIN>) {
		if (/^(\w+)\t([^\n]+)/) { my $a = '_' . $1 ; $self->{$a} = $2 }
		elsif (/#/) { next } # comment
		else { if (/^(\w+)\t([^\n]+)/) {warn "unknown experiment parameter ($a|$2) "} }
	}
	close(FIN);
	foreach my $attribute ($self->_all_attributes()) {
		my($argument) = ($attribute);# =~ /^_(.*)/);
		if (exists $self->{$argument})                      { }# print STDERR "($argument|".$self->{$argument}.")\n" } # argument defined
		elsif ($self->_permissions($attribute, 'required')) { croak("Missing mandatory argument ($argument)") }
		else {
			unless ($attribute =~ /(mappingid|annotationid)/) {
				warn("WARNING: Using default argument [".substr($argument,1,length($argument)-1)."|".$self->_attribute_default($attribute)."]\n");
				$self->{$attribute} = $self->_attribute_default($attribute);
			}
		}
	}
	#mappingID($self); # called to set initial mapping id
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

#returns strand specific suffixes (fwd,rev)
sub strands {
	my ($self, $arg) = @_;
	$arg = ($arg) ? $arg : $self->raparm_coverage;
	
	my $sqldb    = $self->get_db;
	my $sqlhost  = $self->get_sqlhost;
	my $sqlport  = $self->get_sqlport;
	my $sqluser  = $self->get_sqluser;
	
	my $cl = "\"" . $arg . '%' . "\"";
	my @suff;
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my $sql = qq(SHOW TABLES LIKE $cl;);
	my $sth = $dbh->prepare( $sql );
	my $rc = $sth->execute();
	if ($rc) {
		my ( $aa );
		$sth->bind_columns( \$aa );
		while( $sth->fetch() ) {
			if ($aa =~ /$arg/) {
				my $past = $';
				if ($past =~ /\_(\S+)/) {
					push @suff, $1;
				}
			}
		}
	}
	$dbh->disconnect();
	return (scalar(@suff) > 1) ? \@suff : [ 0 ];
}

#######################################################################################################################################
## GLOBALS ############################################################################################################################
#######################################################################################################################################

sub protocol_fields {
	my ($e, $arg) = @_;
	
	my ($sec,$min,$hour,$mday,$mon,$yr,$wday,$yday,$isdst)=localtime(time);
	my $timestamp = sprintf "%4d-%02d-%02d %02d:%02d:%02d", $yr+1900,$mon+1,$mday,$hour,$min,$sec;
	my $pin = "[" . $e->{_mindist} . ',' . $e->{_maxdist} . "]";
	
	my $fields = [
		["id CHAR(64)","species CHAR(64)","gender CHAR(64)","tissue CHAR(64)","type CHAR(64)","sequence_file CHAR(64)","read_length INT","description TEXT","timestamp CHAR(64)","duplication_model INT","exon_tag INT","annotationID CHAR(8)", "mappingID CHAR(8)", 
			"bwt_k INT","bwt_e INT","bwt_m INT","bwt_v INT","bwt_b INT","bwt_l INT","bwt_s INT","overhang INT","exp_date CHAR(64)", "add_retros INT", "strandspecific INT", "best_mapping INT",
			"paired INT", "paired_insert CHAR(32)", "infobase INT", "infoqual INT", "infodrop INT"],
		["\"$e->{_exid}\"","\"$e->{_spec}\"","\"$e->{_gend}\"","\"$e->{_tiss}\"","\"$e->{_type}\"","\"$e->{_read}\"",$e->{_readl},"\"$e->{_desc}\"","\"$timestamp\"",$e->{_emod},$e->{_atag},"\"". annotationID($e) ."\"", "\"". mappingID($e, (-s $e->{_read})) ."\"",
		$e->{_bwtk},$e->{_bwte},$e->{_bwtm},$e->{_bwtv},$e->{_bwtb},$e->{_bwtl},$e->{_bwts},$e->{_ovrh},"\"$e->{_date}\"", $e->{_retros}, $e->{_strands}, $e->{_bestmap}, 
		$e->{_paired}, "\"$pin\"", $e->{_infobase}, $e->{_infoqual}, $e->{_infodrop}]
	];
	return $fields;
}

sub restore_fields {
	my ($e, $arg) = @_;
	my $fields = [
		["id CHAR(64)","species CHAR(64)","gender CHAR(64)","tissue CHAR(64)","type CHAR(64)","sequence_file CHAR(64)","read_length INT","description TEXT","timestamp CHAR(64)","duplication_model INT","exon_tag INT","annotationID CHAR(8)", "mappingID CHAR(8)", 
			"bwt_k INT","bwt_e INT","bwt_m INT","bwt_v INT","bwt_b INT","bwt_l INT","bwt_s INT","overhang INT","exp_date CHAR(64)", "add_retros INT", "strandspecific INT", "best_mapping INT",
			"paired INT", "paired_insert CHAR(32)", "infobase INT", "infoqual INT", "infodrop INT"],
		["_exid","_spec","_gend","_tiss","_type","_read","_readl","_desc","","_emod","_atag","_annotationid", "_mappingid",
		"_bwtk","_bwte","_bwtm","_bwtv","_bwtb","_bwtl","_bwts","_ovrh","_date", "_retros", "_strands", "_bestmap", 
		"_paired", "", "_infobase", "_infoqual", "_infodrop"],
		["edb CHAR(64)","assembly CHAR(64)","exons CHAR(64)","junctions CHAR(64)","min_intron_length INT","splice_exon_stringency INT", "retros CHAR(64)","read_length INT"],
		["_edb","_refs","_cexo","_splc","_mintron","_estr","_rtab","_readl"]
		
	];
	return $fields;
}

sub table_list {
	my ($self, @args) = @_;
	my $tables = [];
	push @{$tables}, $self->rawmap;
	push @{$tables}, $self->splicedmap;
	warn "\nDEBUG: SQL table list still incomplete!\n";
	# raparm tables
	
	# uniqmap
	
	# geco
	
	
	
	return $tables;
}

sub file_list {
	my ($self, @args) = @_;
	my $files = [];
	push @{$files}, $self->annot("splicigarseq");
	push @{$files}, $self->annot("exonseqs");
	push @{$files}, $self->mrdb;
	push @{$files}, $self->mockgenomebwtout;
	push @{$files}, $self->mocksplicebwtout;
	push @{$files}, $self->mapping;
	push @{$files}, $self->splicedbwtout;
	push @{$files}, $self->rfexon;
	push @{$files}, $self->rcexon;
	push @{$files}, $self->rzexon;
	push @{$files}, $self->rdexon;
	return $files;
	
}

#######################################################################################################################################
## DATABASE PROTOCOLS #################################################################################################################
#######################################################################################################################################

sub check_protocol {
	my ($e, $arg) = @_;
	
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};
	my $protocol = $e->{_proto};
	my $expid    = $e->{_exid};
	my $aa;
	my $fields = $e->protocol_fields;
	my $f1 = join ',', @{$fields->[0]};
	
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	$dbh->prepare( qq{ CREATE TABLE IF NOT EXISTS $protocol ($f1); } )->execute();
	
	my $sth = $dbh->prepare( qq{ SELECT id from $protocol WHERE id = \"$expid\"; } );
	$sth->execute();
	$sth->bind_columns( \$aa );
	while ( $sth->fetch() ) {
		if ($aa eq $expid) {
			return 1;
		}
	}
	return 0;
}

sub store_protocol {
	my ($e, $arg) = @_;
	
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};
	my $protocol = $e->{_proto};
	
	my $fields = $e->protocol_fields;
	my @fieldnames;
	foreach my $f (@{$fields->[0]}) {
		if ($f =~ /^(\S+)/) { push @fieldnames, $1 }
		else { die }
	}
	my $f1 = join ',', @{$fields->[0]};
	my $f2 = join ',', @fieldnames;
	my $f3 = join ',', @{$fields->[1]};

	my $sqlstack = [];
	# create table if not exists
	push @{$sqlstack}, "CREATE TABLE IF NOT EXISTS $protocol ($f1);";
	# remove previous entry based on id
	push @{$sqlstack}, "DELETE IGNORE FROM $protocol WHERE ".$fieldnames[0]." = \"$e->{_exid}\";";
	# insert into table
	push @{$sqlstack}, "INSERT INTO $protocol ($f2) VALUES ($f3);";

	# execute sql statements
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	foreach my $sql (@{$sqlstack}) { $dbh->prepare( $sql )->execute() }
	print STDERR "\n -> PROTOCOL UPDATED\n";
	return;
}

### replaces restore_protocol
sub buildFromID {
	my ($self, $arg) = @_;
	my $sqldb    = $self->{_db};
	my $sqlhost  = $self->{_sqlhost};
	my $sqlport  = $self->{_sqlport};
	my $sqluser  = $self->{_sqluser};
	my $protocol = $self->{_proto};
	my $annot    = $self->{_annotid};
	
	### read from protocol
	my $fields = $self->restore_fields;
	my @fieldnames;
	foreach my $f (@{$fields->[0]}) {
		if ($f =~ /^(\S+)/) { push @fieldnames, $1 }
		else { die }
	}
	my $fieldstring = join ", ", @fieldnames;
	my $sql = "SELECT $fieldstring FROM $protocol WHERE id = \"$arg\"";
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my $sth = $dbh->prepare( $sql );
	$sth->execute();
	my $ref = $sth->fetchall_arrayref;
	foreach my $row ( @{$ref} ) {
    	#print STDERR "@$row\n";
		for(my $field = 0; $field < scalar(@{$row}); ++$field) {
			print STDERR $field . "\t". $row->[$field] . " <=> ";
			if ($fields->[1]->[$field]) {
				$self->{$fields->[1]->[$field]} = $row->[$field];
				print STDERR $fields->[1]->[$field];
			}
			print STDERR "\n";
		}
	}
	# read from annotation
	@fieldnames = ();
	foreach my $f (@{$fields->[2]}) {
		if ($f =~ /^(\S+)/) { push @fieldnames, $1 }
		else { die }
	}
	$fieldstring = join ", ", @fieldnames;
	$sql = "SELECT $fieldstring FROM $annot WHERE annotationid = \"".$self->{_annotationid}."\"";
	$sth = $dbh->prepare( $sql );
	$sth->execute();
	$ref = $sth->fetchall_arrayref;
	foreach my $row ( @{$ref} ) {
		for(my $field = 0; $field < scalar(@{$row}); ++$field) {
			print STDERR $field . "\t". $row->[$field] . " <=> ";
			if ($fields->[3]->[$field]) {
				$self->{$fields->[3]->[$field]} = $row->[$field];
				print STDERR $fields->[3]->[$field];
			}
			print STDERR "\n";
		}
	}
	$dbh->disconnect();
	
	print STDERR "\n -> BUILT EXPERIMENT FROM PROTOCOL (still beta)\n\n";
	return $self;
}

sub writeCTL {
	my ($e, $file) = @_;
	open(FOUT, ">$file") or die "FATAL: Cannot open $file for writing.";

	print FOUT "## => FOR RERUN CHANGE (exid) AND (read)\n";
	print FOUT "exid\t" . $e->{_exid} . "\n"; 
	print FOUT "spec\t" . $e->{_spec} . "\n";
	print FOUT "gend\t" . $e->{_gend} . "\n";
	print FOUT "tiss\t" . $e->{_tiss} . "\n";
	print FOUT "type\t" . $e->{_type} . "\n";
	print FOUT "desc\t" . $e->{_desc} . "\n";
	print FOUT "date\t" . $e->{_date} . "\n";
	print FOUT "read\t" . $e->{_read} . "\n";
	print FOUT "infobase\t" . $e->{_infobase} . "\n";
	print FOUT "infoqual\t" . $e->{_infoqual} . "\n";
	print FOUT "infodrop\t" . $e->{_infodrop} . "\n";
	print FOUT "bwtk\t" . $e->{_bwtk} . "\n";
	print FOUT "bwtm\t" . $e->{_bwtm} . "\n";
	print FOUT "bwtv\t" . $e->{_bwtv} . "\n";
	print FOUT "bwtb\t" . $e->{_bwtb} . "\n";
	print FOUT "bwtl\t" . $e->{_bwtl} . "\n";
	print FOUT "bwts\t" . $e->{_bwts} . "\n";
	print FOUT "bwte\t" . $e->{_bwte} . "\n";
	print FOUT "strands\t" . $e->{_strands} . "\n";
	print FOUT "retros\t" . $e->{_retros} . "\n";
	print FOUT "bestmap\t" . $e->{_bestmap} . "\n";
	print FOUT "emod\t" . $e->{_emod} . "\n";
	print FOUT "atag\t" . $e->{_atag} . "\n";
	print FOUT "paired\t" . $e->{_paired} . "\n";
	print FOUT "mindist\t" . $e->{_mindist} . "\n";
	print FOUT "maxdist\t" . $e->{_maxdist} . "\n";
	print FOUT "refs\t" . $e->{_refs} . "\n";
	print FOUT "# DEFAULT\n";
	print FOUT "db\t" . $e->{_db} . "\n"; #DEFAULT
	print FOUT "edb\t" . $e->{_edb} . "\n";
	print FOUT "# DEFAULT\n";
	print FOUT "sqlport\t" . $e->{_sqlport} . "\n"; #DEFAULT
	print FOUT "# DEFAULT\n";
	print FOUT "sqlhost\t" . $e->{_sqlhost} . "\n"; #DEFAULT
	print FOUT "# DEFAULT\n";
	print FOUT "sqluser\t" . $e->{_sqluser} . "\n"; #DEFAULT
	print FOUT "# DEFAULT\n";
	print FOUT "proto\t" . $e->{_proto} . "\n"; #DEFAULT
	print FOUT "# DEFAULT\n";
	print FOUT "annotid\t" . $e->{_annotid} . "\n"; #DEFAULT
	print FOUT "# DEFAULT\n";
	print FOUT "mapid\t" . $e->{_mapid} . "\n"; #DEFAULT
	print FOUT "# DEFAULT\n";
	print FOUT "quba\t" . $e->{_quba} . "\n"; #DEFAULT
	print FOUT "readl\t" . $e->{_readl} . "\n";
	print FOUT "splc\t" . $e->{_splc}. "\n";
	print FOUT "cexo\t" . $e->{_cexo} . "\n";
	print FOUT "rtab\t" . $e->{_rtab} . "\n";
	print FOUT "estr\t" . $e->{_estr} . "\n";
	print FOUT "ovrh\t" . $e->{_ovrh} . "\n";
	print FOUT "mintron\t" . $e->{_mintron} . "\n";
	print FOUT "# DEFAULT\n";
	print FOUT "threads\t" . $e->{_threads} . "\n";
	print FOUT "ctag\t1\n";
	print FOUT "# DEFAULT\n";
	print FOUT "unmp\t" . $e->{_unmp} . "\n";
	print FOUT "# DEFAULT\n";
	print FOUT "unrp\t" . $e->{_unrp} . "\n";

	print FOUT "## BUILT FROM #########\n";
	print FOUT "# ANNOTATIONID ".$e->{_annotationid}."\n";
	print FOUT "# MAPPINGID    ".$e->{_mappingid}."\n";
	print FOUT "## END ################\n";

	close(FOUT);
	return;
}

sub remove_protocol {
	my ($e, $arg) = @_;
	
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};
	my $protocol = $e->{_proto};
	my $expid    = $e->{_exid};
	my $fields = $e->protocol_fields;
	my @fieldnames;
	foreach my $f (@{$fields->[0]}) {
		if ($f =~ /^(\S+)/) { push @fieldnames, $1 }
		else { die }
	}

	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	$dbh->prepare( "DELETE IGNORE FROM $protocol WHERE ".$fieldnames[0]." = \"$expid\";" )->execute();
	
	print STDERR "\n -> PROTOCOL UPDATED (removed experiment $expid)\n";
	return;
}

sub drop_table {
	my ($e, $arg) = @_;
	
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};

	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	$dbh->prepare( "DROP TABLE IF EXISTS $arg;" )->execute();	
	return;
}

sub getSpliceRedundancy {
	my ($e, $arg) = @_;
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};
	my $protocol = $e->{_annotid};
	my $aid = annotationID($e);
	# get sample name by sequence file
	my $sql = "SELECT DISTINCT splice_seq_redundancy FROM $protocol WHERE annotationID  = \"$aid\"";
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my $sth = $dbh->prepare( $sql );
	$sth->execute();
	my ($aa);
	$sth->bind_columns( \$aa );
	while( $sth->fetch() ) {
		return $aa;
	}
	# fallback
	warn "\nWARNING: redundancy not calculated yet (this may take a few seconds)\n";
	return rexARE::spliceSeqRedu($e);
	
}

sub updateSpliceRedundancy {
	my ($e, $arg) = @_;
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};
	my $protocol = $e->{_annotid};
	my $aid = annotationID($e);
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	$dbh->prepare( qq{ UPDATE $protocol SET splice_seq_redundancy = $arg WHERE annotationID = \"$aid\" } )->execute();
	$dbh->disconnect();
	print STDERR "\n -> UPDATED ANNOTATION INFO\n";
}

sub store_annotationID {
	my ($e, $arg) = @_;
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};
	my $protocol = $e->{_annotid};
	
	my $aid = annotationID($e);
	my $fields = [
		["annotationID CHAR(8)","edb CHAR(64)","assembly CHAR(64)","exons CHAR(64)","junctions CHAR(64)","min_intron_length INT","splice_exon_stringency INT", "retros CHAR(64)","read_length INT", "splice_seq_redundancy INT"],
		["\"$aid\"","\"$e->{_edb}\"","\"$e->{_refs}\"","\"$e->{_cexo}\"","\"$e->{_splc}\"","\"$e->{_mintron}\"","\"$e->{_estr}\"","\"$e->{_rtab}\"","\"$e->{_readl}\"", "NULL"]
	];
	my @fieldnames;
	foreach my $f (@{$fields->[0]}) {
		if ($f =~ /^(\S+)/) { push @fieldnames, $1 }
		else { die }
	}
	my $sqlstack = [];
	my $f1 = (join ',', @{$fields->[0]}) . ", PRIMARY KEY (annotationID)";
	my $f2 = join ',', @fieldnames;
	my $f3 = join ',', @{$fields->[1]};
	
	
	# create table if not exists
	push @{$sqlstack}, "CREATE TABLE IF NOT EXISTS $protocol ($f1);";
	# insert into table
	push @{$sqlstack}, "INSERT IGNORE INTO $protocol ($f2) VALUES ($f3);";

	# execute sql statements
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	foreach my $sql (@{$sqlstack}) { $dbh->prepare( $sql )->execute() }
	$dbh->disconnect();
	print STDERR "\n -> ANNOTATION TABLE UPDATED\n";
	return;
}

sub store_mappingID {
	my ($e, $arg) = @_;
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};
	my $protocol = $e->{_mapid};
	
	unless ($arg) { $arg = headstring($e->{_read}) }
	my $aid = mappingID($e);
	
	my $fields = [
		["mappingID CHAR(8)","edb CHAR(64)","tissue CHAR(64)","gender CHAR(64)","readfile CHAR(64)","read_header CHAR(64)", "annotationID CHAR(8)", "postfilter_reads INT", "mapped_reads INT"],
		["\"$aid\"","\"$e->{_edb}\"","\"$e->{_tiss}\"","\"$e->{_gend}\"","\"$e->{_read}\"","\"$arg\"", "\"". annotationID($e) ."\"", "NULL", "NULL"]
	];
	my @fieldnames;
	foreach my $f (@{$fields->[0]}) {
		if ($f =~ /^(\S+)/) { push @fieldnames, $1 }
		else { die }
	}
	my $sqlstack = [];
	my $f1 = (join ',', @{$fields->[0]}) . ", PRIMARY KEY (mappingID)";
	my $f2 = join ',', @fieldnames;
	my $f3 = join ',', @{$fields->[1]};
	
	
	# create table if not exists
	push @{$sqlstack}, "CREATE TABLE IF NOT EXISTS $protocol ($f1);";
	# insert into table
	push @{$sqlstack}, "INSERT IGNORE INTO $protocol ($f2) VALUES ($f3);";

	# execute sql statements
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	foreach my $sql (@{$sqlstack}) { $dbh->prepare( $sql )->execute() }
	$dbh->disconnect();
	return;
}

sub countReads {
	my ($e, $arg) = @_;
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};
	my $protocol = $e->{_mapid};
	my $readfile = $e->{_read};
	my $m = $e->mappingID;
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	if (-s $readfile) {
		my $readcount = 0;
		print STDERR " -> counting reads...\n";
		open(RD, "grep -c \"^+\$\" $readfile |") or warn "READ FILE ($readfile) NOT FOUND\n";
		while(<RD>) {
			if (/(\d+)/) { $readcount = $1;last }
		}
		close(RD);
		if ($readcount == 0) {
			print STDERR "\n -> COULD NOT UPDATE READ STATISTICS FOR SET $m (parsing error)\n";
		}
		else {
			$dbh->prepare( qq{ UPDATE $protocol SET postfilter_reads = $readcount WHERE mappingID = \"$m\" } )->execute();
			print STDERR "\n -> UPDATED READ STATISTICS INFO FOR SET $m\n";
		}
	}
	else {
		print STDERR "\n -> COULD NOT UPDATE READ STATISTICS FOR SET $m (readfile $readfile not found)\n";
	}
	$dbh->disconnect();
	return;
}

sub updateReadCount {
	my ($e, $arg) = @_;
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};
	my $protocol = $e->{_mapid};
	my $m = $e->mappingID;
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my $readtable = reads($e);
	my $sth = $dbh->prepare( qq{ select count(id) from \`$readtable\` } );
	$sth->execute();
	my ($rc, $readcount);
	$sth->bind_columns( \$rc );
	while( $sth->fetch() ) { $readcount = $rc }
	print STDERR "\t$readcount\n";
	$dbh->prepare( qq{ UPDATE $protocol SET postfilter_reads = $readcount WHERE mappingID = \"$m\" } )->execute();
	print STDERR "\n -> UPDATED READ STATISTICS INFO FOR SET $m\n";
	$dbh->disconnect();
	return;
}

sub updateMapCount {
	my ($e, $arg) = @_;
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};
	my $protocol = $e->{_mapid};
	my $m = $e->mappingID;
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my $rawmap = $e->rawmap();
	my $splmap = $e->splicedmap();
	print STDERR "\nCounting Mapped Reads for mappingID $m\n";
	# count mappings
	my ($rid, $count); #the read read ids
	my %readids;
	## from rawmap
	my $sth = $dbh->prepare( qq{ select readid from \`$rawmap\` } );
	$sth->execute();
	$sth->bind_columns( \$rid );
	while( $sth->fetch() ) { $readids{$rid} = 1; if (++$count % 1000 == 0) { print STDERR "\r\t$count" } }
	print STDERR "\r\t$count\n";
	## from splicedmap
	$sth = $dbh->prepare( qq{ select readid from \`$splmap\` } );
	$sth->execute();
	$sth->bind_columns( \$rid );
	while( $sth->fetch() ) { $readids{$rid} = 1; if (++$count % 1000 == 0) { print STDERR "\r\t$count" } }
	print STDERR "\r\t$count\n";
	my $mappedreads = scalar(keys(%readids));
	# update mapping table
	$dbh->prepare( qq{ UPDATE $protocol SET mapped_reads = $mappedreads WHERE mappingID = \"$m\" } )->execute();
	print STDERR "\n -> UPDATED MAPPING STATISTICS INFO FOR SET $m\n";
	$dbh->disconnect();
	return;
}


sub updateReadStats {
	my ($e, $arg) = @_;
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};
	my $protocol = $e->{_mapid};
	
	my $nobuffer = 1;
	
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	
	# get mapping ids without stats
	my $sth = $dbh->prepare( ($arg == 0) ? qq{select mappingID from $protocol } : qq{ select mappingID from $protocol where postfilter_reads IS NULL OR mapped_reads IS NULL OR mapped_reads < 10 });
	$sth->execute();
	my ($mapid, @ids);
	$sth->bind_columns( \$mapid );
	while( $sth->fetch() ) { push @ids, $mapid }
	print STDERR "\nWill update read statistics for " . scalar(@ids) . " mappings\n";
	if (scalar(@ids) == 0) { return }
	# get complete table list
	my ($tab, %tables);
	$sth = $dbh->prepare( qq{ SHOW TABLES } );
	$sth->execute();
	$sth->bind_columns(\$tab);
	while($sth->fetch()) { $tables{$tab}++ }
	$dbh->disconnect();
	foreach my $m (@ids) {
		$dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport:mysql_server_prepare=1","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
		$e->{_mappingid} = $m;
		my $readtable = reads($e);
		my $rawmap = $e->rawmap();
		my $splmap = $e->splicedmap();
		print STDERR "\n$m\n";
		# if forced reset all
		if ($arg == 0) {
			# set to 0
			$dbh->prepare( qq{ UPDATE $protocol SET postfilter_reads = NULL WHERE mappingID = \"$m\" } )->execute();
			$dbh->prepare( qq{ UPDATE $protocol SET mapped_reads     = NULL WHERE mappingID = \"$m\" } )->execute();
		}
		if ($arg =~ /(0|1|3)/ && $tables{$readtable}) {
			# count reads
			my $readcountquery = qq{ select count(id) from \`$readtable\` };
			$sth = $dbh->prepare( $readcountquery );
			$sth->execute();
			my ($rc, $readcount);
			$sth->bind_columns( \$rc );
			while( $sth->fetch() ) { $readcount = $rc }
			print STDERR "\t$readcount\n";
			$dbh->prepare( qq{ UPDATE $protocol SET postfilter_reads = $readcount WHERE mappingID = \"$m\" } )->execute();
		}
		else {
			#DANGEROUS $dbh->prepare( qq{ DELETE FROM $protocol WHERE mappingID = \"$m\" } )->execute(); # remove entry
		}
		
		if ($arg =~ /(0|2|3)/ && $tables{$rawmap} && $tables{$splmap}) {
			# count mappings
			my ($rid, $count); #the read read ids
			my %readids;
			if ($nobuffer) {
				my $localexecute = "mysql --skip-column-names --quick -u $sqluser -h $sqlhost $sqldb";
				## from rawmap
				my $rawqry = qq{ SELECT readid from \\\`$rawmap\\\` };
				open(EIN, "$localexecute --execute=\"$rawqry\" |") or die;
				while(<EIN>) {
					if (/^(\S+)/) { $readids{$1} = 1; if (++$count % 1000 == 0) { print STDERR "\r\t$count" } }
					else {die}
				}
				close(EIN);
				print STDERR "\r\t$count\n";
				## from splicedmap
				my $splqry = qq{ SELECT readid from \\\`$splmap\\\` };
				open(EIN, "$localexecute --execute=\"$splqry\" |") or die;
				while(<EIN>) {
					if (/^(\S+)/) { $readids{$1} = 1; if (++$count % 1000 == 0) { print STDERR "\r\t$count" } }
					else {die}
				}
				close(EIN);
				print STDERR "\r\t$count\n";
			}
			else {
				## from rawmap
				$sth = $dbh->prepare( qq{ select readid from \`$rawmap\` } );
				$sth->execute();
				$sth->bind_columns( \$rid );
				while( $sth->fetch() ) { $readids{$rid} = 1; if (++$count % 1000 == 0) { print STDERR "\r\t$count" } }
				print STDERR "\r\t$count\n";
				## from splicedmap
				$sth = $dbh->prepare( qq{ select readid from \`$splmap\` } );
				$sth->execute();
				$sth->bind_columns( \$rid );
				while( $sth->fetch() ) { $readids{$rid} = 1; if (++$count % 1000 == 0) { print STDERR "\r\t$count" } }
				print STDERR "\r\t$count\n";
			}
			my $mappedreads = scalar(keys(%readids));

			# update mapping table
			$dbh->prepare( qq{ UPDATE $protocol SET mapped_reads     = $mappedreads WHERE mappingID = \"$m\" } )->execute();
		}
		print STDERR "\n -> UPDATED READ STATISTICS INFO FOR SET $m\n";
		$dbh->disconnect();
	}
	print STDERR "\nALL DONE!\n\n";
	return;
}

sub trygetRIN {
	# tries to get the rin value from samples table (fuzzy search)
	my ($e, $arg) = @_;
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};
	my $samples  = 'samples_sequencing_status';
	my $rin = "unknown RIN";
	# get sample name by sequence file
	my $samplename = ($e->{_exid} =~ /(\d+)/) ? $1 : "";
	if (($samplename !~ /\_/) && ($samplename =~ /^(\d+)/)) { $samplename = $1 } # take obly numeric if there is some
	if ($samplename) {
		my $sql = "SELECT DISTINCT rin FROM $samples WHERE sample_id LIKE \"%".$samplename."%\"";
		my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
		my $sth = $dbh->prepare( $sql );
		$sth->execute();
		my ($aa, $r);
		$sth->bind_columns( \$aa );
		while( $sth->fetch() ) {
			if (++$r > 1 && $aa ne $rin) {
				$rin = "ambigous";
				last;
			}
			$rin = $aa;
		}
	}
	return $rin;
}

sub database_busy {
	my ($e, $arg) = @_;
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};
	my ($busycount1, $busycount2, $busycount3) = (0,0,0);
	my $sql = "SHOW PROCESSLIST;";
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my $sth = $dbh->prepare( $sql );
	$sth->execute();
	my ( $id, $user, $host, $db, $command, $time, $state, $info );
	$sth->bind_columns( \$id, \$user, \$host, \$db, \$command, \$time, \$state, \$info );
	while( $sth->fetch() ) {
		next unless $info;
		#print STDERR $id . "\t";
		#print STDERR $user . "\t";
		#print STDERR $host . "\t";
		#print STDERR $db . "\t";
		#print STDERR $command . "\t";
		#print STDERR $time . "\t";
		#print STDERR substr($info,0,50) . "\n";
		if ($db eq $e->{_db}) {
			if ($info =~ /index/i) { 
				$busycount1++;
			}
			elsif ($info =~ /load/i || $info =~ /insert/i) {
				$busycount2++;
			}
			elsif ($info =~ /select/i) {
				$busycount3++;
			}
		}
	}
	if ($busycount1 || ($busycount3 + $busycount2) > 1) { # busycount3 is not dramatic
		return [$busycount1, $busycount2, $busycount3];
	}
	return 0;
}

sub wait4database {
	my ($e, $arg) = @_;
	$arg = ($arg) ? $arg : 20;
	$arg /= 10;
	sleep 2;
	my $waiter = time();
	while (my $business = $e->database_busy()) {
		print STDERR "\rSQL DATABASE BUSY ($business->[0]|$business->[1]|$business->[2]) ";
		for (my $i=0;$i<$arg;$i++) {
			sleep 10;
			print STDERR ".";
		}
		print STDERR "\r                              " . (" " x $arg) ;
	}
	$waiter = time() - $waiter;
	print STDERR "\rDATABASE AVAILABLE (waited $waiter seconds)\n" if ($waiter > 5);
	return;
}

##################
# LEGACY ########
###############


sub restore_protocol {
	my ($e, $arg) = @_;
	my $sqldb    = $e->{_db};
	my $sqlhost  = $e->{_sqlhost};
	my $sqlport  = $e->{_sqlport};
	my $sqluser  = $e->{_sqluser};
	my $protocol = $e->{_proto};
	my $sql = "SELECT id, species, gender, tissue, type, sequence_file, read_length, description, exon_tag, ensembl_annotation, retro_annotation, reference_dna, bwt_k, bwt_e, bwt_m, bwt_v, overhang, timestamp FROM $protocol WHERE id = \"$e->{_exid}\"";
	my $dbh = DBI->connect("DBI:mysql:$sqldb:$sqlhost:$sqlport","$sqluser","") or die "Database connection not made: $DBI::errstr"; 
	my $sth = $dbh->prepare( $sql );
	$sth->execute();
	my        ( $rcount, $aa,  $bb,  $cc,  $dd,  $ee,  $ff,  $gg,  $hh,  $ii,  $jj,  $kk,  $ll,  $mm,  $nn,  $oo,  $pp,  $qq,  $rr );
	$sth->bind_columns( \$aa, \$bb, \$cc, \$dd, \$ee, \$ff, \$gg, \$hh, \$ii, \$jj, \$kk, \$ll, \$mm, \$nn, \$oo, \$pp, \$qq, \$rr );
	while( $sth->fetch() ) {
		$e->{_exid}  = $aa;
		$e->{_spec}  = $bb;
		$e->{_gend}  = $cc;
		$e->{_tiss}  = $dd;
		$e->{_type}  = $ee;
		$e->{_read}  = $ff;
		$e->{_readl} = $gg;
		$e->{_desc}  = $hh;
		$e->{_atag}  = $ii;
		$e->{_edb}   = $jj;
		$e->{_rtab}  = $kk;
		$e->{_refs}  = $ll;
		$e->{_bwtk}  = $mm;
		$e->{_bwte}  = $nn;
		$e->{_bwtm}  = $oo;
		$e->{_bwtv}  = $pp;
		$e->{_ovrh}  = $qq;
		$e->{_date}  = $rr;
		die "ERROR: Cannot restore experiment object from SQL, ambiguity!" if (++$rcount > 1);
	}
	die "restore failed" unless ($rcount);
	return $e;
}


1;

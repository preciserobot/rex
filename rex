#!/usr/bin/perl -w


## Modified version from syntenyslices not using synteny but orthology chains as scoring for synteny confidence on unannotated genomes

#use lib "/Library/Perl/5.10.0/darwin-thread-multi-2level/";

# common Libs
use DBI;
use DBD::mysql;
use strict;
use Getopt::Long;
use Digest::MD5 qw(md5 md5_hex md5_base64);

# own libs
use Usu;
use rexExperiment;
use rexDB;
use rexMAP;
use rexJUNCTION;
use rexRAPARM;
use rexARE;
use rexSPLIT;
use rexFILTER;
use rexENSEMBL;

######
my $prog_title = "REX \"SmallRex's Big Brother\"";
######

Usu::header( $0, $prog_title,["RNA-seq analysis pipeline",
								"with read distribution,",
								"paired-end filtering,",
								"customized annotation interface,",
								"multisplice awareness,",
								"retro detection,",
								"read swampling,",
								"and many other things...",
								"(including lack of free time, 2 weeks in one...)"], "DC Brawand" );
my ($sec,$min,$hour,$mday,$mon,$yr,$wday,$yday,$isdst)=localtime(time);
printf STDERR "\n\tPID $$ started at %4d-%02d-%02d %02d:%02d:%02d\n\n\n", $yr+1900,$mon+1,$mday,$hour,$min,$sec;

# getting options
my ($expfile, $slimpipe, $are, $jump, $sandbox, $ens, $forcelocal, $cleanup, $compress, $killexp, $readstats, $override, $rungeco, $live, $swow, $countreads, $swample, $getdetails, $indexmappings, $uniquemappings, $writectl, $fast) = ('','','','','','','','','','','','','','2','','','','','','','','');
GetOptions(	'e|experiment:s'   => \$expfile,        #
			'z|run'            => \$slimpipe,       #
			'a|are'            => \$are,            #
			'j|jump:s'         => \$jump,           #
			'b|sandbox'        => \$sandbox,        #
			'n|ensembl'        => \$ens,            #
			'l|forcelocal'     => \$forcelocal,     #
			'c|clean'          => \$cleanup,        #
			'k|kill'           => \$killexp,        #
			'p|compress'       => \$compress,       #
			'r|readstats'      => \$readstats,      # Calculates read stats (counts read table and distinct mappings foreach mappingID)
			'o|override:s'     => \$override,       # will use a previous annotation mockmap to run ames-ms
			'g|geco'           => \$rungeco,        # will use a previous annotation mockmap to run ames-ms
			'v|live:s'         => \$live,           # use live read id translation
			'w|swow'           => \$swow,           # SUPER SWOW FUNCTION
			'm|swample:s'      => \$swample,        # runs raparm with reswampled reads
			's|readupdate'     => \$countreads,     # uses quick routine to count reads (directly on readfile)
			'd|details'        => \$getdetails,     # get mapping details from RAPARM-MS
			'i|indexmappings'  => \$indexmappings,  # index mapping and read tables on readid to allow for faster extraction for (rapa|u)rm-ms
			'u|uniquemappings' => \$uniquemappings, # analysis only for unque mappings () 
			'f|writectl:s'     => \$writectl,
			'x|expresslane'    => \$fast
			);


my $experiment;
if ($expfile) {
	print STDERR "\n";
	if (-s $expfile) {
		$experiment = rexExperiment->load($expfile);
		print STDERR "\tControl File: " . $expfile . "\n";
		if (-s $experiment->get_read && $experiment->get_exid =~ /_/) { $experiment->store_protocol() }
	}
	else {
		# try to restore protocol from id
		$experiment = rexExperiment->new(exid => $expfile);
		$experiment->buildFromID($expfile);
		if ($writectl) { $experiment->writeCTL($writectl); exit(0) }
	}
	print STDERR "\tAnnotationID: " . $experiment->annotationID . "\n";
	print STDERR "\t   MappingID: " . $experiment->mappingID . "\n" if (-s $experiment->get_read || !(-s $expfile));
	print STDERR "\t  AnalysisID: " . $experiment->get_exid . "\n";
	print STDERR "\n";
	
} else { Usage() }


if ($killexp) {
	Usu::confirm("THIS WILL DROP ALL ANALYSIS TABLES FOR GIVEN CONTROLFILE\n");
	foreach my $t (@{$experiment->table_list}) { $experiment->drop_table($t) }
	foreach my $f (@{$experiment->file_list})  { Usu::countdownRemove(5, $f) }
	$experiment->remove_protocol();
	exit(0);
}

if ($cleanup) {
	Usu::confirm("THIS REMOVES ALL LOCAL FILES\n");
	foreach my $f (@{$experiment->file_list}) { Usu::countdownRemove(5, $f) }
	exit(0);
}

if ($compress) {
	Usu::confirm("THIS COMPRESSES ALL (OBSOLETE) LOCAL FILES\n");
	my $files = $experiment->file_list();
	foreach my $f (@{$files}) {
		Usu::compress('pbzip2', $f);
	}
	exit(0);
}

if ($readstats) { #r
	print STDERR "\nEnter method to rebuild stats\n\t0:reset and rebuild\n\t1:count postfilter reads\n\t2:count mapped reads\n\t3:count postfilter and mapped reads\n? ";
	my $method = <STDIN>;
	$experiment->updateReadStats($method); # argument forces update and cleans annotation table (removes old non-mapped entries)
	exit(0);
}

if ($countreads) { #s
	rexRAPARM::rebuildreadindex($experiment);
	$experiment->updateReadCount;
	$experiment->updateMapCount;
	#$experiment->countReads;
	exit(0);
}


if ($ens) {
	my $newmethod = 1;
	if (Usu::TimeAndSwitch(1, "Selecting Coordinate System", $jump))                          { rexENSEMBL::selectCoordinateSystem($experiment) } # refetches seq_region_id table to have a full copy
	if (Usu::TimeAndSwitch(2, "Fetching Annotation", $jump))                                  { rexENSEMBL::fetchTables($experiment) }
	if (Usu::TimeAndSwitch(3, "Building Junction File", $jump))                               { rexENSEMBL::buildJunctionFile($experiment) }
	if (Usu::TimeAndSwitch(4, "Building Core exon File", $jump))                              { rexENSEMBL::buildCoreExons($experiment) }  ##\
	exit(0);
}

if ($are) {
	if ($override) {
		Usu::confirm("\nWARNING: The new exon file has to be compatible (strictly overlapping)!\nAll other parameters MUST be unchanged!\n");
		$jump = "2,10,11";
		$experiment->store_annotationID();
		rexDB::cloneLibrary($experiment,$override); # should clone seq_region index as well!
	}

	# fetch ensembl tables
	if ((!($jump) || (($jump =~ /(\d+)/) ? $1 < 2 : 0)) && !($override)) { Usu::confirm("\nWARNING: Updating annotation may render previous analyses to this reference inconsistent!\nshortID is ".$experiment->shortchk."\n") }
	if (Usu::TimeAndSwitch(0, "Building Integer Index", $jump))                               { rexDB::buildIntegerIndex($experiment); $experiment->store_annotationID() }
	if (Usu::TimeAndSwitch(1, "Preparing Reference", $jump))                                  { rexMAP::prepref($experiment) }   
	if (Usu::TimeAndSwitch(2, "Loading Exons (and merge with retros)", $jump))                { rexDB::writeCombinedExonFile($experiment) }
	if (Usu::TimeAndSwitch(3, "Loading junction file", $jump))                                { rexDB::loadJunctions($experiment, $experiment->get_splc) }
	if (Usu::TimeAndSwitch(4, "Calculating Splice library coordinates", $jump))               { rexJUNCTION::generateMultisplice($experiment); rexARE::spliceSeqRedu($experiment) }
	if (Usu::TimeAndSwitch(5, "Building Splice Library", $jump))                              { rexJUNCTION::buildMultispliceSeqs($experiment) }
	if (Usu::TimeAndSwitch(6, "Building Exon Sequences", $jump))                              { rexDB::buildExonSeqs($experiment) } # builds mockreads from transcripts (add them to sql with ignore)
	if (Usu::TimeAndSwitch(7, "Building Mockreads from junctions and exons", $jump))          { rexARE::build_reads_from_junctions_and_exons($experiment) } # builds mockreads from transcripts (add them to sql with ignore) and calciulated splice redundancy
	if (Usu::TimeAndSwitch(8, "Writing Mockreads", $jump))                                    { rexARE::WriteMockReads($experiment) }  # samples max 200000000 reads from database
	if (Usu::TimeAndSwitch(9, "Mapping Mockreads", $jump))                                    { rexARE::MockMap($experiment) } # map on reference and splice sites (second arguments permits skipping of mysql tasks) (QUICK)
	if (Usu::TimeAndSwitch(10, "Running AmEs-ms", $jump))                                     { rexARE::AmEsMs($experiment, $override) } # => RCexons (from mockmap)
	if (Usu::TimeAndSwitch(11, "Running PseudoAmEs", $jump))                                  { rexARE::PseudoAmEs($experiment) } # => RDexons / RZexon
	$experiment->store_annotationID();
	exit(0);
}

if ($slimpipe) {
	# EDIT OVERRIDE TO USE PREVIOUS MAPPING TO RUN RAPARM-MS
	if ($override) {
		Usu::confirm("\nWARNING: The given mapping MUST be compatible with the new referenece (exons overlap / identical splices)!\n");
		$jump = "7,8";
	}
	if ($swample) {
		print STDERR "\nNOTICE: Will run RAPARM after swampling $swample reads\n";
		$jump = "5,7";
		$live = 0;
	}
	if ($uniquemappings) {
		print STDERR "\nNOTICE: Will run reanalysis for unqiue reads (running RM-MS with local results only)\n";
		print STDERR "\n        Make sure mapping and read tables are indexed by readid or use -i option\n" unless ($indexmappings);
		$jump = (-e $experiment->rawmap && -e $experiment->splicedmap) ? "7" : "5-7";
		$indexmappings = (-e $experiment->rawmap && -e $experiment->splicedmap) ? 0 : $indexmappings;
		$live = 0;
	}
	
	# update mappingID (no longer done on control file load)
	$experiment->mappingID(); ################################################################################=====>> Be careful to store appropriate ids when override is effective!
	##### ALL IN ONE #####                                                                            
	if ($jump !~ /(7|5)/ && $experiment->check_protocol && (-s $experiment->mapping))                      { Usu::confirm("\nWARNING: Found previous (partial) analysis. Tables might be overwritten.\n") }
	if (Usu::TimeAndSwitch(0, "Filter", $jump))                                                            { rexFILTER::run($experiment) }
	if (!($live) && Usu::TimeAndSwitch(1, "Storing Sequences", $jump))                                     { rexDB::loadseq($experiment) } # just stores id and creates a read index
	if (Usu::TimeAndSwitch(2, "Mapping to genome", $jump))                                                 { rexMAP::genomemap($experiment, 1) } # seqmap is always forced
	if (Usu::TimeAndSwitch(3, "Mapping to Splice Sites", $jump))                                           { rexMAP::splicemap($experiment, 1) } # seqmap is always forced
	if (!($live) && Usu::TimeAndSwitch(4, "Storing Mappings", $jump))                                      { rexDB::loadmappings($experiment, !($indexmappings)) } # load both mappings
	if ($indexmappings)                                                                                    { rexDB::indexmappings($experiment) } # indexes reads and mappings (when rerunning without files to accelerate SQL extraction)
	if (!($live) && Usu::TimeAndSwitch(5, "Fetching encoded mappings", $jump))                             { rexRAPARM::quickFetchMappings($experiment,$swample) } # uses SQL
	if ($live && Usu::TimeAndSwitch(5, "Preparing mapping files for RAPARM", $jump))                       { rexRAPARM::liveFetch($experiment,$live) } # also stores mappings and seqeunces if requested with live option
	if (Usu::TimeAndSwitch(6, "Filtering AllMap by Pairs", ($experiment->get_paired) ? $jump : -1))        { rexSPLIT::mapfilter($experiment,$forcelocal) }
	if (Usu::TimeAndSwitch(7, "Running (RAPA)RM-MS", $jump))                                               { rexRAPARM::runRaparm($experiment,$forcelocal,$override,$swample,$getdetails,$uniquemappings,$fast) }
	if (!($uniquemappings) && Usu::TimeAndSwitch(8, "Storing Results", $jump))                             { rexRAPARM::storeRaparmOutput($experiment) }
	if (($jump =~ /9/ || $rungeco) && Usu::TimeAndSwitch(9, "Run GECO run!", $jump))                       { rexENSEMBL::geco($experiment) }
	
	#$experiment->store_mappingID();
	#$experiment->store_protocol();
	exit(0);
}

if ($sandbox) {
	Usu::confirm("YOU KNOW WHAT YOU DO!\n");
	exit(0);
}



######################################################################################################################
## LEGACY #############################################################################################################
########################################################################################################################

#############################################################
### SUBS ####################################################
#############################################################
sub Usage {
	Usu::usage($prog_title,["",
							"-e <experiment control file> or <ExperimentID>",
							"-f write control file",
							"",
							"\t-n Get Exon/Junctions from ensembl (required for genome browser)",
							"",
							"\t-a Prepare database (INITIAL)",
							"\t\t-o Exon Override <ANNOTATIONID> (will use mockmap from previous run)",
							"\t\t-j <STEP> Jumps to given step number",
							"",
							"\t-z Analyze sample (generates MAPFILE, RPKM, ISLANDS, COVERAGE, SPLICEFREQ)",
							"\t\t-o Exon Override <MAPPINGID> (will use mappings from previous run, single-end only!)",
							"\t\t-j <STEP> Jumps to given step number",
							"\t\t-v LiveFetch (0:sql|1:live|2:live,store|3:live,store,index) [2]",
							"\t\t-l force local",
							"\t\t-d write and fetch mapping details (HUGE!)",
							"\t\t-g run also Geco (will be skipped otherwise)",
							"\t\t-m <INT> just run raparm after swampling reads",
							"\t\t-i indexes SQL mapping\/read tables",
							"\t\t-u analyze unique reads only (local results only, consider using -i flag)",
							"\t\t-x Expresslane: (RAPA)RM-MS without coverage, islands and splicefreq (Expression only)",
							"",
							"\t-c Local directory cleanup",
							"\t-p Local directory compression",
							"\t-k Kill Experiment (local & sql)",
							"",
							"\t-r Manual read statistics update (whole database)",
							"\t-s Rebuild read index and update counts for current experiment",
							"",
							"\t-b SANDBOX",
							""]);
	exit(0);
}


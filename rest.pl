#!c:/perl/bin/perl.exe
##
##
## rest.pl
## Copyright (C) 2012-2013 Ronald C Beavis, all rights reserved
## The Global Proteome Machine 
## This software is a component of the X! proteomics software
## development project
##
## Use of this software governed by the Artistic license,
## as reproduced at http://www.opensource.org/licenses/artistic-license.php
##
## Version 2013.05.04
## The interface dispatcher for the GPMDB REST 1 interface



use strict;
use CGI qw(:all);
use CGI::Carp qw(fatalsToBrowser);
use HTTP::Request::Common qw(POST GET);
use LWP::UserAgent; 

require "./gpmdb_rest.pl";
require "./gpm_rest.pl";

my $user = "USER_NAME";
my $password = "PASSWORD";

my $cgi = CGI->new();

#print qq(Content-type: text/plain\n\n);
print qq(Content-type: application/json\n\n);

my $type = $cgi->param('type');
my $acc = $cgi->param('acc');
my $seq = $cgi->param('seq');
my $key = $cgi->param('key');
my $w =  $cgi->param('w');
my $filter = $cgi->param('filter');
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$mon++;
$year += 1900;
my $log = "logs/$year\_$mon\_$mday.log";
if(not -d "logs")	{
	system("mkdir logs");
}
my %record = ('SCRIPT_URI' => 1, 'REMOTE_ADDR' => 1,'QUERY_STRING' => 1);
open(LOG,">>$log");
print LOG qq($year-$mon-$mday $hour:$min:$sec\n);
foreach my $key (keys(%ENV)) {
	if($record{$key})	{ 
	    print LOG "\t$key = $ENV{$key}\n";
	}
}
close(LOG);

if($type eq 'interface_version')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbVersion();
}
elsif($type eq 'interface_help')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbHelp();
}
elsif($type eq 'peptide_accessions')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPeptideAccessions($cgi->param('seq'));
}
elsif($type eq 'peptide_count')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPeptideTotal($cgi->param('seq'));
}
elsif($type eq 'peptide_count_z')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPeptideTotalZ($cgi->param('seq'));
}
elsif($type eq 'protein_keyword')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinKeyword($key,$filter);
}
elsif($type eq 'protein_species')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinSpecies($acc);
}
elsif($type eq 'protein_modifications')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinModifications($acc,$cgi->param('mod'),$cgi->param('res'),$cgi->param('max_e'));
}
elsif($type eq 'protein_peptides_z')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinPeptidesCharges($acc);
}
elsif($type eq 'protein_peptides_total')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinPeptidesTotal($acc);
}
elsif($type eq 'protein_peptide_count')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinPeptideCount($acc,$seq);
}
elsif($type eq 'protein_peptide_sequences')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinPeptideSequences($acc);
}
elsif($type eq 'protein_count')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinCount($acc);
}
elsif($type eq 'protein_best_e')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinBestExpect($acc);
}
elsif($type eq 'protein_evidence')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinEvidence($acc);
}
elsif($type eq 'protein_omega')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinOmega($acc,$seq);
}
elsif($type eq 'psyt_motif')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPsytMotif($seq,$w);
}
elsif($type eq 'psyt_predict' or $type eq 'psyt_model')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPsytPredict($seq,$w);
}
elsif($type eq 'protein_polymorphisms' or $type eq 'protein_savs')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinPolymorphisms($acc);
}
elsif($type eq 'protein_sequence')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinSequence($acc);
}
elsif($type eq 'protein_description')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinDescription($acc);
}
elsif($type eq 'protein_model')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbProteinModel($acc);
}
elsif($type eq 'model_metadata')	{
	my $gpm = GPM->new($cgi->param('gpm'),$user,$password);
	print $gpm->gpmModelMetadata();
}
elsif($type eq 'db_metadata')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmDbMetadata($cgi->param('gpm'));
}
elsif($type eq 'db_nssnv' || $type eq 'db_saav')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbNssnv($acc);
}
elsif($type eq 'peptide_pf')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPeptidePf($acc,$cgi->param('pos'),$w,'p');
}
elsif($type eq 'peptide_af')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPeptidePf($acc,$cgi->param('pos'),$w,'a');
}
elsif($type eq 'peptide_uf')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPeptidePf($acc,$cgi->param('pos'),$w,'u');
}
elsif($type eq 'peptide_di')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPeptidePf($acc,$cgi->param('pos'),$w,'d');
}
elsif($type eq 'peptide_ox')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPeptidePf($acc,$cgi->param('pos'),$w,'o');
}
elsif($type eq 'peptide_su')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPeptidePf($acc,$cgi->param('pos'),$w,'s');
}
elsif($type eq 'peptide_nq')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPeptidePf($acc,$cgi->param('pos'),$w,'nq');
}
elsif($type eq 'peptide_ol')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPeptidePf($acc,$cgi->param('pos'),$w,'ol');
}
elsif($type eq 'peptide_ct')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPeptidePf($acc,$cgi->param('pos'),$w,'ct');
}
elsif($type eq 'peptide_qc')	{
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbPeptidePf($acc,$cgi->param('pos'),$w,'qc');
}
elsif($type eq 'model_proteins')	{
	my $gpm = GPM->new($cgi->param('gpm'),$user,$password);
	print $gpm->gpmModelProteins();
}
elsif($type eq 'model_protein_peptides')	{
	my $gpm = GPM->new($cgi->param('gpm'),$user,$password);
	print $gpm->gpmModelProteinPeptides($acc);
}
elsif($type eq 'model_protein_modifications')	{
	my $gpm = GPM->new($cgi->param('gpm'),$user,$password);
	print $gpm->gpmModelProteinModifications($acc);
}
elsif($type eq 'model_protein_polymorphisms' or $type eq 'model_protein_savs')	{
	my $gpm = GPM->new($cgi->param('gpm'),$user,$password);
	print $gpm->gpmModelProteinAps($acc);
}
elsif($type eq 'model_protein_sequence')	{
	my $gpm = GPM->new($cgi->param('gpm'),$user,$password);
	print $gpm->gpmModelProteinSequence($acc);
}
elsif($type eq 'db_nssnv')	{
	my $gpm = GPM->new($cgi->param('gpm'),$user,$password);
	print $gpm->gpmdbProteinNssnv($acc);
}
else	{
	print "Error 404 - interface $type not found\n\nAvailable interfaces are as follows:\n\n";
	my $gpmdb = GPMDB->new($user,$password);
	print $gpmdb->gpmdbHelp();
}

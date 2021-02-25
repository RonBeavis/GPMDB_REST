#!c:/perl/bin/perl.exe
##
# Copyright (C) 2012-2014 Ronald C. Beavis, All rights reserved
#
# Use of this software governed by the Artistic license,
# as reproduced at http://www.opensource.org/licenses/artistic-license.php
# Version 2014.10.18.4:22

package GPMDB;

use DBI;
use JSON;
use List::Util qw[min max];
use strict;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [ qw() ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw( );

our $VERSION = '1.0';

sub new
{
	my ($self,$_u,$_p,$_h) = @_;
	my $h = $_h;
	if(not $h)	{
		$h = "192.168.1.5";
	}
	my $self = {	_user=>$_u,
			_pass=>$_p,
			_host=>$h,
			_mess=>"",
			_dbh_g=>0,
			_dbh_e=>0,
			_dbh_p=>0};

	bless $self, 'GPMDB';
	return $self;
}

sub DESTROY {
	my ($class) = @_;
	my $g = $class->{_dbh_g};
	if($g)	{
		$g->disconnect();
	}
	my $e = $class->{_dbh_e};
	if($e)	{
		e->disconnect();
	}
	my $p = $class->{_dbh_p};
	if($p)	{
		$p->disconnect();
	}
}

sub login	{
	my ($class,$type) = @_;
	my $dsn = "";
	if($type eq "g" and $class->{_dbh_g})	{
		return $class->{_dbh_g};
	}
	elsif($type eq "e" and $class->{dbh_e})	{
		return $class->{_dbh_e};
	}
	elsif($type eq "p" and $class->{dbh_p})	{
		return $class->{_dbh_p};
	}
	elsif($type eq "g" and not $class->{_dbh_g})	{
		$dsn = "DBI:mysql:host=$class->{_host};database=gpmdb";
		$class->{_dbh_g} = DBI->connect($dsn,$class->{_user},$class->{_pass},
				{PrintError => 0, RaiseError => 1}) or die "could not connect to gpmdb";
		return $class->{_dbh_g};
	}
	elsif($type eq "e" and not $class->{dbh_e})	{
		$dsn = "DBI:mysql:host=$class->{_host};database=enspmapdb";
		$class->{_dbh_e} = DBI->connect($dsn,$class->{_user},$class->{_pass},
				{PrintError => 0, RaiseError => 1}) or die "could not connect to enspmapdb";
		return $class->{_dbh_e};
	}
	elsif($type eq "p" and not $class->{dbh_p})	{
		$dsn = "DBI:mysql:host=$class->{_host};database=peakdb";
		$class->{_dbh_p} = DBI->connect($dsn,$class->{_user},$class->{_pass},
				{PrintError => 0, RaiseError => 1}) or die "could not connect to peakdb";
		return $class->{_dbh_p};
	}
	return 0;
}

#input: protein key word
#return: list of protein accession numbers and descriptions
sub gpmdbProteinKeyword {
	my ($class,$key,$type) = @_;
	if(not $type)	{
		$type = "human";
	}           
	my $dbh = login($class,"e");
	my $sql;
	my $sth;
	if ((length($key) <  3 && (length($key) > 0)) || ($key =~ /[a-z]{3}\-\w+/) or $key =~ /\,$/) {
		$sql = "SELECT DISTINCT pid,label,description,db from map where description REGEXP ?";
		$sth = $dbh->prepare($sql) or warn "prepare failed: $DBI::errstr($DBI::err)\n";
		$key =~ s/\,$//;		
		$sth->execute("^$key\,") or warn "execute failed: $DBI::errstr($DBI::err)\n";
	}
	else	{
		$sql = "SELECT DISTINCT pid,label,description,db, match(description) against (? in boolean mode) as score FROM Map WHERE match(description) against(? in boolean mode)";
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute("\"$key\"","\"$key\"") or die "could not execute \"$sql\"";
	}
	my @vs;
	my $out;
	my $ok;
	while(@vs = $sth->fetchrow_array())	{
		$ok = 0;
		if($type eq "human" and @vs[1] =~ /^ENSP\d/)	{
			$ok = 1;
		}
		elsif($type eq "mouse" and @vs[1] =~ /^ENSMUSP\d/)	{
			$ok = 1;
		}
		elsif($type eq "rat" and @vs[1] =~ /^ENSRNOP\d/)	{
			$ok = 1;
		}
		elsif($type eq "dog" and @vs[1] =~ /^ENSCAFP\d/)	{
			$ok = 1;
		}
		elsif($type eq "cow" and @vs[1] =~ /^ENSBTAP\d/)	{
			$ok = 1;
		}
		elsif($type eq "zebrafish" and @vs[1] =~ /^ENSDARP\d/)	{
			$ok = 1;
		}
		elsif($type eq "fruitfly" and @vs[1] =~ /^FBpp\d/)	{
			$ok = 1;
		}
		elsif($type eq "yeast" and (@vs[1] =~ /^Y[A-Z][LR]\d+/ or @vs[1] =~ /^[QR]0\d{3}/))	{
			$ok = 1;
		}
		elsif($type eq "all")	{
			$ok = 1;
		}
		elsif(length($type) > 0 and @vs[1] =~ /$type/)	{
			$ok = 1;
		}
		if($ok)	{
			@vs[2] =~ s/\"/ /g;
			@vs[2] =~ s/\[[^\[\]]+\]$//;
			@vs[2] =~ s/\s+$//;
			my $obs = gpmdbProteinCount($class,@vs[1]);
			$obs =~ s/[\[\]\"]//g;
			if(not $obs)	{
				$obs = 0;
			}
			my $spec = gpmdbProteinPeptidesTotal($class,@vs[1]);
			$spec =~ s/[\[\]\"]//g;
			if(not $spec)	{
				$spec = 0;
			}
			if(1 or ($obs and $spec))	{
				$out .= "{\t\"accession\":\"@vs[1]\",\n\t\"description\":\"@vs[2]\",\n\t\"observed\":\"$obs\",\n\t\"spectra\"\:\"$spec\"\n}\n,";
			}
		}
	}
	$out =~ s/\,$//;
	$sth->finish();
	return "[$out]";
}

#input: protein accession number, modification mass, modified residues, log10(maximum E)
#return: tuple of modified residue:#observations (r1:#1,r2:#2,...,rN:#N);
sub gpmdbProteinModifications {
	my ($class,$acc,$modmass,$mod_res,$max_e) = @_;   
	$mod_res =~ tr/[a-z]/[A-Z]/;
	$mod_res =~ s/\s+//g;                  
	my $maxexpect = -2.0;
	my $sequence = gpmdbProteinSequence($class,$acc);
	$sequence =~ s/[\"\[\]]//g;
	my @seqs = split //,$sequence;
	if(length($max_e))	{
		$maxexpect = $max_e;
	}
	$maxexpect = exp(log(10)*$maxexpect);
	my ($minmass, $maxmass) = ();
	my $delta = 0.5;

	if ($modmass =~ /\.\d+/) {  # change the delta based on the number of sig figs
		my ($mantissa) = $modmass =~ /\.(\d+)/;  # extract mantissa
		$delta /= min(exp(log(10)*length($mantissa)), 1000);  # limit to no smaller than .0005
	}  # end if
	if($modmass >= 0)	{
		$minmass = $modmass - $delta;
		$maxmass = $modmass + $delta;
	}
	else	{
		$minmass = $modmass + $delta;
		$maxmass = $modmass - $delta;
	}
	my $dbh = login($class,"g");
	my $sql = 'select p.pepid, p.start, p.seq from peptide p, protein pro, proseq ps, result r where ps.label = ? and p.proid = pro.proid and 
			pro.proseqid = ps.proseqid and r.active = 1 and pro.resultid = r.resultid and p.expect <= ? order by p.start, p.expect';
	my $sth;
	if ($acc =~ /HIT[0-9]+/ || $acc  =~ /IPI[0-9]+/ || $acc  =~  /At\dg\d+/ || $acc  =~  /^ATCG\d+/i) {
		$sql = 'select p.pepid, p.start, p.seq from peptide p, protein pro, proseq ps, result r where ps.label like ? and p.proid = pro.proid and 
			pro.proseqid = ps.proseqid and r.active = 1 and pro.resultid = r.resultid and p.expect <= ? order by p.start, p.expect';
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute("$acc%",$maxexpect) or die "could not execute \"$sql\"";
	}
	else {
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute($acc,$maxexpect) or die "could not execute \"$sql\"";
	}
	my ($one_pepid, $one_start, $one_seq) = ();
	my %positions;
	while (($one_pepid, $one_start, $one_seq) = $sth->fetchrow_array()) {
		my $sql2 = "";
		my $sth2 = "";
		my $rows2 = 0;
		if($modmass !~ /\.\d/) {
			$sql2 = 'select * from aa where pepid = ? and ROUND(modified) = ? and pm is null'; # sql to get modifications
			$sth2 = $dbh->prepare($sql2) or warn "Error preparing modification fetch: $DBI::errstr($DBI::err)\n";
			$rows2 = $sth2->execute($one_pepid, int($maxmass)) or warn "Error executing modification fetch: $DBI::errstr($DBI::err)\n";
		}
		else	{
			$sql2 = 'select * from aa where pepid = ? and modified between ? and ? and pm is null'; # sql to get modifications
			$sth2 = $dbh->prepare($sql2) or warn "Error preparing modification fetch: $DBI::errstr($DBI::err)\n";
			$rows2 = $sth2->execute($one_pepid, $minmass, $maxmass) or warn "Error executing modification fetch: $DBI::errstr($DBI::err)\n";
		}
		if ($rows2 > 0) { # store modification/mass by start/pepid
			my @res = (); # holder for encountered modified residues
			my @modlist = (); # holder for modifications
			my $flag = 0;
			my @result2;
			while (@result2 = $sth2->fetchrow_array()) {
				push(@res, $result2[2]);
				push(@modlist, $result2[3]);
			}
			my $p = 0;
			my $v;
			foreach $v(@res) {  # for each modified residue
				if (not $mod_res or ($mod_res =~ /$v/ and @seqs[@modlist[$p]-1] eq $v)) {
					$positions{@modlist[$p]} = $positions{@modlist[$p]} + 1;
				}
				$p++;
			}
		}
		$sth2->finish();
	}
	$sth->finish();

	my $c = "";
	my @keys = sort {$a <=> $b} keys(%positions); 
	my $v;
	foreach $v(@keys)	{
		$c .= "\t\"$v\":$positions{$v},\n";
	}
	$c =~ s/\,(\s*)$/$1/;
	return "{$c}";
}

#input: protein accession number, peptide sequences
#return: omega frequency values for each sequence
sub gpmdbProteinOmega {
	my ($class,$acc,$seq) = @_;                     
	my $dbh = login($class,"p");
	my $sql = "select sum(z1total),sum(z2total),sum(z3total),sum(z4total) from protein_omega_count where label=?";
	my $sql_p = "select z1total,z2total,z3total,z4total from protein_omega_count where label=? and seq=?";
	my $sth;
	if ($acc =~ /HIT[0-9]+/ || $acc  =~ /IPI[0-9]+/ || $acc  =~  /At\dg\d+/ || $acc  =~  /^ATCG\d+/i) {
		$sql = "select sum(z1total),sum(z2total),sum(z3total),sum(z4total) from protein_omega_count where label like ?";
		$sql_p = "select z1total,z2total,z3total,z4total from protein_omega_count where label like ? and seq=?";
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute("$acc%") or die "could not execute \"$sql\"";
	}
	else {
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute($acc) or die "could not execute \"$sql\"";
	}
	my @totals = $sth->fetchrow_array();
	$sth->finish();
	$seq =~ s/\s+//g;
	my @peptides = split /\,/,$seq;
	my $pep;
	my @values;
	my $return = "{\"$acc\":[";
	foreach $pep(@peptides)	{
		if(not $pep)	{
			$return .= sprintf("\n{\"%s\":[undef,undef,undef,undef]},",$pep);
		}
		else	{
			if ($acc =~ /HIT[0-9]+/ || $acc  =~ /IPI[0-9]+/ || $acc  =~  /At\dg\d+/ || $acc  =~  /^ATCG\d+/i) {
				$sth = $dbh->prepare($sql_p) or die "could not prepare \"$sql_p\"";
				$sth->execute("$acc%",$pep) or die "could not execute \"$sql_p\"";
			}
			else {
				$sth = $dbh->prepare($sql_p) or die "could not prepare \"$sql_p\"";
				$sth->execute($acc,$pep) or die "could not execute \"$sql\"";
			}
			my @cs = $sth->fetchrow_array();
			my $total;
			$return .= "{\"$pep\":[";
			my $c = 0;
			foreach $total(@totals)	{
				if($total)	{
					$return .= sprintf("%f,",@cs[$c]/$total);
				}
				else	{
					$return .= "undef,";
				}
				$c++;
			}
			$return =~ s/\,$//;
			$return .= "]},\n";
			$sth->finish();
		}
	}
	$return =~ s/\,(\s*)$/$1/;
	return "$return]}";
}


#input: protein accession number
#return: tuple of the individual sums of all peptide observations with (z1,z2,z3)
sub gpmdbProteinPeptidesCharges {
	my ($class,$acc) = @_;                     
	my $dbh = login($class,"p");
	my $sql = "select sum(z1total),sum(z2total),sum(z3total),sum(z4total) from protein_omega_count where label=?";
	my $sth;
	if ($acc =~ /HIT[0-9]+/ || $acc  =~ /IPI[0-9]+/ || $acc  =~  /At\dg\d+/ || $acc  =~  /^ATCG\d+/i) {
		$sql = "select sum(z1total),sum(z2total),sum(z3total),sum(z4total) from protein_omega_count where label like ?";
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute("$acc%") or die "could not execute \"$sql\"";
	}
	else {
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute($acc) or die "could not execute \"$sql\"";
	}
	my @cs = $sth->fetchrow_array();
	my $c;
	my $return;
	my $a = 1;
	foreach $c(@cs)	{
		$return .= sprintf("\"%i\":%i,",$a,$c);
		$a++;
	}
	$return =~ s/\,$//;
	$sth->finish();
	return "{$return}";
}

#input: nsSNV asscession rsid
#return: information about nsSNV consequences
sub gpmdbNssnv {
	my ($class,$acc) = @_;                     
	my $dbh = login($class,"e");
	my $sql = "select * from rsid where rsid=?";
	my $sth;
	my $nacc = $acc;
	my %aas = ("A" => 71.037114,"R" => 156.101111,"N" => 114.042927,"D" => 115.026943,"C" => 103.009185,"E" => 129.042593,"Q" => 128.058578,"G" => 57.021464,"H" => 137.058912,"I" => 113.084064,"L" => 113.084064,"K" => 128.094963,"M" => 131.040485,"F" => 147.068414,"P" => 97.052764,"S" => 87.032028,"T" => 101.047679,"U" => 150.95363,"W" => 186.079313,"Y" => 163.06332,"V" => 99.068414);
	$nacc =~ s/^rs//;
	$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
	$sth->execute($nacc) or die "could not execute \"$sql\"";
	my @cs = $sth->fetchrow_array();
	my $c;
	my $return;
	my $a = 1;
	my $hgvs = "";
	my $hgvs_p = "";
	my $dm = 0.0;
	while(@cs)	{
		$hgvs = sprintf("%s:g.%i%s>%s",@cs[10],@cs[8],@cs[7],@cs[9]);
		$hgvs_p = sprintf("%s:p.%s%i%s",@cs[3],@cs[4],@cs[5],@cs[6]);
		$dm = $aas{@cs[6]} - $aas{@cs[4]};
		$return .= sprintf("{\"rsid\":\"rs%i\",\"af\":%.1e,\"label\":\"%s\",\"ref\":\"%s\",\"at\":%i,\"var\":\"%s\",\"hgvs:g\":\"%s\",\"hgvs:p\":\"%s\",\"dm\":%.4f},",@cs[1],@cs[2],@cs[3],@cs[4],@cs[5],@cs[6],$hgvs,$hgvs_p,$dm);
		@cs = $sth->fetchrow_array();
	}
	$return =~ s/\,$//;
	$sth->finish();
	$return =~ s/\"\:/\":    \t/g;
	$return =~ s/\,/,\n\t/g;
	$return =~ s/\{/\n{\t/g;
	return "[\t$return\n]";
}

#input: protein accession number
#return: total number of peptide observations
sub gpmdbProteinPeptidesTotal {
	my ($class,$acc) = @_;                     
	my $dbh = login($class,"p");
	my $sql = "select sum(z1total+z2total+z3total+z4total) from protein_omega_count where label=?";
	my $sth;
	if ($acc =~ /HIT[0-9]+/ || $acc  =~ /IPI[0-9]+/ || $acc  =~  /At\dg\d+/ || $acc  =~  /^ATCG\d+/i) {
		$sql = "select sum(z1total+z2total+z3total+z4total) from protein_omega_count where label like ?";
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute("$acc%") or die "could not execute \"$sql\"";
	}
	else {
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute($acc) or die "could not execute \"$sql\"";
	}
	my $c = $sth->fetchrow_array();
	$sth->finish();
	return "[$c]";
}

#input: protein accession number,peptide sequence
#return: tuple of the observations of a peptide with (z1,z2,z3)
sub gpmdbProteinPeptideCount {
	my ($class,$acc,$seq) = @_;                     

	my $dbh = login($class,"p");
	my $sql = "select z1total,z2total,z3total,z4total from protein_omega_count where label=? and seq=?";
	my $sth;
	if ($acc =~ /HIT[0-9]+/ || $acc  =~ /IPI[0-9]+/ || $acc  =~  /At\dg\d+/ || $acc  =~  /^ATCG\d+/i) {
		$sql = "select z1total,z2total,z3total,z4total from protein_omega_count where label like ? and seq=?";
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute("$acc%",$seq) or die "could not execute \"$sql\"";
	}
	else {
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute($acc,$seq) or die "could not execute \"$sql\"";
	}
	my @cs = $sth->fetchrow_array();
	my $c;
	my $return;
	my $a = 1;
	if(scalar(@cs))	{
		foreach $c(@cs)	{
			$return .= sprintf("\"%i\":%i,",$a,$c);
			$a++;
		}
	}
	else	{
		$return = "\"1\":0,\"2\":0,\"3\":0,\"4\":0";
	}
	$return =~ s/\,$//;
	$sth->finish();
	return "{$return}";
}

#input: protein accession number
#return: tuple of all peptide sequences observed (seq1,seq2,...,seqN)
sub gpmdbProteinPeptideSequences {
	my ($class,$acc) = @_;                     
	my $dbh = login($class,"p");
	my $sql = "select seq from protein_omega_count where label=? and z1total+z2total+z3total+z4total > 0";
	my $sth;
	if ($acc =~ /HIT[0-9]+/ || $acc  =~ /IPI[0-9]+/ || $acc  =~  /At\dg\d+/ || $acc  =~  /^ATCG\d+/i) {
		$sql = "select seq from protein_omega_count where label like ? and z1total+z2total+z3total+z4total > 0";
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute("$acc%") or die "could not execute \"$sql\"";
	}
	else {
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute($acc) or die "could not execute \"$sql\"";
	}
	my $c = $sth->fetchrow_array();
	my $return;
	while($c)	{
		$return .= "\t\"$c\",\n";
		$c = $sth->fetchrow_array();
	}
	$return =~ s/\,(\s*)$/$1/;
	$sth->finish();
	return "[$return]";
}

#input: protein accession number
#return: total number of protein observations
sub gpmdbProteinCount {
	my ($class,$acc) = @_;                     
	my $dbh = login($class,"g");
	my $sql = "select count from best_expect where label=?";
	my $sth;
	if ($acc =~ /HIT[0-9]+/ || $acc  =~ /IPI[0-9]+/ || $acc  =~  /At\dg\d+/ || $acc  =~  /^ATCG\d+/i) {
		$sql = "select count from best_expect where label like ?";
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute("$acc%") or die "could not execute \"$sql\"";
	}
	else {
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute("$acc") or die "could not execute \"$sql\"";
	}
	my ($c) = $sth->fetchrow_array();
	$sth->finish();
	return "[$c]";
}

#input: protein accession number
#return: lowest log(E) value for the protein specified
sub gpmdbProteinBestExpect {
	my ($class,$acc) = @_;                     
	my $dbh = login($class,"g");
	my $sql = "select expect from best_expect where label=?";
	my $sth;
	if ($acc =~ /HIT[0-9]+/ || $acc  =~ /IPI[0-9]+/ || $acc  =~  /At\dg\d+/ || $acc  =~  /^ATCG\d+/i) {
		$sql = "select expect from best_expect where label like ?";
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute("$acc%") or die "could not execute \"$sql\"";
	}
	else {
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute($acc) or die "could not execute \"$sql\"";
	}
	my ($c) = $sth->fetchrow_array();
	$sth->finish();
	return "[$c]";
}

#input: protein accession number
#return: amino acid polymorphisms
sub gpmdbProteinPolymorphisms {
	my ($class,$acc) = @_;
	my $dbh = login($class,"g");
	my $sql = "select hgvs from peptide_mut where hgvs regexp ?";
	my $sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
	$sth->execute("^$acc\:") or die "could not execute \"$sql\"";
	my %hgvs;
	my @cs = $sth->fetchrow_array();
	while(@cs)	{
		$hgvs{@cs[0]} = $hgvs{@cs[0]} + 1;
		@cs = $sth->fetchrow_array();
	}
	my @keys = keys(%hgvs);
	my $k;
	my $return = "[";
	foreach $k(@keys)	{
		if(gpmdbSnapCheck($class,$k))	{
			$return .= "\t{\"$k\":$hgvs{$k}},\n";
		}
	}
	$return =~ s/\,(\s*)$/$1/;
	$sth->finish();
	$return .= "]";
	return $return;
}

sub gpmdbSnapCheck {
	my ($class,$snap) = @_;
	my ($one,$two) = $snap =~ /\:p\.([A-Z])\d+([A-Z])/;
	if($one eq 'M' and $two eq 'F')	{
		return 0;
	}
	elsif($one eq 'D' and $two eq 'H')	{
		return 0;
	}
	elsif($one =~ /DN/ and $two =~  /DN/)	{
		return 0;
	}
	elsif($one =~ /EQK/ and $two =~ /EQK/)	{
		return 0;
	}
	return 1;
}

#input: protein accession number
#return: sequence of the protein specified
sub gpmdbProteinSequence {
	my ($class,$acc) = @_;
	my $dbh = login($class,"g");
	my $sql = "select seq from ProSeq where";
	my $sth;
	if($acc =~ /HIT[0-9]+/ || $acc =~ /IPI[0-9]+/ || $acc =~  /At[\dC]g\d+/i)	{
		$sql .= ' label like ?';
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute("$acc%") or die "could not execute \"$sql\"";
	} 
	else	{
		$sql .= ' label = ?';
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute($acc) or die "could not execute \"$sql\"";
	}
	my ($c) = $sth->fetchrow_array();
	$sth->finish();
	return "[\"$c\"]";
}

#input: protein accession number
#return: description of the protein specified
sub gpmdbProteinDescription {
	my ($class,$acc) = @_;                     
	my $dbh = login($class,"e");
	my $sql = "select description from Map where";
	my $sth;
	if($acc =~ /HIT[0-9]+/ || $acc =~ /IPI[0-9]+/ || $acc =~  /At[\dC]g\d+/i)	{
		$sql .= ' label like ?';
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute("$acc%") or die "could not execute \"$sql\"";
	} 
	else	{
		$sql .= ' label = ?';
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute($acc) or die "could not execute \"$sql\"";
	}
	my $sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
	$sth->execute($acc) or die "could not execute \"$sql\"";
	my ($c) = $sth->fetchrow_array();
	$sth->finish();
	return "[\"$c\"]";
}

#input: protein accession number
#return: best guess as the species associated with that accession number
sub gpmdbProteinSpecies {
	my ($class,$acc) = @_;
	my $species = "";
	my @accs = split /\,/,$acc;
	my $a;
	my %key_values = (	"AGAP"	=> 	"Anopheles gambiae",
		"AT"	=>	"Arabidopsis thaliana",
		"Bra"	=>	"Brassica rapa",
		"ATCG"	=>	"Arabidopsis thaliana",
		"ATMG"	=>	"Arabidopsis thaliana",
		"LOC_Os"	=>	"Oryza sativa",
		"LOC_Osp"	=>	"Oryza sativa",
		"LOC_Osm"	=>	"Oryza sativa",
		"FBpp" => "Drosophila melanogaster",
		"DappuP" => "Daphnia pulex",
		"PPA"  =>	"Pristionchus pacificus",
		"BRADI"	=>	"Brachypodium distachyon",
		"AAZ"	=>	"Trypanosoma brucei",
		"AAQ"	=>	"Trypanosoma brucei",
		"CAJ"	=>	"Trypanosoma brucei",
		"EAN"	=>	"Trypanosoma brucei",
		"GRMZM"	=>	"Zea mays",
		"SPU_"	=>	"Strongylocentrotus purpuratus",
		"ACYPI"	=>	"Acyrthosiphon pisum",
		"CADAFUAP" =>	"Aspergillus fumigatus",
		"ENSP" => "Homo sapiens",
		"ENSSSCP" => "Sus scrofa",
		"ENSMMUP" => "Macaca mulatta",
		"ENSMODP" => "Monodelphis domestica",
		"ENSCAFP" => "Canis familiaris",
		"ENSMUSP" => "Mus musculus",
		"ENSXETP" => "Xenopus tropicalis",
		"ENSGALP" => "Gallus gallus",
		"ENSRNOP" => "Rattus norvegicus",
		"ENSANGP" => "Anopheles gambiae",
		"AGAP" => "Anopheles gambiae",
		"DappuP" => "Daphnia pulex",
		"ENSDARP" => "Danio rerio",
		"ENSBTAP" => "Bos taurus",
		"ENSAPMP" => "Apis mellifera",
		"ENSOCUP" => "Oryctolagus cuniculus",
		"ENSCPOG" => "Cavia porcellus",
		"ENSCPOP" => "Cavia porcellus",
		"NEWSINFRUP" => "Fugu rubripes",
		"ENSPTRP" => "Pan troglodytes",
		"ENSFCAP" => "Felis catus",
		"GSTENP"  => "Tetraodon nigroviridis",
		"ENSTNIP" => "Tetraodon nigroviridis",
		"ENSTRUP" => "Takifugu rubripes",
		"ENSECAP" => "Equus caballus",
		"ENSMEUP" => "Macropus eugenii",
		"ENSCINP" => "Ciona intestinalis",
		"FBpp" => "Drosophila melanogaster",
		"ENSACAP" => "Anolis carolinensis",
		"ENSTGUP" => "Taeniopygia guttata",
		"ENSLAFP" => "Loxodonta africana",
		"ENSMGAP" => "Meleagris gallopavo",
		"ENSGACP"=> "Gasterosteus aculeatus" );
	foreach $a(@accs)	{
		my $v = "unknown";
		my ($r) = $a =~ /(.+?)\d/;
		if($key_values{$r})	{
			$v = $key_values{$r};
		}
		elsif($a =~ /^AT[\dC]G/i)	{
			$v = "Arabidopsis thaliana";
		}
		elsif($r =~ /^Y[A-Z][LR]/)	{
			$v = "Saccharomyces cerevisiae";
		}
		elsif($r =~ /^HIT/)	{
			$v = "Homo sapiens";
		}
		elsif($r =~ /^HIP/)	{
			$v = "Homo sapiens";
		}
		# deal with SwissProt accessions
		elsif($a =~ /^sp\|/)	{
			$v = "";
			if($a =~ /\_HUMAN/)	{
				$v = "Homo sapiens";
			}
			elsif($a =~ /\_MOUSE/)	{
				$v = "Mus musculus";
			}
			elsif($a =~ /\_RAT/)	{
				$v = "Rattus norvegicus";
			}
			elsif($a =~ /\_CHICK/)	{
				$v = "Gallus gallus";
			}
			elsif($a =~ /\_BOVIN/)	{
				$v = "Bos taurus";
			}
			elsif($a =~ /\_SHEEP/)	{
				$v = "Ovis aries";
			}
			elsif($a =~ /\_RABIT/)	{
				$v = "Oryctolagus cuniculus";
			}
		}
		elsif($a =~ /^gi\|/)	{
			my $d = gpmdbProteinDescription($class,$a);
			($v) = $d =~ /.+\[(.+?)\]/;		
		}
		$species .= "\t\"$a\":\"$v\",\n";
	}
	$species =~ s/\,(\s)$/$1/s;   
	return "{$species}";
}

#input: peptide sequence
#return: tuple of the protein accessions and individual sums of all peptide observations with (z1,z2,z3,z4)
sub gpmdbPeptideAccessions {
	my ($class,$seq) = @_;                     
	my $dbh = login($class,"p");
	my $sql = "select label,z1total+z2total+z3total+z4total from protein_omega_count where seq=? and z1total+z2total+z3total+z4total>0";
	my $sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
	$sth->execute($seq) or die "could not execute \"$sql\"";
	my $return;
	my @cs;
	while(@cs = $sth->fetchrow_array())	{
		$return .= sprintf("\t\"%s\"\t:\t%i,\n",@cs[0],@cs[1]);
	}
	$return =~ s/\,$//;
	$sth->finish();
	return "{\n$return}";
}

#input: peptide sequence
#return: the total number of observations of the sequence
sub gpmdbPeptideTotal {
	my ($class,$seq) = @_;                     
	my $dbh = login($class,"p");
	my $sql = "select sum(z1total+z2total+z3total+z4total) from protein_omega_count where seq=? and z1total+z2total+z3total+z4total > 0";
	my $sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
	$sth->execute($seq) or die "could not execute \"$sql\"";
	my $return = 0;
	my @cs;
	while(@cs = $sth->fetchrow_array())	{
		$return += @cs[0];
	}
	$return = sprintf("%i",$return);
	$sth->finish();
	return "[$return]";
}

#input: peptide sequence
#return: a tuple of the number of observations of the sequence for (z1,z2,z3,z4)
sub gpmdbPeptideTotalZ {
	my ($class,$seq) = @_;                     

	my $dbh = login($class,"p");
	my $sql = "select z1total,z2total,z3total,z4total from protein_omega_count where seq=? and z1total+z2total+z3total+z4total>0";
	my $sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
	$sth->execute($seq) or die "could not execute \"$sql\"";
	my $return;
	my %z;
	my @cs;
	while(@cs = $sth->fetchrow_array())	{
		$z{1} = $z{1} + @cs[0];
		$z{2} = $z{2} + @cs[1];
		$z{3} = $z{3} + @cs[2];
		$z{4} = $z{4} + @cs[3];
	}
	$return = sprintf("\"1\":%i,\"2\":%i,\"3\":%i,\"4\":%i",$z{1},$z{2},$z{3},$z{4});
	$sth->finish();
	return "{$return}";
}

#input: protein accession, list of potentially phosphorylated residues
#return: frequency of observation
sub gpmdbPeptidePf {
	my ($class,$acc,$pos,$w,$t) = @_;                  
	my @ps = ();
	my $a;
	if($pos =~ /\d\-\d/)	{
		my @as = split(/-/,$pos);
		if($as[1] < $as[0])	{
			@as = reverse(@as);
		}
		$a = $as[0];
		while($a < $as[1])	{
			push(@ps,$a);
			$a++;
		}
		push(@ps,$a);
	}
	else	{
		@ps = split(/,/,$pos);
	}
	my $dbh = login($class,"g");
	my %freq = {};
	my $total = 0;
	my $return = "";
	my $sql = "select proseqid from proseq where label=?";
	my $sth;
	my @accs = get_accs($class,$acc);
	my @ids;
	my $id;
	push(@accs,$acc);
	for $a(@accs)	{
		$sth = $dbh->prepare($sql);
		$sth->execute($a);
		while(($id) = $sth->fetchrow_array())	{
			push(@ids,$id);
		}
	}
	if(scalar(@ids) == 0)	{
		$return = "";
		for $a(@ps)	{
			$return .= "\"$a\":null,";
		}
		$return =~ s/\,$//;
		$sth->finish();
		return "{$return}";
	}
	$sth->finish();
	$sql = "select freq from phos_freq where proseqid=? and at=?";
	if($t eq 'a')	{
		$sql = "select freq from acetyl_freq where proseqid=? and at=?";
	}
	elsif($t eq 'u')	{
		$sql = "select freq from ubi_freq where proseqid=? and at=?";
	}
	elsif($t eq 'd')	{
		$sql = "select freq from dim_freq where proseqid=? and at=?";
	}
	elsif($t eq 'o')	{
		$sql = "select freq from oxy_freq where proseqid=? and at=?";
	}
	elsif($t eq 's')	{
		$sql = "select freq from sumo_freq where proseqid=? and at=?";
	}
	elsif($t eq 'nq')	{
		$sql = "select freq from deam_freq where proseqid=? and at=?";
	}
	elsif($t eq 'ol')	{
		$sql = "select freq from ogly_freq where proseqid=? and at=?";
	}
	elsif($t eq 'qc')	{
		$sql = "select freq from qcyc_freq where proseqid=? and at=?";
	}
	elsif($t eq 'ct')	{
		$sql = "select freq from citr_freq where proseqid=? and at=?";
	}
	$sth = $dbh->prepare($sql);
	my $f;
	my $b_total = 0;
	my %b_freq = ();
	for $id(@ids)	{
		$total = 0;
		%freq = ();
		for $a(@ps)	{
			$sth->execute($id,$a);
			($f) = $sth->fetchrow_array();
			if($f > 0)	{
				$total += $f;
				$freq{$a} = $f;
			}
			else	{
				$freq{$a} = 0;
			}
		}
		if($total > $b_total)	{
			$b_total = $total;
			%b_freq = %freq;
		}
	}
	if($b_total == 0)	{
		$return = "";
		for $a(@ps)	{
			if($w eq 'n')	{
				$return .= "\"$a\":0,";
			}
			else	{
				$return .= "\"$a\":null,";
			}
		}
		$return =~ s/\,$//;
		$sth->finish();
		return "{$return}";
	}
	$return = "";
	for $a(@ps)	{
		if($w eq 'n')	{
			$return .= sprintf("\"%i\":%i,",$a,$b_freq{$a});
		}
		else	{
			$return .= sprintf("\"%i\":%.3f,",$a,$b_freq{$a}/$b_total);
		}
	}
	$return =~ s/\,$//;
	$sth->finish();
	return "{$return}";
}

sub get_accs	{
	my ($class,$acc) = @_;
	my $dbh = login($class,"e");
	my $sql = "select ens from up_to_ens where up=?";
	my $sth = $dbh->prepare($sql);
	$sth->execute($acc);
	my @accs = ();
	my $a;
	while(($a) = $sth->fetchrow_array())	{
		push(@accs,$a);
	}
	if(scalar(@accs) == 0)	{
		$acc =~ s/\-\d+//;
		$sth->execute($acc);
		while(($a) = $sth->fetchrow_array())	{
			push(@accs,$a);
		}
	}
	if(scalar(@accs) == 0)	{
		push(@accs,$acc);
	}
	$sth->finish();
	return @accs;         
}

#input: protein accession number
#return: has containing the relavent protein evidence code information
sub gpmdbProteinEvidence {
	my ($class,$acc) = @_;                     
	my ($_acc) = @_;
	my $count = gpmdbProteinCount($class,$acc);
	$count =~ s/[\[\"\]]//g;
	my $best = gpmdbProteinBestExpect($class,$acc);
	my $rating = "black";
	my %code = ("black" => 1, "red" => 2, "yellow" => 3, "green" => 4);
	my %hues = ("black" => "#CCCCCC", "red" => "#FFB5B5", "yellow" => "#FFFF66", "green"  => "#B3FF99");
	my %means = ("black" => "no credible evidence of translation", "red" => "poor evidence of translation", "yellow" => "modest quality evidence of translation", "green"  => "good evidence of translation");
	
	$best =~ s/[\[\"\]]//g;
	my $min_e = -20.5;
	my $min_count = 500;
	if(not $count)	{
		$rating = 'black';
	}
	elsif($best < $min_e and $count > $min_count)	{
		$rating = 'green';
	}
	else	{
		my $ret = gpmdbProteinModel($class,$acc,1);
		my $ref = decode_json($ret);
		if($$ref{'green'})	{
			$rating = 'green';
		}
		elsif($$ref{'yellow'})	{
			$rating = 'yellow';
		}
		elsif($$ref{'red'})	{
			$rating = 'red';
		}
		else	{
			$rating = 'black';
		}
	}
	my $return = qq({"code":"$rating","level":"$code{$rating}","text":"$means{$rating}"});
	return $return;
}

sub gpmdbPsytPredict
{
	my ($class,$_seq,$_w) = @_;
	my $seq = $_seq;
	$seq =~ tr/[a-z]/[A-Z]/;
	$seq =~ s/[^A-Z]//g;
	if(not $seq =~ /[STY]/)	{
		return "{}";  
	}
	my (@s) = split('',$seq);
	my $length = scalar(@s);
	my $a = 0;
	my $return = "{\n";
	while($a < $length)	{
		if(@s[$a] =~ /[STY]/)	{
			my $t;
			my $b = 0;
			while($b < $length)	{
				if($b == $a)	{
					if(@s[$b] eq 'S')	{
						$t .= 's';
					}
					elsif(@s[$b] eq 'T')	{
						$t .= 't';
					}
					elsif(@s[$b] eq 'Y')	{
						$t .= 'y';
					}
				}
				else	{
					$t .= @s[$b];
				}
				$b++;
			}
			my $ret = gpmdbPsytMotif($class,$t,$_w);
			$ret =~ s/[\[\]]//g;
			$return .= sprintf("%i:%i,\n",$a+1,$ret);
		}
		$a++;
	}
	$return =~ s/\,$//;
	$return .= '}';
	return $return;			
}

sub gpmdbPsytMotif
{
	my ($class,$_seq,$_w) = @_; 
	if(not $_seq =~ /[sty]/)	{
		return "[\"error: no site in sequence $_seq\"]";  
	}
	my ($f,$res,$l) = $_seq =~ /(\w+)([sty])(\w+)/;
	$f =~ tr/[a-z]/[A-Z/;
	$l =~ tr/[a-z]/[A-Z/;
	my @first = reverse(split('',$f));
	my @last = split('',$l);
	my $proline = @last[0];
	if($res eq 'y')	{
		$proline = "";
	}
	my $file = "pY_all.tsv";
	if($res eq 's' and $proline eq 'P')	{
		$file = "pS_+P.tsv";
	}
	elsif($res eq 's' and $proline ne 'P')	{
		$file = "pS_-P.tsv";
	}
	elsif($res eq 't' and $proline eq 'P')	{
		$file = "pT_+P.tsv";
	}
	elsif($res eq 't' and $proline ne 'P')	{
		$file = "pT_-P.tsv";
	}
	open(IN,"<psyt/$file") or return "[\"error: could not open parameter file psyt/$file\"]";
	my %factors;
	$_ = <IN>;
	while(<IN>)	{
		chomp($_);
		my @vs = split /\t/,$_;
		$factors{@vs[0]} = \@vs;
		if(@vs[0] eq 'Y')	{
			last;
		}
	}
	close(IN);
	my $score = 0;
	$a = 0;
	my $ref;
	my $count = 0;
	my $width = $_w;
	if(not $width)	{
		$width = 3;
	}
	if($width > 9)	{
		$width = 9;
	}
	if($width < 1)	{
		$width = 1;
	}
	while($a < $width)	{
		if($a == 0 and $proline eq 'P')	{
			$score += 250;
			$count++;
			if(@first[$a])	{
				$ref = $factors{@first[$a]};
				$score += @$ref[10-$a];
				$count++;
			}
		}
		else	{
			if(@last[$a])	{
				$ref = $factors{@last[$a]};
				$score += @$ref[11+$a];
				$count++;
			}
			if(@first[$a])	{
				$ref = $factors{@first[$a]};
				$score += @$ref[10-$a];
				$count++;
			}
		}
		$a++;
	}
	if($count == 0)	{
		$count = 1;
	}
	my $return = sprintf("[%.0f]",$score);
	return $return;
}

#input: protein accession number
#return: rating for existing protein observations
sub gpmdbProteinModel {
	my ($class,$acc,$type) = @_;                     
	my $dbh = login($class,"g");
	my $sql = "SELECT r.file,proid,pathid,uid,expect,seq,gnotes.notes FROM Protein p, Result r, ProSeq ps, enspmapdb.gpmnotes gnotes WHERE r.resultid = p.resultid AND p.proseqid = ps.proseqid AND r.file = gnotes.file and active = 1";
	my $sth;
	if ($acc =~ /HIT[0-9]+/ || $acc  =~ /IPI[0-9]+/ || $acc  =~  /At\dg\d+/ || $acc  =~  /^ATCG\d+/i) {
		$sql .= ' and label like ? order by expect';
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute("$acc%") or die "could not execute \"$sql\"";
	}
	else	{
		$sql .= ' and label = ? order by expect';
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute($acc) or die "could not execute \"$sql\"";
	}
	my @clean_results;
	my @rs;
	my $rows = 0;
	while(@rs=$sth->fetchrow_array())	{
		push(@clean_results,[@rs]);
		$rows++;
	}
	$sth->finish();
	$sql = "SELECT protocol,server,localpath,relpath FROM Paths WHERE pathid = ?";
	my $sql3 = "SELECT start,end,expect from Peptide WHERE proid = ?";
	my $sth3 = $dbh->prepare($sql3);

	my $protein_seq = gpmdbProteinSequence($class,$acc);
	my $protein_length = length($protein_seq)+1;

	my %peps;
	my $result_counter = 0;
	while($result_counter < $rows) {
		my $r = @clean_results[$result_counter];
		my @result1 = @$r;
		$result_counter++;
	
		my $sequence = @result1[5];
		my $sth2 = $dbh->prepare($sql);
		$sth2->execute(@result1[2]);
		my @result2=$sth2->fetchrow_array();
		$sth3->execute(@result1[1]);

		my $loge;
		my @result3;
		while(@result3=$sth3->fetchrow_array()){
			if($result3[2] == 0.0)	{
				next;
			}
			$loge = int(log($result3[2])/2.3026);
			if($peps{"@result3[0] @result3[1]"})	{
				my $ref = $peps{"@result3[0] @result3[1]"};
				$$ref{$loge} = $$ref{$loge} + 1;
			}
			else	{
				my %es;
				$es{$loge} = 1;
				$peps{"@result3[0] @result3[1]"} = \%es;
			}
		}
		$sth2->finish();
	}
	my %score;
	my @keys = keys(%peps);
	my %skews;
	my %kertosi;
	my %wags;
	my %scores;
	my %types;
	my $kmax = 1.5;
	my $smax = 1.5;
	my @vs;
	my $limit = 4;
	my $min_length = 5;
#open(LOG,">a.log");
	foreach $_(@keys)	{ 
		if($type and $types{'green'})	{
			last;
		}
		my @ls = split / /,$_;
		if(@ls[0] > $protein_length or @ls[1] > $protein_length)	{
			next;
		}
		my $r = $peps{$_};
		my $o = -2;
		my @vs;
		my $t = 0;
		my $l = @ls[1]-@ls[0]+1;
		while($o >= -14)	{
			if($$r{$o})	{
				$score{$_} = $score{$_} + -1*$o;
			}
			push(@vs,$$r{$o});
			$t += $$r{$o};
			$o -= 1;
		}
		if($t)	{
			my $s = 0;
			my $sum = 0;
			foreach $s(@vs)	{
				$sum += $s;
			}
			$skews{$_} = skew(\@vs);
			$kertosi{$_} = kertosis(\@vs);
			$wags{$_} = wag(\@vs,-2.0);
			if($l > $min_length and $sum > $limit and (($skews{$_} <= $smax and $kertosi{$_} <= $kmax) or $wags{$_} < -5.5))	{
				$scores{$_} = 1;
				$types{'green'} = $types{'green'} + 1;
#print LOG "$_ green $l and $sum and $skews{$_} and $kertosi{$_}\n";
			}
			elsif($l > $min_length and $sum > $limit and ($skews{$_} <= $smax or $kertosi{$_} <= $kmax or $wags{$_} < -3.5))	{
				$scores{$_} = 2;
				$types{'yellow'} = $types{'yellow'} + 1;
#print LOG "$_ yellow $l and $sum and $skews{$_} and $kertosi{$_}\n";
			}
			else	{
				$scores{$_} = 3;
				$types{'red'} = $types{'red'} + 1;
#print LOG "$_ red $l and $sum and $skews{$_} and $kertosi{$_}\n";
			}
		}
	}
#close(LOG);
	$sth->finish();
	$sth3->finish();
	return sprintf("{\"green\":%i,\"yellow\":%i,\"red\":%i}",$types{'green'},$types{'yellow'},$types{'red'});
}

#input: gpmdb accession number
#return: textblock containing for the metadata
sub gpmDbMetadata {
	my ($class,$acc) = @_;                     
	my $dbh = login($class,"e");
	my $sql = "select notes from gpmnotes where file = ?";
	my $sth;
	if ($acc =~ /^GPM\d+/) {
		if($acc !~ /\.xml/i)	{
			$acc = $acc . ".xml";
		}
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute("$acc") or die "could not execute \"$sql\"";
	}
	else {
		$sth = $dbh->prepare($sql) or die "could not prepare \"$sql\"";
		$sth->execute($acc) or die "could not execute \"$sql\"";
	}
	my $c = $sth->fetchrow_array();
	my $return;
	while($c)	{
		$return .= $c;
		$c = $sth->fetchrow_array();
	}
	$sth->finish();
	my (@lines) = $return =~ /(.+?\n)/sg;
	my $l;
	my %entries;
	my $key;
	my $value;
	foreach $l(@lines)	{
		chomp($l);
		if($l =~ /.+\:/)	{
			($key) = $l =~ /(.+?)\:/;
			($value) = $l =~ /.+?\:(.+)/;
			$key =~ s/^ +//;
			$value =~ s/^ +//;
			$value =~ s/[\r\n]//sg;
			$value =~ s/\"/\\\"/sg;
			if($key)	{
				$entries{$key} = $value;
			}
			$value = "";
		}
		else	{
			if($key and $entries{$key})	{
				$entries{$key} .= " $value";
			}
		}
	}
	my @keys = keys(%entries);
	$return = "";
	foreach $key(@keys)	{
		$return .= "\"$key\":\"$entries{$key}\",\n";
	}
	$return =~ s/\,$//;
	return "{$return}";
}

sub stdev
{
	my ($_v) = @_;
	my @vs = @$_v;
	my $ave = 0;
	my $v;
	foreach $v(@vs)	{
		$ave += $v;
	}
	my $n = scalar(@vs);
	$ave /= $n;
	my $std = 0;
	foreach $v(@vs)	{
		$std += ($v - $ave)*($v-$ave);
	}
	$std /= $n - 1;
	return ($ave,sqrt($std));
}

sub skew
{
	my ($_v) = @_;
	my @vs = @$_v;
	my ($ave,$std) = stdev($_v);
	my $t = 0;
	my $z;
	my $v;
	foreach $v(@vs)	{
		$z = ($v - $ave)/$std;
		$t += $z*$z*$z;
	}
	my $n = scalar(@vs);
	return $t * $n/(($n-1)*($n-2));
#	return $t / ($n-1);
}

sub wag
{
	my ($_v,$_s) = @_;
	my @vs = @$_v;
	my $t = 0;
	my $z = $_s;
	my $v = 0;
	my $sum = 0;
	foreach $v(@vs)	{
		$t += $z*$v;
		$sum += $v;
		$z -= 1;
	}
	if($sum > 0)	{
		$t = $t/$sum;
	}
	return $t;
}

sub kertosis
{
	my ($_v) = @_;
	my @vs = @$_v;
	my ($ave,$std) = stdev($_v);
	my $t = 0;
	my $z;
	my $v;
	foreach $v(@vs)	{
		$z = ($v - $ave)/$std;
		$t += $z*$z*$z*$z;
	}
	my $n = scalar(@vs);
	$t = $t * $n*($n+1)/(($n-1)*($n-2)*($n-3));
	return $t - (3*($n-1)*($n-1)/(($n-2)*($n-3)));
#	return $t/($n-1);
}

#input: none
#output: list of methods in GPMDB REST 1 interface
sub gpmdbHelp
{
	my ($class) = @_;
	my $return = qq(["
GET /1/interface/version => [ARRAY - string]
GET /1/interface/help => [ARRAY - string]

GET /1/model/metadata/gpm=GPM => {OBJECT}
GET /1/model/proteins/gpm=GPM => [ARRAY]
GET /1/model/protein_modifications/gpm=GPM&acc=ACC => [ARRAY]
GET /1/model/protein_peptides/gpm=GPM&acc=ACC => {OBJECT}
GET /1/model/protein_savs/gpm=GPM&acc=ACC => [ARRAY]
GET /1/model/protein_sequence/gpm=GPM&acc=ACC => [ARRAY]

GET /1/peptide/accessions/seq=SEQ => {OBJECT - string:int}
GET /1/peptide/count/seq=SEQ => [ARRAY - int]
GET /1/peptide/count_z/seq=SEQ => {OBJECT - string:int}

GET /1/protein/best_e/acc=ACC => [ARRAY - float]
GET /1/protein/count/acc=ACC => [ARRAY - int]
GET /1/protein/description/acc=ACC => [ARRAY - string]
GET /1/protein/evidence/acc=ACC => {OBJECT - string:string}
GET /1/protein/keyword/key=KEY => [ARRAY - OBJECT {string:string}]
GET /1/protein/modifications/acc=ACC&mod=MOD&res=RES&maxe=MAXE => {OBJECT - string:int}
GET /1/protein/omega/acc=ACC&seq=SEQ => {OBJECT - string:[{string:[float]}]}
GET /1/protein/peptide_count/acc=ACC&seq=SEQ => {OBJECT - string:int}
GET /1/protein/peptide_sequences/acc=ACC => [ARRAY - string]
GET /1/protein/peptides_z/acc=ACC => {OBJECT - string:int}
GET /1/protein/peptides_total/acc=ACC => [ARRAY - int]
GET /1/protein/savs/acc=ACC => [ARRAY - OBJECT - string:int]
GET /1/protein/sequence/acc=ACC => [ARRAY - string]
GET /1/everything/else => 404

acc=accession number
maxe=log10(maximum allowed expectation value)
mod=modification mass (Da)
res=list of residues
seq=peptide sequence
key=keyword

documentation at http://wiki.thegpm.org/wiki/GPMDB_REST
"]);
	return $return;
}

#input: none
#return: current interface version
sub gpmdbVersion	{
	my ($class) = @_;   
	return "[\"2015.02.19\"]";
}

1;

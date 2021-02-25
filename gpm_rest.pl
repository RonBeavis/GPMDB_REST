#
# GPM REST interface
# Copyright (C) 2012 Ronald C. Beavis, All rights reserved
#
# Use of this software governed by the Artistic license,
# as reproduced at http://www.opensource.org/licenses/artistic-license.php
#  

package GPM;

use DBI;

use strict;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [ qw() ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw( );

our $VERSION = '1.0';

sub new
{
	my ($self,$_g,$_u,$_p) = @_;
	my $dsn = "DBI:mysql:host=192.168.1.5;database=gpmdb";
	my $dbh_g = DBI->connect($dsn,$_u,$_p,
			{PrintError => 0, RaiseError => 1}) or die "could not connect to gpmdb";

	my $self = {	_gpm=>$_g,
			_path=>"",
			_mess=>"",
			_user=>$_u,
			_pass=>$_p,
			_dbh_g=>$dbh_g};

	bless $self, 'GPM';
	return $self;
}

sub DESTROY
{
	my ($class) = @_;
	my $g = $class->{_dbh_g};
	$g->disconnect();
}

sub getPath
{
	my ($class) = @_;
	if($class->{_path})	{
		return;
	}
	my $gpm = $class->{_gpm};
	my ($dir) = $gpm =~ /GPM(\d\d\d)/;
	my $path = "..\\gpm\\archive\\$dir";
	my $mess;
	if(-e "$path\\$gpm.xml")	{
		$path = "$path\\$gpm.xml";
		$mess = "OK";
	}
	elsif(-e "$path\\$gpm.xml.gz")	{
		$mess = `pigz -d -f -q $path\\$gpm.xml.gz`;
		if(-e "$path\\$gpm.xml")	{
			$path = "$path\\$gpm.xml";
			$mess = "OK";
		}
		else	{
			$mess = "Error: file not found";
			$path = "";
		}
	}
	else	{
		$mess = "Error: file not found";
		$path = "";
	}
	$class->{_path} = $path;
	$class->{_mess} = $mess;
}

sub gpmModelProteins
{
	my ($class) = @_; 
	my $gpm = $class->{_gpm} . ".xml";                 

	my $dbh = $class->{_dbh_g};
	my $sql = "SELECT DISTINCT file, p.proid, pathid, uid, p.expect, label, ps.seq FROM Protein p, Result r, Peptide pep, ProSeq ps WHERE r.resultid = p.resultid AND p.proid = pep.proid AND ps.proseqid = p.proseqid AND didb = 1 AND file = ? AND active = 1 ORDER BY expect";
	my $sth = $dbh->prepare($sql) or warn "Prepare failed: $DBI::errstr($DBI::err)\n";
	my $rows = $sth->execute($gpm) or warn "Execute failed: $DBI::errstr($DBI::err)\n";
	my @r;
	my $out = "[";
	while(@r = $sth->fetchrow_array())	{
		$out .= "\t\"@r[5]\",\n";
	}
	$out =~ s/\,$//;
	$out .= "]";
	$sth->finish();
	return $out;
}

sub gpmModelProteinPeptides
{
	my ($class,$acc) = @_; 
	my $gpm = $class->{_gpm} . ".xml";                 
	my $dbh = $class->{_dbh_g};
	my $sql = "select p.proid from Protein p,Result r,ProSeq ps where r.resultid=p.resultid and ps.proseqid=p.proseqid and r.file= ? and ps.label= ?";
	my $sth = $dbh->prepare($sql) or warn "Prepare failed: $DBI::errstr($DBI::err)\n";
	$sth->execute($gpm,$acc) or warn "Execute failed: $DBI::errstr($DBI::err)\n";
	my @r;
	@r = $sth->fetchrow_array();
	if(not @r)	{
		$sth->finish();
		return "{ }";
	}
	$sql = "select pep.seq from Peptide pep where proid = ?";
	$sth = $dbh->prepare($sql) or warn "Prepare failed: $DBI::errstr($DBI::err)\n";
	$sth->execute(@r[0]) or warn "Execute failed: $DBI::errstr($DBI::err)\n";
	my %sequences;
	while(@r = $sth->fetchrow_array())	{
		$sequences{@r[0]} = $sequences{@r[0]} + 1;
	}
	my @keys = keys(%sequences);
	my $out = "{";
	foreach $_(@keys)	{
		$out .= "\t\"$_\" : $sequences{$_},\n";
	}
	$out =~ s/\,$//;
	$out .= "}";
	$sth->finish();
	return $out;
}

sub gpmModelProteinModifications
{
	my ($class,$acc) = @_; 
	my $gpm = $class->{_gpm} . ".xml";                 
	my $dbh = $class->{_dbh_g};
	my $sql = "select p.proid from Protein p,Result r,ProSeq ps where r.resultid=p.resultid and ps.proseqid=p.proseqid and r.file= ? and ps.label= ?";
	my $sth = $dbh->prepare($sql) or warn "Prepare failed: $DBI::errstr($DBI::err)\n";
	$sth->execute($gpm,$acc) or warn "Execute failed: $DBI::errstr($DBI::err)\n";
	my @r;
	@r = $sth->fetchrow_array();
	if(not @r)	{
		$sth->finish();
		return "[ ]";
	}
	$sql = "select a.type,a.at,a.modified,a.pm from Aa a, Peptide pep where a.pepid=pep.pepid and pep.proid = ?";
	$sth = $dbh->prepare($sql) or warn "Prepare failed: $DBI::errstr($DBI::err)\n";
	$sth->execute(@r[0]) or warn "Execute failed: $DBI::errstr($DBI::err)\n";
	my @observed;
	while(@r = $sth->fetchrow_array())	{
		if(not @r[3])	{
			push(@observed,"$acc:pm.@r[0]@r[1]#@r[2];");
		}
	}
	my $out = "[";
	foreach $_(@observed)	{
		$out .= "\t\"$_\",\n";
	}
	$out =~ s/\,$//;
	$out .= "]";
	$sth->finish();
	return $out;
}

sub gpmModelProteinAps
{
	my ($class,$acc) = @_; 
	my $gpm = $class->{_gpm} . ".xml";                 
	my $dbh = $class->{_dbh_g};
	my $sql = "select p.proid from Protein p,Result r,ProSeq ps where r.resultid=p.resultid and ps.proseqid=p.proseqid and r.file= ? and ps.label= ?";
	my $sth = $dbh->prepare($sql) or warn "Prepare failed: $DBI::errstr($DBI::err)\n";
	$sth->execute($gpm,$acc) or warn "Execute failed: $DBI::errstr($DBI::err)\n";
	my @r;
	@r = $sth->fetchrow_array();
	if(not @r)	{
		$sth->finish();
		return "[ ]";
	}
	$sql = "select a.type,a.at,a.modified,a.pm from Aa a, Peptide pep where a.pepid=pep.pepid and pep.proid = ?";
	$sth = $dbh->prepare($sql) or warn "Prepare failed: $DBI::errstr($DBI::err)\n";
	$sth->execute(@r[0]) or warn "Execute failed: $DBI::errstr($DBI::err)\n";
	my @observed;
	while(@r = $sth->fetchrow_array())	{
		if(@r[3])	{
			push(@observed,"$acc:p.@r[0]@r[1]@r[3];");
		}
	}
	my $out = "[";
	foreach $_(@observed)	{
		$out .= "\t\"$_\",\n";
	}
	$out =~ s/\,$//;
	$out .= "]";
	$sth->finish();
	return $out;
}

sub gpmModelMetadata
{
	my ($class) = @_;
	getPath($class);
	if(not $class->{_path})	{
		return $class->{_mess};
	}
	open(IN,"<$class->{_path}") or die "could not open $class->{_gpm}";
	my %input;
	my %unused;
	my %perform;
	my $r;
	my $v;
	while(<IN>)	{
		if(/\<group\s+label=\"input\sparameters\"\s+type=\"parameters\"/)	{
			$_ = <IN>;
			while($_ and not /\<\/group>/)	{
				($r) = $_ =~ /label=\"(.+?)\"/;
				if($r)	{
					($v) = $_ =~ /\>(.+?)\</;
					$r =~ s/\W+/\_/g;
					$input{$r} = $v;
				}
				$_ = <IN>;
			}
			if(not $_)	{
				last;
			}
		}
		elsif(/\<group\s+label=\"unused\sinput\sparameters\"\s+type=\"parameters\"/)	{
			$_ = <IN>;
			while($_ and not /\<\/group>/)	{
				($r) = $_ =~ /label=\"(.+?)\"/;
				if($r and $r =~ /^gpmdb\,/)	{
					($v) = $_ =~ /\>(.+?)\</;
					$r =~ s/\W+/\_/g;
					$unused{$r} = $v;
				}
				$_ = <IN>;
			}
			if(not $_)	{
				last;
			}
		}
		elsif(/\<group\s+label=\"performance\sparameters\"\s+type=\"parameters\"/)	{
			$_ = <IN>;
			while($_ and not /\<\/group>/)	{
				($r) = $_ =~ /label=\"(.+?)\"/;
				if($r)	{
					($v) = $_ =~ /\>(.+?)\</;
					$r =~ s/\W+/\_/g;
					$perform{$r} = $v;
				}
				$_ = <IN>;
			}
			if(not $_)	{
				last;
			}
		}
	}
	close(IN);
	my @keys = sort {$a cmp $b} keys(%input); 
	my $output = "{\n\"parameter\":{\n";
	foreach $r(@keys)	{
		if($r eq 'output_message'){
			next;
		}
		$output .= qq(\t"$r" : "$input{$r}",\n);
	}
	$output =~ s/\,$//;
	$output .= "},\n\"sample\":{\n";
	@keys = sort {$a cmp $b} keys(%unused); 
	foreach $r(@keys)	{
		$output .= qq(\t"$r" : "$unused{$r}",\n);
	}
	$output =~ s/\,$//;
	$output .= "},\n\"result\":{\n";
	@keys = sort {$a cmp $b} keys(%perform); 
	foreach $r(@keys)	{
		$output .= qq(\t"$r" : "$perform{$r}",\n);
	}
	$output =~ s/\,$//;
	$output .= "}\n}";
	return $output;	
}

sub gpmModelProteinSequence
{
	my ($class,$acc) = @_;
	getPath($class);
	if(not $class->{_path})	{
		return $class->{_mess};
	}
	open(IN,"<$class->{_path}") or die "could not open $class->{_gpm}";
	my %input;
	my %unused;
	my %perform;
	my $r;
	my $v;
	while(<IN>)	{
		if(/\<protein .+ label\=\"$acc\"/)	{
			$_ = <IN>;
			$v = "";
			while($_ and not /\<peptide/)	{
				$_ = <IN>;
			}
			$_ = <IN>;
			while($_ and not /\</)	{
				$v .= $_;
				$_ = <IN>;
			}
			last;
		}
	}
	close(IN);
	$v =~ s/[^A-Z]//g;
	return "[\"$v\"]";	
}

1;

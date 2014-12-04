#
#===============================================================================
#
#         FILE:  Affy.pm
#
#  DESCRIPTION:  Package to handle Affymetrix files
#         BUGS:  ---
#        NOTES:  The package is spun off from the AnnotationIO.pm
#       AUTHOR:  Chaolin Zhang (cz), czhang@rockefeller.edu
#      COMPANY:  Rockefeller University
#      VERSION:  1.0
#      CREATED:  12/17/2010
#     REVISION:  ---
#===============================================================================

package Affy;

our $VERSION = 1.01;


=head1 NAME

Affy - read and write Affymetrix files

subroutines starting with a hyphen should not be called outside

=cut

use strict;
use warnings;
use Data::Dumper;
use Carp;


=head2 read1LQFile

read .1lq file that contains probe coordinates and sequences

my $ret = read1LQFile ($inFile)

$ret is an array reference. each element contains:

	X=>
	Y=>
	SEQUENCE=>
	...

=cut
sub read1LQFile
{
	my $in = $_[0];
	my $fin;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	while (my $line=<$fin>)
	{
		chomp $line;
		last if $line=~/^X\tY/;
	}

	my @ret;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		my @cols = split ("\t", $line);
		push @ret, {
			X=>$cols[0],
			Y=>$cols[1],
			SEQUENCE=>$cols[2],
			DESTYPE=>$cols[3],
			FEATURE=>$cols[4],
			QUALIFIER=>$cols[5],
			EXPOS=>$cols[6],
			PLEN=>$cols[7],
			POS=>$cols[8],
			CBASE=>$cols[9],
			PBASE=>$cols[10],
			TBASE=>$cols[11],
			IPBASE=>$cols[12],
			UNIT=>$cols[13],
			BLOCK=>$cols[14],
			ATOM=>$cols[15]
		};
	}
	close ($fin);
	return \@ret;
}

=head2 readPgfFile

read .pgf file that contains probe coordinates and sequences

my $ret = readPgfFile ($inFile)

$ret is an array reference. each element contains:

	X=>
	Y=>
	SEQUENCE=>
	...

=cut
sub readPgfFile
{
	my $in = $_[0];
	my $fin;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	
	while (my $line=<$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		last if $line =~/^\#\%header2/;
	}

	my @ret;
	my ($probesetId, $probesetType, $probesetName) = ("", "", "");
	my ($atom_id, $exon_position);

	while (my $line=<$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;

		#print $line, "\n";
		if ($line=~/^\S/)
		{
			($probesetId, $probesetType, $probesetName) = split ("\t", $line);
			next;
		}

		#next unless $line =~/pm\:st/;
	
		$line =~/^\s*(\S.*)$/;
		$line = $1;	

		#$line = reverse ($line);

		my @cols = split ("\t", $line);

		if (@cols <=2)
		{
			($atom_id, $exon_position) = @cols;
		}
		else
		{
			push @ret, {
				probeset_id=>$probesetId,
				probeset_type=>$probesetType,
				probeset_name=>$probesetName,
				atom_id=>$atom_id,
				exon_position=>$exon_position,
				probe_id=>$cols[0],
				type=>$cols[1],
				gc_count=>$cols[2],
				probe_length=>$cols[3],
				interrogation_position=>$cols[4],
				probe_sequence=>$cols[5],
			};
			#Carp::croak Dumper (\@ret), "\n";
		}

		#Carp::croak Dumper (pop @ret), "\n";
	}
	close ($fin);
	return \@ret;
}




=head2 readCelFile

Usage:
my $ret = readCelFile ($inFile);

	$ret = {
		header=>
		intensity=>
		masks=>
		outliers=>
		modified=>
	}

=cut

sub readCelFile
{
	my $in = $_[0];
	my %ret;
   	
	my $fin;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	my $line = <$fin>;
	$line = <$fin>;
	chomp $line;
	my ($name, $value) = split ("=", $line);
	Carp::croak "only Cel file version 3 is recognized\n" if $value != 3;
	my $iter = 0;
	while ($line =<$fin>)
	{
		$iter++;
		#last if $iter > 10;
		chomp $line;
		next if $line=~/^\s*$/;
		#print $line, "\n";
		if ($line eq '[HEADER]')
		{
			#print "parsing header...\n";
			$ret{"header"} = parseCelHeader ($fin);
			Carp::croak "missing number of columns\n" unless exists $ret{"header"}->{"Cols"};
			Carp::croak "missing number of rows\n" unless exists $ret{"header"}->{"Rows"};
		}
		elsif ($line eq '[INTENSITY]')
		{
			#print "parsing intensity...\n";
			$ret{"intensity"} = parseCelIntensity ($fin);
			Carp::croak "incorrect number of columns\n" if @{$ret{"intensity"}} != $ret{"header"}->{"Cols"};
			Carp::croak "incorrect number of rows\n" if @{$ret{"intensity"}->[0]} != $ret{"header"}->{"Rows"};
		}
		elsif ($line eq '[MASKS]')
		{
			#print "parsing masks...\n";
			$ret{"masks"} = parseCelMasksOrOutliers ($fin);
		}
		elsif ($line eq '[OUTLIERS]')
		{
			$ret{"outliers"} = parseCelMasksOrOutliers ($fin);
		}
		elsif ($line eq '[MODIFIED]')
		{
			$ret{"modified"} = parseCelModified ($fin);
		}
	}
	close ($fin);
	return \%ret;	
}


sub parseCelHeader
{
	my $fin = $_[0];
	my %ret;
	while (my $line =<$fin>)
	{
		chomp $line;
		last if $line=~/^\s*$/;
		my @cols = split ("=", $line);
		Carp::croak "incorrect name=value pair: $line\n" if @cols < 2;
		my $name = shift @cols;
		my $value = join ("=", @cols);
		$ret{$name} = $value;
	}
	return \%ret;
}

sub parseCelIntensity
{
	my $fin = $_[0];
	my $ret;
	my $line = <$fin>;
	my ($name, $numCells) = split ("=", $line);
	
	Carp::croak "NumberCells not found:$line\n" if $name ne 'NumberCells';
	$line = <$fin>;
	for (my $i = 0; $i < $numCells; $i++)
	{
		my $line = <$fin>;
		chomp $line;
		my @cols = split (/\s+/, $line);
		shift @cols if @cols > 5;
		Carp::croak "incorrect number of columns at line: $line\n" if @cols != 5;
		my ($x, $y, $mean, $stdev, $npixels) = @cols;
		$ret->[$x][$y] = {MEAN=>$mean, STDV=>$stdev};
	}
	return $ret;
}

sub parseCelMasksOrOutliers
{
	my $fin = $_[0];
	my @ret;
	my $line = <$fin>;
	my ($name, $numCells) = split ("=", $line);
	
	Carp::croak "NumberCells not found:$line\n" if $name ne 'NumberCells';
	$line = <$fin>;
	for (my $i = 0; $i < $numCells; $i++)
	{
		my $line = <$fin>;
		chomp $line;
		my @cols = split (/\s+/, $line);
		shift @cols if @cols > 2;
		Carp::croak "incorrect number of columns at line: $line\n" if @cols != 2;
		my ($x, $y) = @cols;
		push @ret, {X=>$x, Y=>$y};
	}
	return \@ret;
}

sub parseCelModified
{
	my $fin = $_[0];
	my @ret;
	my $line = <$fin>;
	my ($name, $numCells) = split ("=", $line);
	
	Carp::croak "NumberCells not found:$line\n" if $name ne 'NumberCells';
	$line = <$fin>;
	for (my $i = 0; $i < $numCells; $i++)
	{
		my $line = <$fin>;
		chomp $line;
		my @cols = split (/\s+/, $line);
		shift @cols if @cols > 3;
		Carp::croak "incorrect number of columns at line: $line\n" if @cols != 3;
		my ($x, $y, $origMean) = @cols;
		push @ret, {X=>$x, Y=>$y, ORIGMEAN=>$origMean};
	}
	return \@ret;
}


=head2 readSgrFile

my $ret = readSgrFile ($inFile)

$ret is an array reference, with each row containing: 

	chrom=>
	pos=>
	score=>

=cut


sub readSgrFile
{
	my $in = $_[0];
	my $fin;
	my @ret;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	while (my $line =<$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		my ($chr, $pos, $score) = split ("\t", $line);
		push @ret, {chrom=>$chr, pos=>$pos, score=>$score};
	}
	close ($fin);
	return \@ret;
}


1;



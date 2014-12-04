#
#===============================================================================
#
#         FILE:  ScaleFree.pm
#
#  DESCRIPTION:  Package to handle scale free distribution
#         BUGS:  ---
#        NOTES:  
#       AUTHOR:  Chaolin Zhang (cz), czhang@rockefeller.edu
#      COMPANY:  Rockefeller University
#      VERSION:  1.0
#      CREATED:  11/20/2010
#     REVISION:  ---
#===============================================================================


package ScaleFree;

require Exporter;

@ISA = qw (Exporter);

our $VERSION = 1.01;


@EXPORT = qw (
	buildHarmonicNumberTable
	encode
	getHarmonicNumber
	getParams
	initParams
	loadHarmonicNumberTable
	saveHarmonicNumberTable
);



=head1 NAME

ScaleFree - handle scale free distribution


=cut


use strict;
use warnings;


use Data::Dumper;
use Carp;

use Common;

my $debug = 0;


#default parameters

my $step = 0.001;
my $start = 0;
my $end = 8;
my $eps = 1e-6;
my $verbose = 0;


=head2 initParams


=cut
sub initParams
{
	my %param = @_;
	foreach my $name (keys %param)
	{
		#print "$name=", $param{$name}, "\n";

		if ($name eq 'step')
		{
			$step = $param{$name};
		}
		elsif ($name eq 'start')
		{
			$start = $param{$name};
		}
		elsif ($name eq 'end')
		{
			$end = $param{$name};
		}
		elsif ($name eq 'eps')
		{
			$eps = $param{$name};
		}
		elsif ($name eq 'v')
		{
			$verbose = $param{$name};
		}
		else
		{
			Carp::croak "unknown parameter: $name\n";
		}
	}
}

sub getParams
{
	my $name=$_[0];
	if ($name eq 'step')
	{
		return $step;
	}
	elsif ($name eq 'start')
	{
		return $start;
	}
	elsif ($name eq 'end')
	{
		return $end;
	}
	elsif ($name eq 'eps')
	{
		return $eps;
	}
	elsif ($name eq 'v')
	{
		return $verbose;
	}
	else
	{
		Carp::croak "unknown parameter: $name\n";
	}
}

sub encode
{
	my ($s) = @_;
	my $decimal_point = int(0.5-log($step) / log(10));
	return sprintf ("%.$decimal_point"."f", $s);
}


sub buildHarmonicNumberTable
{
	my $maxN = $_[0];
	my %table;
	for (my $s = $start; $s <= $end; $s+= $step)
	{
		my $s2 = encode ($s);
		print "s = $s2\n" if $verbose;
		for (my $N = 1; $N <= $maxN; $N++)
		{
			my $HN = 0;
			map {$HN += 1/Common::pow ($_, $s)} (1..$N);
			$table{$s2}{$N} = $HN;
		}
	}
	return \%table;
}


sub saveHarmonicNumberTable
{
	my ($table, $outFile) = @_;

	my $fout;
	open ($fout, ">$outFile") || Carp::croak "cannot open $outFile to write\n";
	
	foreach my $s (sort {$a <=> $b} keys %$table)
	{
		my $s2 = encode ($s);

		print "s=$s2\n" if $verbose;
		foreach my $N (sort {$a <=> $b} keys %{$table->{$s2}})
		{
			print $fout join ("\t", $s2, $N, $table->{$s2}{$N}), "\n";
		}
	}
	close ($fout);
}


sub loadHarmonicNumberTable
{
	my ($inFile) = @_;

	my %table;
	my $fin;
	my $i = 0;
	open ($fin, "<$inFile") || Carp::croak "cannot open $inFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;

		print "$i ...\n" if $verbose && $i % 100000 == 0;
		$i++;

		my ($s, $N, $z) = split (/\t/, $line);
		my $s2 = encode ($s);
		$table{$s2}{$N} = $z;
	}
	close ($fin);
	return \%table;
}

sub getHarmonicNumber
{
	my ($table, $s, $N) = @_;
	my $s2 = encode ($s);
	return exists $table->{$s2} && exists $table->{$s2}{$N} ? $table->{$s2}{$N} : -1;
}

1;


=end
sub buildZetaTable
{
	my %zeta;
	for (my $s = $start; $s <= $end; $s+=$step)
	{
		#print "s=$s\n";
		my $i;
		my $s2 = encode ($s, $step);
		for ($i = 1; $i < 1e7; $i++)
		{
			my $p = Common::pow ($i, -$s);
			if (exists $zeta{$s2})
			{
				last if $p / $zeta{$s2} < $eps;
			}
			$zeta{$s2} += $p;
		}
	
		print "zeta($s2)=", $zeta{$s2}, "\n" if $verbose;
		#print ", mean=", $zeta{encode($s-1, $step)} / $zeta{encode ($s, $step)} if exists $zeta{encode ($s-1, $step)};
		#print ", inaccurate" if $i == 1e5;
		#print "\n";	
	}
	return \%zeta;
}

sub saveZetaTable
{
	my ($zeta, $outFile) = @_;

	my $fout;
	open ($fout, ">$outFile") || Carp::croak "cannot open $outFile to write\n";
	
	foreach my $s (sort {$a <=> $b} keys %$zeta)
	{
		my $s2 = encode ($s, $step);
		print $fout join ("\t", $s2, $zeta->{$s2}), "\n";
	}
	close ($fout);
}

sub loadZetaTable
{
	my ($inFile) = @_;

	my %zeta;
	my $fin;
	open ($fin, "<$inFile") || Carp::croak "cannot open $inFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;

		my ($s, $z) = split (/\t/, $line);
		my $s2 = encode ($s, $step);
		$zeta{$s2} = $z;
	}
	close ($fin);
	return \%zeta;
}

sub getZetaNumber
{
	my ($zetaTable, $s) = @_;
	my $s2 = encode ($s, $step);
	return exists $zetaTable->{$s2} ? $zetaTable->{$s2} : -1;
}

1;

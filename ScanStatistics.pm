#
#===============================================================================
#
#         FILE:  ScanStatistics.pm
#
#  DESCRIPTION:  Package to calculate scan statistics
#         BUGS:  ---
#        NOTES:  
#       AUTHOR:  Chaolin Zhang (cz), czhang@rockefeller.edu
#      COMPANY:  Rockefeller University
#      VERSION:  1.0
#      CREATED:  10/21/2010
#     REVISION:  ---
#===============================================================================


package ScanStatistics;

require Exporter;

@ISA = qw (Exporter);

our $VERSION = 1.01;


@EXPORT = qw (
	calcScanStatistic
);


=head1 NAME

ScaleFree - calculate scan statistics

Reference

J Glaz, J Naus, and S Wallenstein 2001 Scan Statistics, Springer, P.28

usage:

my $pvalue = calcScanStatistic ($k, $phi, $N);

$k: observed peak height
$phi: the expected peak height
$N: number of nonoverlapping window

=cut


use strict;
use warnings;


#dictionary stores some intermediate results to save time

my %dictionary;


#  lamda_g=a_g * T /sum(a_g * (L_g-w)) is the expected number of reads at a single position, a_g and L_g are the abundance of gene g, w is the read length
#  phi_g = lamda_g * w
#  N_g = L_g / w is the number of nonoverlapping window
#

sub calcScanStatistic 
{
	my ($k, $phi, $N) = @_;

	my $p = 1 - Q2($k, $phi) * (Q3 ($k, $phi) / Q2($k, $phi))** ($N - 2);
	$p = 0 if $p < 0; #rounding error when p is really small
	%dictionary = ();

	return $p;
}

sub Q2
{
	my ($k, $phi) = @_;

	$dictionary{'Fp'}{$k-1}{$phi} = Fp ($k -1, $phi) unless exists $dictionary{'Fp'}{$k-1}{$phi};
	$dictionary{'Fp'}{$k-3}{$phi} = Fp ($k -3, $phi) unless exists $dictionary{'Fp'}{$k-3}{$phi};
	$dictionary{'p'}{$k}{$phi} = p($k, $phi) unless exists $dictionary{'p'}{$k}{$phi};
	$dictionary{'p'}{$k-2}{$phi} = p($k-2, $phi) unless exists $dictionary{'p'}{$k-2}{$phi};

	#my $ret = Fp ($k -1, $phi)**2 - ($k-1) * p($k, $phi) * p ($k-2, $phi) - ($k-1 - $phi) * p ($k, $phi) * Fp ($k-3, $phi);
	my $ret = $dictionary{'Fp'}{$k-1}{$phi} ** 2 
		- ($k-1) * $dictionary{'p'}{$k}{$phi} * $dictionary{'p'}{$k-2}{$phi} 
		- ($k-1 - $phi) * $dictionary{'p'}{$k}{$phi} * $dictionary{'Fp'}{$k-3}{$phi};
	
	#print "Q2($k, $phi) = $ret\n";
	return $ret;
}

sub p
{
	my ($j, $phi) = @_;

	my $s = 0;

	for (my $i = 1; $i <= $j; $i++)
	{
		$s+= log($i);
	}

	my $ret = -$phi + $j * log($phi) - $s;

	#print "p($j, $phi)=", exp($ret), "\n";
	return exp($ret);
}

sub Fp
{
	my ($k, $phi) = @_;

	return 0 if $k < 0;
	my $ret = 0;

	for (my $i = 0; $i <= $k; $i++)
	{
		$dictionary{'p'}{$i}{$phi} = p($i, $phi) unless exists $dictionary{'p'}{$i}{$phi};
		#$ret += p ($i, $phi);
		$ret += $dictionary{'p'}{$i}{$phi};
	}
	#print "Fp ($k, $phi)=$ret\n";
	return $ret;
}


sub Q3
{
	my ($k, $phi) = @_;
	
	$dictionary{'Fp'}{$k-1}{$phi} = Fp ($k -1, $phi) unless exists $dictionary{'Fp'}{$k-1}{$phi};
	#my $ret = Fp ($k-1, $phi)**3 - A1 ($k, $phi) + A2 ($k, $phi) + A3 ($k, $phi) - A4 ($k, $phi);
	my $ret = $dictionary{'Fp'}{$k-1}{$phi} ** 3 - A1 ($k, $phi) + A2 ($k, $phi) + A3 ($k, $phi) - A4 ($k, $phi);

	#print "Q3($k, $phi)=$ret\n";
	return $ret;
}


sub A1
{
	my ($k, $phi) = @_;
	
	$dictionary{'Fp'}{$k-1}{$phi} = Fp ($k -1, $phi) unless exists $dictionary{'Fp'}{$k-1}{$phi};
	$dictionary{'Fp'}{$k-2}{$phi} = Fp ($k -2, $phi) unless exists $dictionary{'Fp'}{$k-2}{$phi};
	$dictionary{'Fp'}{$k-3}{$phi} = Fp ($k -3, $phi) unless exists $dictionary{'Fp'}{$k-3}{$phi};
	$dictionary{'p'}{$k}{$phi} = p($k, $phi) unless exists $dictionary{'p'}{$k}{$phi};

	
	#my $ret = 2 * p($k, $phi) * Fp ($k-1, $phi) * (($k-1) * Fp ($k-2, $phi) - $phi * Fp ($k-3, $phi));
	my $ret = 2 * $dictionary{'p'}{$k}{$phi} * $dictionary{'Fp'}{$k-1}{$phi} * 
			(($k-1) * $dictionary{'Fp'}{$k-2}{$phi} - $phi * $dictionary{'Fp'}{$k-3}{$phi});
	#print "A1=$ret\n";
	return $ret;
}

sub A2
{
	my ($k, $phi) = @_;
	
	$dictionary{'Fp'}{$k-3}{$phi} = Fp ($k -3, $phi) unless exists $dictionary{'Fp'}{$k-3}{$phi};
	$dictionary{'Fp'}{$k-4}{$phi} = Fp ($k -4, $phi) unless exists $dictionary{'Fp'}{$k-4}{$phi};
	$dictionary{'Fp'}{$k-5}{$phi} = Fp ($k -5, $phi) unless exists $dictionary{'Fp'}{$k-5}{$phi};
	$dictionary{'p'}{$k}{$phi} = p($k, $phi) unless exists $dictionary{'p'}{$k}{$phi};

	#my $ret = 0.5 * p($k, $phi)**2 * (($k-1)* ($k-2) * Fp ($k-3, $phi) - 2 * ($k-2) * $phi * Fp ($k-4, $phi)+ $phi**2 * Fp ($k-5, $phi));
	my $ret = 0.5 * $dictionary{'p'}{$k}{$phi}**2 * 
		(($k-1)* ($k-2) * $dictionary{'Fp'}{$k-3}{$phi} - 2 * ($k-2) * $phi * $dictionary{'Fp'}{$k-4}{$phi} + $phi**2 * $dictionary{'Fp'}{$k-5}{$phi});
	#print "A2=$ret\n";
	return $ret;
}

sub A3
{
	my ($k, $phi) = @_;
	my $ret = 0;

	for (my $i=1; $i <= $k-1; $i++)
	{
		$dictionary{'Fp'}{$i-1}{$phi} = Fp ($i-1, $phi) unless exists $dictionary{'Fp'}{$i-1}{$phi};
		$dictionary{'p'}{2*$k-$i}{$phi} = p(2*$k-$i, $phi) unless exists $dictionary{'p'}{2*$k-$i}{$phi};

		#$ret += p(2*$k - $i, $phi) * Fp($i-1, $phi)**2;
		$ret += $dictionary{'p'}{2*$k-$i}{$phi} * $dictionary{'Fp'}{$i-1}{$phi} **2;
	}
	#print "A3=$ret\n";
	return $ret;
}

sub A4
{
	my ($k, $phi) = @_;
	my $ret = 0;
	for (my $i = 2; $i <= $k-1; $i++)
	{
		$dictionary{'Fp'}{$i-2}{$phi} = Fp ($i-2, $phi) unless exists $dictionary{'Fp'}{$i-2}{$phi};
		$dictionary{'Fp'}{$i-3}{$phi} = Fp ($i-3, $phi) unless exists $dictionary{'Fp'}{$i-3}{$phi};

		$dictionary{'p'}{2*$k-$i}{$phi} = p(2*$k-$i, $phi) unless exists $dictionary{'p'}{2*$k-$i}{$phi};
		$dictionary{'p'}{$i}{$phi} = p($i, $phi) unless exists $dictionary{'p'}{$i}{$phi};

		#$ret += p(2*$k-$i, $phi) * p($i, $phi) * (($i-1)*Fp ($i-2, $phi) - $phi * Fp ($i-3, $phi));
		$ret += $dictionary{'p'}{2*$k-$i}{$phi} * $dictionary{'p'}{$i}{$phi} * (($i-1)* $dictionary{'Fp'}{$i-2}{$phi} - $phi * $dictionary{'Fp'}{$i-3}{$phi});
	}
	#print "A4=$ret\n";
	return $ret;
}


1;

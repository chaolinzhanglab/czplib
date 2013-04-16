#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  test.pl
#
#        USAGE:  ./test.pl  
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Chaolin Zhang (cz), czhang@rockefeller.edu
#      COMPANY:  Rockefeller University
#      VERSION:  1.0
#      CREATED:  04/05/11 10:08:23
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

use Bed;
use Common;

use Carp;
die "script <breakdown.bed> <contig.bed> <out.bed>\n" if (@ARGV != 3);


my ($breakdownBedFile, $contigBedFile, $outBedFile) = @ARGV;

print "read breakdown ...\n";
my $breakdown = readBedFile ($breakdownBedFile, 1);

print "read contig ...\n";
my $contigs = readBedFile ($contigBedFile, 1);


print "overlap ...\n";

getGenomicBreakdown ($breakdown, $contigs, 1);


print "dump output\n";
my $fout;
open ($fout, ">$outBedFile");
for my $c (@$contigs)
{
	my $b = $c->{'breakdown'};
	foreach my $r (@$b)
	{
		print $fout bedToLine ($r), "\n";
	}
}
close ($fout);

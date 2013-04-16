#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  test_Maf.pl
#
#        USAGE:  ./test_Maf.pl  
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
#      CREATED:  03/30/11 10:35:20
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

use Maf;


die "script <in.maf> <out.maf>\n" if @ARGV != 2;

my ($inFile, $outFile) = @ARGV;

my $fin;
my $fout;
open ($fin, "<$inFile");
open ($fout, ">$outFile");

my $iter = 0;
while (my $block = readNextMafBlock ($fin))
{
	print "$iter ...\n"; $iter++;
	writeMafBlock ($fout, $block);
}


close ($fin);
close ($fout);


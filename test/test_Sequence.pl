#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  test2.pl
#
#        USAGE:  ./test2.pl  
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
#      CREATED:  12/29/10 15:41:52
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

use Sequence;

my $seq = Sequence::enumerateSeq (4);

print join ("\n", @$seq), "\n";



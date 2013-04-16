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
#      CREATED:  12/27/10 20:25:11
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

use Bio::SeqIO;
use BWT;
use Sequence;

my $seqIO = Bio::SeqIO-> new (-file=>'chr1.fa', format=>'fasta');

print "reading chr1\n";
my $seq = $seqIO->next_seq ();

print $seq->length(), " nt loaded\n";

my $seqStr = uc ($seq->seq());
$seqStr .= " ";

print "performing bwt ...\n";

my $str = "abracadabra";
#my $seqStr = substr ($seq->seq(), 70000*50, 10000);

my $W = BWT::bwt2 ($seqStr, undef, 1);
#or
#my $W = BWT2::bwt_naive ($str, 1);

my $fout;
open ($fout, ">out.txt");
print $fout segmentStr ($W, 60), "\n";
close ($fout);


=end
#my $str = "abracadabra";
##my $str = "acaaccg";
#my $str = 'g';



my %ret = BWT2::bwt2 ($str, undef, 1);


print "bwt=", $ret{'W'}, "\n";

my $SA = $ret{'SA'};
print "SA= ", join ("\t", @$SA), "\n";

my $SAinv = BWT2::invert ($SA);

print "SAinv=", join ("\t", @$SAinv), "\n";

my $psi = BWT2::SA2Psi ($SA);

print "psi =", join ("\t", @$psi), "\n";


$SAinv = BWT2::psi2SAinv ($psi);

print "psi2SAinv = ", join ("\t", @$SAinv), "\n";


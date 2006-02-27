#!/usr/bin/perl -w
use strict;
use Common;

my @seq = (0..9);
my $bootstrap = Common::bootstrapArray (\@seq);
my $diff = Common::diffArray (\@seq, $bootstrap);

print "bootstrap:", join ("\t", @$bootstrap), "\n";
print "left out:", join ("\t", @$diff), "\n";

#!/usr/bin/perl -w

use lib "/data/zhangc/perl_lib2/";

use strict;
use AnnotationIO;

die "script <in> <out>\n" if @ARGV != 2;

my ($in, $out) = @ARGV;

print "reading motifs...\n";
my $motifs = AnnotationIO::readMotifFile ($in);

print "writing motifs ...\n";
AnnotationIO::writeMotifFile ($motifs, $out);



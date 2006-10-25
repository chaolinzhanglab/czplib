#!/usr/bin/perl -w
#
use lib "/mnt/raygun/zhangc/perl_lib2";
use Data::Dumper;
use Common;

die "script <db.nib> <query.fa> <out.psl>\n" if @ARGV != 3;

my ($db, $query, $out) = @ARGV;

my $results = Common::blat ("blat", $db, $query);
AnnotationIO::writePslFile ($results, $out);


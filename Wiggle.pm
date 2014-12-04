#
#===============================================================================
#
#         FILE:  Wiggle.pm
#
#  DESCRIPTION:  Package to handle Wiggle file
#         BUGS:  ---
#        NOTES:  The package is spun off from the AnnotationIO.pm
#       AUTHOR:  Chaolin Zhang (cz), czhang@rockefeller.edu
#      COMPANY:  Rockefeller University
#      VERSION:  1.0
#      CREATED:  12/17/10
#     REVISION:  ---
#===============================================================================


package Wiggle;


require Exporter;

our $VERSION = 1.01;

@ISA = qw (Exporter);

@EXPORT = qw (
	readBedGraphFile
	writeBedGraphFile
);


=head1 NAME

Wiggle - read and write wiggle like files
subroutines starting with a hyphen should not be called outside

=cut

use strict;
use warnings;
use Data::Dumper;
use Carp ();




=head2 readBedGraphFile

Note: original name: readWigFile

my $ret = readBedGraphFile ($inFile)

=cut

sub readWigFile
{
	Carp::croak "obsolete subroutine, call readBedGraphFile\n";
}

sub readBedGraphFile
{
	my $in = $_[0];
	my $fin;
	my @ret;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		next if $line =~/^\#/;
		my ($chrom, $chromStart, $chromEnd, $score) = split (/\s/, $line);
		push @ret, {chrom=>$chrom, chromStart=>$chromStart, chromEnd=>$chromEnd-1, score=>$score};
	}
	close ($fin);
	return \@ret;
}

=head2 writeBedGraphFile

Note: original name: writeWigFile

writeBedGraphFile (\@regions,  $outFile, $headerLine);

=cut

sub writeWigFile
{
	Carp::croak "obsolete subroutine, call writeBedGraphFile\n";
}
sub writeBedGraphFile
{
	my ($regions, $out, $header) = @_;
	my $fout;
	open ($fout, ">$out") || Carp::croak "cannot open file $out to write\n";
	if ($header ne '')
	{
		print $fout $header, "\n";
	}

	foreach my $r (@$regions)
	{
		print $fout join ("\t", $r->{"chrom"},
				$r->{"chromStart"},
				$r->{"chromEnd"} + 1,
				$r->{"score"}), "\n";
	}
	close ($fout);
}



1;



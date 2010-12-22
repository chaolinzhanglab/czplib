package ViennaRNA;

=head1 NAME

ViennaRNA - interface of ViennaRNA package
subroutines starting with a hyphen should not be called outside

=cut

use strict;
use warnings;
use Data::Dumper;
use Carp;

=head2 readRNAplfoldFile

my $ret = readRNAplfoldFile

=cut


sub readRNAplfoldFile
{
	my $inFile = $_[0];

	my $fin = new FileHandle;
	open ($fin, "<$inFile") || Carp::croak "can not open file $inFile to read\n";

	my @ret;
	my $seqLen = -1;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;

		if ($line =~/^\/sequence \{/)
		{
			my $seq = <$fin>;
			chomp $seq;
			$seqLen = length ($seq) - 1;

			#print "seqlen = $seqLen\n";
			for (my $i = 0; $i < $seqLen; $i++)
			{
				$ret[$i] = 0;
			}
		}
		elsif($line=~/^(\d+)\s+(\d+)\s+(\S+)\s+ubox$/)
		{
			my $i = $1 - 1;
			my $j = $2 - 1;
			my $p = $3 * $3;

			$ret[$i] += $p;
			$ret[$j] += $p;
		}
	}
	close ($fin);
	return \@ret;
}


1;



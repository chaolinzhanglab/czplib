package ViennaRNA;

require Exporter;

our $VERSION = 1.01;

@ISA = qw (Exporter);

@EXPORT = qw (
	readRNAplfoldFile
	readRNAduplexFile
);



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

	my $fin;
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

=head2 readRNAduplexFile

we assume target sequence goes first

=cut

sub readRNAduplexFile
{
    my ($inFile, $engine, $sort) = @_;

    Carp::croak "engine=$engine not recognized\n" unless $engine eq 'RNAplex' || $engine eq 'RNAduplex';

    my @ret;

    my $fin;
    open ($fin, "<$inFile") || Carp::croak "cannot open file $inFile to read\n";

    while (my $line = <$fin>)
    {
        next if $line=~/^\s*$/;
		next if $line=~/^\>/;
        $line =~/^(\S*?)\s+(\d+)\,(\d+)\s+\:\s+(\d+)\,(\d+)\s+\(\s*(\S*?)\).*$/;
        my %entry = (struct=>$1,
					targetStart=>$2,
                    targetEnd=>$3,
  					queryStart=>$4,
                    queryEnd=>$5,
                    score=>$6);

        push @ret, \%entry;
    }
    close ($fin);

	@ret = sort {$a->{'score'} <=> $b->{'score'}} @ret if $sort;

    return \@ret;
}

1;



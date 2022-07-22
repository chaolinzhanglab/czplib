package ASO;


require Exporter;


our $VERSION = 1.00;


@ISA = qw (Exporter);

@EXPORT = qw (
    scoreASO
	formatASO
);


=head1 NAME

ASO - subroutines to handle antisense oligos

subroutines starting with a hyphen should not be called outside

=cut

use strict;
use Data::Dumper;
use Carp;

use Sequence;


=head2 scoreASO

characterize ASO such as GC content, RNA secondary structure etc

=cut


sub scoreASO
{
    my $seq = $_[0];
    my $bc = baseComp ($seq);

    my $gc = $bc->{'G'} + $bc->{'C'};

    return {'GCContent'=>$gc};
}


=head2 formatASO

=cut

sub formatASO
{
    my ($seq, $format) = @_;

    #/52MOErT/*/i2MOErC/*/i2MOErT/*/i2MOErG/*/i2MOErA/*/i2MOErA/*/i2MOErG/*/i2MOErG/*/i2MOErG/*/i2MOErG/*/i2MOErG/*/i2MOErG/*/i2MOErC/*/i2MOErC/*/i2MOErC/*/i2MOErA/*/i2MOErC/*/i2MOErC/*/i2MOErT/*/32MOErT/    TCTGAAGGGGGGCCCACCTT

    my $formattedSeq = "";
    if ($format eq '2MOE_PS')
    {
        my @bases =split (//, $seq);
        my @bases2;
        for (my $i = 0; $i < @bases; $i++)
        {
            my $b = $bases[$i];

            Carp::croak "ambiguous base: $b\n" unless $b eq 'A' || $b eq 'C' || $b eq 'G' || $b eq 'T';

            if ($i == 0)
            {
                $b = "/52MOEr" . $b . "/";
            }
            elsif ($i == $#bases)
            {
                $b = "/32MOEr" . $b . "/";
            }
            else
            {
                $b = "/i2MOEr" . $b . "/";
            }
            push @bases2, $b;

            $formattedSeq = join ("*", @bases2);
        }
    }
    else
    {
        Carp::croak "wrong format = $format\n";
    }

    return $formattedSeq;
}


1;

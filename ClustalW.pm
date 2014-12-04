package ClustalW;

require Exporter;

our $VERSION = 1.01;

@ISA = qw (Exporter);

@EXPORT = qw (
    readClustalWFile
	writeClustalWFile
);



=head1 NAME

ClustalW - interface to interact Clustalw
subroutines starting with a hyphen should not be called outside

=cut

use strict;
use warnings;
use Data::Dumper;
use Carp;

use Common;


=head2

my $aln = readClustalWFile ($inFile)
read clustal multiple alignment file

=cut


sub readClustalWFile
{
	my $inFile = $_[0];
	my $fin;
	open ($fin, "<$inFile") || Carp::croak "cannot open file $inFile to read\n";

	my $status = "";
	
	my %sequences;
	my $identity = "";
	my $leadingSpace = 0; #for the identity line
	my $iter = 0;

	while (my $line=<$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		
		if ($line =~/^\#/ || $line=~/^\%/)
		{
			next;
		}
		elsif ($line =~/^\S/)
		{
			$line =~/^(\S*?)(\s+)(\S*?)$/;
			my ($id, $space, $seqStr) = ($1, $2, $3);
			$leadingSpace = length($id . $space);
			
			if (exists $sequences{$id})
			{
				$sequences{$id}{'seq'} .= $seqStr;
			}
			else
			{
				$sequences{$id} = {id=>$id, iter=>$iter++, seq=>$seqStr};
			}
		}
		elsif ($line =~/^\s/)
		{
			#identity line
			$identity .= substr ($line, $leadingSpace);
		}
	}	
	close ($fin);
	
	my @sequences = sort {$a->{'iter'} <=> $b->{'iter'}} values %sequences;
	return {leadingSpace=>$leadingSpace, sequences=>\@sequences, identity=>$identity};
}

sub writeClustalWFile
{
	my ($aln, $outFile) = @_;
	my $fout;

	open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";
	
	print $fout "CLUSTAL 2.1 multiple sequence alignment\n\n";	

	my $sequences = $aln->{'sequences'};
	my $leadingSpace = 0;
	if (exists $aln->{'leadingSpace'})
	{
		$leadingSpace = $aln->{'leadingSpace'};
	}
	else
	{
		foreach my $s (@$sequences)
		{
			$leadingSpace = max($leadingSpace, length($s->{'id'}));
		}
		$leadingSpace += 4;
	}
	
	foreach my $s (@$sequences)
	{
		my $id = $s->{'id'};
		my $seq = $s->{'seq'};

		$id = sprintf ("%-*s", $leadingSpace, $id);
		print $fout join ("", $id, $seq), "\n";
	}
	my $space = sprintf("%-*s", $leadingSpace, "");
	print $fout join ("", $space, $aln->{'identity'}), "\n";

	close ($fout);	
}


1;



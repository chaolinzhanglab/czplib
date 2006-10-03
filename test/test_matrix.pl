#!/usr/bin/perl -w

use strict;
use AnnotationIO;
use Bio::SeqIO;
use Common;

die "script <motif-file> <bgSeq> <fgseq>\n" if @ARGV != 3;

my ($motifFile, $bgSeqFile, $fgSeqFile) = @ARGV;

my $motifs = AnnotationIO::readMotifFile ($motifFile);

my $seqIO = Bio::SeqIO->new (-file=>$bgSeqFile, -format=>'Fasta');
my @seqs;
while (my $seq=$seqIO->next_seq())
{
	push @seqs, $seq->seq();
}

my $baseComp = Common::baseComp (\@seqs);
@seqs =();

print "Base Composition:\n";
print join ("\t", "A=".$baseComp->{"A"},
				  "C=".$baseComp->{"C"},
				  "G=".$baseComp->{"G"},
				  "T=".$baseComp->{"T"}), "\n";

$seqIO = Bio::SeqIO->new (-file=>$fgSeqFile, -format=>'Fasta');
while (my $seq=$seqIO->next_seq())
{
	push @seqs, $seq;
}


foreach my $motif (@$motifs)
{
	my $matrix = $motif->{"MT"};
	Common::ToStormoMatrix ($matrix, $baseComp);
	my $maxScore = Common::getMaxMatrixScore ($matrix);
	my $minScore = Common::getMinMatrixScore ($matrix);

	print "Matrix=", $motif->{"AC"}, "max score=$maxScore, min score=$minScore\n";
	foreach my $seq (@seqs)
	{
		my $score = Common::getMatrixScore ($matrix, $seq->seq());
		print join ("\t", $seq->id(), $seq->seq(), $score), "\n";
	}
}


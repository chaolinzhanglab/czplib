#
#===============================================================================
#
#         FILE:  Sam.pm
#
#  DESCRIPTION:  Package to handle sam files
#         BUGS:  ---
#        NOTES:  
#       AUTHOR:  Chaolin Zhang (cz), cz2294@columbia.edu
#      COMPANY:  Columbia University
#      VERSION:  1.0
#      CREATED:  12/24/2013
#     REVISION:  ---
#===============================================================================

package Sam;


require Exporter;

our $VERSION = 1.01;

@ISA = qw (Exporter);

@EXPORT = qw (
	lineToSam
	decodeSamFlag
	readSamFile
	samToLine
	samToBed
	writeSamFile
);



=head1 NAME

Sam - subroutines that deal with sequence alignment

subroutines starting with a hyphen should not be called outside

=cut

use strict;
use warnings;
use Data::Dumper;
use Carp;


use Common;

sub readSamFile
{
	my ($inFile, $verbose) = @_;
	my $fin;
	open ($fin, "<$inFile") || Carp::croak "cannot open file $inFile to read\n";
	my @ret;
	my $iter = 0;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		next if $line =~/^\@/;

		print "$iter ...\n" if $verbose && $iter % 100000 == 0;
		$iter++;
		my $s = lineToSam ($line);
		push @ret, $s;
	}
	close ($fin);
	return \@ret;
}

sub writeSamFile
{
	my ($sam, $outFile) = @_;
	my $fout;
	open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";
	foreach my $s (@$sam)
	{
		print $fout samToLine ($s), "\n";
	}
	close ($fout);
}



=head2 lineToSam

very light-weight function, with essentially no decoding

=cut

sub lineToSam
{
	my $line = $_[0];
	my ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, $MRNM, $MPOS, $ISIZE, $SEQ, $QUAL, $TAGS) = split (/\s+/, $line, 12);
	return {
	QNAME => $QNAME,
	FLAG=> $FLAG,
	RNAME=>$RNAME,
	POS=>$POS,
	MAPQ=>$MAPQ,
	CIGAR=>$CIGAR,
	MRNM=>$MRNM,
	MPOS=>$MPOS,
	ISIZE=>$ISIZE,
	SEQ=>$SEQ,
	QUAL=>$QUAL,
	TAGS=>$TAGS
	};
}

sub samToLine
{
	my $sam = $_[0];

	my $line = join ("\t", 
		$sam->{'QNAME'},
		$sam->{'FLAG'},
		$sam->{'RNAME'},
		$sam->{'POS'},
		$sam->{'MAPQ'},
		$sam->{'CIGAR'},
		$sam->{'MRNM'},
		$sam->{'MPOS'},
		$sam->{'ISIZE'},
		$sam->{'SEQ'},
		$sam->{'QUAL'},
		$sam->{'TAGS'});
	return $line;
}

=head2

my $bed = samToBed ($sam)
or 
my $bed = samToBed ($sam, $useRNAStrand)

return 0 if no alignment
besides all columns for bed format, an additional key is flagInfo, which is decoded flags from SAM format

#TODO: deal with soft clipping

#NOTE: sam2bed.pl in OLego has a separate copy of this function.  Any changes here need to be incoporated in that program
=cut

sub samToBed
{
	my ($sam, $useRNAStrand) = @_;
	$useRNAStrand = 0 unless $useRNAStrand;

	return 0 if $sam->{"CIGAR"} eq '*'; #no alignment
	
	my $flagInfo = decodeSamFlag ($sam->{"FLAG"});
	
	my $strand = $flagInfo->{'query_strand'};
	
	my $TAGS = "";
	$TAGS = $sam->{"TAGS"} if $sam->{"TAGS"};

	if ($useRNAStrand)
	{
		if ($TAGS=~/XS\:\S*\:([-+\.])/)
		{
			$strand = $1;
			$strand = '+' if $strand eq '.';
		}
	}
	my $read1_or_2 = $flagInfo->{'read_1_or_2'};
	
	my $name = $sam->{"QNAME"};
	my $chrom = $sam->{"RNAME"};
	my $chromStart = $sam->{"POS"} - 1;

	my $score = 0;
	if ($TAGS=~/NM\:\S*\:(\d+)/)
	{
		$score = $1;
	}

	my $CIGAR = $sam->{"CIGAR"};
	my $QNAME = $sam->{"QNAME"};
	my $SEQ = $sam->{"SEQ"};

    #remove soft cliped nucleotides
    if ($CIGAR =~/^\d+S(.*?)$/)
    {
        $CIGAR = $1;
    }

    if ($CIGAR =~/^(.*?)\d+S$/)
    {
        $CIGAR = $1;
    }

	#deal with the rest
	if ($sam->{"CIGAR"}=~/[^\d+|M|N|I|D]/g)
	{
		Carp::croak "unexpected CIGAR string: $CIGAR in $QNAME: $SEQ\n";
	}

	my (@blockSizes, @blockStarts);
	
	my $currLen = 0;
	my $extendBlock = 0;

	while ($CIGAR =~/(\d+)([M|N|I|D])/g)
	{
		my ($size, $type) = ($1, $2);
		if ($type eq 'I' || $type eq 'D')
		{
			#insertion in reads
			$extendBlock = 1;
			if ($type eq 'D')
			{
				my $n = @blockSizes;
				if ($n < 1)
				{
					$chromStart += $size;	
				}
				else
				{
					$blockSizes[$#blockSizes] += $size;
					$currLen += $size;
				}
			}
			next;
		}

		if ($type eq 'M')
		{
			if ($extendBlock && @blockSizes > 0)
			{
				#extend the previous block
				my $n = @blockSizes;
				$blockSizes[$n-1] += $size;
			}
			else
			{
				push @blockSizes, $size;
				push @blockStarts, $currLen;
			}
			$extendBlock = 0;
		}
		$currLen += $size;
	}
	
	my $blockCount = @blockSizes;
	my $chromEnd = $chromStart + $blockStarts[$blockCount-1] + $blockSizes[$blockCount-1] - 1;
	
	my $bed = {
		chrom=>$chrom,
		chromStart=>$chromStart,
		chromEnd=>$chromEnd,
		name=>$name,
		score=>$score,
		strand=>$strand,
		thickStart=>$chromStart,
		thickEnd=>$chromEnd,
		itemRgb=>0,
		blockCount=>$blockCount,
		blockSizes=>\@blockSizes,
		blockStarts=>\@blockStarts,
		flagInfo=>$flagInfo
	};
}	



sub decodeSamFlag
{
	my $flag = $_[0];
	$flag = sprintf ("%012b", $flag);
	my @flags = split (//, $flag);

	my $flagInfo = {
		PE=>$flags[11],					#1 means paired-end data
		PE_map=>$flags[10],				#1 means each end is properly aligned according to the aligner
		query_nomap=>$flags[9],				#1 means this read is unmapped
		mate_nomap=>$flags[8],				#1 means its mate is unmapped
		query_strand=>$flags[7] == 0 ? '+' : '-',	#1 means the strand of this read is on the negative strand
		mate_strand=>$flags[6] == 0 ? '+' : '-',	#1 means its mate is on the negative strand
		read_1_or_2=> $flags[5] == 1 ? 1 : 2 		#1 means this is read1
	};
	return $flagInfo;
}	




1;



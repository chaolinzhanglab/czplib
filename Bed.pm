#
#===============================================================================
#
#         FILE:  Bed.pm
#
#  DESCRIPTION:  Package to handle Bed file
#         BUGS:  ---
#        NOTES:  The package is spun off from the AnnotationIO.pm
#       AUTHOR:  Chaolin Zhang (cz), czhang@rockefeller.edu
#      COMPANY:  Rockefeller University
#      VERSION:  1.0
#      CREATED:  12/16/10 20:58:53
#     REVISION:  ---
#===============================================================================


package Bed;

require Exporter;

@ISA = qw (Exporter);

@EXPORT = qw (
	bedToFull
	bedToLine
	clusterRegions
	collapseReads
	combineRegions
	contigToGenome
	copyBedRegion
	geneToExon
	genomeToContig
	getUniqPaths
	lineToBed
	printBedRegion
	readBedFile
	readNextBedBlock
	readNextBedLine
	removeSmallGap
	segmentRegion
	sortBedFile
	splitBedFileByChrom
	writeBedFile
);
	
=head1 Bed

Subroutines to handle BED file

=cut

use strict;
use warnings;

use Data::Dumper;
use Carp;



#
#object-oriented interface will be added later
#

my $debug = 0;


=head2 readBedFile

Read a Bed File. Only one track is allowed in the file now and header will be ignored

my $regions = readBedFile ($inFile, $verbose);

inFile [string]: input Bed file
verbose        : verbose mode

return         : reference to an array
=cut

sub readBedFile
{
	my ($inBedFile, $verbose) = @_;

	my $fin;
	my @ret;

	open ($fin, "<$inBedFile")||Carp::croak "cannot open file $inBedFile to read\n";
	
	my $i = 0;
	while (my $bedLine = readNextBedLine ($fin))
	{
		print "$i ...\n" if $verbose && $i % 100000 == 0;
		$i++;
		push @ret, $bedLine;
	}
	close ($fin);
	return \@ret;
}

=head2 readNextBedBlock

Read the next "block" of regions. 
*Regions in each block are separated by no more than $maxGap
*Regions of different blocks are separated by more than $maxGap

if $minBlockSize is specified
neighboring blocks might be merged if the size is too small


Usage:
my $block = readNextBedBlock ($fin, $maxGap, $separateStrand, %params);

Important notes: the input file has to be sorted, first by strand (if $separateStrand ==1),
                 then by chrom, chromStart, and chromEnd

$fin                  : File handle
$maxGap         [int] : maxGap allowed between regions in a block
$separateStrand [bool]: regions on different strands were considered separately
%params               : minBlockSize=> N, min block size 

=cut

sub readNextBedBlock
{
		my ($fin, $maxGap, $separateStrand, %params) = @_;

		my $currBlockEnd = -1;
		my @bedBlock;

		my $minBlockSize = exists $params{'minBlockSize'} ? $params{'minBlockSize'} : 0;


		Carp::croak "gap cannot be negative\n" if $maxGap < 0;

		my $currPointer = tell ($fin);

		while (my $bedLine = readNextBedLine ($fin))
		{
			if ($separateStrand)
			{
				Carp::croak "no strand information in ", Dumper ($bedLine), "\n" unless exists $bedLine->{'strand'};
			}
			
			Carp::croak "negative coordinates in ", Dumper ($bedLine), "\n" if $bedLine->{'chromEnd'} < 0;

			my $expand = 0;

			if ($currBlockEnd < 0 || @bedBlock < $minBlockSize) #this is the first line of the block
			{
				$expand = 1;
			}
			else
			{
				if ($separateStrand)
				{
					$expand = 1 if $bedBlock[$#bedBlock]->{'strand'} eq $bedLine->{'strand'} 
						&& $bedBlock[$#bedBlock]->{'chrom'} eq $bedLine->{'chrom'} 
						&& $bedLine->{'chromStart'} - $currBlockEnd  - 1 <= $maxGap;
				}
				else
				{
					$expand = 1 if $bedBlock[$#bedBlock]->{'chrom'} eq $bedLine->{'chrom'} 
						&& $bedLine->{'chromStart'} - $currBlockEnd  - 1 <= $maxGap;
				}
			}

			if ($expand)
			{	#expand the block
				push @bedBlock, $bedLine;
				$currBlockEnd = $bedLine->{'chromEnd'} if $bedLine->{'chromEnd'} > $currBlockEnd;
				$currPointer = tell ($fin);
			}
			else
			{
				#reached the end of the current block
				seek ($fin, $currPointer, 0);
				return \@bedBlock;
			}
		}
		return @bedBlock > 0 ? \@bedBlock : 0;#the last block
}


=head2 readNextBedLine

read the next bed line
return zero when it reaches the end

Usage: 

my $bedLine = readNextBedLine ($fin)

$fin: file handle

=cut

sub readNextBedLine
{
		my $fin = $_[0];
		while (my $line = <$fin>)
		{
				chomp $line;
				next if $line =~/^\s*$/;
				next if $line =~/^\#/;
				next if $line =~/^track/;
				next if $line =~/^browser/;

				return lineToBed ($line);
		}
		return "";
}


=head2 lineToBed

parse a line to a BED region
any column is fine, but it has to be in correct format

=cut

sub line2bed
{
	Carp::croak "obsolete function, call lineToBed instead\n";
}

sub lineToBed
{
	my $line = $_[0];

	my @cols = split (/\s+/, $line);
	Carp::croak "less than three columns in line: $line\n" if @cols < 3;
	
	#	print join ("\t", @cols), "\n";
	my @colNames = qw (chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts);
	my $i;
	my $entry = {};
	for ($i = 0; $i < @colNames; $i++)
	{
		last if ($#cols < $i);
		$entry->{$colNames[$i]} = $cols[$i];
	}

	if (exists $entry->{"blockSizes"})
	{
		my $blockSizes = $entry->{"blockSizes"};
		my @bs = split (/\,/, $blockSizes);
		Carp::croak "in correct number of blocks at line $line\n" if $#bs != $entry->{"blockCount"} - 1;
		$entry->{"blockSizes"} = \@bs;
	}

	if (exists $entry->{"blockStarts"})
	{
		my $blockStarts = $entry->{"blockStarts"};
		my @bs = split (/\,/, $blockStarts);
		Carp::croak "in correct number of blocks at line $line\n" if $#bs != $entry->{"blockCount"} - 1;
		$entry->{"blockStarts"} = \@bs;
	}
		
	$entry->{"chromEnd"} -= 1;
	$entry->{"thickEnd"} -= 1 if (exists $entry->{"thickEnd"});
		
	Carp::croak "chromStart (" . $entry->{"chromStart"} . ")  > chromEnd (" . $entry->{"chromEnd"} . ")\n" 
	if ($entry->{"chromStart"} > $entry->{"chromEnd"});
		
		#print join ("\t", $entry->{"chromStart"}, $entry->{"chromEnd"}), "\n";
	return $entry;
}


=head2 writeBedFile

Usage:

writeBedFile ($regions, $out, $header = "", $append = 0);
append  [string] :  'a' means append


Note: the order of parameters changed on 12/17/2010

=cut

sub writeBedFile
{
	
	my ($regions, $out, $header, $append) = @_;
	my $fout;
	
	Carp::croak "empty output file name\n" unless $out;
	
	if ($append && $append eq 'a')
	{
		open ($fout, ">>$out") || Carp::croak "cannot open file $out to append\n";
	}
	else
	{
		open ($fout, ">$out") || Carp::croak "cannot open file $out to write\n";
	}

	if (@$regions <1)
	{
		close ($fout);
		return;
	}

	print $fout $header, "\n" if $header && $header =~/^track/;
	
	map {print $fout bedToLine ($_), "\n";} @$regions;

	close ($fout);
}


=head2 sortBedFile

sort regions by strand (if $separateStrand == 1), then by
chrom, chromStart, and chromEnd

Note:  1. it depends on grep and sort
       2. the input file is assumed to have NO header (to be improved in the future)
       3. input and output files cannot be the same file
=cut

sub sortBedFile
{
		my ($inFile, $outFile, $separateStrand) = @_;
		my $ncols = `head -n 1 $inFile | awk '{print NF}'`; chomp $ncols;
		Carp::croak "less than 6 columns\n" if $separateStrand && $ncols < 6;

		my $cmd = "grep -v \"^track\"  $inFile | grep -v \"^#\" | sort";
		$cmd .= " -k 6,6" if $separateStrand;
		$cmd .= " -k 1,1 -k 2,2n -k 3,3n > $outFile";

		#print $cmd, "\n" if $verbose;
		my $ret = system ($cmd);
		Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;
}

=head2 splitBedFileByChrom

Usage:

my $info = splitBedFileByChrom ($inBedFile, $outDir, %params)

 inBedFile [string]:
 outDir    [string]:
 %params           :
	=> v        : verbose mode (0|1)
    => sort     : sort mode
                  sort == 0: no sort
                  sort == 1: sort, but with the two strands together
                  sort == 2: sort, with the two strands separate

Return: 
 $info->{$chrom}->{n=>$tagNum, f=>$chrFile}
 	

=cut


sub splitBedFileByChrom
{
	#my ($inBedFile, $outDir, $verbose) = @_;

	my $inBedFile = shift @_;
	my $outDir = shift @_;

	#optional parameters
	my $verbose = 0; #verbose mode
	my $sort = 0;	#sort by chrom, then by chromStart, and then by chromEnd, if sort ==2, sort by strand first
	my %params;

	if (@_ == 1)
	{
		$verbose = $_[0]; #this is for back compatibility
	}
	else
	{
		%params = @_;
		$verbose = $params{'v'} if exists $params{'v'};
		$sort = $params {'sort'} if exists $params {'sort'};
	}

	Carp::croak "$inBedFile does not exist\n" unless -f $inBedFile;
	Carp::croak "$outDir does not exist\n" unless -d $outDir;

	my %tagCount;

	my $fin = new FileHandle;
	#print "reading tags from $inBedFile ...\n" if $verbose;
	open ($fin, "<$inBedFile") || Carp::croak "can not open file $inBedFile to read\n";

	#my %tagCount;
	#write tags according to chromosomes
	my $i = 0;
	my %fhHash;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		next if $line=~/^track name/;
		next if $line=~/^browser/;

		$i++;

		print "$i ...\n" if $i % 100000 == 0 && $verbose;
	
		$line =~/(\S+)\s/;

		my $chrom = $1;
		#print "chrom = $chrom\n";
		my $tmpFile = "$outDir/$chrom.bed";
	
		if (not exists $fhHash{$chrom})
		{
			my $fout = new FileHandle;
			open ($fout, ">>$tmpFile") || Carp::croak "can not open file $tmpFile to write\n";
			$fhHash{$chrom} = $fout;
		}
	
		#open ($fout, ">>$tmpFile") || Carp::croak "can not open file $tmpFile to write\n";
	
		my $fout = $fhHash{$chrom};
		print $fout $line, "\n";
		#close ($fout);

		if (exists $tagCount{$chrom})
		{
			$tagCount{$chrom}->{"n"} ++;
		}
		else
		{
			$tagCount{$chrom} = {f=>"$outDir/$chrom.bed", n=> 1};
		}
	}
	close ($fin);
	#close all file handles
	foreach my $chrom (keys %fhHash)
	{
		close ($fhHash{$chrom});

		if ($sort)
		{
				print "sorting $chrom ...\n" if $verbose;
				my $f = $tagCount{$chrom}->{'f'};
				my $f2 = "$f.sort";

				my $separateStrand = $sort > 1 ? 1 : 0;
				sortBedFile ($f, $f2, $separateStrand);
				system ("mv $f2 $f");
		}
	}
	return \%tagCount;
}


=head2 bed2Full

expand to 12 column format

Usage: 
my $r = bed2Full ($r)

=cut

sub bed2Full
{
	Carp::croak "obsolete function, call bedToFull instead\n";
}

sub bedToFull
{
	my $region = $_[0];
	Carp::croak "at least three columns should exist\n" 
	unless exists $region->{"chrom"} && exists $region->{"chromStart"} && exists $region->{"chromEnd"};
	
	$region->{"name"} = $region->{"chrom"} . ":" . $region->{"chromStart"} . "-" . ($region->{"chromEnd"} + 1)
	unless exists $region->{"name"};

	$region->{"score"} = 0 unless exists $region->{"score"};
	$region->{"strand"} = "+" unless exists $region->{"strand"};

	$region->{"thickStart"} = $region->{"chromStart"} unless exists $region->{"thickStart"};
	$region->{"thickEnd"} = $region->{"chromEnd"} unless exists $region->{"thickEnd"};
	$region->{"itemRgb"} = "0,0,0" unless exists $region->{"itemRgb"};

	$region->{"blockCount"} = 1 unless exists $region->{"blockCount"};
	my @blockSizes = ($region->{"chromEnd"} - $region->{"chromStart"} + 1);
	$region->{"blockSizes"} = \@blockSizes unless exists $region->{"blockSizes"};
	
	my @blockStarts = (0);
	$region->{"blockStarts"} = \@blockStarts unless exists $region->{"blockStarts"};

	return $region;
}


=head2 copyBedRegion

make an independent copy of a region
my $to = copyRegion ($from);

=cut

sub copyBedRegion
{
	my $region = $_[0];
	my %regionCopy = %$region;
	if (exists $regionCopy{'blockStarts'})
	{
		my @blockStarts = %{$regionCopy{'blockStarts'}};
		$regionCopy{'blockStarts'}= \@blockStarts;
	
		my @blockSizes = %{$regionCopy{'blockSizes'}};
		$regionCopy{'blockSizes'} = \@blockSizes;
	}
	return \%regionCopy;
}


=head2 printBedRegion

printBedRegion ($region);

=cut

sub printBedRegion
{
	my $region = $_[0];
	my $str = bedToLine ($region);
	print $str, "\n";
}


=head2 bedToLine

The same as printBedRegion
bed2line ($region);

=cut

sub bed2line
{
	Carp::croak "obsolete function, call bedToLine instead\n";
	
}

sub printBedRegionToString
{
	Carp::croak "obsolete function, call bedToLine instead\n";
}

sub bedToLine
{
	my $region = $_[0];
	
	my @colNames = qw (chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts);
	
	my $colNum = 12; #keys %$region;
	
	my %rCopy = %$region;
	$rCopy{"chromEnd"} += 1;
	if (exists $rCopy{'thickEnd'})
	{
		$rCopy{'thickEnd'} += 1;
	}
	
	if (exists $rCopy{'blockCount'})
	{
		Carp::croak "no blockSizes\n" unless exists $rCopy {'blockSizes'};
		Carp::croak "no blockStarts\n" unless exists $rCopy {'blockStarts'};
			
		$rCopy{'blockSizes'} = join (",", @{$rCopy{'blockSizes'}});
		$rCopy{'blockStarts'} = join (",", @{$rCopy{'blockStarts'}});
	}
	
	my $ret = join ("\t", $rCopy{"chrom"}, $rCopy{"chromStart"}, $rCopy{"chromEnd"});
	for (my $i = 3; $i < $colNum; $i++)
	{
		my $col = $colNames[$i];
		if (exists $rCopy{$col})
		{
			$ret .= "\t" . $rCopy{$col};
		}
		else
		{
			last;
			#Carp::croak "col=$col is not defined\n"; 
		}
	}
	return $ret;
}



=head2 removeSmallGap

remove small gaps in a BED region, which are typically introduced by insertions
 or deletions during alignemnt

my $r2 = removeSmallGap ($r, $gapSize);

=cut

sub removeSmallGap
{
	my ($r, $gapSize) = @_;
	return $r unless exists $r->{"blockCount"} && $r->{"blockCount"} > 1;

	my @blockStarts;
	my @blockSizes;

	my $currBlockStart = 0; 
	my $currBlockSize = $r->{"blockSizes"}->[0];

	for (my $i = 1; $i < @{$r->{"blockStarts"}}; $i++)
	{
		my $currBlockEnd = $currBlockStart + $currBlockSize - 1;
		if ($r->{"blockStarts"}->[$i] > $currBlockEnd + $gapSize + 1) # a new block
		{
			push @blockStarts, $currBlockStart;
			push @blockSizes, $currBlockSize;
			$currBlockStart = $r->{"blockStarts"}->[$i];
			$currBlockSize = $r->{"blockSizes"}->[$i];
		}
		else
		{
			#extend the current block
			$currBlockSize = $r->{"blockStarts"}->[$i] + $r->{"blockSizes"}->[$i] - $currBlockStart;
		}
	}
	push @blockStarts, $currBlockStart;
	push @blockSizes, $currBlockSize;
	
	$r->{"blockCount"} = @blockStarts;
	$r->{"blockSizes"} = \@blockSizes;
	$r->{"blockStarts"} = \@blockStarts;

	return $r;
}

=head2 geneToExon

=cut


sub gene2exon
{
	Carp::croak "obsolete function, call geneToExon instead\n";
}

sub geneToExon
{
    my $g = $_[0];
    my $nexon = $g->{"blockCount"};
    my @exons;
    for (my $i = 0; $i < $nexon; $i++)
    {
        my $exonStart = $g->{"chromStart"} + $g->{"blockStarts"}->[$i];
        my $exonEnd = $exonStart + $g->{"blockSizes"}->[$i] - 1;
        my $e = {
            chrom => $g->{"chrom"},
            chromStart=> $exonStart,
            chromEnd => $exonEnd,
            name => join (":", $g->{"name"}, $i),
            score => $g->{"score"},
            strand => $g->{"strand"}
        };
        push @exons, $e;
    }
    return \@exons;
}


=head2 contigToGenome

map coordinates on a contig to the genome


my $genomeCoord = contigToGenome ($genomeCoordofContig, $startOnContig, $endOnContig)

$genomeCoordofConfig is a simple interval
$contig->{chrom=>chr1, chromStart=>2, chromEnd=>3, strand=>4}

$zero based coordinates


=cut

sub contig2genome
{
	Carp::croak "obsolete function, call contigToGenome instead\n";
}

sub contigToGenome
{
	my ($contig, $start, $end) = @_;
	my $contigChrom = $contig->{"chrom"};
	my $contigStart = $contig->{"chromStart"};
	my $contigEnd = $contig->{"chromEnd"};
	my $contigStrand = $contig->{"strand"};

	my ($genomeStart, $genomeEnd);
	if ($contigStrand eq '+')
	{
		$genomeStart = $contigStart + $start;
		$genomeEnd = $contigStart + $end;
	}
	else
	{
		$genomeEnd = $contigEnd - $start;
		$genomeStart = $contigEnd - $end;
	}
	return {chrom=> $contigChrom, 
			strand=>$contigStrand, 
			chromStart=>$genomeStart, 
			chromEnd=>$genomeEnd};
}

=head2 genomeToContig

genome coordinates to contig coordinates
$contig->{chrom=>chr1, chromStart=>2, chromEnd=>3, strand=>4}
$zero based coordinates

=cut

sub genome2contig
{
	Carp::croak "obsolete function, call genomeToContig instead\n";	
}

sub genomeToContig
{
	my ($contig, $tsOnGenomeStart, $tsOnGenomeEnd) = @_;
	my $contigChrom = $contig->{"chrom"};
	my $contigStart = $contig->{"chromStart"};
	my $contigEnd = $contig->{"chromEnd"};
	my $contigStrand = $contig->{"strand"};

	my ($tsOnContigStart, $tsOnContigEnd);
	if ($contigStrand eq '+')
	{
		$tsOnContigStart = $tsOnGenomeStart - $contigStart;
		$tsOnContigEnd = $tsOnGenomeEnd - $contigStart;
	}
	else
	{
		$tsOnContigStart = $contigEnd - $tsOnGenomeEnd;
		$tsOnContigEnd = $contigEnd - $tsOnGenomeStart;
	}
	return {chrom=> $contig->{"name"}, 
			strand=>$contigStrand, 
			chromStart=>$tsOnContigStart, 
			chromEnd=>$tsOnContigEnd};
}

=head2 getUniqPaths

get unique paths(transcripts segments) between two coordinates from a set of paths(transcripts)

my $uniqPath = getUniqPaths ($paths, $chromStart, $chromEnd)

$paths: a reference to an array, each element is a bed row
$chromStart and $chromEnd: zero-based coordinates of the query interval 

Return:  an array ref. to the unique paths (transcript segments) between 
$chromStart and $chromEnd

=cut
sub getUniqPaths
{
	my ($paths, $chromStart, $chromEnd) = @_;
	
	my %uniqPaths;
	foreach my $p (@$paths)
	{
		my $segment = segmentRegion ($p, $chromStart, $chromEnd);
		next unless ref($segment) && $segment->{"chromStart"} == $chromStart && $segment->{"chromEnd"} == $chromEnd;

		if (not exists ($segment->{"blockCount"}))
		{
			$segment = bed2Full($segment);
		}

		#with the same chromStart and end, block sizes and starts are sufficient to distinguish diffrerent clusters

		my $key = join ("-", @{$segment->{"blockSizes"}},"//", @{$segment->{"blockStarts"}});
		$uniqPaths{$key} = $segment;
	}
	my @uniqPaths = values %uniqPaths;
	return \@uniqPaths;	
}


=head2 segmentRegion
#get exon/intron structure between two coordinates specifying an interval of interest

my $segment = segmentRegion ($ts, $chromStart, $chromEnd)

if both $chromStart and $chromEnd are on exons, the returned region should have 
the same $chromStart and $chromEnd as the input.  Otherwise, the intronic part 
at the end will be truncated.

=cut

sub segmentRegion
{ 
	my ($ts, $chromStart, $chromEnd) = @_;

	Carp::croak "region beyond the start and end of transcript:", Dumper ($ts), "\n" 
	unless $ts->{"chromStart"} <= $chromStart && $ts->{"chromEnd"} >= $chromEnd;

	my $blockCount = exists $ts->{'blockCount'} ? $ts->{'blockCount'} : 1;
	if ($blockCount == 1)
	{
		my %tsNew = %$ts;

		$tsNew{"chromStart"} = $chromStart;
		$tsNew{"chromEnd"} = $chromEnd;
		$tsNew{"thickStart"} = $chromStart if exists $ts->{"thickStart"};
		$tsNew{"thickEnd"} = $chromEnd if exists $ts->{"thickEnd"};
		$tsNew{"blockSizes"} = [$chromEnd - $chromStart + 1] if exists $ts->{"blockSizes"};
		return \%tsNew;
	}

	#there are multiple blocks
	
	my $firstBlockIdx = -1;
	my $lastBlockIdx = -1;
	for (my $i = 0; $i < $ts->{"blockCount"}; $i++)
	{
		my $exonStart = $ts->{"chromStart"} + $ts->{"blockStarts"}->[$i];
		my $exonEnd = $ts->{"chromStart"} + $ts->{"blockStarts"}->[$i] + $ts->{"blockSizes"}->[$i] - 1;
		
		if ($exonEnd >= $chromStart)
		{
			$firstBlockIdx = $i;
			last;
		}
	}

	for (my $i = $ts->{"blockCount"} - 1; $i >= 0; $i--)
	{
		my $exonStart = $ts->{"chromStart"} + $ts->{"blockStarts"}->[$i];
		my $exonEnd = $ts->{"chromStart"} + $ts->{"blockStarts"}->[$i] + $ts->{"blockSizes"}->[$i] - 1;
		if ($exonStart <= $chromEnd)
		{
			$lastBlockIdx = $i;
			last;
		}
	}
	
	return 0 unless $lastBlockIdx >= $firstBlockIdx; 
	#otherwise, $chromStart and $chromEnd will specify a region located in an intron and spanning no exons
	
	my %tsNew = (chrom=>$ts->{"chrom"},
				chromStart=>$chromStart,
				chromEnd=>$chromEnd,
				name=>$ts->{"name"},
				score=>$ts->{"score"},
				strand=>$ts->{"strand"});
	
	$blockCount = $lastBlockIdx - $firstBlockIdx + 1;

	#get absolute start
	my @blockStartsNew = map {$_ + $ts->{"chromStart"}} @{$ts->{"blockStarts"}}[$firstBlockIdx .. $lastBlockIdx];
	my @blockSizesNew = @{$ts->{"blockSizes"}}[$firstBlockIdx .. $lastBlockIdx];

	if ($chromStart < $blockStartsNew[0]) 
	{
		#start is in intronic region
		$tsNew{'chromStart'} = $blockStartsNew[0];
	}
	else
	{
		#start is in exonic region
		$blockSizesNew[0] = $blockStartsNew[0] + $blockSizesNew[0] - $chromStart;
		$blockStartsNew[0] = $chromStart;
	}

	if ($chromEnd > $blockStartsNew[$blockCount - 1] + $blockSizesNew[$blockCount - 1] - 1)
	{
		#end is in intronic region
		$tsNew{'chromEnd'} = $blockStartsNew[$blockCount - 1] + $blockSizesNew[$blockCount - 1] - 1;
	}
	else
	{
		$blockSizesNew[$blockCount - 1] = $chromEnd - $blockStartsNew[$blockCount -1] + 1;
	}
	@blockStartsNew = map {$_ - $tsNew{'chromStart'}} @blockStartsNew;
	
	$tsNew{'itemRgb'} = $ts->{'itemRgb'};
	$tsNew{'thickStart'} = $tsNew{'chromStart'};
	$tsNew{'thickEnd'} = $tsNew{'chromEnd'};
	$tsNew{'blockCount'} = $lastBlockIdx - $firstBlockIdx + 1;
	$tsNew{'blockStarts'} = \@blockStartsNew;
	$tsNew{'blockSizes'} = \@blockSizesNew;

	return \%tsNew;
}


=head2 combineRegion

Note: change of interface, pass by reference now
combine multiple regions (on the same strand) to get the union (in terms of exons)

my $ts = combineRegions (\@regions)

=cut

sub combineRegion
{
	Carp::croak "obsolete function, call combineRegions\n";
}

sub combineRegions
{
	my $regions = $_[0];
	my $nregions = @$regions;
	
	my @blocks;

	my $region1 = $regions->[0];

	#my @regions = ($region1, $region2);

	foreach my $region (@$regions)
	{
		if (exists $region->{"blockStarts"})
		{
			for (my $i = 0; $i < @{$region->{"blockStarts"}}; $i++)
			{
				push @blocks, {start=>$region->{"chromStart"} + $region->{"blockStarts"}->[$i], 
						end=> $region->{"chromStart"} + $region->{"blockStarts"}->[$i] + $region->{"blockSizes"}->[$i] - 1};
			}
		}
		else
		{
			push @blocks, {start=>$region->{"chromStart"}, end=>$region->{"chromEnd"}};
		}
	}

	@blocks = sort {$a->{"end"} <=> $b->{"end"}} @blocks;
	@blocks = sort {$a->{"start"} <=> $b->{"start"}} @blocks;


	my @clusters;
	
	my $currClusterStart = -1;
	my $currClusterEnd = - 1;
	
	#my $n = @blocks;
	my $i = 0;

	foreach my $b (@blocks)
	{
		$i++;
		my $chromStart = $b->{"start"};
		my $chromEnd = $b->{"end"};

		my $openNewCluster = 0;

		$openNewCluster = $chromStart > $currClusterEnd;
		if ($openNewCluster)
		{
			if ($currClusterEnd >= 0)
			{
				push @clusters, {start=>$currClusterStart, end=>$currClusterEnd};
			}

			$currClusterStart = $chromStart;
			$currClusterEnd = $chromEnd;
		}
		else
		{
			$currClusterEnd = $chromEnd if $currClusterEnd < $chromEnd;
		}

	}
	
	if ($currClusterEnd >= 0)
	{
		push @clusters, {start=>$currClusterStart, end=>$currClusterEnd};
	}

	my $n = @clusters;

	my $chromStart =  $clusters[0]->{"start"};
	my $chromEnd = $clusters[$n-1]->{"end"};

	my $blockCount = @clusters;

	my @blockStarts;
	my @blockSizes;
	for (my $i = 0; $i < $n; $i++)
	{
		$blockStarts[$i] = $clusters[$i]->{"start"} - $chromStart;
		$blockSizes[$i] = $clusters[$i]->{"end"} - $clusters[$i]->{"start"} + 1;
	}


	my $region = {
		chrom => $region1->{"chrom"},
		chromStart => $chromStart,
		chromEnd => $chromEnd, 
		name => $region1->{"name"},
		score => $nregions, #$region1->{"score"} + $region2->{"score"},
		strand=> $region1->{"strand"},
		thickStart=> $chromStart,
		thickEnd=> $chromEnd,
		itemRgb=> "0,0,0",
		blockCount=> $blockCount,
		blockSizes=> \@blockSizes,
		blockStarts=>\@blockStarts
	};

	return $region;
}


=head2 clusterRegions
cluster regions on the same chromosome

my $clusters = clusterRegions ($regionsOnChrom, $strand, $maxGap, $overlapFraction, $collapse); 

regionsOnChrom : regions on a chrom, bed format
strand         : consider only tags on the give strand, could be +, - b (both)
maxGap         : max gap allowed for regions in a cluster (minimum overlap if maxGap < 0)
overlapFraction: if minGap < 0, minimum fraction of overlap to consider a match
collapse       : 0 (no collapse), 1 (exact match), 2 (if one read are contained by other)

return         : indices in the same clusters are put in an array


=cut

sub clusterRegions
{
	my ($regionsOnChrom, $strand, $maxGap, $overlapFraction, $collapse) = @_;
	
	my $minOverlap = -1;#no requirement on overlap, when maxGap > 0

	if ($maxGap < 0)
	{
		$minOverlap = -$maxGap; #with requirement on overlap
		$maxGap = 0;
	}

	my @regionsOnChrom = @$regionsOnChrom; #make a copy

	#my $chrom = $regionsOnChrom[0]->{"chrom"};
	for (my $i = 0; $i < @regionsOnChrom; $i++)
	{
		my $r = $regionsOnChrom[$i];
		#Carp::croak "regions are on different chromosomes\n" unless $chrom eq $r->{"chrom"};
		$r->{"idx"} = $i; #add indices in the copy, not the original
	}

	#sort ascendingly according to the start point
	#perl used a stable sort, so we want if the start point the same, the data is sorted according to chromEnd
	#this is important for the collpase mode
	
	my @regionsOnChromSorted = sort {$a->{"chromEnd"} <=> $b->{"chromEnd"}} @regionsOnChrom;
	@regionsOnChromSorted = sort {$a->{"chromStart"} <=> $b->{"chromStart"}} @regionsOnChromSorted;

	my $n = @regionsOnChromSorted;

	print "\n\nclustering begins with $n regions...\n" if $debug;
	
	my @clusters;
	
	my @currCluster;
	my $currClusterEnd = -$maxGap - 1;
	my $currClusterStart = -1;
	
	my $i = 0;
	foreach my $r (@regionsOnChromSorted)
	{
		#next if ($strand ne 'b' && $strand ne $r->{"strand"});
		$i++;

		next if $strand ne 'b' && $strand ne $r->{"strand"};

		Carp::croak "negative coordinates\n", Dumper ($r), "\n" if $r->{'chromStart'} < 0 || $r->{'chromEnd'} < 0;
		my $chromStart = $r->{"chromStart"};
		my $chromEnd = $r->{"chromEnd"};

		if ($maxGap > 0)
		{
			#we extend each region, so that an overlap in the extended region means a distance < maxgap in the original region
			$chromStart = $r->{"chromStart"} - $maxGap / 2;
			$chromEnd = $r->{"chromEnd"} + $maxGap / 2;
			#extension will not affect the collapse mode
		}

		my $idx = $r->{"idx"};
		
		print "$i: Previous clusterEnd = $currClusterEnd, curr read chromStart = $chromStart, chromEnd = $chromEnd, idx=$idx\n" if $debug;
		
		my $openNewCluster = 0;
		
		if ($collapse == 0)
		{
			my $overlapLen = Common::min ($currClusterEnd, $chromEnd) - Common::max($currClusterStart, $chromStart) + 1;

			if ($minOverlap < 0)
			{
				#does not require overlap, but a max gap
				#print "chromStart = $chromStart, curr cluster end = $currClusterEnd\n";
				$openNewCluster = ($chromStart - $currClusterEnd -1 > 0);
			}
			else
			{
				#requirement on overlap 
				$openNewCluster = ($overlapLen < $minOverlap) || ($overlapLen / ($chromEnd - $chromStart + 1) < $overlapFraction);
			}
		}
		elsif ($collapse == 1)
		{
			$openNewCluster = ($chromStart > $currClusterStart || $chromEnd > $currClusterEnd);
		}
		elsif ($collapse == 2)
		{
			$openNewCluster = ($chromStart > $currClusterStart && $chromEnd > $currClusterEnd);
		}
		else
		{
			Carp::croak "invalid value of collapse mode\n";
		}
		
		#begin a new cluster
		if ($openNewCluster)
		{
			my $n = @currCluster;
			print "close the old cluster with $n regions ...\n" if $debug;
			print join ("\t", @currCluster), "\n" if $debug;

			if (@currCluster > 0)
			{
				my @currClusterCpy = @currCluster;
				push @clusters, \@currClusterCpy;
			}
			print "begin a new cluster...\n" if $debug;
			@currCluster = ();
			push @currCluster, $r->{"idx"};
			$currClusterStart = $chromStart;
			$currClusterEnd = $chromEnd;
		}
		else
		{
			my $n = @currCluster;
			print "expand the current cluster with $n regions ...\n" if $debug;
			#expand the old cluster
			#
			push @currCluster, $r->{"idx"};
			#$currClusterStart = $chromStart if $currClusterStart > $chromStart;
			$currClusterEnd = $chromEnd if $currClusterEnd < $chromEnd;
		}
	}
	
	if (@currCluster >0)
	{
		my @currClusterCpy = @currCluster;
		push @clusters, \@currClusterCpy;
	}
	
	$n = @clusters;
	print "\n\n $n clusters found\n\n" if $debug;
	return \@clusters;

}

=head2 collapseReads

collapseReads by the 5'end

my $clusters = collapseReads ($readsOnChrom, $strand);

Reads are supposed to be on the same chromosome.
only reads on $strand are considered

=cut

sub collapseReads
{
	my ($readsOnChrom, $strand) = @_;
	
	my $key = "chromStart";
	$key = "chromEnd" if $strand eq '-';

	my %readsHash;

   	foreach my $r (@$readsOnChrom)
	{
		next unless $r->{"strand"} eq $strand;
		push @{$readsHash{$r->{$key}}}, $r;
	}

	my @clusters;
	foreach my $pos (sort {$a <=> $b} keys %readsHash)
	{
		my $readsInCluster = $readsHash{$pos};
		my @readsInClusterSorted = sort { $a->{"score"} <=> $b->{"score"} } @$readsInCluster; #sort according the number of mismatches
		#sort according to length
		@readsInClusterSorted = sort { $b->{"chromEnd"} - $b->{"chromStart"} <=> $a->{"chromEnd"} - $a->{"chromStart"} } @readsInClusterSorted;
		push @clusters, \@readsInClusterSorted;
	}
	return \@clusters;
}

1;

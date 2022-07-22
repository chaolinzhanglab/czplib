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

our $VERSION = 1.02;

@ISA = qw (Exporter);

@EXPORT = qw (
	bedToFull
	bedToLine
	clusterRegions
	collapseReads
	combineRegions
	contigToGenome
	copyBedRegion
	determineExonReadingFrame
	geneToSize
	geneToExon
	genomeToContig
	genomeToTranscript
	getUniqPaths
	getGenomicBreakdown
	isUpstreamRegion
	lineToBed
	overlapRegions
	printBedRegion
	readBedFile
	readNextBedBlock
	readNextBedLine
	removeSmallGap
	segmentRegion
	sortBedFile
	splitBedFileByChrom
	transcriptToGenome
	writeBedFile
);
	
=head1 Bed

Subroutines to handle BED file

=cut

use strict;
use warnings;

use Data::Dumper;
use Carp;

use sort 'stable'; #added by Chaolin Zhang on Feb 14, 2018

use Common;


#
#object-oriented interface will be added later
#

my $debug = 0;


=head2 readBedFile

Read a Bed File. Only one track is allowed in the file now and header will be ignored

my $regions = readBedFile ($inFile, $verbose, $msgio);

inFile [string]: input Bed file, gzip compressed files with .gz extension is allowed
               : use '-' for stdin
verbose        : verbose mode
msgio          : file handle for message output

return         : reference to an array
=cut

sub readBedFile
{
	my ($inBedFile, $verbose, $msgio) = @_;

	my $fin;
	my @ret;

	$msgio = *STDOUT unless $msgio;

	if ($inBedFile eq '-')
	{
		$fin = *STDIN;
	}
	else
	{
		if ($inBedFile =~/\.gz$/)
		{
			open ($fin, "gunzip -c $inBedFile | ")||Carp::croak "cannot open file $inBedFile to read\n";
		}
		elsif ($inBedFile =~/\.bz2$/)
        {
            open ($fin, "bunzip2 -c $inBedFile | ")||Carp::croak "cannot open file $inBedFile to read\n";
        }
		else
		{
			open ($fin, "<$inBedFile")||Carp::croak "cannot open file $inBedFile to read\n";
		}
	}

	my $i = 0;
	while (my $bedLine = readNextBedLine ($fin))
	{
		print $msgio "$i ...\n" if $verbose && $i % 100000 == 0;
		$i++;
		push @ret, $bedLine;
	}
	close ($fin) if $inBedFile ne '-';
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
		Carp::croak "incorrect number of blocks at line $line\n" if $#bs != $entry->{"blockCount"} - 1;
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
	
	if ($out eq '-')
	{
		$fout = *STDOUT;
	}
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

	close ($fout) if $out ne '-';
}


=head2 isUpstreamRegion

if $r1 is upstream(5') of $r2;

return:
1 - yes,
-1 - no,
0 - overlapping or cannot compare (e.g. not on the same chromosome or strand)

=cut


sub isUpstreamRegion
{
	my ($r1, $r2) = @_;
	
	Carp::croak "no strand information\n" unless exists $r1->{'strand'} && exists $r2->{'strand'};
	Carp::croak "illegal strand information\n" unless $r1->{'strand'} eq '+' || $r1->{'strand'} eq '-';

	return 0 if $r1->{'chrom'} ne $r2->{'chrom'} || $r1->{'strand'} ne $r2->{'strand'};
	
	my $ret = 0;

	if ($r1->{'chromEnd'} < $r2->{'chromStart'})
	{
		$ret = 1;
	}
	elsif ($r1->{'chromStart'} > $r2->{'chromEnd'})
	{
		$ret = -1;
	}

	$ret *= -1 if $r1->{'strand'} eq '-';
	
	return $ret;
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
		my ($inFile, $outFile, $separateStrand, $tmpDir) = @_;
		my $ncols = `head -n 1 $inFile | awk '{print NF}'`; chomp $ncols;
		Carp::croak "less than 6 columns\n" if $separateStrand && $ncols < 6;

		my $cmd = "grep -v \"^track\"  $inFile | grep -v \"^#\" | sort";
		$cmd .= " -T $tmpDir" if $tmpDir && (-d $tmpDir);  #this could be very larger than /tmp, so we might need to specify a separate address
		$cmd .= " -s -k 6,6" if $separateStrand;
		$cmd .= " -k 1,1 -k 2,2n -k 3,3n > $outFile";

		#print $cmd, "\n" if $verbose;
		my $ret = system ($cmd);
		Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;
}

=head2 splitBedFileByChrom

Usage:

my $info = splitBedFileByChrom ($inBedFile, $outDir, %params)

 inBedFile [string]: gzip compressed file with .gz input is allowed
                   : use '-' for stdin	 
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
	my $msgio = *STDOUT;

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
		$msgio = $params {'msgio'} if exists $params {'msgio'};
	}

	if ($inBedFile ne '-')
	{
		Carp::croak "$inBedFile does not exist\n" unless -f $inBedFile;
	}

	Carp::croak "$outDir does not exist\n" unless -d $outDir;

	my %tagCount;

	my $fin;
	#print "reading tags from $inBedFile ...\n" if $verbose;
	
	if ($inBedFile eq '-')
	{
		$fin = *STDIN;
	}
	else
	{
		if ($inBedFile =~/\.gz$/)
		{
			open ($fin, "gunzip -c $inBedFile |") || Carp::croak "can not open file $inBedFile to read\n";
		}
		elsif ($inBedFile =~/\.bz2$/)
		{
			open ($fin, "bunzip2 -c $inBedFile |") || Carp::croak "can not open file $inBedFile to read\n";
		}
		else
		{
			open ($fin, "<$inBedFile") || Carp::croak "can not open file $inBedFile to read\n";
		}
	}

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

		print $msgio "$i ...\n" if $i % 100000 == 0 && $verbose;
	
		$line =~/(\S+)\s/;

		my $chrom = $1;
		#print "chrom = $chrom\n";
		my $tmpFile = "$outDir/$chrom.bed";
	
		if (not exists $fhHash{$chrom})
		{
			my $fout;
			open ($fout, ">>$tmpFile") || Carp::croak "can not open file $tmpFile to write: $?\n";
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
	close ($fin) if $inBedFile ne '-';

	#close all file handles
	foreach my $chrom (sort keys %fhHash)
	{
		close ($fhHash{$chrom});

		if ($sort)
		{
				print $msgio "sorting $chrom ...\n" if $verbose;
				my $f = $tagCount{$chrom}->{'f'};
				my $f2 = "$f.sort";

				my $separateStrand = $sort > 1 ? 1 : 0;
				sortBedFile ($f, $f2, $separateStrand, $outDir);
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

sub geneToSize
{
	my $ts = $_[0];
	my $tsSize = $ts->{'chromEnd'} - $ts->{'chromStart'} + 1;

    if (exists $ts->{'blockSizes'})
    {
        $tsSize = sum ($ts->{'blockSizes'});
    }
	return $tsSize;
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

	my @exons;
	if (not exists $g->{'blockCount'})
	{
		my $e = {
            chrom => $g->{"chrom"},
            chromStart=> $g->{'chromStart'},
            chromEnd => $g->{'chromEnd'},
            name => join (":", $g->{"name"}, 0),
            score => $g->{"score"},
            strand => $g->{"strand"}
        };
        push @exons, $e;
    }
	else
	{
	    my $nexon = $g->{"blockCount"};
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
    }
    return \@exons;
}

=head2	hasPTC

determine if the stop codon (specified by thickStart, or thickEnd is PTC)
according to the 50 nt rule

return: 1 if yes, 0 if no

Note:This subroutine will not check the consistency of the reading frame


05/07/2011

=cut


sub hasPTC
{
	my $g = $_[0];
	Carp::croak "No strand defined: ", Dumper ($g), "\n" unless exists $g->{'strand'};
	Carp::croak "No blockStarts:", Dumper ($g), "\n" unless exists $g->{'blockStarts'};
	
	return 0 if $g->{'blockCount'} < 2;

	my $stopCodon = $g->{'strand'} eq '+' ? $g->{'thickEnd'} : $g->{'thickStart'};
	my $stopCodonOnTs = genomeToTranscript ($g, $stopCodon);

	my $lastExonJunctionStart = $g->{'strand'} eq '+' ?
	        $g->{'chromStart'} + $g->{'blockStarts'}->[$g->{'blockCount'}-2] + $g->{'blockSizes'}->[$g->{'blockCount'}-2] - 1 :
			$g->{'chromStart'} + $g->{'blockStarts'}->[1];
	
	my $lastExonJunctionStartOnTs = genomeToTranscript ($g, $lastExonJunctionStart);

	return $lastExonJunctionStartOnTs - $stopCodonOnTs >= 50 ? 1 : 0;
}


=head2 determineExonReadingFrame

determine the reading frame of each exon in a transcript
each exon is assigned two numbers: chopStart and chopEnd, whose value can be 0, 1, 2, or -1

chopStart: the number of nucleotide to be removed to start a complete codon
chopEnd: the number of nucleotide after the last complete codon

its easy to see that prevExon.chopStart + currExon.chopEnd = 0 or 3

return 0 if successful, return -1 if the transcript is problematic

We assume the transcript is a proper protein coding transcript, without PTC

=cut

sub determineExonReadingFrame
{
	my $g = $_[0];
	Carp::croak "No strand defined: ", Dumper ($g), "\n" unless exists $g->{'strand'};
	$g = bedToFull ($g);

	#return -1 unless $g->{'name'} eq 'NM_010828';

	#return -1 if hasPTC ($g);

	my $chrom = $g->{"chrom"};
	my $chromStart = $g->{"chromStart"};
	my $chromEnd = $g->{"chromEnd"};
	my $name = $g->{"name"};
	
	my $strand = $g->{"strand"};
	my $blockSizes = $g->{"blockSizes"};
	my $blockStarts = $g->{"blockStarts"};
	my $blockNum = @$blockSizes;
	
	my $thickStart = $g->{"thickStart"};
	my $thickEnd = $g->{"thickEnd"};


	for (my $i = 0; $i < $blockNum; $i++)
	{
		$g->{"chopStart"}->[$i] = -1;
		$g->{"chopEnd"}->[$i] = -1;
	}

	if ($thickStart == $chromStart || $thickEnd == $chromEnd)
	{
		#Carp::croak "incomplete ORF: ", Dumper ($g), "\n";
		#return -1;
	}
	
	#intronless gene
	if ($blockNum == 1)
	{
		my $s = $g->{"blockSizes"}->[0] - ($thickStart - $chromStart) - ($chromEnd - $thickEnd);
		if ($s % 3 != 0)
		{
			Carp::croak "the length of ORF is not the multiple of three: ", Dumper ($g), "\n", bedToLine ($g), "\n";
			return -1;
		}
		
		$g->{"chopStart"}->[0] = $thickStart - $chromStart;
		$g->{"chopEnd"}->[0] = $chromEnd - $thickEnd;
		return 0;
	}


	#print Dumper ($g), "\n";
	#determine the length of cds
	my $cdsLen = 0;

	for (my $i = 0; $i < $blockNum; $i++)
	{
		$cdsLen += $blockSizes->[$i];

		if($chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 < $thickStart)
		{
			#5'/3' UTR exon
			$cdsLen -= $blockSizes->[$i];
			#next;
		}
		elsif ($chromStart + $blockStarts->[$i] <= $thickStart 
				&& $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 >= $thickStart)
		{
			#the (partial) coding exon
			my $utr = $thickStart -($chromStart + $blockStarts->[$i]);
			$cdsLen -= $utr;
		}
		
		if ($chromStart + $blockStarts->[$i] > $thickEnd)
		{
			#5'/3' exon
			$cdsLen -= $blockSizes->[$i];
		}
		elsif ($chromStart + $blockStarts->[$i] <= $thickEnd
			&& $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 >= $thickEnd)
		{
			my $utr = $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 - $thickEnd;
			$cdsLen -= $utr;
		}
	}

	#print "cds Len = $cdsLen\n";
	if ($cdsLen % 3 != 0)
	{
		Carp::croak "the length of ORF ($cdsLen) is not the multiple of three: ", Dumper ($g), "\n", bedToLine ($g), "\n";
		return -1;
	}

	#determine the complete codons
	my $currCDSLen = 0;
	
	for (my $i = 0; $i < $blockNum; $i++)
	{
		#5'/3'utr exon on the left of $thickStart
		if ($chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 < $thickStart)
		{
			$g->{"chopStart"}->[$i] = -1;
			$g->{"chopEnd"}->[$i] = -1;
		}
		elsif (    $chromStart + $blockStarts->[$i] <= $thickStart 
				&& $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 >= $thickEnd)
		
		{
			#both thickStart and thickEnd is in the exon
			$g->{"chopStart"}->[$i] = $thickStart - ($chromStart + $blockStarts->[$i]);
			$currCDSLen = $thickEnd - $thickStart + 1;
			$g->{"chopEnd"}->[$i] = $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 - $thickEnd;
		}
		elsif ($chromStart + $blockStarts->[$i] <= $thickStart 
				&& $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 >= $thickStart
				&& $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 < $thickEnd)
		{
			#the exon overlap with thickStart, but not thickEnd

			$g->{"chopStart"}->[$i] = $thickStart - ($chromStart + $blockStarts->[$i]);
			$currCDSLen = $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - $thickStart;
			$g->{"chopEnd"}->[$i] = $currCDSLen - int(($currCDSLen+0.5)/3) * 3;
			
			my $s = $g->{"blockSizes"}->[$i] - $g->{"chopStart"}->[$i] - $g->{"chopEnd"}->[$i];

			my $exonStart = $chromStart + $blockStarts->[$i]; # + $g->{"chopStart"}->[$i];
			my $exonEnd = $chromStart + $blockStarts->[$i] + $blockSizes->[$i]; # - $g->{"chopEnd"}->[$i];
		
			Carp::croak "$chrom:$exonStart-$exonEnd: $name, block = $i, CDS len = $s, can not be divided by 3\n" 
			unless $s % 3 == 0;
		}
		elsif ($chromStart + $blockStarts->[$i] >= $thickStart && 
				$chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 <= $thickEnd)
		{
			#complete coding exon
			if ($g->{"chopEnd"}->[$i-1] >= 0)
			{
				$g->{"chopStart"}->[$i] = (3 - $g->{"chopEnd"}->[$i-1]) % 3;
			}
			else
			{
				$g->{"chopStart"}->[$i] = 0;
			}
			
			$currCDSLen += $blockSizes->[$i];
			$g->{"chopEnd"}->[$i] = $currCDSLen - int(($currCDSLen+0.5)/3) * 3;
			my $s = $g->{"blockSizes"}->[$i] - $g->{"chopStart"}->[$i] - $g->{"chopEnd"}->[$i];
	
			my $exonStart = $chromStart + $blockStarts->[$i]; # + $g->{"chopStart"}->[$i];
			my $exonEnd = $chromStart + $blockStarts->[$i] + $blockSizes->[$i]; # - $g->{"chopEnd"}->[$i];
			Carp::croak "$chrom:$exonStart-$exonEnd: $name, block = $i, CDS len = $s, can not be divided by 3\n" 
			unless $s % 3 == 0;
	
		}
		elsif ($chromStart + $blockStarts->[$i] >= $thickStart &&
				$chromStart + $blockStarts->[$i] <= $thickEnd &&
				$chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 > $thickEnd)
		{
			#the exon overlaps with thickEnd, but not thickStart

			if ($g->{"chopEnd"}->[$i-1] >= 0)
			{
				$g->{"chopStart"}->[$i] = (3 - $g->{"chopEnd"}->[$i-1]) % 3;
			}
			else
			{
				$g->{"chopStart"}->[$i] = 0;
			}
			
			$g->{"chopEnd"}->[$i] = $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 - $thickEnd;

			my $s = $g->{"blockSizes"}->[$i] - $g->{"chopStart"}->[$i] - $g->{"chopEnd"}->[$i];
			
			my $exonStart = $chromStart + $blockStarts->[$i]; # + $g->{"chopStart"}->[$i];
			my $exonEnd = $chromStart + $blockStarts->[$i] + $blockSizes->[$i]; # - $g->{"chopEnd"}->[$i];
		
			Carp::croak "$chrom:$exonStart-$exonEnd: $name, block = $i, CDS len = $s, can not be divided by 3\n" 
			unless $s % 3 == 0;
		}
		elsif ($chromStart + $blockStarts->[$i] > $thickEnd)
		{
			#3'/5' UTR exon
			$g->{"chopStart"}->[$i] = -1;
			$g->{"chopEnd"}->[$i] = -1;
		}
		else
		{
			Carp::croak "$name: block=$i, something wrong not handled properly...\n";
		}
	}
	return 0;
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

=head2 genomeToTranscript

#map a genomic position onto the transcript coordinate system
#return -1 if no overlap or pos is in intron

=cut
sub genomeToTranscript
{
	my ($ts, $pos) = @_;

	return -1 unless $pos >= $ts->{'chromStart'} && $pos <= $ts->{'chromEnd'};

	bedToFull ($ts) unless exists $ts->{'blockStarts'};
	my $chromStart = $ts->{'chromStart'};

	my $posOnTs = 0;
	for (my $i = 0; $i < $ts->{'blockCount'}; $i++)
	{
		my $blockStart = $chromStart + $ts->{'blockStarts'}->[$i];
		my $blockEnd = $blockStart + $ts->{'blockSizes'}->[$i] - 1;

		if ($pos > $blockEnd)
		{
			$posOnTs += $ts->{'blockSizes'}->[$i];
		}
		else
		{
			return -1 unless $pos >= $blockStart; #in intron
			$posOnTs += ($pos - $blockStart);
			last;
		}
	}

	if ($ts->{'strand'} eq '-')
	{
		my $tsLen = sum ($ts->{'blockSizes'});
		$posOnTs = $tsLen - 1 - $posOnTs;
	}
	return $posOnTs;
}

=head2 transcriptToGenome
#map a transcript position onto the genomic coordinate system
#return -1 if no overlap

=cut

sub transcriptToGenome
{
	my ($ts, $pos) = @_;
	
	#print "pos=$pos, ts=", Dumper ($ts), "\n";
	bedToFull ($ts) unless exists $ts->{'blockStarts'};

	my $blockSizes = $ts->{'blockSizes'};
	my $tsLen = sum ($blockSizes);
	return -1 if $pos < 0 || $pos >= $tsLen;

	if ($ts->{'strand'} eq '-')
	{
		$pos = $tsLen - 1 - $pos;
	}

	my $chromStart = $ts->{'chromStart'};

	my $tsLen2 = 0;

	for (my $i = 0; $i < $ts->{'blockCount'}; $i++)
	{
		if ($pos < $tsLen2 + $ts->{'blockSizes'}->[$i])
		{
			return $ts->{'chromStart'} + $ts->{'blockStarts'}->[$i] + $pos - $tsLen2;
		}
		$tsLen2 += $ts->{'blockSizes'}->[$i];
	}
	return -1;
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
	my ($regions, $checkStrand) = @_;
	my $nregions = @$regions;
	
	my @blocks;

	my $region1 = $regions->[0];

	#my @regions = ($region1, $region2);

	my $chrom = $region1->{'chrom'};
	my $strand = exists $region1->{'strand'} ? $region1->{'strand'} : "";
	foreach my $region (@$regions)
	{
		#add integrity check, make sure the regions to be combined are on the same chromosome, and on the same strand (if requested)
		#Chaolin Zhang, Dec 20, 2013
		my $chrom2 = $region->{'chrom'};
		return {} unless $chrom eq $chrom2;

		if ($checkStrand)
		{
			Carp::croak "no strand information in ", Dumper ($region), "\n" unless exists $region->{'strand'};

			my $strand2 = $region->{'strand'};
			return {} unless $strand eq $strand2;
		}

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


=head2 overlapRegions

overlapRegions ($regions1UnSorted, $regions2UnSorted)

assume the regions are all on the proper strand of the same chromosomes

For each region $r1 in $regions1, find all regions in $regions2 that overlap with $r1
return: add an array at $r->{'overlap'} and each entry is the reference to $r2

=cut
sub overlapRegions
{
	my ($regions1UnSorted, $regions2UnSorted) = @_;

	my @regions2 = sort {$a->{"chromStart"} <=> $b->{"chromStart"}} @$regions2UnSorted;
	my $regions2 = \@regions2;

	my $chromEndExt = 0;
	foreach my $r2 (@$regions2)
	{
		$chromEndExt = $r2->{"chromEnd"} if $r2->{"chromEnd"} > $chromEndExt;
		$r2->{"chromEndExt"} = $chromEndExt;
	}
	
	my @regions1 = sort {$a->{"chromStart"} <=> $b->{"chromStart"}} @$regions1UnSorted;
	my $regions1 = \@regions1;
	
	my $firstRegion2Idx = 0; #the first $r2 that overlaps with or on the right of the current window

	foreach my $r1 (@$regions1)
	{
		
		my $chromStart = $r1->{"chromStart"};
		my $chromEnd = $r1->{"chromEnd"};
		
		while ($firstRegion2Idx < @$regions2 && $regions2->[$firstRegion2Idx]->{"chromEndExt"} < $chromStart)
		{
			$firstRegion2Idx++;
		}
		
		my $i = $firstRegion2Idx;
		while ($i < @$regions2 && $regions2->[$i]->{"chromStart"} <= $chromEnd)
		{
			push @{$r1->{'overlap'}}, $regions2->[$i] if $regions2->[$i]->{"chromEnd"} >= $chromStart;
			$i++;
		}
	}
}


=head2 getGenomicBreakdown

get genomic breakdown annotation for $contigs

getGenomicBreakdown ($breakdown, $contigs, $verbose);

$beakdown: nonoverlapping bed intervals with genomic breakdown (generated by gene2breakdown.pl
$contigs: bed intervals to be annotated
$verbose:

return: 

add 'breakdown' to each entry $c in $contigs
each entry in $c->{'breakdown'} is an interval of CDS (score=1), 5'UTR (score=5) and 3'UTR (score=3)
they are contig-based coordinates and sorted in ascending order

April 5, 2011

=cut
sub getGenomicBreakdown
{
	my ($breakdown, $contigs, $verbose) = @_;
	$verbose = 0 unless $verbose;

	my %breakdownHash;
	map {push @{$breakdownHash{$_->{'strand'}}->{$_->{'chrom'}}}, $_} @$breakdown;

	my %contigHash;
	map {push @{$contigHash{$_->{'strand'}}->{$_->{'chrom'}}}, $_} @$contigs;

	for my $strand (qw(+ -))
	{
		next unless exists $contigHash{$strand};

		print "processing strand $strand ...\n" if $verbose;
		foreach my $chrom (sort keys %{$contigHash{$strand}})
		{
			next unless exists $breakdownHash{$strand} && exists $breakdownHash{$strand}->{$chrom};

			print "processing chrom $chrom ...\n" if $verbose;
			overlapRegions ($contigHash{$strand}->{$chrom}, $breakdownHash{$strand}->{$chrom});

			foreach my $c (@{$contigHash{$strand}->{$chrom}})
			{
				next unless exists $c->{'overlap'};

				my $overlap = $c->{'overlap'};
				foreach my $r (@$overlap)
				{
					my $chromStart = max ($c->{'chromStart'}, $r->{'chromStart'});
					my $chromEnd = min ($c->{'chromEnd'}, $r->{'chromEnd'});
					Carp::croak "no overlap between", Dumper ($c), " and ", Dumper ($r), "\n" if $chromStart > $chromEnd;
					
					my $r2c = genomeToContig ($c, $chromStart, $chromEnd);
					$r2c->{'name'} = $r->{'name'};
					$r2c->{'score'} = $r->{'score'};
					$r2c->{'strand'} = $r->{'strand'} eq $c->{'strand'} ? '+' : '-';
					push @{$c->{'breakdown'}}, $r2c;
					
					my @r2c = sort {$a->{'chromStart'} <=> $b->{'chromStart'}} @{$c->{'breakdown'}};
					$c->{'breakdown'} = \@r2c;
				}
				delete $c->{'overlap'};
			}
		}
	}
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

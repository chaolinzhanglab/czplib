#
#===============================================================================
#
#         FILE:  BedFile.pm
#
#  DESCRIPTION:  Package to handle Bed file
#        FILES:  BedFile.pm
#         BUGS:  ---
#        NOTES:  The package is spun off from the AnnotationIO.pm
#       AUTHOR:  Chaolin Zhang (cz), czhang@rockefeller.edu
#      COMPANY:  Rockefeller University
#      VERSION:  1.0
#      CREATED:  12/16/10 20:58:53
#     REVISION:  ---
#===============================================================================


package BedFile;

=head1 BedFile

Subroutines to handle BED file

=cut

use strict;
use warnings;

use FileHandle;
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
	while (my $bedLine = readNextBedFline ($fin))
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

				return line2bed ($line);
		}
		return "";
}


=head2 line2bed

parse a line to a BED region
any column is fine, but it has to be in correct format

=cut

sub line2bed
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

writeBedFile ($regions, $header, $out, $append);
append  [string] :  'a' means append

=cut

sub writeBedFile
{
	
	my ($regions, $header, $out, $append) = @_;
	my $fout;
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

	print $fout $header, "\n" if ($header =~/^track/);
	
	map {print printBedRegionToString ($_), "\n";} @$regions;

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
	my $str = printBedRegionToString ($region);
	print $str, "\n";
}


=head2 bed2line

The same as printBedRegion
bed2line ($region);

=cut

sub bed2line
{
	my $region = $_[0];
	return printBedRegionToString ($region);
}


=head2 printBedRegionToString

generate a line in BED format
my $line = printBedRegionToString ($region);

=cut

sub printBedRegionToString
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

1;

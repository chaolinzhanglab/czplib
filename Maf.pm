package Maf;

require Exporter;

@ISA = qw (Exporter);

@EXPORT = qw (
	readMafFile
	readNextMafBlock
	writeMafFile
	writeMafBlock
	buildGapMap
	genomeToMafCoord
	getSpeciesWithSite
	siteToTotalBranchLength
	speciesToTotalBranchLength
);


use strict;
use warnings;

use Carp;


use Sequence;
use PhyloTree;


my $debug = 0;
sub readMafFile
{
	my ($inMafFile, %options) = @_;

	my $verbose = exists $options {'verbose'} ? $options {'verbose'} : 0;
	my $species = exists $options {'species'} ? $options {'species'} : 0;

	my $fin;
	open ($fin, "<$inMafFile") || Carp::croak "can not open file $inMafFile to read\n";
	my $iter = 0;

	my @blocks;
	while (my $mafBlock = readNextMafBlock ($fin, $species))
	{
		print "$iter ...\n" if $iter % 5000 == 0 && $verbose;
		$iter++;
		push @blocks, $mafBlock;
	}
	close ($fin);
	return \@blocks;
}


sub readNextMafBlock
{
	my ($fin, $species) = @_;

	my %speciesHash;
	%speciesHash = map {$_=>1} @$species if $species;
	
	my $blockFound = 0;
	my %block;
	while (my $line = <$fin>)
	{
		chomp $line;

		if ($line =~/^\#BLOCK_NAME=(.*?)$/)
		{
			#the name is nonstandard
			$block{'name'} = $1;
		}
		elsif ($line =~/^a/)
		{
			$blockFound = 1;
			$line=~s/^a\s+//;
			$block {'annotation'} = $line;
		}
		elsif ($line =~/^s/)
		{
			my ($tmp, $src, $start, $size, $strand, $srcSize, $text) = split (/\s+/, $line);
			if ($species)
			{
				my $db = $src;
				$db = $1 if $db=~/^(\w+)\./;
				if (not exists $speciesHash{$db})
				{
					print "$db does not exist in the given list of species\n";
					next;
				}
				#next unless exists $speciesHash{$db};
			}
			push @{$block{'seq'}}, {src=>$src, start=>$start, size=>$size, strand=>$strand, srcSize=>$srcSize, text=>$text};
		}
		elsif ($line =~/^\s*$/)
		{
			if ($blockFound)
			{
				#end of the block
				last;
			}
			else
			{
				next;
			}
		}
	}

	if ($blockFound)
	{	
		return \%block;
	}
	else
	{	return 0;
	}
}

sub buildGapMap
{
	my $block = $_[0];

	foreach my $seq (@{$block->{'seq'}})
	{
		my $seqStr = $seq->{'text'};
		my @gapMap;
		my $gapLen = 0;

		if ($seqStr !~/^\-/)
		{
			push @gapMap, [0, 0];
		}

		while ($seqStr =~/(\-+)/g)
		{
			$gapLen += length ($1);
			my $gapEnd = pos($seqStr);
			my $gaplessEnd = $gapEnd - $gapLen;

			push @gapMap, [$gaplessEnd, $gapEnd];
		}
		$seq->{'gapMap'} = \@gapMap;
	}
}

#map genomic coordinates to the maf coordinates
#one site per call, but make sure 
#1. the sites are sorted before calling this function
#2. update $mapStartIdx after each call
#
sub genomeToMafCoord 
{
	my ($chromStart, $chromEnd, $map, $mapStartIdx) = @_;

	for (; $mapStartIdx < @$map; $mapStartIdx++)
	{
		my $m = $map->[$mapStartIdx];
		last if $m->[0] > $chromStart; #from
	}

	my $j = $mapStartIdx;
	if ($mapStartIdx < @$map)
	{
		print "j =$j, chromStart = $chromStart, last map entry=", $map->[$j-1]->[0], ", current map entry=", $map->[$j]->[0], "\n" if $debug;
	}
	else
	{
		print "j = $j, chromStart = $chromStart, last map entry=", $map->[$j-1]->[0], ", last block\n" if $debug;
	}

	my $mStart = $map->[$j-1];
	Carp::croak "chromStart = $chromStart, blockStart = ", $mStart->[0], ", inconsistent\n" unless $chromStart >= $mStart->[0];

	my $chromStartNew = $mStart->[1] + $chromStart - $mStart->[0];

	for (; $j < @$map; $j++)
	{
		my $m = $map->[$j];
		last if $m->[0] > $chromEnd;
	}
    
	if ($j < @$map)
	{
		print "j = $j, chromEnd = $chromEnd, last map entry=", $map->[$j-1]->[0], ", current map entry=", $map->[$j]->[0], "\n" if $debug;
	}
	else
	{
		print "j = $j, chromEnd = $chromEnd, last map entry=", $map->[$j-1]->[0], ", last block\n" if $debug;
	}
	
	my $mEnd = $map->[$j-1];
	Carp::croak "chromStart = $chromEnd, blockStart = ", $mEnd->[0], ", inconsistent\n" unless $chromEnd >= $mEnd->[0];

	my $chromEndNew = $mEnd->[1] + $chromEnd - $mEnd->[0];
    $mapStartIdx--;

	return {chromStart=>$chromStartNew, chromEnd=>$chromEndNew, mapStartIdx=>$mapStartIdx};
}


sub getSpeciesWithSite
{
	my ($start, $sequences, $word) = @_;
	
	my $nspec = @$sequences;
	
	my $wordSize = length($word);
	my @speciesToKeep;

	for (my $j = 0; $j < $nspec; $j++)
	{
		my $seq = $sequences->[$j];
		my $w = substr ($seq->{"text"}, $start, $wordSize);
		
		if ($w eq $word) 
		{
			my $db = $seq->{'src'};
			$db = $1 if $db=~/^(\w+)\./;
			push @speciesToKeep, $db;
		}
	}
	return \@speciesToKeep;
}

sub speciesToTotalBranchLength
{
	my ($speciesToKeep, $tree) = @_;
	return 0 if @$speciesToKeep < 2;

	my $treeCpy = copyTree ($tree);
	my $siteTree = subTree ($treeCpy, $speciesToKeep);
		
	my $tbl = totalBranchLength ($siteTree);
	releaseTree ($siteTree);
	return $tbl;
}


#
sub siteToTotalBranchLength
{
	my ($start, $end, $sequences, $tree) = @_;

	my $nspec = @$sequences;
	return 0 if $nspec < 2;
	
	my $wordSize = $end - $start + 1;
	my $ref = $sequences->[0];
	my $word = substr ($ref->{'text'}, $start, $wordSize);
	
	my $speciesToKeep = getSpeciesWithSite ($start, $sequences, $word);
	return speciesToTotalBranchLength ($speciesToKeep, $tree);
}



sub writeMafBlock
{
	my ($fout, $block) = @_;
	print $fout "#BLOCK_NAME=", $block->{'name'}, "\n" if exists  $block->{'name'};
	print $fout join ("\t", "a", $block->{'annotation'}), "\n";

	foreach my $seq (@{$block->{'seq'}})
	{
		print $fout join ("\t", 's', $seq->{'src'}, $seq->{'start'}, $seq->{'size'}, $seq->{'strand'}, $seq->{'srcSize'}, $seq->{'text'}), "\n";
	}
	print $fout "\n";
}




sub writeMafFile
{
	my ($outFile, $blocks, $verbose) = @_;
	
	$verbose = 0 unless $verbose;

	my $fout;
	open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";
	my $iter = 0;
	foreach my $block (@$blocks)
	{
		print "$iter ...\n" if $iter % 5000 == 0 && $verbose;
		writeMafBlock ($fout, $block);
	}
	close ($fout);
}

1;


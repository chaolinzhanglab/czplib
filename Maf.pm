package Maf;

require Exporter;


our $VERSION = 1.01;


@ISA = qw (Exporter);

@EXPORT = qw (
	buildGapMap
	genomeToMafCoord
	getSpeciesInMafBlock
	getSpeciesWithSite
	indexBigMafFile
	mafBlockToClustalW
	maskMafBlock
	readBigMafFile
	readMafFile
	readNextMafBlock
	siteToTotalBranchLength
	speciesToTotalBranchLength
	writeMafFile
	writeMafBlock
);


use strict;
use warnings;

use Carp;
use Data::Dumper;


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

=head2 readNextMafBlock

my $block = readNextMafBlock ($fin, $species);

The second parameter, an array of species to extract, is optional

return a hash of a block

=cut
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


sub indexBigMafFile
{
    my ($in, $verbose, $msgio) = @_;
    my @ret;

	$msgio = *STDOUT unless $msgio;

    my $fin;
    open ($fin, "<$in") || Carp::croak "can not open file $fin to read\n";

    #seek ($fin, 0, 0);
    my $currPos = tell ($fin);
    my $line = <$fin>;

    my $n = 0;

	my $blockName = "";
	my $iter = 0;

    while (1)
    {
        chomp $line;
        if ($line =~/^\s*$/)
        {
            $currPos = tell($fin);
            my $more = ($line =<$fin>);
            last unless $more;
            next;
        }

        if ($line =~/^\#BLOCK_NAME=(\S*?)$/)
        {
			#this is the nonstandard field
			#if we find a block name, we record it; otherwise, the $blockName will be empty
			#the pointer will point to the 'a score=XXX' line
			$blockName = $1;
            $currPos = tell($fin);
            my $more = ($line =<$fin>);
            last unless $more;
            next;
        }

        if ($line=~/^a/)
        {
            my $entry = {id=>$blockName, pointer=>$currPos};
            push @ret, $entry;

			print $msgio "$iter ...\n" if $verbose && $iter % 500 == 0;
			$iter++;
        }
        $currPos = tell($fin);

        my $more = ($line =<$fin>);
        last unless $more;
    }
    close ($fin);
    return {file=>Common::getFullPath($in), index=>\@ret};
}

sub readBigMafFile
{
	my ($inFile, $blockInfo, $species) = @_;

	#Carp::croak Dumper ($blockInfo), "\n";
	my $fin;
	open ($fin, "<$inFile") || Carp::croak "cannot open file $inFile to read\n";
	my $pointer = $blockInfo->{'pointer'};
	my $id = $blockInfo->{'id'};

	seek ($fin, $pointer, 0);
	my $block = readNextMafBlock ($fin, $species);
	$block->{'name'} = $id if $id ne '';
	close ($fin);
	return $block;
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


sub getSpeciesInMafBlock
{
	my $block = $_[0];
	my @species = map {$_->{'src'}} @{$block->{'seq'}};
    my $refSpecies = $species[0];
    if ($refSpecies=~/^(\S*?)\./)
    {
		#ref species have chromosome attached, remove that
        $species[0] = $refSpecies;
    }
	return \@species;
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

=head2 mafBlockToClustalW

my $aln = mafBlockToClustalW ($block, $species);
convert a block to clustalw alignment format

$species is an array of species to extract, and is optional

=cut

sub mafBlockToClustalW
{
	my ($block, $species) = @_;

	my %speciesHash;
	%speciesHash = map {$_=>1} @$species if $species;
	
	my $blockName = $block->{'name'};

	my @sequences;

	my @identity = split(//, $block->{'seq'}[0]->{'text'});
	foreach my $seq (@{$block->{'seq'}})
	{
		my $src = $seq->{'src'};
		
		if ($src=~/^(\S*?)\./)
		{
			$src = $1;
		}
		
		if ($species)
		{
			next unless exists $speciesHash{$src};
		}

		my $id = $blockName ne '' ? "$blockName.$src" : $src;
		my $seqStr = $seq->{'text'};
		$seqStr =~s/\./N/g;

		#dot is not a legal character in clustalw, replace with N
		#this will cause some ambiguity, but this is probably the best we can do
		push @sequences, {id=>$id, seq=>$seqStr};
		
		my @bases = split(//, $seqStr);
		for (my $i = 0; $i < @identity; $i++)
		{
			$identity[$i] = " " if $identity[$i] ne $bases[$i];
		}
	}
	
	my @identity2 = map{($_=~/[ACGTacgt]/) ? "*" : " "} @identity;
	
	return {sequences=>\@sequences, identity=>join("", @identity2)};
}

=head2 maskMafBlock

maskMafBlock ($block);


mask small letters (repetitive region) by n

=cut

sub maskMafBlock
{
	my $block = $_[0];
	my $seq = $block->{'seq'};
    foreach my $s (@$seq)
    {
    	$s->{'text'} =~tr/a-z/n/;
    }
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


# create on Nov 22, 2004
# by Chaolin Zhang
#
# Common subroutines
#!/usr/bin/perl -w

package Common;
use strict;
use Carp;
use AnnotationIO;
use Data::Dumper;
use Bio::SeqIO;
use File::Temp qw(:mktemp);
no warnings 'recursion';


=head1 NAME

Common - commonly used subroutines to manipulate numbers,
arrays, strings, etc

=head1 AUTHOR

Chaolin Zhang (zhangc@cshl.edu)
Aug, 2005

=cut

my $debug = 0;



sub chop_carriage
{
	my $str = $_[0];
	$str =~ /^(.*?)\s*?$/;
	$str = $1; 
	return $str;
}
sub clean_rep
{
	my $str = $_[0];
	#$str =~ s/\./\\\./g;
	return quotemeta($str);
}	

sub list_to_rep
{
	my ($symbols, $boundary) = @_;
	my @tmp_list;
	foreach my $s (@$symbols)
	{
		push @tmp_list, clean_rep ($s);
	}

	my $rep= "";
	if ($boundary)
	{
		$rep = join ("\\b|\\b", @tmp_list);
		$rep = "\\b". $rep. "\\b";
	}
	else
	{
		$rep = join ("|", @tmp_list);
	}
	return $rep;
}	


#get full path from relative path
use Cwd;
sub getFullPath
{
	my $path = $_[0];
	return $path if $path=~/^\//;
	my $pwd = cwd ();
	return "$pwd/$path";
}

#get taxonomic id from common names
sub getTaxId
{
	my $org = $_[0];
	Carp::croak "The name of organism $org is wrong\n" if (length($org) < 2);
	$org = substr($org, 0, 2);
	my $map = {hg=>9606,    #human
			mm=>10090,  #mouse
			rn=>10116,  #rat
			dm=>7227,   #fly
			ce=>6239,	#worm
			dr=>7955	#zebrafish
	};
    Carp::croak "The name of organism $org can not be found\n" unless exists $map->{$org};
    return exists $map->{$org} ? $map->{$org} : 0;
}

#check enum param type
sub checkParamType
{
	my ($type, $enum) = @_;
	foreach my $t (@$enum)
	{
		return $type if $type eq $t;
	}
	return '';
}


sub ls
{
	Carp::croak "need one or two arguments (dir, suffix)\n" if @_ != 1 && @_ != 2;
	my $dir = $_[0];
	my $suffix = $_[1] if @_ > 1;
	
	Carp::croak "$dir does not exists\n" unless -d $dir;

	opendir (DIR, "$dir") || Carp::croak "can not open $dir to read\n";
	my @files;

	while (my $f = readdir(DIR))
	{
		next if $f eq '.';
		next if $f eq '..';

		if ($suffix)
		{
			next unless $f=~/\.$suffix$/;
		}
		push @files, $f;
	}
	close (DIR);
	return \@files;
}

#
sub segmentStr
{
	my $str = $_[0];
	my $colNum = 60;
	$colNum = $_[1] if @_ > 1;
	

	#my ($str, $colNum) = @_;
	my $strNew = "";
	my $start = 0;

	while ($start < length ($str))
	{
		my $length = min ($colNum, length($str) - $start);
		$strNew .= substr ($str, $start, $length) . "\n";
		#$str = substr ($str, $colNum);
		$start += $colNum;
	}
	chomp $strNew; #remove the last "\n"
	#$strNew .= $str;
	return $strNew;
}

#inefficient, should not be used
sub segmentStr2
{
	my $str = $_[0];
	my $colNum = 60;
	$colNum = $_[1] if @_ > 1;
	
	#my ($str, $colNum) = @_;
	my $strNew = "";
	while (length ($str) > $colNum)
	{
		$strNew .= substr ($str, 0, $colNum) . "\n";
		$str = substr ($str, $colNum);
	}
	$strNew .= $str;
	return $strNew;
}



#//////////////////////////Array Manipulation/////////////////////////
sub bootstrapArray
{
	my $arrayRef = $_[0];
	my $len = @$arrayRef;
	my @idx;
	my $i;
	for ($i = 0; $i < $len; $i++)
	{
		push @idx, int(rand($len));
	}
	my @arrayNew = @$arrayRef;
	@arrayNew = @arrayNew[@idx];
	return \@arrayNew;
}


=head2 diffArray
diffArray - find set difference of two arrays A\B
	input: $array1_ref, $array2_ref
	output: the reference of the set diff array
	status: tested
	date: 08/27/2005
=cut


sub diffArray
{
	die "arrayDiff: incorrect number of parameters.\n" if @_!= 2;
	my ($array1, $array2) = @_;
	my @arrayNew;
	my $elem;
	foreach $elem (@$array1)
	{
		if (locateArrayElem ($array2, $elem) < 0)
		{
			push @arrayNew, $elem;
		}	
	}	
	return \@arrayNew;
}	

=head2 diffArray
diffArray - find intersect of two arrays A\B
	input: $array1_ref, $array2_ref
	output: the reference of the set diff array
	status: tested
	date: 08/27/2005
=cut
sub intersectArray
{
	die "arrayDiff: incorrect number of parameters.\n" if @_!= 2;
	my ($array1, $array2) = @_;
	my @arrayNew;
	my $elem;
	foreach $elem (@$array1)
	{
		if (locateArrayElem ($array2, $elem) >= 0)
		{
			push @arrayNew, $elem;
		}	
	}	
	return \@arrayNew;
}	


# return the index of the first matched element
sub locateArrayElem
{
	die "arrayDiff: incorrect number of parameters.\n" if @_!= 2;
	my ($array, $elem) = @_;
	my $i;
	for ($i = 0; $i <@$array; $i++)
	{
		if ($array->[$i] == $elem)
		{
			return $i;
		}	
	}	
	return -1;
}	

#randomly shuffle an array
sub shuffleArray
{
	die "shuffleArray: incorrect number of paramters.\n" if @_!= 1;
	my $array = $_[0];
	my $len = @$array;
	my $ret = [];
	return $ret if $len <= 0;

	my $randIdx = randSeq (0, $len);
	my @randIdx = @$randIdx;
	my @arrayNew = @$array;
	@arrayNew = @arrayNew[@randIdx];
	return (\@arrayNew);
}	

# generate a random sequence from $start to $start +$len -1
sub randSeq
{
	die "randSeq: incorrect number of paramters.\n" if @_!= 2;
	my ($start, $len) = @_;
	my $ret = [];
	return $ret if $len <= 0;
	my %randHash;
	my $i;
	for ($i = $start; $i < $len + $start; $i++)
	{
		$randHash{$i} = rand (1);
	}	  
	my @ret = sort{$randHash{$b} <=> $randHash{$a}} keys %randHash; 
	#print join ("\t", @ret), "\n";
	return (\@ret);
}

sub sampleSeq
{
	die "sampleSeq: incorrect number of parameters\n" if @_!= 3;
	my ($start, $lenTotal, $lenSample) = @_;
	my $ret = [];
	return $ret if $lenSample <= 0;

	my $end = $start + $lenTotal - 1;
	my @seq = ($start..$end);
	my $seqRand = randSeq ($start, $lenTotal);
	my @seqRand = @$seqRand;
	my @seqSample = @seq[@seqRand[(0..($lenSample-1))]];
	return \@seqSample;
}


#/////////////////////////////Math//////////////////////////////////
# pow
sub pow
{
	die "pow: incorrect number of parameters.\n" if @_!= 2;
	my ($base, $expr) = @_;
	my $result = exp (log ($base) * $expr);
	return $result;
}	

# norm of a vector
sub norm
{
	die "norm: incorrect number of parameters.\n" if @_ != 1;
	my $array = $_[0];
	my $sum = 0;
	my $elem;
	foreach $elem (@$array)
	{
		$sum += $elem * $elem;
	}
	$sum = sqrt ($sum);
	return $sum;
}


sub max
{
	my @array = @_;

	@array = @{$array[0]} if ref ($array[0]);

	my $m = $array[0];
	foreach my $elem (@array)
	{
		$m = ($elem > $m)? $elem : $m;
	}
	return $m;
}

sub min
{
	my @array = @_;

	@array = @{$array[0]} if ref ($array[0]);

	my $m = $array[0];
	foreach my $elem (@array)
	{
		$m = ($elem < $m)? $elem : $m;
	}
	return $m;
}


sub ABS
{
	Carp::croak __PACKAGE__ . "::ABS: incorrect number of arguments\n" if @_ != 1;
	my $in = $_[0];
	return ($in >= 0)? $in : -$in;
}

sub mod
{
	Carp::croak __PACKAGE__ . "::mod: incorrect number of arguments\n" if @_ != 2;
	my ($in, $div) = @_;

	Carp::croak "only integer is accepted in mod\n" unless $in - int ($in) == 0 && $div - int($div) == 0;
	
	return $in - int(($in+0.5)/$div) * $div;
}


#sum
sub sum
{
	my $array = $_[0];
	my $sum = 0;
	my $elem;
	foreach $elem (@$array)
	{
		$sum += $elem;
	}
	return $sum;
}

#
sub mean
{
	my $array = $_[0];
	my $len = @$array;
	die "empty array\n" if $len == 0;
	return sum ($array) / $len;
}

#combination of choosing k elements from a total of n
sub cnk
{
	die "cnk: incorrect number of paramters.\n" if @_!= 2;
	my ($n, $k) = @_;
	print "n=$n, k=$k\n";
	
	my $result = 1.0;
	my $i;
	for ($i = 0; $i < $k; $i++)
	{
		$result *= ($n - $i) / ($i + 1);
	}
	$result = int ($result + 0.5);
	return $result;
}

#enumerate k elements out of n
sub enumerateKofN
{
	die "enumerateKofN: incorrect number of parameters.\n" if @_ > 2;
	my ($n, $k) = @_;
	my $array;

	my $i;

	#initiate
	for ($i = 0; $i < $n - $k +1; $i++)
	{
		my $set;
		$set->[0] = $i;
		push @$array, $set;
	}

	#extending
	my $iter;
	for ($iter = 1; $iter < $k; $iter++)
	{
		my $len = $iter;	
		#my $len = length (@{$array->[0]});
		my @arrayNew;
		my ($set, $elem);
		foreach $set (@$array)	#each enumeration of length $len
		{
			$elem = $set->[$len-1];	#the last element of each enumeration
			my $i;
			for ($i = $elem+1; $i < $n - $k +$iter + 1; $i++)	#append one more element
			{
				my @setNew = @$set;
				push @setNew, $i;
				push @arrayNew, \@setNew;	#push extended enumeration into the container
			}
		}
		$array = \@arrayNew;
	}
	return $array;
}

#/////////////////////////Clustering////////////////////////////////
sub matrix2clusters
{
	my $matrix = $_[0];
	my $n = @$matrix;
	my @doneFlag;
	my @clusters;
	#indicate whether a vertex has been clustered
	for (my $i = 0; $i < $n; $i++)
	{
		$doneFlag[$i] = 0;
	}
	
	for (my $i = 0; $i < $n; $i++)
	{
		next if $doneFlag[$i];

		$doneFlag [$i] = 1;
		my @vertices = ($i);
		_addNeighbor ($matrix, \@doneFlag, $i, \@vertices);
		push @clusters, \@vertices;
	}
	return \@clusters;
}

sub _addNeighbor 
{
	my ($matrix, $doneFlag, $i, $vertices) =@_;
	my $n = @$matrix;
	for (my $j = 0; $j < $n; $j++)
	{
		if ($doneFlag->[$j] == 0 && $matrix->[$i]->[$j] > 0 && $i != $j)
		{
			$doneFlag->[$j] = 1;
			push @$vertices, $j;
			_addNeighbor ($matrix, $doneFlag, $j, $vertices);
		}
	}
}

#///////////////////////Sequence manipulation//////////////////////////
#reverse complementary nucleotide sequence
sub revcom
{
	my $str = $_[0];
	$str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	return CORE::reverse ($str);
}



sub nibFrag
{
	my ($nibFrag, $chromNib, $chromStart, $chromEnd, $strand, $name, $cacheDir) = @_;

	Carp::croak "$cacheDir does not exists\n" unless -d $cacheDir;

	$strand = 'm' if $strand eq '-';

	my $tmpFile = mktemp ("XXXXXX");

	my $nibFragTmp = "$cacheDir/$tmpFile.fa";
	my $cmd = "$nibFrag -masked -name=\"$name\" $chromNib $chromStart $chromEnd $strand $nibFragTmp";
	system ($cmd);
	my $seqIO = Bio::SeqIO->new (-file=>$nibFragTmp, format=>'fasta');
	my $seq=$seqIO->next_seq();
	unlink $nibFragTmp;
	
	return $seq;
}

#contig coordinates to genome coordinates
#$contig->{chrom=>chr1, chromStart=>2, chromEnd=>3, strand=>4}
#$zero based coordinates
#
sub contig2genome
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

#genome coordinates to contig coordinates
#$contig->{chrom=>chr1, chromStart=>2, chromEnd=>3, strand=>4}
#$zero based coordinates
#
sub genome2contig
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


#get unique paths(transcripts segments) between two coordinates from a set of paths(transcripts)
#$paths: a refeence to an array, each element is a bed row
#$chromStart and $chromEnd: zero-based coordinates of the query interval 
#return an array ref. to the unique paths (transcript segments)

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



#get exon/intron structure between two coordinates specifying an interval of interest

sub segmentRegion
{ 
	my ($ts, $chromStart, $chromEnd) = @_;

	Carp::croak "region beyound the start and end of transcript:", Dumper ($ts), "\n" 
	unless $ts->{"chromStart"} <= $chromStart && $ts->{"chromEnd"} >= $chromEnd;

	if (not exists $ts->{"blockCount"} || $ts->{"blockCount"} == 1)
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
	
	my %tsNew = (chrom=>$ts->{"chrom"},
				chromStart=>$chromStart,
				chromEnd=>$chromEnd,
				name=>$ts->{"name"},
				score=>$ts->{"score"},
				strand=>$ts->{"strand"});
	
	my $blockCount = $lastBlockIdx - $firstBlockIdx + 1;

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




#cluster regions on the same chromosome
#arg:
#regionsOnChrom: regions on a chrom, bed format
#strand filter : consider only tags on the give strand, could be +, - b (both)
#maxGap        : max gap allowed to consider an overlap
#overlapFraction: minimum fraction of overlap to consider a match
#collapse      : 0 (no collapse), 1 (exact match), 2 (if one read are contained by other)
#
#return        : indices in the same clusters are put in an array


sub clusterRegions
{
	my ($regionsOnChrom, $strand, $maxGap, $overlapFraction, $collapse) = @_;
	

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

		my $chromStart = $r->{"chromStart"} - $maxGap / 2;
		my $chromEnd = $r->{"chromEnd"} + $maxGap / 2;
		my $idx = $r->{"idx"};
		
		print "$i: Previous clusterEnd = $currClusterEnd, curr read chromStart = $chromStart, chromEnd = $chromEnd, idx=$idx\n" if $debug;
		
		my $openNewCluster = 0;
		
		if ($collapse == 0)
		{
			my $overlapLen = Common::min ($currClusterEnd, $chromEnd) - Common::max($currClusterStart, $chromStart) + 1;
			$openNewCluster = ($chromStart > $currClusterEnd) || ($overlapLen / ($chromEnd - $chromStart + 1) < $overlapFraction);
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

#base composition
sub baseComp
{
	my $seqs = $_[0];
	my ($a, $c, $g, $t) = (0, 0, 0, 0);
	foreach my $seqStr (@$seqs)
	{
		my $a2 = ($seqStr=~tr/aA//);
		my $c2 = ($seqStr=~tr/cC//);
		my $g2 = ($seqStr=~tr/gG//);
		my $t2 = ($seqStr=~tr/tTuU//);
		$a += $a2;
		$c += $c2;
		$g += $g2;
		$t += $t2;
	}
	my $n = $a + $c +$g + $t;
	$a /= $n;
	$c /= $n;
	$g /= $n;
	$t /= $n;
	return {A=>$a, C=>$c, G=>$g, T=>$t, N=>$n};
}

#base composition
sub baseCount
{
	my $seqs = $_[0];
	my ($a, $c, $g, $t) = (0, 0, 0, 0);
	foreach my $seqStr (@$seqs)
	{
		my $a2 = ($seqStr=~tr/aA//);
		my $c2 = ($seqStr=~tr/cC//);
		my $g2 = ($seqStr=~tr/gG//);
		my $t2 = ($seqStr=~tr/tTuU//);
		$a += $a2;
		$c += $c2;
		$g += $g2;
		$t += $t2;
	}
	my $n = $a + $c +$g + $t;
	return {A=>$a, C=>$c, G=>$g, T=>$t, N=>$n};
}



#calculate percent identity from psl line
#according to the pslCalcMilliBad function
#http://genome.ucsc.edu/FAQ/FAQblat.html

sub calcPslPercentIdentity
{
	my ($psl, $isMrna) = @_;
	my $sizeMul = 1; #3 for protein
	my $milliBad = 0;
	
	my $qAliSize = $sizeMul * ($psl->{"qEnd"} - $psl->{"qStart"} + 1);
	my $tAliSize = $psl->{"tEnd"} - $psl->{"tStart"} + 1;
	my $aliSize = Common::min($qAliSize, $tAliSize);
	return 0 if $aliSize <= 0;

	my $sizeDif = $qAliSize - $tAliSize;
	if ($sizeDif < 0)
	{
		$sizeDif = 0 - $sizeDif;
		$sizeDif = 0 if $isMrna;
	}
	my $insertFactor = $psl->{"qNumInsert"};
	$insertFactor += $psl->{"tNumInsert"} if ($isMrna == 0);

	my $total = $sizeMul * ($psl->{"matches"} + $psl->{"repMatches"} + $psl->{"misMatches"});
	
	if ($total != 0)
	{
		$milliBad = $psl->{"misMatches"} * $sizeMul + $insertFactor + int(3*log(1+$sizeDif)+0.5);
		$milliBad /= $total;
	}
	return (1- $milliBad);
}

sub blat
{
	Carp::croak "three or four arguments expected\n" unless (@_ >= 4 && @_ <= 5);
	
	my ($blat, $db, $query, $cache) = @_;
	
	#my $cache = "/tmp";
	my $options = $_[4] if (@_ == 5); # && (-d $_[3]));
	my $optPair = {
		t=>1,
		q=>1,
		occ=>1,
		tileSize=>1,
		stepSize=>1,
		oneOff=>1,
		minMatch=>1,
		minScore=>1,
		minIdentity=>1,
		maxGap=>1,
		makeOcc=>1,
		repMatch=>1,
		mask=>1,
		qMask=>1,
		repeats=>1,
		minRepDivergence=>1,
		dots=>1,
		out=>1,
		maxIntron=>1
	};
	my $optLeft = {
		prot=>1,
		noHead=>1,
		trimT=>1,
		noTrimA=>1,
		trimHardA=>1,
		fastMap=>1,
		fine=>1,
		extendThroughN=>1
	};
	
	my $optStr = "";
	foreach my $opt (keys %$options)
	{
		if (exists $optPair->{$opt})
		{
			$optStr .= " -$opt=" . $options->{$opt};
		}
		elsif (exists $optLeft->{$opt})
		{
			$optStr .= " -$opt";
		}
		else
		{
			Carp::croak "bad options for blat: $opt\n";
		}
	}
	
	my $out = "$cache/blatout.".time ().rand();
	Carp::croak "the output $out already exists\n" if -f $out;

	my $cmd = "$blat $db $query $optStr $out >& /dev/null";
	
	#print $cmd, "\n";
	my $ret = system ($cmd);
	Carp::croak "Command has been crashed ($cmd): $?\n" unless $ret == 0;

	my $result = AnnotationIO::readPslFile ($out);
	unlink $out if -f $out;
	return $result;
}

#need to add code to handle alignment text
sub checkSim4AlignQuality
{
	my ($aligns, $gene, $coverage, $identity, $verbose) = @_;
	my @goodAligns;
	
	foreach my $a (@$aligns)
	{
		#print "align = ", Dumper ($a), "\n";

		my @seq1blocks = @{$a->{"seq1blocks"}};#genomic
		my @seq2blocks = @{$a->{"seq2blocks"}};#transcript
		my @blockscore = @{$a->{"blockscore"}};
		
		my $nblocks = @{$a->{"seq2blocks"}};
		if ($nblocks < 2)
		{
			print "gene=$gene vs transcript=", $a->{"seq2id"}, ": no splice junction spaned, failed\n" if $verbose;
			next;
		}
		
		#OK now we have at least two blocks, so there should be no error in the following
		#code to check the first and last blocks
	
		#remove the first block if less than 25 nt
		my $firstBlockSize = $seq1blocks[0]->{"end"} - $seq1blocks[0]->{"start"} + 1;
		
		#############################################################
		#TO BE FIXIED
		#use the while loop in the next version
		#
		#############################################################
		#
		#while ($nblocks >= 2 && $firstBlockSize < 25)
		if ($firstBlockSize < 25)
		{
			$a->{"match"} -= $firstBlockSize;
			$a->{"identity"} -= $firstBlockSize * $blockscore[0] / 100;
			shift @seq1blocks;
			shift @seq2blocks;
			shift @blockscore;
			$nblocks -= 1;
			#$firstBlockSize = $seq1blocks[0]->{"end"} - $seq1blocks[0]->{"start"} + 1;
		}
		
		#remove the last block if less than 25 nt 
		my $lastBlockSize = $seq1blocks[$nblocks-1]->{"end"} - $seq1blocks[$nblocks-1]->{"start"} + 1;
		

		#############################################################
		#TO BE FIXIED
		#use the while loop in the next version
		#
		#############################################################
		#
		
		#while ($nblocks >= 2 && $lastBlockSize < 25)
		if ($lastBlockSize < 25)
		{
			$a->{"match"} -= $lastBlockSize;
			$a->{"identity"} -= $lastBlockSize * $blockscore[$nblocks - 1] / 100;
			pop @seq1blocks;
			pop @seq2blocks;
			pop @blockscore;
			$nblocks -= 1;
			#$lastBlockSize = $seq1blocks[$nblocks-1]->{"end"} - $seq1blocks[$nblocks-1]->{"start"} + 1;
		}
		
		if ($nblocks < @{$a->{"seq2blocks"}})
		{
			$a->{"seq1blocks"} = \@seq1blocks;
			$a->{"seq2blocks"} = \@seq2blocks;
			$a->{"blockscore"} = \@blockscore;
		}
		
		#make sure we still have at least two blocks
		$nblocks = @{$a->{"seq2blocks"}};
		if ($nblocks < 2)
		{
			print "gene=$gene vs transcript=", $a->{"seq2id"}, ": no splice junction spaned, failed\n" if $verbose;
			next;
		}
		
		#check coverage and identity
		my $cov = $a->{"match"} / $a->{'seq2len'}; 	#percent coverage
		my $id = $a->{"identity"} / $a->{"match"};	#percent identity
		
		if ($cov < $coverage || $id < $identity)
		{
			print "gene=$gene vs transcript=", $a->{"seq2id"}, ": coverage=$cov, identity=$id, failed\n" if $verbose;
			next;
		}

		my $contineous = 1;
		my $tsblocks = $a->{"seq2blocks"};
		for (my $i = 0; $i < $nblocks - 1; $i++)
		{
			if ($tsblocks->[$i]->{"end"} + 1 != $tsblocks->[$i+1]->{"start"})
			{
				$contineous = 0;
				last;
			}
		}
		
		if ($contineous != 1)
		{
			print "gene=$gene vs transcript=", $a->{"seq2id"}, ": gaps in transcript, failed\n" if $verbose;
			next;
		}

		push @goodAligns, $a;
	}

	my $n = @goodAligns;
	my $ntotal = @$aligns;
	
	print "$n of $ntotal aligned transcripts passed quality control\n" if $verbose;
	
	return \@goodAligns;
}



#search a nucleotide word in a sequence
#return the start position of each hit
#
#support IUB code


sub allowIUB
{
	my $word = $_[0];
	$word=~s/R/[A|G]/ig;
	$word=~s/Y/[C|T]/ig;
	$word=~s/M/[A|C]/ig;
	$word=~s/K/[G|T]/ig;
	$word=~s/S/[C|G]/ig;
	$word=~s/W/[A|T]/ig;
	$word=~s/H/[A|C|T]/ig;
	$word=~s/B/[C|G|T]/ig;
	$word=~s/V/[A|C|G]/ig;
	$word=~s/D/[A|G|T]/ig;
	$word=~s/N/[A|C|G|T]/ig;
	return $word;
}

sub searchWord
{
	my ($str, $word) = @_;
	my @hits;

	$word=allowIUB ($word);
	
	while ($str =~/($word)/ig)
	{
		my $start = pos ($str) - length ($1);
		my $end = $start + length ($1) - 1;
		push @hits, [$start, $end];
		pos ($str) = $start + 1;
	}
	return \@hits;	
}


sub countMismatch
{
	my ($w1, $w2, $ignoreCase) = @_;

	$w1 =~tr/a-z/A-Z/ if $ignoreCase;
	$w2 =~tr/a-z/A-Z/ if $ignoreCase;

	Carp::croak "unequal length of $w1 and $w2\n" unless length ($w1) == length ($w2);

	my @w1 = split (//, $w1);
	my @w2 = split (//, $w2);

	my $n = 0;
	for(my $i = 0; $i < @w1; $i++)
	{
		$n++ if $w1[$i] ne $w2[$i];
	}
	return $n;
}

#Convert count matrix to stormo matrix
#see AnnotationIO for the data structure of matrix
#
#status: tested
#Date: 09/29/2006
#
sub ToStormoMatrix
{
	my ($matrix, $baseComp) = @_;
	my @bases = ('A', 'C', 'G', 'T');
	
	my $width = @$matrix;
	my $nseq = $matrix->[0]->{'A'} + $matrix->[0]->{'C'} + $matrix->[0]->{'G'} + $matrix->[0]->{'T'};
	$nseq = int ($nseq + 0.5);
	my $tmp = $baseComp->{'A'} + $baseComp->{'C'} + $baseComp->{'G'} + $baseComp->{'T'};
	Carp::croak "incorrect base composition? sum=$tmp\n" if (ABS ($tmp - 1) > 1e-2);
	Carp::croak "number of sequences is $nseq, unlikely a count matrix\n" if ($nseq < 2);
	
	for (my $p = 0; $p < $width; $p++)
	{
		#every position may have a different number of sequences
		my $nseq = $matrix->[$p]->{'A'} + $matrix->[$p]->{'C'} + $matrix->[$p]->{'G'} + $matrix->[$p]->{'T'};
		foreach my $b (@bases)
		{
			Croak::carp ("negative matrix elements (p=$p, base=$b, value=" . $matrix->[$p]->{$b} . ")\n") if $matrix->[$p]->{$b} < 0;
			$matrix->[$p]->{$b} = log (($matrix->[$p]->{$b} + $baseComp->{$b})/($nseq+1)/$baseComp->{$b})/log(2);
		}
	}
	return $matrix;
}

#status: tested
#Date: 09/29/2006
#
sub getMaxMatrixScore
{
	my $matrix = $_[0];
	my $width = @$matrix;
	my $score = 0;
	for (my $p = 0; $p < $width; $p++)
	{
		$score += max ($matrix->[$p]->{'A'}, $matrix->[$p]->{'C'}, $matrix->[$p]->{'G'}, $matrix->[$p]->{'T'});
	}
	return $score;
}

#status: tested
#Date: 09/29/2006
#
sub getMinMatrixScore
{
	my $matrix = $_[0];
	my $width = @$matrix;
	my $score = 0;
	for (my $p = 0; $p < $width; $p++)
	{
		$score += min ($matrix->[$p]->{'A'}, $matrix->[$p]->{'C'}, $matrix->[$p]->{'G'}, $matrix->[$p]->{'T'});
	}
	return $score;
}

#status: tested
#Date: 09/29/2006
sub getMatrixScore
{
	my ($matrix, $site) = @_;
	my $width = @$matrix;
	Carp::croak "incorrect site length (site=$site, matrix width=$width)\n" if length($site) != $width;
	
	$site =~tr/acgut/ACGUT/;
	my $score = 0;
	for (my $p = 0; $p < $width; $p++)
	{
		my $c = substr ($site, $p, 1);
		if (exists $matrix->[$p]->{$c})
		{
			$score += $matrix->[$p]->{$c};
		}
		else
		{
			$score += min ($matrix->[$p]->{'A'}, $matrix->[$p]->{'C'}, $matrix->[$p]->{'G'}, $matrix->[$p]->{'T'});
		}
	}
	return $score;
}

#Date: 09/29/2006
sub ToFrequencyMatrix
{
	my $matrix = $_[0];
	my @bases = ('A', 'C', 'G', 'T');
	my $width = @$matrix;
	my $nseq = $matrix->[0]->{'A'} + $matrix->[0]->{'C'} + $matrix->[0]->{'G'} + $matrix->[0]->{'T'};
	$nseq = int ($nseq + 0.5);

	for (my $p = 0; $p < $width; $p++)
	{
		my $nseq = $matrix->[$p]->{'A'} + $matrix->[$p]->{'C'} + $matrix->[$p]->{'G'} + $matrix->[$p]->{'T'};
		foreach my $b (@bases)
		{
			Croak::carp ("negative matrix elements (p=$p, base=$b, value=" . $matrix->[$p]->{$b} . ")\n") if $matrix->[$p]->{$b} < 0;
			$matrix->[$p]->{$b} /= $nseq; 
		}
	}
}

#Date: 09/29/2006
sub ToBayesianMatrix
{
	my ($matrix, $baseComp) = @_;
	my $prior = 0.5;
	$prior = $_[2] if (@_ > 2);

	my @bases = ('A', 'C', 'G', 'T');
	my $width = @$matrix;
	my $nseq = $matrix->[0]->{'A'} + $matrix->[0]->{'C'} + $matrix->[0]->{'G'} + $matrix->[0]->{'T'};
	$nseq = int ($nseq + 0.5);

	for (my $p = 0; $p < $width; $p++)
	{
		my $nseq = $matrix->[$p]->{'A'} + $matrix->[$p]->{'C'} + $matrix->[$p]->{'G'} + $matrix->[$p]->{'T'};
		foreach my $b (@bases)
		{
			Croak::carp ("negative matrix elements (p=$p, base=$b, value=" . $matrix->[$p]->{$b} . ")\n") if $matrix->[$p]->{$b} < 0;
			$matrix->[$p]->{$b} /= $nseq; 
			$matrix->[$p]->{$b} = log (($matrix->[$p]->{$b} + $prior * $baseComp->{$b}) / ($baseComp->{$b} * (1 + $prior))) / log (2);
		}
	}
}

sub waitUntilQsubDone
{
	my ($user, $jobName, $nSplit, $verbose) = @_;
	my $secondSlept = 0;
	while (1)
	{
		my @qstat = `qstat -u $user`;
		#remove title lines
		if (@qstat)
		{
			shift @qstat;
			shift @qstat;
		}
		my $jobNotFinished = 0;
	
		foreach my $line (@qstat)
		{
			chomp $line;
			my @cols = split (/\s+/, $line);
			shift @cols if $cols[0] eq '';
			my $jname = $cols[2];
			#print $jname, "\n";

			$jobNotFinished++ if ($jname=~/^$jobName/);
		}

		if ($jobNotFinished > 0)
		{
			$secondSlept += 10;
			if ($verbose)
			{
				my $date = `date`;
				chomp $date;
				print "$date: $jobNotFinished of $nSplit jobs are still running...\n" if ($secondSlept - int($secondSlept / 60) * 60 == 0)
			}
			sleep (10);	#10 second
		}
		else
		{
			print "done\n" if $verbose;
			last;
		}
	}
}




#phylogenetic tree manipulation
#

sub nodeInfo
{
	my $node = $_[0];
	return $node->{"iter"} . "($node)";
}

sub copyTree
{
	my $rootFrom = $_[0];
	my $rootTo = {};

	$rootTo = $_[1] if (@_ > 1);

	$rootTo->{"iter"} = $rootFrom->{"iter"};
	$rootTo->{"id"} = $rootFrom->{"id"};
	$rootTo->{"blen"} = $rootFrom->{"blen"};

	if (exists $rootFrom->{"left"})
	{
		
		my $leftChildTo = {};
		my $rightChildTo = {};
		$rootTo->{"left"} = $leftChildTo;
		$leftChildTo->{"parent"} = $rootTo;
		copyTree ($rootFrom->{"left"}, $rootTo->{"left"});
				
		$rootTo->{"right"} = $rightChildTo;
		$rightChildTo->{"parent"} = $rootTo;
		copyTree ($rootFrom->{"right"}, $rootTo->{"right"});
	}
	
	return $rootTo;
}

sub codeTree
{
	my $root = $_[0];
	my $str = "";
	$str = $_[1] if (@_ > 1);

	if ($root == 0)
	{
		return ""; #empty
	}
	if (exists $root->{"left"})
	{
		$str .= "(";
		$str = codeTree ($root->{"left"}, $str);
		$str .= ",";
		$str = codeTree ($root->{"right"}, $str);
		$str .= ")";
		$str .= ":" . $root->{"blen"} if exists $root->{"parent"};
	}
	else #leaf
	{
		$str .= $root->{"id"};
	    $str .=	":" . $root->{"blen"} if exists $root->{"parent"};
	}
	return $str;
}


sub getNodes
{
	my $root = $_[0];
	my $nodes = [];
	$nodes = $_[1] if @_ > 1;

	push @$nodes, $root;
	if ($root->{"left"})
	{
		getNodes ($root->{"left"}, $nodes);
		getNodes ($root->{"right"}, $nodes);
	}
	return $nodes;
}

#keep only the minimal subtree that covers all the given species 
sub subTree
{
	my ($tree, $species) = @_;
	#my $treeCpy = copyTree ($tree);
	my $leaves = getLeafNodes ($tree);

	my %speciesToKeep;

	foreach my $s (@$species)
	{
		Carp::croak "species $s does not exist in the tree\n" unless exists $leaves->{$s};
		$speciesToKeep{$s} = 1;
	}
	
	my @leavesToRm;
	foreach my $s (keys %$leaves)
	{
		push @leavesToRm, $leaves->{$s} unless exists $speciesToKeep{$s};
	}
	return removeLeafNodes ($tree, \@leavesToRm);
}

sub removeSpecies
{
	my ($tree, $species) = @_;
	#my $treeCpy = copyTree ($tree);
	my $leaves = getLeafNodes ($tree);

	my @leavesToRm;
	foreach my $s (@$species)
	{
		Carp::croak "species $s does not exist in the three\n" unless exists $leaves->{$s};
		push @leavesToRm, $leaves->{$s}; # if exists $leaves->{$s};
	}
	return removeLeafNodes ($tree, \@leavesToRm);
}

#remove nodes
sub removeLeafNodes
{
	my ($root, $nodes) = @_;
	
	my $i = 0;
	foreach my $n (@$nodes)
	{
		Carp::croak "The node is not a leaf node: ", Dumper ($n), "\n" if exists $n->{"left"};
		
		print $i++,",  node to remove: ", $n->{"id"}, "\n" if $debug;
		
		if (not exists $n->{"parent"})
		{#single node tree
			return 0;
		}
		elsif (not exists $n->{"parent"}->{"parent"})
		{
			#parent is the root
			my $branch = ($n->{"parent"}->{"left"}->{"iter"} == $n->{"iter"})? 'right' : 'left';
			
			$root = $root->{$branch};
			$root->{"blen"} = 0;
			delete $root->{"parent"};
		}
		elsif ($n->{"parent"}->{"left"}->{"iter"} == $n->{"iter"})
		{
			#delete left leaf
			my $rightSibling = $n->{"parent"}->{"right"};
			$rightSibling->{"blen"} += $n->{"parent"}->{"blen"};
			my $grandParent = $n->{"parent"}->{"parent"};
			
			my $branch = ($grandParent->{"left"}->{"iter"} == $n->{"parent"}->{"iter"})? 'left' : 'right';
			$grandParent->{$branch} = $rightSibling;
			$rightSibling->{"parent"} = $grandParent;
			
		}
		else
		{
			#delete right leaf
			print "delete right leaf\n" if $debug;
			print "left sibling: ", $n->{"parent"}->{"left"}->{"id"}, "\n" if $debug;
			my $leftSibling = $n->{"parent"}->{"left"};
			$leftSibling->{"blen"} += $n->{"parent"}->{"blen"};
			my $grandParent = $n->{"parent"}->{"parent"};
			my $branch = ($grandParent->{"left"}->{"iter"} == $n->{"parent"}->{"iter"})? 'left' : 'right';
			$grandParent->{$branch} = $leftSibling;
			$leftSibling->{"parent"} = $grandParent;
		}

		releaseNode ($n->{"parent"}) if exists $n->{"parent"};
		releaseNode ($n);

		print codeTree ($root), "\n\n" if $debug;
	}
	return $root;
}

sub releaseNode
{
	my $node = $_[0];
	foreach my $k (keys %$node)
	{
		delete $node->{$k};
	}
	return;
}

sub releaseTree
{
	my $root = $_[0];
	if (exists $root->{"left"})
	{
		releaseTree ($root->{"left"});
		releaseTree ($root->{"right"});
	}
	
	releaseNode ($root);
}

sub totalBranchLength
{
	my $root = $_[0];
	print "current node: ", nodeInfo ($root), "id=", $root->{"id"}, ", blen=", $root->{"blen"}, "\n" if $debug;
	my $tbn = 0;
	#$root->{"blen"} if exists $root->{"parent"};
	
	if (exists $root->{"left"})
	{
		$tbn += $root->{"left"}->{"blen"} + totalBranchLength ($root->{"left"});
		$tbn += $root->{"right"}->{"blen"} + totalBranchLength ($root->{"right"});
	}
	
	#$tbn = $root->{"blen"} if exists $root->{"parent"};
	return $tbn;
}

sub getLeafNodes
{
	my $root = $_[0];
	my $nodes = getNodes ($root); 
	#print Dumper ($nodes), "\n";
	my %leaves;
	foreach my $n (@$nodes)
	{
		if (not exists $n->{"left"})
		{
			Carp::croak "no id for leaf node ", Dumper ($n), "\n" unless exists $n->{"id"} && $n->{"id"} ne '';
			$leaves{$n->{"id"}} = $n;
		}
	}
	return \%leaves;
}




1;


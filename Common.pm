# create on Nov 22, 2004
# by Chaolin Zhang
#
# Common subroutines
#!/usr/bin/perl -w

package Common;
use strict;
use Carp;
use AnnotationIO;

=head1 NAME

Common - commonly used subroutines to manipulate numbers,
arrays, strings, etc

=head1 AUTHOR

Chaolin Zhang (zhangc@cshl.edu)
Aug, 2005

=cut


use Cwd;
sub getFullPath
{
	my $path = $_[0];
	return $path if $path=~/^\//;
	my $pwd = cwd ();
	return "$pwd/$path";
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


my $debug = 0;
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

sub blat
{
	Carp::croak "three or four arguments expected\n" unless (@_ >= 3 && @_ <= 4);
	
	my ($blat, $db, $query) = @_;
	
	my $cache = "/tmp";
	$cache = $_[3] if (@_ == 4 && (-d $_[3]));
	
	my $out = "$cache/blatout.".time ().rand();
	Carp::croak "the output $out already exists\n" if -f $out;
	
	system ("$blat $db $query $out >& /dev/null");
	my $result = AnnotationIO::readPslFile ($out);
	unlink $out if -f $out;
	return $result;
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

1;



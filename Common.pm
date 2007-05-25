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

sub sampleSeq
{
	die "sampleSeq: incorrect number of parameters\n" if @_!= 3;
	my ($start, $lenTotal, $lenSample) = @_;
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
			my $jname = $cols[2];

			$jobNotFinished++ if ($jname eq $jobName);
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



1;



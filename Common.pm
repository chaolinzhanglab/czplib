#
#===============================================================================
#
#         FILE:  Common.pm
#
#  DESCRIPTION:  common utility subroutines
#         BUGS:  ---
#        NOTES:  create on Nov 22, 2004
#       AUTHOR:  Chaolin Zhang (cz), czhang@rockefeller.edu
#      COMPANY:  Rockefeller University
#      VERSION:  1.0
#      CREATED:  12/17/10
#     REVISION:  ---
#===============================================================================



package Common;

require Exporter;


our $VERSION = 1.02;


@ISA = qw (Exporter);

@EXPORT = qw (
	ABS
	binomTest
	bootstrapArray
	bSearchSupLessThan
	checkParamType
	chisq
	chop_carriage
	clean_rep
	cnk
	cnkln
	count2prob
	dbinom
	diffArray
	entropy
	enumerateKofN
	getFullPath
	getTaxId
	gscore
	hypergeoTest
	hypergeo
	intersectArray
	list_to_rep
	locateArrayElem
	ls
	matrix2clusters
	max
	mean
	min
	norm
	pow
	stdev
	randSeq
	readConfig
	sampleSeq
	shuffleArray
	splitFileByRow
	sum
);

=head1 NAME

Common - commonly used subroutines to manipulate numbers,
arrays, strings, etc

=head1 AUTHOR

Chaolin Zhang (czhang@rockefeller.edu)
Created on Nov 22, 2004
Last revision on Dec 17, 2010

=cut


use strict;
use warnings;

use Carp;
use Data::Dumper;

no warnings 'recursion';





my $debug = 0;


#read simple configuration files
#each line is name/value pairs
#lines that begin with # are comments
#name = value

sub readConfig
{
    my $in = $_[0];
    my %conf;
    my $fin;
    open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
    while (my $line=<$fin>)
    {
        chomp $line;
        next if $line=~/^\s*$/;
        next if $line=~/^\#/;
        my ($name, $value) = split (/\s*\=\s*/, $line);
        $conf{$name} = $value;
    }
    close ($fin);
    return \%conf;
}


sub splitFileByRow
{
	my ($inFile, $outFileStem, $n, $verbose, $append) = @_;

	my $fin;
	my $fout;
	open ($fin, "<$inFile") || Carp::croak "cannot open $inFile to read\n";
	
	my $i = 0;
	my $iSplit = -1;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;

		if ($i % $n == 0)
		{
			if ($iSplit >=0)
			{
				close ($fout);
			}

			$iSplit++;

			print "starting a new split $iSplit ...\n" if $verbose;

			my $outFile = "$outFileStem.$iSplit";
			print "outFile = $outFile\n" if $verbose;
			
			my $appendFlag = $append ? ">>" : ">";
			open ($fout, $appendFlag . $outFile) || Carp::croak "cannot open $outFile to write\n";
		}
		$i++;
		print $fout $line, "\n";
	}

	close ($fin);
	close ($fout) if $iSplit > -1;
}




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





#//////////////////////////Array Manipulation/////////////////////////


=head2 bootstrapArray

sample with replicates

my $newArray = bootstrapArray (\@array)

=cut

sub bootstrapArray
{
	my $arrayRef = $_[0];
	my $len = @$arrayRef;
	my $i;
	
	#my @idx;
	#for ($i = 0; $i < $len; $i++)
	#{
	#	push @idx, int(rand($len));
	#}
		
	#my @arrayNew = @$arrayRef;
	#@arrayNew = @arrayNew[@idx];
	
	my @idx = map {int(rand($len))} @$arrayRef;
	my @arrayNew = @{$arrayRef}[@idx];
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

	my %newArray2 = map {$_=>1} @$array2;
	my @ret;
	foreach my $elem (@$array1)
	{
		push @ret, $elem if exists $newArray2{$elem};
	}	
	return \@ret;
}	


=head2 locateArrayElem
return the index of the first matched element

Note: inefficient implementation, need to be fixed in the future

=cut

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

=head bSearchSupLessThan
binary search to find the max i in [$lowerBound, $upperBound],
so that $psi->[$i] < $elem

my $i = bSearchSupLessThan ($psi, $elem, $lowerBound = 0, $upperBound = length(@$psi)-1);


=cut

sub bSearchSupLessThan
{
	my ($psi, $elem, $lowerBound, $upperBound) = @_;

	$lowerBound = 0 unless defined $lowerBound;
	$upperBound = @$psi - 1 unless defined $upperBound;

	#print "entering bsearch\n";
	#print "lb=$lowerBound, ub=$upperBound\n";

	return -1 if $lowerBound > $upperBound || $psi->[$lowerBound] >= $elem;
	return $upperBound if $psi->[$upperBound] < $elem;

	#i is somewhere in between
	my ($l, $u) = ($lowerBound, $upperBound);  # lower, upper end of search interval
	
	my $i;                       # index of probe
	while ($l <= $u) 
	{
		#print "l= $l, u= $u\n";
		$i = int(($l + $u)/2);
		if ($psi->[$i+1] < $elem)
		{
			#go up
			$l = $i+1;
		}
		elsif ($psi->[$i] >= $elem) 
		{
			#go down
			$u = $i-1;
		} 
		elsif ($psi->[$i] < $elem && $psi->[$i+1] >= $elem)
		{
			return $i; # found
		}
	}
	return -1;     # not found, it will actually never reaches here
}



=head2 shuffleArray
#randomly shuffle an array


=cut
sub shuffleArray
{
	die "shuffleArray: incorrect number of paramters.\n" if @_!= 1;
	my $array = $_[0];
	my $len = @$array;
	my $ret = [];
	return $ret if $len <= 0;

	my $randIdx = randSeq (0, $len);
	my @arrayNew = @{$array}[@$randIdx];
	return (\@arrayNew);
}	


=head2 randSeq

generate an array with $len numbers from $start in random order

my $array = readSeq (0, 10);


=cut

sub randSeq
{
	die "randSeq: incorrect number of paramters.\n" if @_!= 2;
	my ($start, $len) = @_;
	my $ret = [];
	return $ret if $len <= 0;
	
	my %randHash = map {$_=> rand (1)} ($start .. ($start + $len -1));
	
	#my $i;
	#for ($i = $start; $i < $len + $start; $i++)
	#{
	#	$randHash{$i} = rand (1);
	#}	  
	
	my @ret = sort{$randHash{$b} <=> $randHash{$a}} keys %randHash; 
	#print join ("\t", @ret), "\n";
	return (\@ret);
}


=head2 sampleSeq

a subsample of a random sequence from $start to $start+$len -1

=cut

sub sampleSeq
{
	die "sampleSeq: incorrect number of parameters\n" if @_!= 3;
	
	my ($start, $lenTotal, $lenSample) = @_;
	my $ret = [];
	return $ret if $lenSample <= 0;

	my $seqRand = randSeq ($start, $lenTotal);
	
	my @seqSample = @{$seqRand}[(0..($lenSample-1))];
	return \@seqSample;
}


=obsolete
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
=cut




#/////////////////////////////Math//////////////////////////////////



# stdev of a vector
sub stdev
{
	die "stdev: incorrect number of parameters.\n" if @_ != 1;
	my $array = $_[0];
	
	Carp::croak "array length is 1\n" if @$array == 1;

	my $sum = 0;
	
	my $m = mean ($array);	

	map {$sum+= ($_-$m) **2} @$array;

	my $n = @$array;
	return sqrt ($sum/($n-1));
}




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
	
	map {$sum+= $_ **2} @$array;

	#my $elem;
	#foreach $elem (@$array)
	#{
	#	$sum += $elem * $elem;
	#}
	$sum = sqrt ($sum);
	return $sum;
}

sub count2prob
{
	die "count2prob: incorrect number of parameters.\n" if @_ != 1 && @_!= 2;
    my ($array, $pseudoCount) = @_;
	$pseudoCount = 0 unless $pseudoCount;

    my $s = sum ($array);
	my $n = @$array;

	my @ret = map{($_+$pseudoCount)/($s+$pseudoCount * $n)} @$array;
	return \@ret;
}



sub max
{
	my @array = @_;

	@array = @{$array[0]} if ref ($array[0]);

	my $m = $array[0];
	map {$m = $_ if $m < $_} @array;

	#foreach my $elem (@array)
	#{
	#	$m = ($elem > $m)? $elem : $m;
	#}
	return $m;
}

sub min
{
	my @array = @_;

	@array = @{$array[0]} if ref ($array[0]);

	my $m = $array[0];
	map {$m = $_ if $m > $_} @array;

	#foreach my $elem (@array)
	#{
	#	$m = ($elem < $m)? $elem : $m;
	#}
	return $m;
}


sub ABS
{
	Carp::croak __PACKAGE__ . "::ABS: incorrect number of arguments\n" if @_ != 1;
	my $in = $_[0];
	return ($in >= 0)? $in : -$in;
}


=obsolete

sub mod
{
	Carp::croak __PACKAGE__ . "::mod: incorrect number of arguments\n" if @_ != 2;
	my ($in, $div) = @_;

	Carp::croak "only integer is accepted in mod\n" unless $in - int ($in) == 0 && $div - int($div) == 0;
	
	return $in - int(($in+0.5)/$div) * $div;
}
=cut

#sum
sub sum
{
	my $array = $_[0];
	my $sum = 0;
	map {$sum+= $_} @$array;

	#my $elem;
	#foreach $elem (@$array)
	#{
	#	$sum += $elem;
	#}
	return $sum;
}


=obsolete
sub sum2
{
	my $array = $_[0];
	my $sum = 0;
	map {$sum+=$_} @$array;
	return $sum;
}
=cut


#
sub mean
{
	my $array = $_[0];
	my $len = @$array;
	Carp::croak "empty array\n" if $len == 0;
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


sub entropy
{
	my $dist = $_[0];
	my $entropy = 0;

	my @d = @$dist;
	my $s = sum($dist);
	if ($s != 1)
	{
		@d = map {$_/$s} @d;	
	}

	foreach my $p (@d)
	{
		$entropy -= $p > 0 ? $p * log($p) / log(2) : 0;
	}
	return $entropy;
}


sub chisq
{
    my $dat = $_[0];

#    print Dumper ($dat), "\n";
    my $nrow = @$dat;
    my $ncol = @{$dat->[0]};

    my @a; #rowSum / total
    my @b; #colSum / total

    for (my $i = 0; $i < $nrow; $i++)
    {
        $a[$i] = sum($dat->[$i]);
    }

    for (my $j = 0; $j < $ncol; $j++)
    {
        my @x = map {$dat->[$_]->[$j]} (0 .. ($nrow -1));
        $b[$j] = sum(\@x);
    }

    my $total = sum (\@a);
    return 0 if $total <= 0;

	#print "a=", join ("\t", @a), "\n";
	#print "b=", join ("\t", @b), "\n";

    @a = map {$_/$total} @a;
    @b = map {$_/$total} @b;

	#print "a=", join ("\t", @a), "\n";
	#print "b=", join ("\t", @b), "\n";

    my $chisq = 0;

	#my @tmp;
    for (my $i = 0; $i < $nrow; $i++)
    {
        for (my $j = 0; $j < $ncol; $j++)
        {
            my $dat_exp = $a[$i] * $b[$j] * $total;
            $chisq += ($dat->[$i][$j] - $dat_exp) * ($dat->[$i][$j] - $dat_exp) / $dat_exp if $dat_exp > 0;
			#$tmp[$i][$j] = $dat_exp;
        }
    }

	#print "tmp=", Dumper (\@tmp), "\n";
    #print "chisq=$chisq\n";
    return $chisq;
}


sub gscore
{
    my $dat = $_[0];

#    print Dumper ($dat), "\n";
    my $nrow = @$dat;
    my $ncol = @{$dat->[0]};

    my @a; #rowSum / total
    my @b; #colSum / total

    for (my $i = 0; $i < $nrow; $i++)
    {
        $a[$i] = sum($dat->[$i]);
    }

    for (my $j = 0; $j < $ncol; $j++)
    {
        my @x = map {$dat->[$_]->[$j]} (0 .. ($nrow -1));
        $b[$j] = sum(\@x);
    }

    my $total = sum (\@a);
    return 0 if $total <= 0;

	#print "a=", join ("\t", @a), "\n";
	#print "b=", join ("\t", @b), "\n";

    @a = map {$_/$total} @a;
    @b = map {$_/$total} @b;

	#print "a=", join ("\t", @a), "\n";
	#print "b=", join ("\t", @b), "\n";

    my $chisq = 0;

	#my @tmp;
    for (my $i = 0; $i < $nrow; $i++)
    {
        for (my $j = 0; $j < $ncol; $j++)
        {
            my $dat_exp = $a[$i] * $b[$j] * $total;
            $chisq += 2 * $dat->[$i][$j] * log($dat->[$i][$j] / $dat_exp) if $dat_exp > 0 && $dat->[$i][$j] > 0;

			#$chisq += ($dat->[$i][$j] - $dat_exp) * ($dat->[$i][$j] - $dat_exp) / $dat_exp if $dat_exp > 0;
			#$tmp[$i][$j] = $dat_exp;
        }
    }

	#print "tmp=", Dumper (\@tmp), "\n";
    #print "chisq=$chisq\n";
    return $chisq;
}







#two tailed binom test
#a re-implementation of the R binom.test
sub binomTest
{
	my ($x, $n, $p) = @_;
	my $ret;

	if ($p == 0)
	{
		return $x == 0 ? 1 : 0;
	}
	elsif ($p == 1)
	{
		return $x == 1 ? 1 : 0;
	}

	#p > 0 &&  p < 1
	
	my $relErr = 1 + 1e-7;
	my $d = dbinom ($x, $n, $p);
	my $m = $n * $p;
	
	if ($x == $m)
	{
		return 1;
	}
	elsif ($x < $m)
	{
		my $i = int ($m); #ceiling
		$i++ if $i < $m; #ceiling $m

		while ($i <= $n)
		{
			last if dbinom ($i, $n, $p) <= $d * $relErr;
			$i++;
		}
		#y=n-i+1
		return pbinom ($x, $n, $p) + 1 - pbinom ($i-1, $n, $p);
	}
	else
	{
		my $i = int($m); #floor $m
		while ($i >= 0)
		{
			last if dbinom ($i, $n, $p) <= $d * $relErr;
			$i--;
		}
		#y=i+1
		my $ret = ($i < 0 ? 0 : pbinom ($i, $n, $p)) + 1 - pbinom ($x -1, $n, $p);

		return $ret;
	}
}


use Math::CDF qw(:all);
#the probability mass function
sub dbinom
{
	my ($x, $n, $p) = @_;
	return $x >= 1 ? pbinom ($x, $n, $p) - pbinom ($x - 1, $n, $p) : pbinom($x, $n, $p);
}


#hypergeometric probability

sub hypergeoTest
{
    my ($m, $k, $N, $n) = @_;

    my $pval = 0;
    for (my $i = $k; $i <= Common::min($n, $m); $i++)
    {
        my $p = hypergeo ($m, $i, $N, $n);
        #print "hypergeo ($m, $i, $N, $n) = $p ...\n";
        $pval += $p;
    }
    return $pval;
}

sub hypergeo
{
    my ($m, $k, $N, $n) = @_;
    my $logP = cnkln($m, $k) + cnkln ($N-$m, $n-$k) - cnkln ($N, $n);
    return exp($logP);
}

{
    my @cache;
    $cache[0] = 0;
    $cache[1] = 0;

    sub factorln
    {
        my $n = shift;
        return $cache[$n] if defined $cache[$n];
        return undef if $n < 0;
        for (my $i = scalar(@cache); $i <= $n; $i++)
        {
            $cache[$i] = $cache[$i-1] + log($i);
        }
        return $cache[$n];
    }
}

sub cnkln
{
    my ($n, $k) = @_;
    return factorln($n) - factorln($k) - factorln($n-$k);
}


#/////////////////////////Clustering////////////////////////////////
sub matrix2clusters
{
	my $matrix = $_[0];
	
	my $n = ref $matrix eq 'ARRAY' ? @$matrix : keys %$matrix;
	my @clusters;
	
	my @doneFlag;
	#indicate whether a vertex has been clustered
	for (my $i = 0; $i < $n; $i++)
	{
		$doneFlag[$i] = 0;
	}
	
	if (ref $matrix eq 'ARRAY')
	{
		for (my $i = 0; $i < $n; $i++)
		{
			next if $doneFlag[$i];
	
			$doneFlag [$i] = 1;
			my @vertices = ($i);
			_addNeighbor ($matrix, \@doneFlag, $i, \@vertices);
			push @clusters, \@vertices;
		}
	}
	else
	{
		my @idx = sort keys %$matrix;
		for (my $i = 0; $i < $n; $i++)
        {
            next if $doneFlag[$i];

            $doneFlag [$i] = 1;
			my $v = $idx[$i];
            my @vertices = ($v);
            _addNeighbor ($matrix, \@doneFlag, $v, \@vertices, \@idx);
            push @clusters, \@vertices;
        }
	}

	return \@clusters;
}

sub _addNeighbor 
{
	my ($matrix, $doneFlag, $i, $vertices, $idx) =@_;

	if (ref $matrix eq 'ARRAY')
	{
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
	else
	{
		my $n = @$idx;

		my $v = $i;
		for (my $j = 0; $j < $n; $j++)
        {
			my $v2 = $idx->[$j];
			#print "v = $v, v2=$v2\n";
            if ($doneFlag->[$j] == 0 && $v ne $v2 && exists $matrix->{$v}{$v2} && $matrix->{$v}{$v2} > 0)
            {
                $doneFlag->[$j] = 1;
                push @$vertices, $v2;
                _addNeighbor ($matrix, $doneFlag, $v2, $vertices, $idx);
            }
        }
	}
}



1;


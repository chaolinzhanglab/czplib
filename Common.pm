# create on Nov 22, 2004
# by Chaolin Zhang
#
# Common subroutines
#!/usr/bin/perl -w

package Common;
use strict;

=head1 NAME

Common - commonly used subroutines to manipulate numbers,
arrays, strings, etc

=head1 AUTHOR

Chaolin Zhang (zhangc@cshl.edu)
Aug, 2005

=cut


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
		
1;

#
#===============================================================================
#
#         FILE:  BWT.pm
#
#  DESCRIPTION:  Burrows-Wheeler transform
#
#        FILES:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Chaolin Zhang (cz), czhang@rockefeller.edu
#      COMPANY:  Rockefeller University
#      VERSION:  1.0
#      CREATED:  12/27/10 12:06:20
#     REVISION:  ---
#===============================================================================

package BWT;

use strict;
use warnings;

use Carp;
use Data::Dumper;
use POSIX qw(ceil floor);

use Common;


my $terminator = ' ';
my $debug = 1;

sub log2
{
	my $in = $_[0];
	return log($in) / log(2);
}


=head2 buildCyclicSuffixTable

my $table = buildCyclicSuffixTable ($str, $l = 0);
when $l is not specified, $table is the full cyclic table

when $l is specified and $l < length ($str), on the first $l
suffixes are built, and these suffixes are truncated to size $l

=cut


sub buildCyclicSuffixTable
{
	my ($str, $l) = @_;
	my $n = length ($str);
	
	$l = $n unless defined $l;

	print "entering buildCyclicSuffixTable ...", `date` if $debug;	
	my $v = $n;
	$v = $l if $l && $l < $n;
	my %cyclicSuffixTable = (0 => $str);

	for (my $i = 1; $i < $v; $i++)
	{
		#print "$i ...\n" if $debug;

		my $firstChar = substr ($cyclicSuffixTable{$i-1}, 0, 1);
		my $theRest = substr ($cyclicSuffixTable{$i-1}, 1, 2 * $l);

		$cyclicSuffixTable{$i} = $theRest . $firstChar;
	}

	#get the first $l characters if necessary
	if ($l && $l < $n && $n < 2 * $l)
	{
		print "truncate to $l characters ...", `date` if $debug;
		foreach my $i (sort {$a <=> $b} keys %cyclicSuffixTable)
		{
			#print "$i ...\n" if $debug;
			$cyclicSuffixTable{$i} = substr ($cyclicSuffixTable{$i}, 0, $l);
		}
	}

	return \%cyclicSuffixTable;
}

=head2 bwt_naive

bwt by cyclic shifting to enumerate all suffixes, which are then sorted

my $ret = bwt_naive ($str); #$ret is the BWT string
my %ret = bwt_naive ($str);

$ret{'W'} is the BWT string and $ret{'SA'} is the suffix array


=cut

sub bwt_naive
{
	my ($str) = @_;

	my $strLen = length ($str);
	my $lastChar = substr ($str, $strLen -1 , 1);

	$str .= $terminator if $lastChar ne $terminator;

	my $cyclicSuffixTable = buildCyclicSuffixTable ($str);
	
	my @SA = sort {$cyclicSuffixTable->{$a} cmp $cyclicSuffixTable->{$b}} keys %$cyclicSuffixTable;
	my $W = join ("", map {substr($cyclicSuffixTable->{$_}, -1, 1)} @SA);
	
	return wantarray ? (W=>$W, SA=>\@SA) : $W;
}

#get the SAinv array from SA array
sub invert
{
	my $orig = $_[0];
	my $n = @$orig;

	my @inv;
	map {$inv[$orig->[$_]] = $_} (0 ..($n-1));
	return \@inv;
}


sub SA2Psi
{
	my $SA = $_[0];
	my $n = @$SA;

	my $SAinv = invert ($SA);

	my @psi;

	map {$psi[$_] = $SA->[$_] == $n-1 ? $SAinv->[0] : $SAinv->[$SA->[$_] + 1]} (0 .. ($n-1));

	return \@psi;
}

=head2 psi2SAinv

convert psi to SAinv

my ($ret = psi2SAinv ($psi, $l=@$psi);

if $l is specified, it returns only the first $l elements

=cut
sub psi2SAinv
{
	my ($psi, $l) = @_;
	my $n = @$psi;
	if ($l)
	{
		$l = $l <= $n ? $l : $n;
	}
	else
	{
		$l = $n;
	}

	my @SAinv;
	$SAinv[0] = $psi->[0];
	
	#for (my $i = 1; $i < $n; $i++)
	#{
	#	$SAinv[$i] = $psi->[$SAinv[$i-1]];
	#}
	map {$SAinv[$_] = $psi->[$SAinv[$_-1]]} (1 .. ($l-1));

	return \@SAinv;
}

sub psi2SA
{
	my $psi = $_[0];
	my $n = @$psi;
	my $p = $psi->[0];

	my @SA;
	map {$SA[$p] = $_; $p=$psi->[$p]} (0 .. ($n-1));
	return \@SA;

	#return invert (psi2SAinv($psi));
}

sub SA2bwt
{
	my ($str, $SA) = @_;

	my $n = length ($str);

	#append terminator character if necessary
	my $lastChar = substr ($str, $n-1, 1);
	$str .= $terminator unless $lastChar eq $terminator;
	$n = length ($str);
	my $m = @$SA;

	Carp::croak "expecting $n elements in SA, observed $m elements\n" unless @$SA == $n;

	my $W = "";
	map {$W = $W . ($SA->[$_] > 0 ? substr($str, $SA->[$_]-1,1) : substr($str, $n-1, 1))} (0 ..($n-1));

	return $W;
}

sub psi2bwt
{
	my ($str, $psi) = @_;
	my $SA = psi2SA ($psi);
	
	return SA2bwt ($str, $SA);
}


=head2 sortLongSuffixByPrefix
calculate the P vector: i.e., the rank of the first l characters of the 
long suffixes of A, or SAinv (A[0..(l-1)])

To be replaced by the suffix sorting algorithm of Larsson and Sadakane
(for ceil[log2(l)] round)

=cut 
sub sortLongSuffixByPrefix
{
	my ($A, $l) = @_;

	print "entering sortLongSuffixByPrefix ...", `date` if $debug;
	Carp::croak "The length of the input string must be larger than $l\n" unless length ($A) > $l;

	#get the cyclic suffix table, with only the first $l suffixes and each suffix is truncated to the size $l
	my $cyclicSuffixTable = buildCyclicSuffixTable ($A, $l);
	my @SA = sort {$cyclicSuffixTable->{$a} cmp $cyclicSuffixTable->{$b}} keys %$cyclicSuffixTable;

	print "leaving sortLongSuffixByPrefix ...", `date` if $debug;
	return invert (\@SA);
}


=head2 sortLongSuffixInA
sort long suffixes

my %ret = sortLongSuffixInA ($A, $l, $psi_B);
#%ret = (M=>$M, Minv=>$Minv);
or 
my $M = sortLongSuffixInA ($A, $l, $psi_B);


=cut

sub sortLongSuffixInA
{
	my ($A, $l, $psi_B) = @_;
	
	print "entering sortLongSuffixInA ...", `date` if $debug;

	my $P = sortLongSuffixByPrefix ($A, $l);
	my $Q = psi2SAinv ($psi_B, $l); #SAinv

	print "got P & Q ...", `date` if $debug;
	$#{$Q} = $l - 1; 
	#if @$B < $l, the rank can be determined from the first $l characters, so 
	#we can simly add zero at the end

	my @Minv = sort {$P->[$a] <=> $P->[$b] || $Q->[$a] <=> $Q->[$b]} (0..($l-1));

	print "leaving sortLongSuffixInA ...", `date` if $debug;
	return wantarray ? (M=>invert (\@Minv), Minv=>\@Minv) : invert (\@Minv);
}


=head2 countChar

count the number of occurrances

my $sharp = countChar ($str, $alphabets = 0);

$alphabets: the reference to the alphabets. if specified, a zero is included for those characters that
do not occur in $str, and croak if $str has characters not in $alphabets

=cut
sub countChar
{
	my ($str, $alphabets) = @_;

	my $n = length ($str);

	my %ret;

	print "split the text ...", `date` if $debug;
	my @array = unpack ("C*", $str);

	print "count characters ...", `date` if $debug;
	map{$ret{chr($_)} +=1} @array;

	@array = ();
	#map {$ret{substr($str, $_, 1)} += 1} (0 .. ($n - 1));

	if ($alphabets)
	{
		foreach my $char (@$alphabets)
		{
			$ret{$char} = 0 unless exists $ret{$char};
		}
		Carp::croak "unknown characters in input text:", Dumper (%ret), "\n" if @$alphabets != keys %ret;
	}

	return \%ret;
}

=head2 countChar
count the number of occurrances of characters lexically smaller than every character in the alphabet

my $alpha = countCharLess ($sharp, $alphabets = undef);

$sharp: the occurrance of each character
$alphabets: if $alphabets is specified (as an array), it has to be sorted lexicographically

return:

$alpha{$c} is the number of characters lexicographically smaller than $c

=cut

sub countCharLess
{
	my ($sharp, $alphabets) = @_;

	#print "entering countCharLess\n";
	my %alpha;

	if ($alphabets)
	{
		foreach my $char (@$alphabets)
		{
			$sharp->{$char} = 0 unless exists $sharp->{$char};
		}
		Carp::croak "input string has characters not in alphabets: \nsharp=", Dumper ($sharp), "\n", 
		"alphabets=", Dumper ($alphabets), "\n" if keys %$sharp != @$alphabets;
	}
	else
	{
		my @alphabets = sort keys %$sharp;
		$alphabets = \@alphabets;
	}
	my $n = @$alphabets;
	
	#print "sharp=", Dumper ($sharp), "\n";
	#print "alphabets=", join ("\t", @$alphabets), "\n";

	$alpha{$alphabets->[0]} = 0;

	map {$alpha{$alphabets->[$_]} = $alpha{$alphabets->[$_-1]} + $sharp->{$alphabets->[$_-1]}} (1 ..($n-1)) ;
	
	#print "alpha=", Dumper (\%alpha), "\n";

	#print "exiting countCharLess\n";
	return \%alpha;

	#for (my $i = 1; $i < @alphabets; $i++)
	#{
	#	my $char0 = $alphabets[$i-1]
	#	my $char = $alphabets[$i];
	#	$alpha{$char} = $alpha{$char0} + $sharp{$char0};
	#}
}


=head2 sortLongSuffixInB

my $L = sortLongSuffixInB ($A, $l, $psi_B, $alphabets);

=cut
sub sortLongSuffixInB
{
	my ($A, $l, $psi_B, $sharp_B, $alpha_B) = @_;
	
	my $B = substr ($A, $l);
	
	#my $sharp = countChar ($B, $alphabets);
	#my $alpha = countCharLess ($sharp, $alphabets); 
	#$alphabets is optional now because it was specified above
	
	#print "sharp_B=", Dumper ($sharp), "\n";
	#print "alpha_B=", Dumper ($alpha), "\n";

	my @L;
	$L[$l] = $psi_B->[0];
	for (my $k = $l -1; $k >= 0; $k--)
	{
		my $char = substr ($A, $k, 1);
		my $rMax = -1;

		my @range = ($alpha_B->{$char} .. ($alpha_B->{$char} + $sharp_B->{$char} - 1));
		my @psi = @$psi_B[@range];
		
		my $lb = $alpha_B->{$char};
		my $ub = $alpha_B->{$char} + $sharp_B->{$char} - 1;

		#print "lb= $lb, ub=$ub\n";
		my $r = bSearchSupLessThan ($psi_B, $L[$k+1], $alpha_B->{$char}, $alpha_B->{$char} + $sharp_B->{$char} - 1);

		$L[$k] = $r >= 0 ? ($r+1) : $alpha_B->{$char};

		#print "c=$char, r=$r, L[$k]=", $L[$k], "L[$k+1]=", $L[$k+1], "\n";
	}
	pop @L;

	return \@L;
}


=head2 buildIndexableDictionary

my $ret = buildIndexableDictionary ($V, $alphabets);

$V: is an array of characters
$alphabets: an array of legitimate alphabets

return:

$ret->{'rank'}->{$char}->{$pos}: is the number of $char in $V preceeding $pos
$ret->{'select'}->{$char}->{$rank}: is the index of $rank's $char in $V

=cut

sub buildIndexableDictionary
{
	my ($V, $alphabets) = @_;

	my $l = @$V;
	
	my %occSoFar; #the occurrance of each character preceeding each position
				  #this support the rank operator
	my %occTotal;
	my $sum;

	map {$occSoFar{$_}->[0] = 0} @$alphabets;
	foreach my $c (@$alphabets)
	{
		map {$occSoFar{$c}->[$_] = $V->[$_-1] eq $c ? $occSoFar{$c}->[$_-1]+1 : $occSoFar{$c}->[$_-1]} (1..($l-1));
		$occTotal{$c} = $V->[$l-1] eq $c ? $occSoFar{$c}->[$l-1] + 1 : $occSoFar{$c}->[$l-1];
		$sum += $occTotal{$c};
	}
	
	Carp::croak "unknown character in input: $V\n" if ($sum != $l);

	#print "rank=", Dumper (\%occSoFar), "\n";

	my %indexFromRank; #support the select operator
	foreach my $c (@$alphabets)
	{
		my @rank = map {$V->[$_] eq $c ? $occSoFar{$c}->[$_]+1 : -1} (0..($l-1));
		#a dummy index -1 if position $i is not $c
		map {$rank[$_] >= 0 ? $indexFromRank{$c}->[$rank[$_]] = $_ : 0} (0..($l-1));
		$indexFromRank{$c}->[0] = -1;
		$#{$indexFromRank{$c}} = $occTotal{$c}; #resize to the correct size, note that the first element is not used
	}
	
	return {rank=> \%occSoFar, "select"=>\%indexFromRank};
}


=head2 bwt2

my $W = bwt2 ($str, $alphabets = undef, $verbose = 0);
or
my %ret = bwt2 ($str, $alphabets = undef, $verbose = 0);

return: 
1) If a scalar is wanted:
	the BWT string
2) if a hash is wanted:
	$ret->{'W'} is the BWT string
	$ret->{'SA'} is the suffix array
	$ret->{'psi'} is the CSA array
=cut

sub bwt2
{
	my ($str, $alphabets, $verbose) = @_;
	
	#append the terminating char if necessary
	my $strLen = length ($str);
	my $lastChar = substr ($str, $strLen -1 , 1);
	$str .= $terminator if $lastChar ne $terminator;

	my $n = length ($str);
	
	my @alphabets;
	if ($alphabets)
	{
		#if alphabets are specified, make sure they are in the correct oder and
		#the teminator is also included
		@alphabets = sort @$alphabets;
		unshift @alphabets, $terminator if $alphabets[0] ne $terminator;
	}
	else
	{
		#if not specified, the alphabets are all characters that occur in the text
		print "determining alphabets ...\n" if $verbose;
		my $sharp = countChar ($str);
		@alphabets = sort keys %$sharp;
	}

	#my $l = floor (log2($n)*10);
	my $l = floor ($n / log2($n)); #length of each segment
	#my $l = floor (sqrt($n));
	my $ns = ceil ($n / $l); #number of segment

	print "The input text has $n characters; broken into $ns segments, each with $l characters\n" if $verbose;

	#the last segment calls for the brute force BWT to initiate
	my $B = substr ($str, $l * ($ns - 1), $l);
	my %ret = bwt_naive ($B);

	my $sharp_B = countChar ($B, \@alphabets);
	my $alpha_B = countCharLess ($sharp_B, \@alphabets);

	#my $W_B = $ret{'W'};
	my $SA_B = $ret{'SA'};
	my $psi_B = SA2Psi ($SA_B);
	
	#the block-by-block incremental BWT
	for (my $i = $ns-2; $i >= 0; $i--)
	{

		print "\n" if $debug;
		print "processing block $i ...\n" if $verbose;
		my $seg = substr ($str, $l * $i, $l);
		my $A = $seg . $B;
		my $lA = length ($A);

		my $sharp_A = countChar ($seg, \@alphabets);
		my $alpha_A = countCharLess ($sharp_A, \@alphabets);

		map {
			$sharp_A->{$_} += $sharp_B->{$_};
			$alpha_A->{$_} += $alpha_B->{$_};
		} keys %$sharp_A;

		#print "\nB=$B|\n";
		#print "psi_B=", join("\t", @$psi_B), "\n";

		#print "A=$A|\n";
		
		#step1: get the rank of long suffixes A0, A1, A{l-1} among themselves
		
		print "sort long suffixes among themselves ...", `date` if $debug;

		my %ret = sortLongSuffixInA ($A, $l, $psi_B);
		
		#my %ret = sortLongSuffixInA ($A, $l, $SA_B);
		my $M = $ret{'M'};

		#my $cyclicSuffixTable = buildCyclicSuffixTable ($A);
		#foreach my $i (0..($l-1))
		#{
		#	print "A[$i] = ", join ("\t", $cyclicSuffixTable->{$i}), "\n";
		#}

		my $Minv = $ret{'Minv'};

		#step 2: get the rank of long suffixes A0, A1, A{l-1} among suffixes of B
		print "sort long suffixes among the suffixes of B ...", `date` if $debug;
		my $L = sortLongSuffixInB ($A, $l, $psi_B, $sharp_B, $alpha_B);

		
		#step 3: #calculate the V bit vector that supports fast rank and select
		#for i=0,...,(l-1), V[SAinv[i]] = 1, corresponding to long suffixes
		#for the rest, V[i] = 0, corresponding to short suffixes
		
		print "building indexable dictionary ...", `date` if $debug;
		my @ML = map {$M->[$_] + $L->[$_]} (0..($l-1));
		my @V = map {0} (0 .. ($lA-1));
		map {$V[$ML[$_]] = 1} (0..($l-1));


		#print "M=", join ("\t", @$M), "\n";
		#print "L=", join ("\t", @$L), "\n";
		#print "ML=", join ("\t", @ML), "\n";	

		my $Vdict = buildIndexableDictionary (\@V, [0,1]);


		#step 4: calculate psi_A
		print "update psi ...", `date` if $debug;
		my @psi_A;
		$psi_A[0] = $M->[0] + $L->[0];
		for (my $r = 1; $r < $lA; $r++)
		{
			if ($V[$r] == 0)
			{
				#The suffix with rank r is a short suffix
				#rank0 ($V, $r);
				my $r2 = $Vdict->{'rank'}->{'0'}->[$r];
				
				#select0 ($V, $psi_B->[$r2]+1);
				$psi_A[$r] = $Vdict->{'select'}->{'0'}->[$psi_B->[$r2]+1];
			}
			else
			{
				#The suffix with rank r is a long suffix
				#my $r2 = rank1 ($V, $r);
				my $r2 = $Vdict->{'rank'}->{'1'}->[$r];
				my $k = $Minv->[$r2];
				$psi_A[$r] = $k < $l -1 ? $M->[$k+1] + $L->[$k+1] : $Vdict->{'select'}->{'0'}->[$psi_B->[0]+1]; #select0 ($V, $psi_B->[0]+1);
			}
		}

		#print "psi_A=", join ("\t", @psi_A), "\n";

		#prepare for the next iteration
		$B = $A;
		$psi_B = \@psi_A;

		#$SA_B = psi2SA ($psi_B);
		$sharp_B = $sharp_A;
		$alpha_B = $alpha_A;
	}
	
	#my $W = SA2bwt ($str, $SA_B);
	my $W = psi2bwt ($str, $psi_B);
	return wantarray ? (W=>$W, psi=>$psi_B) : $W;
}


=head2 occ

$array: array of char

=cut

sub occ
{
	my ($array, $endPos, $char) = @_;
	
	return 0 if $endPos == 0;

	my $ret = 0;
	map {$ret++ if $array->[$_] eq $char} (0 ..($endPos -1));
	return $ret;
}

sub locateChar
{
	my ($array, $char) = @_;
	
	my $i = 0;

	#map {last if $_ eq $char; $i++} @$array;
	for (; $i < @$array; $i++)
	{
		last if $array->[$i] eq $char;
	}
	return $i;
}

=head2 getC
get the number of characters lexigraphically smaller than $char

=cut
sub getC
{
	my ($C, $char) = @_;
	
	my @chars = sort keys %$C;
	my $i = 0;

	#map {last if $_ >= $char; $i++} @chars;
	for (; $i < @chars; $i++)
	{
		last if $chars[$i] ge $char;
	}

	Carp::croak "terminator is not the lexigraphically smallest?\n" if $i == 0;

	my $char2 = $chars[$i - 1];
	return $C->{$char2};
}

=head2 updateC

update %C after insertion of $char

$C{$char} is the number of characters lexigraphically smaller than or equal to $char

=cut
sub updateC
{
	my ($C, $char) = @_;

	my @chars = sort keys %$C;
	my $i = 0;

	#map {last if $_ >= $char; $i++} @chars;
	
	for (; $i < @chars; $i++)
	{
		last if $chars[$i] ge $char;
	}


	Carp::croak "terminator is not the lexigraphically smallest?\n" if $i == 0;
	
	if (not exists $C->{$char})
	{
		my $char1 = $chars[$i-1];
		$C->{$char} = $C->{$char1} + 1;
	}
	
	for (; $i < @chars; $i++)
	{
		my $char2 = $chars[$i];
		$C->{$char2} += 1;	
	}
	return;
}

=head2 bwt

Reference: Lippert, Mobarry and Walenz, 2005, JCB, 12:943-951, p945

=cut
sub bwt
{
	my ($str, $verbose) = @_;

	#print "input = $str\n";
	Carp::croak "the input string cannot have space\n" if $str=~/\s/;

	#print "splitting input string ...\n" if $verbose;
	#my @from = split (//, $str);
	my @to = ($terminator); #start from the empty string with only the terminating char

	my %C;
	$C{$terminator} = 1;
	my $d = 0;

	print "encoding input string ...\n" if $verbose;
	my $iter = 0;
	
	for (my $i = length ($str) -1; $i >= 0; $i--)
	#foreach my $c (reverse @from)
	{
		my $c = substr ($str, $i, 1);

		print "$iter ...\n" if $verbose && $iter % 1000 == 0;
		$iter++;

		#print "\n\nchar = $c\n";
		#print "original to =|", join("", @to), "|\n";
		#my $d = locateChar (\@to, $terminator);
		#print "d=$d\n";
		Carp::croak "cannot locate the teminating char\n" if $d >= @to;
		
		#print "replace terminator at pos $d with $c\n";

		$to[$d] = $c;
		#print "tmp to =|", join('', @to), "|\n";
		#print "C=", Dumper (\%C), "\n";
		#my $a = getC (\%C, $c);
		#my $b =  occ (\@to, $d, $c);
		#print "getC=$a, occ=$b\n";

		#TODO: occ is very inefficent
		$d = getC (\%C, $c) + occ (\@to, $d, $c);
	
		splice @to, $d, 0, $terminator; #insert $terminator at position $d
		
		updateC (\%C, $c);
	}
	print "return encoded string ...\n";
	return join ('', @to);
}


1;


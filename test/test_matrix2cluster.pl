#!/usr/bin/perl -w
use Common;

my $matrix = [
	[0, 1, 0, 0, 0, 0],
	[1, 0, 0, 0, 0, 0],
	[0, 0, 0, 1, 1, 0],
	[0, 0, 1, 0, 1, 0],
	[0, 0, 1, 1, 0, 0],
	[0, 0, 0, 0, 0, 0]
];

my $clusters = Common::matrix2clusters ($matrix);

for (my $i = 0; $i < @$clusters; $i++)
{
	print join ("\t", @{$clusters->[$i]}), "\n";
}

# create on Nov 22, 2004
# by Chaolin Zhang
#
# Manipulation of microarray data
#!/usr/bin/perl -w
package Microarray::arrayOp;

use strict;
use Common;
use Carp;

#@subArrayMat@
#@Input:
#@ $array = $_[0];
#@ $arrayIdx = $_[1];
#@ $geneIdx = $_[2];
#@Output:
#@ $arrayNew
#@ status: tested
sub subArrayMat 
{
#	my $argin = @_;
	die "subArrayMat: incorrect number of parameters.\n" if @_!= 3;
	my $array = $_[0];
	my $arrayIdx = $_[1];
	my $geneIdx = $_[2];
	
	#Carp::croak join ("\t", @$geneIdx), "\n";
	my $arrayNew;
	$arrayNew->{"hasCls"} = $array->{"hasCls"};
	$arrayNew->{"arrayN"} = @$arrayIdx;
	$arrayNew->{"geneN"} = @$geneIdx;
	#print "array=", $arrayNew->{"arrayN"}, "\t gene=", $arrayNew->{"geneN"}, "\n";
	
	die "subArrayMat: incorrect parameters.\n" if ($arrayNew->{"arrayN"} <= 0 || $arrayNew->{"geneN"} <= 0);
	my $arrayNameTmp = $array->{"arrayNames"};
	my @arrayNameTmp = @$arrayNameTmp;
	@arrayNameTmp = @arrayNameTmp[@$arrayIdx];
	
	my $geneNameTmp = $array->{"geneNames"};
	my @geneNameTmp = @$geneNameTmp;	
	@geneNameTmp = @geneNameTmp[@$geneIdx];

	$arrayNew->{"arrayNames"} = \@arrayNameTmp;
	$arrayNew->{"geneNames"} = \@geneNameTmp;

	my ($g, $a);
	
	for ($g = 0; $g < $arrayNew->{"geneN"}; $g++)
	{
		my $dataTmp = $array->{"data"}->[$geneIdx->[$g]];
		my @dataTmp = @$dataTmp;
		@dataTmp = @dataTmp[@$arrayIdx];
		$arrayNew->{"data"}->[$g] = \@dataTmp;
	}
	
	if ($arrayNew->{"hasCls"} == 1)
	{
		my $clsTmp = $array->{"cls"};
		my @clsTmp = @$clsTmp;
		@clsTmp = @clsTmp[@$arrayIdx];
		$arrayNew->{"cls"} = \@clsTmp;
	}	
	return $arrayNew;
}	

#@selArray@
#@Input:
#@ $array = $_[0];
#@ $arrayIdx = $_[1];
#@Output:
#@ $arrayNew
#@ status: tested
sub selArray
{
	die "selArray: incorrect number of parameters.\n" if @_!= 2;
	my $array = $_[0];
	my $arrayIdx = $_[1];
	my @geneIdx = (0..($array->{"geneN"} - 1));
	my $arrayNew = subArrayMat ($array, $arrayIdx, \@geneIdx);
	return $arrayNew;
}	

#@selGene@
#@Input:
#@ $array = $_[0];
#@ $geneIdx = $_[1];
#@Output:
#@ $arrayNew
#@ status: tested
sub selGene
{
	die "selGene: incorrect number of parameters.\n" if @_!= 2;
	my $array = $_[0];
	my $geneIdx = $_[1];
	my @arrayIdx = (0..($array->{"arrayN"} - 1));
	my $arrayNew = subArrayMat ($array, \@arrayIdx, $geneIdx);
	return $arrayNew;
}

sub permuteLabel
{
	die "permuteLabel: incorrect number of parameters.\n" if @_ != 1;
	my $array = $_[0];
	$array->{"cls"} = Common::shuffleArray ($array->{"cls"});
#	return $array;
}	
#@delGene@
#@Input:
#@ $array = $_[0];
#@ $geneIdx = $_[1];
#@Output:
#@ $arrayNew
#@ status: NOT tested
sub delGene
{
	die "delGene: incorrect number of parameters.\n" if @_!= 2;
	my $array = $_[0];
	my $geneIdx = $_[1];
	my @geneIdxAll = (0..($array->{"geneN"} - 1));
	
	$geneIdx = Common::diffArray (\@geneIdxAll, $geneIdx);
	
	my @arrayIdx = (0..($array->{"arrayN"} - 1));
	my $arrayNew = subArrayMat ($array, \@arrayIdx, $geneIdx);
	return $arrayNew; 
}	

#@delArray@
#@Input:
#@ $array = $_[0];
#@ $arrayIdx = $_[1];
#@Output:
#@ $arrayNew
#@ status: NOT tested
sub delArray
{
	die "delArray: incorrect number of parameters.\n" if @_!= 2;
	my $array = $_[0];
	my $arrayIdx = $_[1];
	my @arrayIdxAll = (0..($array->{"arrayN"} - 1));
	
	$arrayIdx = Common::diffArray (\@arrayIdxAll, $arrayIdx);
	
	my @geneIdx = (0..($array->{"geneN"} - 1));
	my $arrayNew = subArrayMat ($array, $arrayIdx, \@geneIdx);
	return $arrayNew; 
}	

1;

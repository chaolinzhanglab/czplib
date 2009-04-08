# create at Nov 22, 2004
# by Chaolin Zhang
#
# IO manipulation of microarray data
#!/usr/bin/perl -w
package Microarray::arrayIO;

use strict;
use FileHandle;


#@readArray@
#@Input:
#@ $dataFile:
#@ $clsFile: class information, optional
#@Output:
#@ $array
sub readArray
{
	my $dataFile;
	my $clsFile;
	my $hasCls = 0;	
	
	my $array;
	my @arrayName;
	my @geneName;
	my @cls;
	my $data;
	my $arrayN = 0;
	my $geneN = 0;
	
	die "readArray: incorrect parameters, die\n" if (@_<1);
	$dataFile = $_[0];
	#print "dataFile=$dataFile\n";
	if (@_ > 1)
	{
		$clsFile = $_[1];
		$hasCls = 1;
	}	
		

	
	my $fd = new FileHandle;
	
	open ($fd, "<$dataFile") || die "readArray: failed to read data from $dataFile\n";
	#read array name
	my $line = <$fd>;
	chomp $line;
	@arrayName = split ("\t", $line);
	shift @arrayName;
	$arrayN = @arrayName;
	
	#read data
	while ($line = <$fd>)
	{
		chomp $line;
		my @tmp = split ("\t", $line);
		$geneName[$geneN] = shift @tmp;
		while (@tmp < $arrayN)
		{
			push @tmp, "";
		}
		if (@tmp != $arrayN)
		{
			die "line number: $geneN, arrayNum=$arrayN, column number in this line:", $#tmp+1, "\n";
		}
		$data->[$geneN] = \@tmp;
		$geneN++;
	}		
	close ($fd);
	if ($hasCls == 1)
	{
		open ($fd, "<$clsFile") || die "readArray: failed to read data from $clsFile\n";
		my $i = 0;
		while ($line = <$fd>)
		{
			chomp $line;
			$line =~/^([\-\d]+)/;
			$cls[$i] = $1;
			$i++;
		}	
		close ($fd);
		die "readArray: data incomplete, the number of array in data file is not equal to that in class file\n" if ($i != $arrayN);
		$array->{"cls"} = \@cls;
	}	
	$array->{"hasCls"} = $hasCls;
	$array->{"data"} = $data;
	$array->{"arrayNames"} = \@arrayName;
	$array->{"geneNames"} = \@geneName;
	$array->{"arrayN"} = $arrayN;
	$array->{"geneN"} = $geneN;
	return $array;
}	

#@writeArray@
#@Input:
#@ $array:
#@ $dataFile:
#@ $classFile: class information, optional
sub writeArray
{
	die "writeArray: incorrect parameters, die\n" if (@_<2);

	my $array = $_[0];
	#my $array = $$arrayref;
	my $dataFile = $_[1];
	my $clsFile;
	my $outCls = 0;	
	
	if (@_ > 2)
	{
		$clsFile = $_[2];
		$outCls = 1;
	}	
	my $fd = new FileHandle;
	
	open ($fd, ">$dataFile") || die "writeArray: failed to write data to $dataFile\n";
	
	#output array names
	print $fd "ProbeID";
	my ($g, $a);
	#print "OK1\n";
	for ($a = 0; $a < $array->{"arrayN"}; $a++)
	{
		print $fd "\t", $array->{"arrayNames"}->[$a];
	}	
	print $fd "\n";
	
	#print "OK2\n";
	#output data
	for ($g = 0; $g < $array->{"geneN"}; $g++)
	{
		print $fd $array->{"geneNames"}->[$g];
		for ($a = 0; $a < $array->{"arrayN"}; $a++)
		{
				#print "g = $g, a = $a\n";
			if (defined $array->{"data"}->[$g]->[$a] && $array->{"data"}->[$g]->[$a] ne '')
			{
				print $fd "\t", $array->{"data"}->[$g]->[$a];
			}
			else
			{
				print $fd "\t", "NA";
			}
		}
		print $fd "\n";
	}	
	close ($fd);
	return if ($outCls == 0);
	
	if ($array->{"hasCls"} == 0)
	{
		warn "writeArray: warining: failed to output class information\n";
		return;
	}	
	
	open ($fd, ">$clsFile") || die "Failed to write data to $clsFile\n";
	#print "OK3\n";
	for ($a = 0; $a < $array->{"arrayN"}; $a++)
	{
		print $fd $array->{"cls"}->[$a], "\n";
	}	
	close ($fd);	
}	

#@writeTorchFormat@
#@Input:
#@ $array:
#@ $dataFile:

sub writeTorchFormat
{
	die "writeTorchFormat: incorrect parameters, die\n" if (@_!=2);
	my $array = $_[0];
	my $dataFile = $_[1];
	die "writeTorchFormat: no class information, die\n" if ($array->{"hasCls"} != 1);
	
	my $fd = new FileHandle;
	open ($fd, ">$dataFile") || die "writeTorchFormat: failed to write to file $dataFile\n";
	print $fd $array->{"arrayN"}, "\t", ($array->{"geneN"} + 1), "\n";
	my ($a, $g);
	for ($a = 0; $a < $array->{"arrayN"}; $a++)
	{
		for ($g = 0; $g < $array->{"geneN"}; $g++)
		{
			print $fd $array->{"data"}->[$g]->[$a], "\t";
		}	
		print $fd $array->{"cls"}->[$a], "\n";
	}	
	close ($fd);
}	

1;

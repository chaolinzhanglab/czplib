##------------------------------------------------------------------------
## ADAM - [A]utomatic [D]etector for [A]bstracts from Pub[M]ed version 1.1
##							 By
## 						Chaolin Zhang 
##
##	  MOE Key Lab for Bioinformatics/Dept of Automation
##					Tsinghua University
##				Beijing, 100084, PR China
##			Email: zcl98@mails.tsinghua.edu.cn
##					Tel: +86-10-62795578
##------------------------------------------------------------------------

##
# File Name: GeneSymbol.pl
# Package Name: GeneSymbol
# Purpose:   Build a hash table of names of gene/protein/other molecules
#			 Note that 
#				(1) The gene list file should be the input.
#				(2)	In the gene list file, each row is the names of a 
#					gene/protein/molecule dilimited by tabs, with the 
#					stardard name at the first column. So the standard 
#					names should be UNIQUE.	
#				(3) If the short descriptions/previous names are not 
#					available the column should be blank. 		
# Create date: 			May 16, 2002
# Last modified:		Dec 2, 2003
# 	
# Example:
# my $hashObj = new GeneSymbol ('/path/to/dir/nomeids.txt');
# my $symbol_hash = $hashObj->build_hash ();
# my $key;
# foreach $key (keys %$symbol_hash)
# {
#		my $symbol = $symbol_hash->{$key}{'symbol'};
#		my $name = $symbol_hash->{$key}{'name'};
#		my $alias = $symbol_hash->{$key}{'alias'};
#		my @alias = @$alias;
# }
##


package GeneSymbol;

use strict;
use FileHandle;
use Common;

our $VERSION = 1.01;

#-------------------------------Globals--------------------------------------# 
my $symbol_file;
my $hash;
#----------------------------------------------------------------------------#

sub new
{ 
	my $class = shift;
	$symbol_file = shift;
	my $self = bless 
				{
					build_hash => undef,
					debug => 0
                    # variable initialization goes here
				};
	return $self;
} # new

sub build_hash
{
	my $self = shift;
	my $fd = new FileHandle;
	open ($fd, "<$symbol_file") || die "Failed to open file $symbol_file to read\n";
	my $line;
	my $iter = 0;
	while ($line = <$fd>)
	{
		$line = Common::chop_carriage ($line);
		next if ($line =~/^\s*$/);
		
		$iter++;
		my @cols = split (/\s*?\t+\s*?/, $line);
		my $boundary = shift @cols;

		my $geneId = shift @cols;
		
		#my $symbol = shift @alias;  #the first column is the symbol or standard name
		my @symbols;
		foreach (@cols)
		{
			push @symbols, $_ if length($_) > 1;
		}	
		$hash->{$geneId}->{"boundary"} = $boundary;
		push @{$hash->{$geneId}->{"keywords"}}, @symbols;
	}
	close ($fd);

	#remove redundant keywords if any
	foreach my $geneId (keys %$hash)
	{
		my $symbols = $hash->{$geneId}->{"keywords"};
		my %symbols;

		foreach my $s (@$symbols)
		{
			$symbols{$s} = 1;
		}
		my @symbols = keys %symbols;
		$hash->{$geneId}->{"keywords"} = \@symbols;
	}
	return $hash;
}

sub dump_hash
{
	my $self = shift;
	my $key; 
	print "id\tsymbols\n";
	foreach $key (sort (keys %$hash))
	{
		my $id = $key; #hash->{$key}{'id'};
		my $symbols = $hash->{$key}->{"keywords"};
		print $id, "\t", join (",", @$symbols), "\n";
	}
}
1;

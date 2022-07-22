#
#===============================================================================
#
#         FILE:  Vcf.pm
#
#  DESCRIPTION:  Package to handle Vcf file
#         BUGS:  ---
#        NOTES:  
#       AUTHOR:  Chaolin Zhang (cz), cz2294@columbia.edu
#      COMPANY:  Columbia University
#      VERSION:  1.0
#      CREATED:  12/07/13
#     REVISION:  ---
#===============================================================================


package Vcf;

require Exporter;

our $VERSION = 1.01;

@ISA = qw (Exporter);

@EXPORT = qw(
	assignAltBase
	generateVcfHeader
	lineToVcf
	readVcfFile
	splitVcfFileByChrom
	vcfToLine
	writeVcfFile
);
	
=head1 Bed

Subroutines to handle BED file

=cut

use strict;
use warnings;

use Data::Dumper;
use Carp;

use Common;


#
#object-oriented interface will be added later
#

my $debug = 0;

my @vcfInfoField = (
		['DP4', 0,'s'],
		['BC', 1, 'i'],
		['snp', 2,'s'],
		['refStrandP', 3,'f'],
		['altStrandP', 4,'f'],
		['altP', 5,'f'],
		['strand', 6,'s'],
		['gene', 7,'s'],
		['region',8,'s']
);

my %vcfInfoFieldHash = (
	'DP4'=>[0,'s'],
	'BC'=>[1,'i'],
	'snp'=> [2,'s'],
	'refStrandP'=>[3,'f'],
	'altStrandP'=>[4,'f'],
	'altP'=>[5,'f'],
	'strand'=>[6,'s'],
	'gene'=>[7,'s'],
	'region'=>[8,'s']
);


sub readVcfFile
{
	my ($inFile, $verbose, $msgio) = @_;

	$msgio = *STDOUT unless $msgio;
	my $validateFields = 0; #prefixed for now.

	my $fin;
	my @ret;

	my $iter = 0;
	open ($fin, "<$inFile") || Carp::croak "cannot open file $inFile to read\n";
	while (my $line =<$fin>)
	{
		chomp $line;
		next if $line =~/^\#/;
		next if $line =~/^\s*$/;
		print $msgio "$iter ...\n" if $verbose && $iter % 500000 == 0;
		$iter++;

		my $snv = lineToVcf ($line, $validateFields);
		push @ret, $snv;
	}	
	close($fin);
	return \@ret;
}


sub generateVcfHeader
{
	
	my $ret = "#" . join ("\t", qw(CHROM POS ID REF ALT QUAL FILTER INFO));
	return $ret;
}

sub writeVcfFile
{
	my ($sites, $outFile, $verbose, $msgio) = @_;

	$msgio = *STDOUT unless $msgio;

	my $fout;
	open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";
	print $fout generateVcfHeader (), "\n";

	for (my $i = 0; $i < @$sites; $i++)
	{
		print $msgio "$i ...\n" if $verbose && $i % 500000 == 0;
		print $fout vcfToLine ($sites->[$i]), "\n";
	}
	close ($fout);
}



#outDir has to be an existing dir

sub splitVcfFileByChrom
{
	my ($inFile, $outDir, $verbose, $msgio) = @_;

	if ($inFile ne '-')
	{
		Carp::croak "$inFile does not exist\n" unless -f $inFile;
	}

	Carp::croak "$outDir does not exist\n" unless -d $outDir;


	$msgio = *STDOUT unless $msgio;

	my %siteCount;

	my $fin;
	#print "reading tags from $inBedFile ...\n" if $verbose;
	
	if ($inFile eq '-')
	{
		$fin = *STDIN;
	}
	else
	{
		if ($inFile =~/\.gz$/)
		{
			open ($fin, "gunzip -c $inFile |") || Carp::croak "can not open file $inFile to read\n";
		}
		elsif ($inFile =~/\.bz2$/)
        {
            open ($fin, "bunzip2 -c $inFile |") || Carp::croak "can not open file $inFile to read\n";
        }
		else
		{
			open ($fin, "<$inFile") || Carp::croak "can not open file $inFile to read\n";
		}
	}

	#my %tagCount;
	#write tags according to chromosomes
	my $i = 0;
	my %fhHash;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		next if $line=~/^\#/;

		$i++;

		print $msgio "$i ...\n" if $i % 500000 == 0 && $verbose;
	
		$line =~/(\S+)\s/;

		my $chrom = $1;
		my $tmpFile = "$outDir/$chrom.vcf";
	
		if (not exists $fhHash{$chrom})
		{
			my $fout;
			open ($fout, ">>$tmpFile") || Carp::croak "can not open file $tmpFile to write: $?\n";
			$fhHash{$chrom} = $fout;
		}
	
		#open ($fout, ">>$tmpFile") || Carp::croak "can not open file $tmpFile to write\n";
	
		my $fout = $fhHash{$chrom};
		print $fout $line, "\n";
		#close ($fout);

		if (exists $siteCount{$chrom})
		{
			$siteCount{$chrom}->{"n"} ++;
		}
		else
		{
			$siteCount{$chrom} = {f=>"$outDir/$chrom.vcf", n=> 1};
		}
	}
	close ($fin) if $inFile ne '-';

	#close all file handles
	foreach my $chrom (sort keys %fhHash)
	{
		close ($fhHash{$chrom});
	}
	return \%siteCount;
}


=head2 lineToVcf

Pass a line of Vcf file

my $vcf = lineToVcf ($line;

=cut

sub lineToVcf
{
	my ($line, $validateFields) = @_;

	my ($chrom, $position, $id, $refBase, $altBase, $qual, $filter, $info, @ignored) = split ("\t", $line);
	my $infoHash = parseVcfInfo ($info, $validateFields);

	my %vcfHash = (chrom=>$chrom, position=>$position-1, id=>$id, refBase=>$refBase, altBase=>$altBase, qual=>$qual, filter=>$filter, info=>$infoHash);
	return \%vcfHash;
}


sub vcfToLine
{
	my $vcf = $_[0];
	my $validateInfoFields = 0;
	my $infoStr = encodeVcfInfo ($vcf->{'info'});
	return join ("\t", $vcf->{'chrom'}, $vcf->{'position'}+1, $vcf->{'id'}, $vcf->{'refBase'}, $vcf->{'altBase'}, $vcf->{'qual'}, $vcf->{'filter'}, $infoStr);
}


sub parseVcfInfo
{
	my ($infoStr, $validateFields) = @_;
	
	my %infoHash;
	my @pairs = split (";", $infoStr);
	foreach my $pair (@pairs)
	{
		next if $pair eq '.';

		my $name = $pair;
		if ($pair=~/^(.*?)\=(.*?)$/)
		{
			#a pair
			$infoHash{$1} = $2;
			$name = $1;
		}
		else
		{
			#not a pair
			$infoHash{$pair} = 1;
		}

		if (defined $validateFields && $validateFields == 1)
		{
			Carp::croak "$name is not a valid VCF info field\n" unless exists $vcfInfoFieldHash{$name};
		}
	}
	return \%infoHash;
}



sub encodeVcfInfo
{
	my ($info, $validateFields) = @_;
	
	if (defined $validateFields && $validateFields == 1)
	{
		#make sure every field is considered
		foreach my $name (keys %$info)
		{
			Carp::croak "$name is not a valid VCF info field\n" unless exists $vcfInfoFieldHash{$name};
		}
	}

	my $ret = "";
	foreach my $f (@vcfInfoField)
	{
		my $name = $f->[0];
		next unless exists $info->{$name};
		if ($f->[2] ne '')
		{
			$ret .= join ("=", $name, $info->{$name}) . ";";
		}
		else
		{
			$ret .= $name . ";";
		}
	}
	chop $ret if $ret ne '';
	return $ret;
}



#
sub assignAltBase
{
	my ($readBaseHash, $refBase) = @_;
	my %readBaseHash = %$readBaseHash;
	
	#add a small number to resolve ties
	map {$readBaseHash{$_} += rand(1) * 0.01} qw(A C G T);

	my @readBaseSort = sort {$readBaseHash{$b} <=> $readBaseHash{$a}} qw(A C G T);
	my $altBase = '';
    
	foreach my $b (@readBaseSort)
    {
    	if ($b ne $refBase)
        {
                $altBase = $b;
                last;
        }
    }
	return $altBase;
}


1;



#
#===============================================================================
#
#         FILE:  Motif.pm
#
#  DESCRIPTION:  Package to handle Motif (matrix or consensus)data
#         BUGS:  ---
#        NOTES:  The package is spun off from the AnnotationIO.pm
#       AUTHOR:  Chaolin Zhang (cz), czhang@rockefeller.edu
#      COMPANY:  Rockefeller University
#      VERSION:  1.0
#      CREATED:  12/16/2010
#     REVISION:  ---
#===============================================================================



package Motif;

require Exporter;

our $VERSION = 1.03;

@ISA = qw (Exporter);

@EXPORT = qw (
	allowIUB
	countMismatch
	countMismatchBase
	countToBayesianMatrix
	countToFrequencyMatrix
	countToStormoMatrix
	getMatrix
	getMatrixScore
	getMaxMatrixScore
	getMinMatrixScore
	maskWord
	printMotif
	readMEMEMotifFile
	readMotifFile
	readMatchFile
	revComMatrix
	searchWord
	writeMotifFile
);

=head1 NAME

Motif - subroutines to handle motifs

=cut

use strict;
use warnings;

use Data::Dumper;
use Carp;

use Common;



#handle matrix based motifs


=head2 readMotifFile

augmented TRANSFAC motif file

=cut


sub _initMotif
{
	my $motif = {
		AC=>''#,	#accession
		#TY=>'',	#'Motif'
		#ID=>'', #ID
		#DT=>[], #date
		#NA=>'',	#name
		#DE=>'',	#description
		#BF=>'', #
		#MT=>[],	#matrix
		#CO=>'',	#copy right
		#BA=>'', #
		#BS=>[],	#site sequence
		#AT=>{},	#attributes
	};
	return $motif;
}


#private
sub readMEMEMotifFile
{
	my ($in) = @_;
	my $prefix = "MEME";

	my $fin;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	my $line;
	my @motif = ();
	my ($BEGIN, $SEQ, $HEAD, $MATRIX) = qw (bg sq hd st mt);
	my $matrix = ();
	my $sites = ();
	my ($width, $nsites, $evalue) = (-1, -1, -1);
	my $status = $BEGIN;
	my $id = 0;
	while ($line =<$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		
		#print $line, "\n";
		if (($line =~/in BLOCKS format/) && $status eq $BEGIN)
		{
			#print "motif hits...\n";
			$line = <$fin>;
			$line = <$fin>;
			while ($line =<$fin>)
		   	{	if ($line!~/^\/\//)
				{
					chomp $line;
					$line =~/^(.*?)\s+\(\s+(\d+)\)\s+(.*?)\s+(.*?)$/;
					my $seq = $1;
					my $pos = $2;
					my $site = $3;
					my $strand = ($4 == 1)?'p':'n';
					push @$sites, {site=>$site, pos=>$pos, seq=>$seq, strand=>$strand};
				}
				else
				{
					$status = $SEQ;
					last;
				}
			}
			next;
		}
	
		if ($line =~/position\-specific probability matrix/ && $status eq $SEQ)
		{
			$status = $HEAD;
			next;
		}

		if ($line =~m@letter-probability matrix: alength= \d+ w= (\d+) nsites= (\d+) E= (.*?)$@ 
				&& $status eq $HEAD)
		{
			$status = $MATRIX;
			$width = $1;
			$nsites = $2;
			$evalue = $3;
			next;
		}
		
		if ($line !~/^\---/ && $status eq $MATRIX)
		{
			my @cols = split (/\s+/, $line);
			shift @cols;
			push @$matrix,
			{A=>$cols[0], C=>$cols[1], G=>$cols[2], T=>$cols[3]};
			next;
		}
		
		if ($line =~/^---/ && $status eq $MATRIX)
		{
			my $motifWidth = @$matrix;
			my $motifId = $prefix . "_W$motifWidth" . "_" . $id;
	
			$id++;
			push @motif, {
					ID=>$motifId,
					AC=>$motifId,
					attr=>{WIDTH=>$width, NSITES=>$nsites, EVALUE=>$evalue}, 
					MT=>$matrix, 
					BS=>$sites
			};
			$matrix = ();
			$sites = ();
			$width = $nsites = $evalue = -1;
			$status = $BEGIN;
		}
	}
	close ($fin);
	return \@motif;
}



=head2 readMotifFile

Transfac Motif file
references are ignored

my $motifs = readMotifFile (fileName, returnHash=0)

if returnHash == 1: $motifs will be a hash, with motif accession used as keys
otherwise: an array

=cut

sub readMotifFile
{
	my $in = $_[0];
	my $returnHash = 0;

	if (@_ > 1)
	{
		$returnHash = $_[1];
	}
	
	my $motifs;

	my $fin;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	my $line;
	
	my $motif;
	while ($line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		next if $line=~/^XX/;

		if ($line =~/^AC/)
		{
			$line =~/^AC\s+(\S+.*?)$/;
			$motif = _initMotif ();
			$motif->{'AC'} = $1;
			next;
		}

		if ($line =~/^TY/)
		{
			$line =~/^TY\s+(\S+.*?)$/;
			$motif->{'TY'} = $1;
			next;
		}

		if ($line =~/^ID/)
		{
			$line =~/^ID\s+(\S+.*?)$/;
			$motif->{'ID'} = $1;
			next;
		}
	
		if ($line =~/^DT/)
		{
			$line =~/^DT\s+(\S+.*?)$/;
			push @{$motif->{'DT'}}, $1;
			next;
		}

		if ($line =~/^NA/)
		{
			$line =~/^NA\s+(\S+.*?)$/;
			$motif->{'NA'} = $1;
			next;
		}
		if ($line =~/^DE/)
		{
			$line =~/^DE\s+(\S+.*?)$/;
			$motif->{'DE'} = $1;
			next;
		}
		if ($line =~/^BF/)
		{
			$line =~/^BF\s+(\S+.*?)$/;
			$motif->{'BF'} = $1;
			next;
		}

		if ($line =~/^IC/)
		{
			$line =~/^IC\s+(\S+.*?)$/;
			my @cols = split (/\s+/, $1);
			$motif->{'IC'} = \@cols;
		}
		
		if ($line =~/^P0/)
		{
			while ($line = <$fin>)
			{
				chomp $line;
				last unless $line =~/\d+/;
				
				#print $line, "\n";
				my @cols = split (/\s+/, $line);
				shift @cols;
				#my $s = $cols[0] + $cols[1] + $cols[2] + $cols[3];
				push @{$motif->{'MT'}},
					{A=>($cols[0]), C=>($cols[1]), G=>($cols[2]), T=>($cols[3])};
			}
			#no next here
		}

		if ($line =~/^XL/)
		{
			while ($line = <$fin>)
            {
                chomp $line;
                last unless $line =~/\d+/;

                #print $line, "\n";
                my @cols = split (/\s+/, $line);
                push @{$motif->{'XL'}}, $cols[1];
            }
		}
		
		if ($line =~/^CO/)
		{
			$line =~/^CO\s+(\S+.*?)$/;
			$motif->{'CO'} = $1;
			next;
		}

		if ($line =~/^BA/)
		{
			$line =~/^BA\s+(\S+.*?)$/;
			$motif->{'BA'} = $1;
			next;
		}

		if ($line =~/^AT/)
		{
			#print $line, "\n";
			my @cols = split (/\s+/, $line);
			#print join ("\n|", @cols), "\n";
			my $at = pop @cols;
			my ($name, $value) = split (/\=/, $at);
			$motif->{'AT'}->{$name} = $value;
			#print "name =$name, value=$value\n";
			next;
		}
		
		if ($line =~/^BS\s/)
		{
			my @cols = split (/\;\s*/, $line);
			#print "seq = \n";
			#print join ("\n|", @cols), "\n";
			#shift @cols;
			my $site = $cols[0];
			$site =~/\s+(\S*?)$/;
			$site = $1;
			
			my $siteInfo = {
				'site'=>$site,
				'seq'=>$cols[1],
				'pos'=>$cols[2], #already zero-based coordinates
				'len'=>$cols[3],
				'len2'=>$cols[4],
				'strand'=>$cols[5]
			};

			if (@cols > 6)
			{
				$siteInfo->{'score'} = $cols[6];
			}

			push @{$motif->{'BS'}}, $siteInfo;
			next;
		}
		
		if ($line =~/^\/\//)
		{
			if ($returnHash)
			{
				$motifs->{$motif->{'AC'}} = $motif;
			}
			else
			{
				push @$motifs, $motif;
			}
		}
	}
	close ($fin);
	return $motifs;
}


=head2 printMotif

printMotif ($motif)

=cut

sub printMotif
{
	my $motif = $_[0];

	print join ("\t", 'AC', $motif->{'AC'}), "\n";
	print "XX\n";

	if (exists $motif->{'TY'} && $motif->{'TY'} ne '')
	{
		print join ("\t", 'TY', $motif->{'TY'}), "\n";
		print "XX\n";
	}
	
	if (exists $motif->{'ID'} && $motif->{'ID'} ne '')
	{
		print join ("\t", 'ID', $motif->{'ID'}), "\n";
		print "XX\n";
	}
	
	my $printSep = 0;
	if (exists $motif->{'DT'})
	{
		for (my $i = 0; $i < @{$motif->{'DT'}}; $i++)
		{
			print join ("\t", 'DT', $motif->{'DT'}->[$i]), "\n";
		}
		$printSep = 1;
	}
	if (exists $motif->{"CO"} && $motif->{'CO'} ne '')
	{
		print join ("\t", "CO", $motif->{'CO'}), "\n";
		$printSep = 1;
	}
	print "XX\n" if $printSep;

	if (exists $motif->{'NA'} && $motif->{'NA'} ne '')
	{
		print join ("\t", 'NA', $motif->{'NA'}), "\n";
		print "XX\n";
	}

	if (exists $motif->{'DE'} && $motif->{'DE'} ne '')
	{
		print join ("\t", 'DE', $motif->{'DE'}), "\n";
		print "XX\n";
	}
	
	if (exists $motif->{'BF'} && $motif->{'BF'} ne '' )
	{
		print join ("\t", 'BF', $motif->{'BF'}), "\n";
		print "XX\n";
	}

	print "P0\t", join ("\t",
		sprintf ("%4s", "A"),
		sprintf ("%4s", "C"),
		sprintf ("%4s", "G"),
		sprintf ("%4s", "T")), "\n";

	my $width = @{$motif->{'MT'}};
	
	for (my $i = 0; $i < $width; $i++)
	{
		my $pos = $motif->{"MT"}->[$i];
		printf ("%02d\t", $i+1);
		my $temp = $motif->{"MT"}->[0]->{"A"} 
			+ $motif->{"MT"}->[0]->{"C"}
			+ $motif->{"MT"}->[0]->{"G"}
			+ $motif->{"MT"}->[0]->{"T"};
		if (Common::ABS ($temp - int($temp + 0.5)) < 1e02 && $temp > 2)
		{	
			print join ("\t", 
				sprintf ("%4d", $pos->{"A"}), 
				sprintf ("%4d", $pos->{"C"}), 
				sprintf ("%4d", $pos->{"G"}), 
				sprintf ("%4d", $pos->{"T"})), "\n";
		}
		else
		{
			print join ("\t", 
				sprintf ("%.4f", $pos->{"A"}), 
				sprintf ("%.4f", $pos->{"C"}), 
				sprintf ("%.4f", $pos->{"G"}), 
				sprintf ("%.4f", $pos->{"T"})), "\n";
		}
	}
	print "XX\n";

	if (exists $motif->{'BA'} && $motif->{'BA'} ne '')
	{
		print join ("\t", "BA", $motif->{'BA'}), "\n";
		print "XX\n";
	}
	

	if (exists $motif->{"AT"})
	{
		my $attrs = $motif->{'AT'};
		foreach (keys %$attrs)
		{
			print "AT\t" . $_ . "=" . $attrs->{$_}, "\n";
		}
		print "XX\n";
	}
	
	if (exists $motif->{"BS"})
	{
		for (my $i = 0; $i < @{$motif->{"BS"}}; $i++)
		{
			my $siteInfo = $motif->{"BS"}->[$i];
			print join ("\t", 
				"BS", 
				join ("; ",
					$siteInfo->{"site"},
					$siteInfo->{"seq"},
					$siteInfo->{"pos"}, # 0-based coordinates
					$siteInfo->{"len"},
					$siteInfo->{"len2"},
					$siteInfo->{"strand"}
				));
			if (exists $siteInfo->{"score"})
			{
				print "; ", $siteInfo->{"score"};
			}
			print "\n";
		}
		print "XX\n";
	}
	print "//\n";
}


=head2 writeMotifFile

writeMotifFile ($motifs, $outFile);

=cut

sub writeMotifFile
{
	my ($motifs, $out) = @_;

	my $fout;
	open ($fout, ">$out") || Carp::croak "can not open file $out to write\n";
	
	my $oldFout = select ($fout);
	my $i;
	for ($i = 0; $i < @$motifs; $i++)
	{
		printMotif ($motifs->[$i]);
	}
	select ($oldFout);
	close ($fout);
}


=head2 readMatchFile

motif files used by the Match program

my $motifs = readMatchFile ($inFile)

=cut
sub readMatchFile
{
	my $in = $_[0];

	my %siteHash;
	my $fin;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	my $seqId;
	while (my $line =<$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		
		if ($line =~/^Inspecting sequence ID\s+(.*?)$/)
		{
			$seqId = $1;
			next;
		}
		if ($line =~/^\s+(\S+)\s+\|\s+(\d+)\s+\((\S)\)\s+\|\s+(\S+)\s+\|\s+(\S+)\s+\|\s+(\S+)$/)
		{
			my $siteInfo = {matrixId=>$1,
					pos=>$2 -1, #convert to 0-based coordinates
					strand=>$3,
					scoreCore=>$4,
					score=>$5,
					site=>$6,
					seq=>$seqId};
			push @{$siteHash{$1}->{$seqId}}, $siteInfo;
			#print Dumper ($siteInfo), "\n";
		}
	}
	close ($fin);
	return \%siteHash;
}


=head2 getMatrix

my $matrix = getMatrix($motif);

=cut
sub getMatrix 
{
	my ($motif) = @_;
	Carp::croak "no matrix found in the motif: ", Dumper ($motif), "\n" unless exists $motif->{'MT'};
	return $motif->{'MT'};
}


=head2 countToStormoMatrix

Note:
original name: ToStormoMatrix

Convert count matrix to stormo matrix
see AnnotationIO for the data structure of matrix

my $stormo = ToStormoMatrix ($matrix)

status: tested
Date: 09/29/2006

=cut


sub countToStormoMatrix
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

=head2 countToFrequencyMatrix

note: original name: ToFrequencyMatrix

my $freqMat = countToFrequencyMatrix ($countMatrix);

#Date: 09/29/2006

=cut
sub countToFrequencyMatrix
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

=head2 countToBayesianMatrix

note: original name: ToBayesianMatrix
my $bMat = countToBayesianMatrix ($matrix)

#Date: 09/29/2006
=cut

sub countToBayesianMatrix
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


=head2 revComMatrix

remove complement a matrix

Note: original name is RevComMatrix

my $matrix_revcomp = revComMatrix($matrix)

=cut

sub revComMatrix
{
	my $matrix = $_[0];
	my $width = @$matrix;
	
	my @matrixRC;

	for (my $i = 0; $i < $width; $i++)
	{
		$matrixRC[$i] = {A=>$matrix->[$i]->{'T'}, C=>  $matrix->[$i]->{'G'}, G=> $matrix->[$i]->{'C'},T=> $matrix->[$i]->{'A'}};
	}
	@matrixRC = reverse (@matrixRC);
	return \@matrixRC;
}


=head2 getMaxMatrixScore

my $maxScore = ($matrix)

status: tested
Date: 09/29/2006

=cut

sub getMaxMatrixScore
{
	my $matrix = $_[0];
	my $width = @$matrix;
	my $score = 0;
	for (my $p = 0; $p < $width; $p++)
	{
		$score += Common::max ($matrix->[$p]->{'A'}, $matrix->[$p]->{'C'}, $matrix->[$p]->{'G'}, $matrix->[$p]->{'T'});
	}
	return $score;
}



=head2 getMinMatrixScore

my $minScore = ($matrix)

status: tested
Date: 09/29/2006

=cut
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

=head2 getMatrixScore

my $score = getMatrixScore ($matrix, $siteSeq);

status: tested
Date: 09/29/2006

=cut
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
			#in case there are ambiguous nucleotides
			$score += Common::min ($matrix->[$p]->{'A'}, $matrix->[$p]->{'C'}, $matrix->[$p]->{'G'}, $matrix->[$p]->{'T'});
		}
	}
	return $score;
}


#handle pattern-based motifs



#search a nucleotide word in a sequence
#return the start position of each hit
#
#support IUB code


sub allowIUB
{
	my $word = $_[0];
	$word=~s/R/[A|G]/ig;
	$word=~s/Y/[C|T]/ig;
	$word=~s/M/[A|C]/ig;
	$word=~s/K/[G|T]/ig;
	$word=~s/S/[C|G]/ig;
	$word=~s/W/[A|T]/ig;
	$word=~s/H/[A|C|T]/ig;
	$word=~s/B/[C|G|T]/ig;
	$word=~s/V/[A|C|G]/ig;
	$word=~s/D/[A|G|T]/ig;
	$word=~s/N/[A|C|G|T]/ig;
	return $word;
}

sub searchWord
{
	my ($str, $word) = @_;
	my @hits;

	$word=allowIUB ($word);
	
	while ($str =~/($word)/ig)
	{
		my $start = pos ($str) - length ($1);
		my $end = $start + length ($1) - 1;
		push @hits, [$start, $end];
		pos ($str) = $start + 1;
	}
	return \@hits;	
}

sub maskWord
{
	my ($str, $word) = @_;
	my @hits;
	while ($str =~/$word/ig)
	{
		my $start = pos ($str) - length ($word);
		push @hits, $start;
		pos ($str) = $start + 1;
	}
	foreach my $h (@hits)
	{
		substr ($str, $h, length ($word)) = 'N' x length ($word);
	}
	return $str;
}


sub countMismatch
{
	my ($w1, $w2, $ignoreCase, $maxDiff) = @_;

	$w1 =~tr/a-z/A-Z/ if $ignoreCase;
	$w2 =~tr/a-z/A-Z/ if $ignoreCase;

	Carp::croak "unequal length of $w1 and $w2\n" unless length ($w1) == length ($w2);

	my @w1 = split (//, $w1);
	my @w2 = split (//, $w2);

	return countMismatchBase (\@w1, \@w2, $maxDiff);
}

#we assume $w1 and $w2 has the same length
#
# my $n = countMismatchBase (\@w1, \@w2, $maxDiff);

sub countMismatchBase
{
	my ($w1, $w2, $maxDiff) = @_;
	my $n = 0;
	my $l = @$w1;
	$maxDiff = $l unless $maxDiff;	

	for (my $i = 0; $i < $l; $i++)
	{
		$n++ if $w1->[$i] ne $w2->[$i];
		last if $n >= $maxDiff;
	}
	return $n;
}

1;



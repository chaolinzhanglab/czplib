#
#===============================================================================
#
#         FILE:  Bed.pm
#
#  DESCRIPTION:  Package to handle sequence alignment
#         BUGS:  ---
#        NOTES:  The package is spun off from the AnnotationIO.pm and Common.pm
#       AUTHOR:  Chaolin Zhang (cz), czhang@rockefeller.edu
#      COMPANY:  Rockefeller University
#      VERSION:  1.0
#      CREATED:  12/17/2010
#     REVISION:  ---
#===============================================================================

package Align;


require Exporter;

@ISA = qw (Exporter);

@EXPORT = qw (
	blat
	skipPslHeader
	calcPslPercentIdentity
	lineToPsl
	readPslFile
	printPsl
	pslToBed
	pslToLine
	writePslFile
	readNextPslLine
	checkSim4AlignQuality
	printSim4Align
	readSim4File
	sim4ToBed
	lineToVulgar
	readVulgarFile
	vulgarToBed
	vulgarToLine
	writeVulgarFile
	lineToSam
	samToBed
);



=head1 NAME

Align - subroutines that deal with sequence alignment

subroutines starting with a hyphen should not be called outside

=cut

use strict;
use warnings;
use Data::Dumper;
use Carp;


use Common;
use Sam;


=head2 blat

my $psl = blat ($blat, $db, $query, $cache, $options);

$blat   : path to the blat binary
$db     : target database for blat
$query  : query sequence for blat
$cache  : cache dir; should be created before this subroutine is called
$options: a hash ref. with options to be passed on to blat

=cut

sub blat
{
	Carp::croak "3-5 arguments expected\n" unless (@_ >= 4 && @_ <= 5);
	
	my ($blat, $db, $query, $cache, $options) = @_;
	
	#my $cache = "/tmp";
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

	my $result = readPslFile ($out);
	unlink $out if -f $out;
	return $result;
}


=head2 skipPslHeader

my $fin = skipPslHeader ($fin);


=cut

sub skipPslHeader
{
	my $fin = $_[0];
	my $line = <$fin>;
	if ($line=~/^psLayout/)	# there is a header
	{
		while ($line = <$fin>)
		{
			chomp $line;
			last if ($line =~/^\-\-\-/);
		}
	}
	elsif ($line=~/^track/)
	{
	}
	elsif ($line=~/^browser/)
	{
	}
	else
	{
		seek ($fin, 0, 0); # no header line, go to the very beginnning
	}
	return $fin;
}


=head2 readNextPslLine

my $psl = readNextPslLine ($fin);


=cut
sub readNextPslLine
{
	my $fin = $_[0];
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		
		return lineToPsl ($line);
	}
	return "";
}


=head2 lineToPsl

parse a Psl line

my $psl = lineToPsl ($line);
 
 
=cut
sub lineToPsl
{
	my $line = $_[0];
	my @cols = split (/\s+/, $line);
		
	shift @cols if ($#cols == 21); #the first column is the bin id
		
	#print join ("|\n", @cols), "|\n";
	Carp::croak "the number of columns is incorrect:", Dumper (\@cols), "\n" if $#cols != 20;

	my @blockSizes = split (/\,/, $cols[18]); #pop @blockSizes;
	my @qStarts = split (/\,/, $cols[19]); #pop @qStarts;
	#print "qStarts=", join ("|\t|", @qStarts), "|\n";
	my @tStarts = split (/\,/, $cols[20]); #pop @tStarts;
	my $ret = {
		matches=>$cols[0],
		misMatches=>$cols[1],
		repMatches=>$cols[2],
		nCount=>$cols[3],
		qNumInsert=>$cols[4],
		qBaseInsert=>$cols[5],
		tNumInsert=>$cols[6],
		tBaseInsert=>$cols[7],
		strand=>$cols[8], 
		qName=>$cols[9],
		qSize=>$cols[10],
		qStart=>$cols[11],	#the convention is the same as BED file, 0 to len
		qEnd=>$cols[12] - 1,
		tName=>$cols[13],
		tSize=>$cols[14],
		tStart=>$cols[15],
		tEnd=>$cols[16] - 1,
		blockCount=>$cols[17],
		blockSizes=>\@blockSizes,
		qStarts=> \@qStarts, 
		tStarts=>\@tStarts};
	return $ret;
}


=head2 readPslFile

Read UCSC PSL file

my $ret = readPslFile ($pslFile)
header line will be ignored

	$ret = {
		matches=>$cols[0],
		misMatches=>$cols[1],
		repMatches=>$cols[2],
		nCount=>$cols[3],
		qNumInsert=>$cols[4],
		qBaseInsert=>$cols[5],
		tNumInsert=>$cols[6],
		tBaseInsert=>$cols[7],
		strand=>$cols[8], 
		qName=>$cols[9],
		qSize=>$cols[10],
		qStart=>$cols[11],	#the convention is the same as BED file, 0 to len
		qEnd=>$cols[12] - 1,
		tName=>$cols[13],
		tSize=>$cols[14],
		tStart=>$cols[15],
		tEnd=>$cols[16] - 1,
		blockCount=>$cols[17],
		blockSizes=>\@blockSizes,
		qStarts=> \@qStarts, 
		tStarts=>\@tStarts};
	}

=cut

sub readPslFile
{
	my ($in, $verbose) = @_;
	
	my $fin;
	my @results;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	$fin = skipPslHeader ($fin);
	
	my $i = 0;
	while (my $psl = readNextPslLine ($fin))
	{
		
		print "$i ...\n" if $verbose && $i % 100000 == 0;		
		$i++;
		push @results, $psl;	
	
	}
	close ($fin);
	
	return \@results;
}



=head2 obsolete
sub readPslFile
{
	my ($in, $verbose) = @_;
	
	my $fin;
	my @results;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	
	$fin = skipPslHeader ($fin);
	
	my $i = 0;
	while ($line =<$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		
		my @cols = split (/\t/, $line);
		
		shift @cols if ($#cols == 21); #the first column is the bin id
		
		#print join ("|\n", @cols), "|\n";
		Carp::croak ("the number of columns is incorrect\n") if $#cols != 20;

		$i++;
		print "$i ...\n" if $verbose && $i % 100000 == 0;		
		my @blockSizes = split (/\,/, $cols[18]); #pop @blockSizes;
		my @qStarts = split (/\,/, $cols[19]); #pop @qStarts;
		#print "qStarts=", join ("|\t|", @qStarts), "|\n";
		my @tStarts = split (/\,/, $cols[20]); #pop @tStarts;
		push @results, {
				matches=>$cols[0],
				misMatches=>$cols[1],
				repMatches=>$cols[2],
				nCount=>$cols[3],
				qNumInsert=>$cols[4],
				qBaseInsert=>$cols[5],
				tNumInsert=>$cols[6],
				tBaseInsert=>$cols[7],
				strand=>$cols[8], 
				qName=>$cols[9],
				qSize=>$cols[10],
				qStart=>$cols[11],	#the convention is the same as BED file, 0 to len
				qEnd=>$cols[12] - 1,
				tName=>$cols[13],
				tSize=>$cols[14],
				tStart=>$cols[15],
				tEnd=>$cols[16] - 1,
				blockCount=>$cols[17],
				blockSizes=>\@blockSizes,
				qStarts=> \@qStarts, 
				tStarts=>\@tStarts};
	}
	close ($fin);
	return \@results;
}
=cut

=head2 pslToBed

my $bed = pslToBed ($psl, $chrom=0, $useQuery = 0)

convert target coordinates to the BED format by assuming the 
query sequence is the genomic sequence. vice versa if $useQuery == 1

if $chrom == 0, its value is determined according to $useQuery
otherwise, use the specified $chrom value


=cut

#note Psl used absolute coordinates

sub pslToBed
{
	my ($align, $chrom, $useQuery) = @_;
	
	my ($chromStart, $chromEnd, $name, @blockStarts);
	
	if ($useQuery)
	{
		$chrom = $align->{'qName'} unless $chrom;
		$chromStart = $align->{"qStart"};
		$chromEnd = $align->{"qEnd"};
		$name = $align->{"tName"};
		@blockStarts = ();
		foreach my $s (@{$align->{"qStarts"}})
		{
			push @blockStarts, $s - $chromStart;
		}
	}
	else
	{
		#default
		$chrom = $align->{'tName'} unless $chrom;
		$chromStart = $align->{"tStart"};
		$chromEnd = $align->{"tEnd"};
		$name = $align->{"qName"};
    	foreach my $s (@{$align->{"tStarts"}})
		{
			push @blockStarts, $s - $chromStart;
		}
	}
	
	
	my $score = calcPslPercentIdentity ($align, 1);
	#my $strand = '+';#always set to '+', the gene strand, although a few transcripts were aligned to the reverse strand
	my $strand = $align->{'strand'};
	$strand=~/(\S)$/; #for translated alignment, we get the genomic strand
	$strand = $1;
	
	my $blockCount = $align->{"blockCount"};
	my @blockSizes = @{$align->{"blockSizes"}};
	
	my $cov = $align->{'matches'} / $align->{'qSize'};

	my $region = {chrom=>$chrom,			#gene id
		chromStart=>$chromStart,	#start on gene contig
		chromEnd=>$chromEnd, 	#end on gene contig
		name=>$name,
		score=>$score, 
		strand=> $strand,
		thickStart=>$chromStart,
		thickEnd=>$chromEnd,
		itemRgb=>"0,0,0",
		blockCount=>$blockCount,
		blockSizes=>\@blockSizes,
		blockStarts=>\@blockStarts,
		cov=>$cov
	};
	return $region;
}


=head2 printPsl

Note: original name printPslAlign

printPsl ($psl);

=cut


sub printPslAlign
{
	Carp::croak "obsolete function, call printPsl instead\n";
}

sub printPsl
{
	my $psl = $_[0];
	print pslToLine ($psl), "\n";
}


=head2 pslToLine

my $str = pslToLine ($psl);

=cut

sub pslToLine
{
	my ($psl) = @_;
	return join ("\t",
		$psl->{"matches"},
		$psl->{"misMatches"},
		$psl->{"repMatches"},
		$psl->{"nCount"},
		$psl->{"qNumInsert"},
		$psl->{"qBaseInsert"},
		$psl->{"tNumInsert"},
		$psl->{"tBaseInsert"},
		$psl->{"strand"},
		$psl->{"qName"},
		$psl->{"qSize"},
		$psl->{"qStart"},
		$psl->{"qEnd"} + 1,
		$psl->{"tName"}, #tName
		$psl->{"tSize"},
		$psl->{"tStart"},#tStart
		$psl->{"tEnd"} + 1, #tEnd
		$psl->{"blockCount"},
		join (",", @{$psl->{"blockSizes"}}) . ",", #blockSizes
		join (",", @{$psl->{"qStarts"}}) . ",", #qStarts
		join (",", @{$psl->{"tStarts"}}) . ",");  #tStarts
}


=head2 writePslFile

writePslFile ($psls, $outFile, $printHeader = 0);


=cut

sub writePslFile
{
	my ($aligns, $out, $printHeader) = @_;
	
	my $header = "psLayout version 3\n";
	$header .= "match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts\n";
	$header .= "     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count\n";
	$header .= "---------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
	my $fout;

	open ($fout, ">$out") || Carp::croak "can not open file $out to write\n";
	print $fout $header if $printHeader;
	
	my $stdout = select ($fout);
	foreach my $align (@$aligns)
	{	
		print pslToLine ($align), "\n";
	}

	close ($fout);
	select ($stdout);
}



=head2 calcPslPercentIdentity

calculate percent identity from psl line
according to the pslCalcMilliBad function at http://genome.ucsc.edu/FAQ/FAQblat.html

my $identity = calcPslPercentIdentity ($psl, $isMrna);

if $isMrna == 0, $psl is treated as amino acids

=cut

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



#assume seq1 is genomic sequence
#seq2 is transcript (from a db)

sub _splitSim4AlignText
{
	my ($seq1Text, $matchText, $seq2Text) = @_;
	
	my @seq1blocksText;
	my @matchBlocksText;
	my @seq2blocksText;

	my $start = 0;
	while ($matchText =~/\>\>\>\.\.\.\>\>\>/g)
	{
		push @seq1blocksText, substr ($seq1Text, $start, pos($matchText) - 9 - $start);
		push @matchBlocksText, substr ($matchText, $start, pos($matchText) - 9 - $start);
		push @seq2blocksText, substr ($seq2Text, $start, pos($matchText) - 9 - $start);
		$start = pos($matchText);
	}
	
	push @seq1blocksText, substr ($seq1Text, $start);
	push @matchBlocksText, substr ($matchText, $start);
	push @seq2blocksText, substr ($seq2Text, $start);

	return {seq1text=>\@seq1blocksText, matchtext=>\@matchBlocksText, seq2text=>\@seq2blocksText};
}

=head2 readSim4File

note: original name: readsim4File

=cut
sub readSim4File
{
	my $in = $_[0];
	my @ret;
	
	my ($seq1, $seq2, $seq2id, $strand) = ("", "", "", '+');
	my ($seq1len, $seq2len) = (0, 0);
	my $seq1blocks = [];
	my $seq2blocks = [];
	my $blockscore = [];
	
	my $alignText = 0;
	my $alignSeq1 = "";
	my $alignMatch = "";
	my $alignSeq2 = "";

	my $match = 0; #total length of matches
	my $identity = 0; #number of correct bases
	
	my $fin;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	while (my $line=<$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		
		if ($line=~/^seq1/)
		{
			if (@$seq1blocks > 0)
			{
				$identity = int ($identity/100 + 0.5);
				Carp::croak "inconsistency found near line $line\n"
				unless length ($alignSeq1) == length ($alignMatch) && length ($alignSeq2) == length ($alignMatch);
				
				my $item = {seqfile1=>$seq1, seqfile2=>$seq2, seq2id=>$seq2id, 
						seq1len=>$seq1len, seq2len=>$seq2len,
						match=> $match, identity=>$identity, strand=>$strand, 
						seq1blocks=>$seq1blocks,
						seq2blocks=>$seq2blocks,
						blockscore=>$blockscore};

				if (length ($alignSeq1) > 0)
				{
					my $align = _splitSim4AlignText ($alignSeq1, $alignMatch, $alignSeq2);
					$item->{"seq1text"} = $align->{"seq1text"};
					$item->{"matchtext"} = $align->{"matchtext"};
					$item->{"seq2text"} = $align->{"seq2text"};

					my $nblock = @{$align->{"seq1text"}};
					if ($nblock != @{$item->{"seq1blocks"}})
					{
						print "inconsistent block number for gene = $seq1, acc = $seq2id\n"; 
						print Dumper ($item), "\n";
						
						$seq2 = "";
						$seq2id = "";
						$seq2len = 0;
						$seq1blocks = [];
						$seq2blocks = [];
						$blockscore = [];
						$match = 0;
						$identity = 0;
						$strand = '+';
						$alignText = 0;
						$alignSeq1 = "";
						$alignSeq2 = "";
						$alignMatch = "";
	
						next;
					}


					for (my $i = 0; $i < $nblock; $i++)
					{
						my $b1Len = $seq1blocks->[$i]->{"end"} - $seq1blocks->[$i]->{"start"} + 1;
						my $b2Len = $seq2blocks->[$i]->{"end"} - $seq2blocks->[$i]->{"start"} + 1;
						my $b1Text = $align->{"seq1text"}->[$i];
						my $b2Text = $align->{"seq2text"}->[$i];

						$b1Text=~s/\s//g;
						$b2Text=~s/\s//g;

						if ($b1Len != length ($b1Text) || $b2Len != length ($b2Text))
						{
							print "inconsistency size of block $i: seq1 = $b1Len, seq1text = ", length ($b1Text), ", seq2=$b2Len, seq2text=", length($b2Text), "\n";
							print Dumper ($item), "\n";
						
							$seq2 = "";
							$seq2id = "";
							$seq2len = 0;
							$seq1blocks = [];
							$seq2blocks = [];
							$blockscore = [];
							$match = 0;
							$identity = 0;
							$strand = '+';
							$alignText = 0;
							$alignSeq1 = "";
							$alignSeq2 = "";
							$alignMatch = "";
	
							next;
						}
					}
				}
					
				push @ret, $item;
			}
			
			$line=~/\=\s+(\S+)\,\s(\d+)\sbp$/;
			$seq1 = $1; $seq1len = $2;

			#initialize
			$seq2 = "";
			$seq2id = "";
			$seq2len = 0;
			$seq1blocks = [];
			$seq2blocks = [];
			$blockscore = [];
			$match = 0;
			$identity = 0;
			$strand = '+';
			$alignText = 0;
			$alignSeq1 = "";
			$alignSeq2 = "";
			$alignMatch = "";
		}
		elsif ($line=~/^seq2/)
		{
			$line=~/\=\s+(\S+)\s\((.*?)\)\,\s(\d+)\sbp$/;
			$seq2 = $1; $seq2id = $2; $seq2len = $3;
		}
		elsif ($line=~/^\(complement\)/)
		{
			$strand = '-';
		}
		elsif ($line=~/^(\d+)\-(\d+)\s+\((\d+)\-(\d+)\)\s+(\d+)\%/)
		{
			#		print join ("\t", $1, $2, $3, $4, $5), "\n";
			my %s1b = (start=>$1-1, end=>$2-1);	#convert to 0-based coordinates
			my %s2b = (start=>$3-1, end=>$4-1);
			
			push @$seq1blocks, \%s1b;
			push @$seq2blocks, \%s2b;
			push @$blockscore, $5;

			#my $c = @$seq1blocks;
			#print "block num = $c\n";
			
			$match += ($2 - $1 + 1);
			$identity += ($2 - $1 + 1) * $5;
		}
		elsif ($line =~/^\s*\d+\s/)
		{
			$alignText = 1;
			
			#do not want to chop space in the end
			$line = <$fin>;
			chop $line;
			$line =~/^\s*\d+\s(.*?)$/;
			$alignSeq1 .= $1;
			my $textLen = length ($1);
			
			$line = <$fin>;
			chop $line;
			$line =~/^\s*(\S.*?)$/;
			
			my $spaceLen = $textLen - length ($1);
			
			$alignMatch .= "" . (" "x $spaceLen) . $1;
		
			$line = <$fin>;
			chop $line;
			$line =~/^\s*\d+\s(.*?)$/;
			$alignSeq2 .= $1;
		}
		else
		{
			Carp::croak "$in has been disrupted\n$line\n";
		}
	}
	if (@$seq2blocks > 0)
	{
		$identity = int ($identity/100 + 0.5);
		Carp::croak "inconsistency found near file end\n"
		unless length ($alignSeq1) == length ($alignMatch) && length ($alignSeq2) == length ($alignMatch);
				
		my $item = {seqfile1=>$seq1, seqfile2=>$seq2, seq2id=>$seq2id, 
			seq1len=>$seq1len, seq2len=>$seq2len,
			match=> $match, identity=>$identity, strand=>$strand, 
			seq1blocks=>$seq1blocks,
			seq2blocks=>$seq2blocks,
			blockscore=>$blockscore};

		if (length ($alignSeq1) > 0)
		{
			my $align = _splitSim4AlignText ($alignSeq1, $alignMatch, $alignSeq2);

			my $nblock = @{$align->{"seq1text"}};
			$item->{"seq1text"} = $align->{"seq1text"};
			$item->{"matchtext"} = $align->{"matchtext"};
			$item->{"seq2text"} = $align->{"seq2text"};

			if ($nblock != @{$item->{"seq1blocks"}})
			{
				print "inconsistent block number for gene = $seq1, acc = $seq2id\n"; 
				print Dumper ($item), "\n";

				$seq2 = "";
				$seq2id = "";
				$seq2len = 0;
				$seq1blocks = [];
				$seq2blocks = [];
				$blockscore = [];
				$match = 0;
				$identity = 0;
				$strand = '+';
				$alignText = 0;
				$alignSeq1 = "";
				$alignSeq2 = "";
				$alignMatch = "";
	
				next;
			}

			for (my $i = 0; $i < $nblock; $i++)
			{
				my $b1Len = $seq1blocks->[$i]->{"end"} - $seq1blocks->[$i]->{"start"} + 1;
				my $b2Len = $seq2blocks->[$i]->{"end"} - $seq2blocks->[$i]->{"start"} + 1;
				my $b1Text = $align->{"seq1text"}->[$i];
				my $b2Text = $align->{"seq2text"}->[$i];

				$b1Text=~s/\s//g;
				$b2Text=~s/\s//g;
				if ($b1Len != length ($b1Text) || $b2Len != length ($b2Text))
				{
					print "inconsistency size of block $i: seq1 = $b1Len, seq1text = ", length ($b1Text), ", seq2=$b2Len, seq2text=", length($b2Text), "\n"; 
					print Dumper ($item), "\n";
					
					$seq2 = "";
					$seq2id = "";
					$seq2len = 0;
					$seq1blocks = [];
					$seq2blocks = [];
					$blockscore = [];
					$match = 0;
					$identity = 0;
					$strand = '+';
					$alignText = 0;
					$alignSeq1 = "";
					$alignSeq2 = "";
					$alignMatch = "";
	
					next;
				}
			}
		}
					
		push @ret, $item;
	}
			
	close ($fin);
	return \@ret;
}


sub printSim4Align
{
	my $align = $_[0];

	print "\n\n";
	print "seq1 = ", $align->{"seqfile1"}, ", ", $align->{"seq1len"}, " bp\n";
	print "seq2 = ", $align->{"seqfile2"}, " (", $align->{"seq2id"}, "), ", $align->{"seq2len"}, " bp\n";

	print "\n";

	my $seq1blocks = $align->{"seq1blocks"};
	my $seq2blocks = $align->{"seq2blocks"};
	
	my $nblocks = @$seq1blocks;
	for (my $i = 0; $i < $nblocks; $i++)
	{
		print	"", ($seq1blocks->[$i]->{"start"} + 1), "-", ($seq1blocks->[$i]->{"end"} + 1), "\t";
		print "(", ($seq2blocks->[$i]->{"start"} + 1), "-", ($seq2blocks->[$i]->{"end"} + 1), ")\t";
		print  $align->{"blockscore"}->[$i], "%", "\n";
	}
}

sub sim4ToBed
{
	my ($align, $gene) = @_;
	
	my @blockSizes;
	my @blockStarts;

	my $score = $align->{"identity"} / $align->{"match"};
	$score = sprintf ("%.2f", $score);
	
	my $geneBlocks = $align->{"seq1blocks"};	
	my $nblocks = @$geneBlocks;
	
	my $chromStart = $geneBlocks->[0]->{"start"};
	my $chromEnd = $geneBlocks->[$nblocks - 1]->{"end"};
		
	for (my $i = 0; $i < $nblocks; $i++)
	{
		$blockSizes[$i] = $geneBlocks->[$i]->{"end"} - $geneBlocks->[$i]->{"start"} + 1;
		$blockStarts[$i] = $geneBlocks->[$i]->{"start"} - $chromStart;
	}
		
	my $region = {chrom=>$gene,			#gene id
			chromStart=>$chromStart,	#start on gene contig
			chromEnd=>$chromEnd, 	#end on gene contig
			name=>$align->{"seq2id"},
			score=>$score, 
			strand=> '+', 		#always set to '+', the gene strand, although a few transcripts were aligned to the reverse strand
			thickStart=>$chromStart,
			thickEnd=>$chromEnd,
			itemRgb=>"0,0,0",
			blockCount=>$nblocks,
			blockSizes=>\@blockSizes,
			blockStarts=>\@blockStarts};
	return $region;
	#removeSmallGap ($region, 0);
}



#need to add code to handle alignment text
sub checkSim4AlignQuality
{
	my ($aligns, $gene, $coverage, $identity, $verbose) = @_;
	my @goodAligns;
	
	foreach my $a (@$aligns)
	{
		#print "align = ", Dumper ($a), "\n";

		my @seq1blocks = @{$a->{"seq1blocks"}};#genomic
		my @seq2blocks = @{$a->{"seq2blocks"}};#transcript
		my @blockscore = @{$a->{"blockscore"}};
		
		my $nblocks = @{$a->{"seq2blocks"}};
		if ($nblocks < 2)
		{
			print "gene=$gene vs transcript=", $a->{"seq2id"}, ": no splice junction spaned, failed\n" if $verbose;
			next;
		}
		
		#OK now we have at least two blocks, so there should be no error in the following
		#code to check the first and last blocks
	
		#remove the first block if less than 25 nt
		my $firstBlockSize = $seq1blocks[0]->{"end"} - $seq1blocks[0]->{"start"} + 1;
		
		#############################################################
		#TO BE FIXIED
		#use the while loop in the next version
		#
		#############################################################
		#
		while ($nblocks >= 2 && $firstBlockSize < 25)
		#if ($firstBlockSize < 25)
		{
			$a->{"match"} -= $firstBlockSize;
			$a->{"identity"} -= $firstBlockSize * $blockscore[0] / 100;
			shift @seq1blocks;
			shift @seq2blocks;
			shift @blockscore;
			$nblocks -= 1;
			#$firstBlockSize = $seq1blocks[0]->{"end"} - $seq1blocks[0]->{"start"} + 1;
		}
		
		#remove the last block if less than 25 nt 
		my $lastBlockSize = $seq1blocks[$nblocks-1]->{"end"} - $seq1blocks[$nblocks-1]->{"start"} + 1;
		

		#############################################################
		#TO BE FIXIED
		#use the while loop in the next version
		#
		#############################################################
		#
		
		while ($nblocks >= 2 && $lastBlockSize < 25)
		#if ($lastBlockSize < 25)
		{
			$a->{"match"} -= $lastBlockSize;
			$a->{"identity"} -= $lastBlockSize * $blockscore[$nblocks - 1] / 100;
			pop @seq1blocks;
			pop @seq2blocks;
			pop @blockscore;
			$nblocks -= 1;
			#$lastBlockSize = $seq1blocks[$nblocks-1]->{"end"} - $seq1blocks[$nblocks-1]->{"start"} + 1;
		}
		
		if ($nblocks < @{$a->{"seq2blocks"}})
		{
			$a->{"seq1blocks"} = \@seq1blocks;
			$a->{"seq2blocks"} = \@seq2blocks;
			$a->{"blockscore"} = \@blockscore;
		}
		
		#make sure we still have at least two blocks
		$nblocks = @{$a->{"seq2blocks"}};
		if ($nblocks < 2)
		{
			print "gene=$gene vs transcript=", $a->{"seq2id"}, ": no splice junction spaned, failed\n" if $verbose;
			next;
		}
		
		#check coverage and identity
		my $cov = $a->{"match"} / $a->{'seq2len'}; 	#percent coverage
		my $id = $a->{"identity"} / $a->{"match"};	#percent identity
		
		if ($cov < $coverage || $id < $identity)
		{
			print "gene=$gene vs transcript=", $a->{"seq2id"}, ": coverage=$cov, identity=$id, failed\n" if $verbose;
			next;
		}

		my $contineous = 1;
		my $tsblocks = $a->{"seq2blocks"};
		for (my $i = 0; $i < $nblocks - 1; $i++)
		{
			if ($tsblocks->[$i]->{"end"} + 1 != $tsblocks->[$i+1]->{"start"})
			{
				$contineous = 0;
				last;
			}
		}
		
		if ($contineous != 1)
		{
			print "gene=$gene vs transcript=", $a->{"seq2id"}, ": gaps in transcript, failed\n" if $verbose;
			next;
		}

		push @goodAligns, $a;
	}

	my $n = @goodAligns;
	my $ntotal = @$aligns;
	
	print "$n of $ntotal aligned transcripts passed quality control\n" if $verbose;
	
	return \@goodAligns;
}

=head2 readVulgarFile

my $regions = readVulgarFile ($inFile, $verbose);

=cut

sub readVulgarFile
{
	my ($inVulgarFile, $verbose) = @_;

	my @regions;
	my $fin;
	open ($fin, "<$inVulgarFile") || Carp::croak "cannot open file $inVulgarFile to read\n";
	my $i = 0;
	while (my $line = <$fin>)
	{
		chomp $line;
		next unless $line=~/^vulgar:/;

		print "$i ...\n" if $verbose && $i % 5000 == 0;
		$i++;

		my $ret = lineToVulgar ($line);
		push @regions, $ret;
	}
	close ($fin);
	return \@regions;

}


=head2 lineToVulgar

my $vulgar = lineToVulgar ($line);

=cut
sub line2vulgar
{
	Carp::croak "obsolete function, call lineToVulgar instead\n";
}

sub lineToVulgar
{
	my $line = $_[0];
	my ($vulgar, $queryName, $queryStart, $queryEnd, $queryStrand, 
			$targetName, $targetStart, $targetEnd, $targetStrand, $score, @cols) = split (/\s+/, $line);
	my $ncol = @cols;

	Carp::croak "match columns must be multiples of three: $line\n" unless $ncol % 3 == 0;
	
	my @blocks;
	while (@cols)
	{
		my $type = shift @cols;
		my $querySize = shift @cols;
		my $targetSize = shift @cols;

		push @blocks, {type=> $type, q=>$querySize, t=>$targetSize};
	}
	#print Dumper (\@blocks), "\n";

	return {queryName=>$queryName, 
			queryStart=>$queryStart, 
			queryEnd=>$queryEnd, 
			queryStrand=>$queryStrand, 
			targetName=>$targetName, 
			targetStart=>$targetStart, 
			targetEnd=>$targetEnd, 
			targetStrand=>$targetStrand, 
			score=>$score, 
			blocks=>\@blocks};
}

=head2 writeVulgarFile

writeVulgarFile ($regions, $outFile, $verbose);


=cut
sub writeVulgarFile
{
	my ($regions, $outFile, $verbose) = @_;

	my $fout;

	open ($fout, ">$outFile") || Carp::croak "can not  open file $outFile to write\n";
	
	my $i = 0;
	foreach my $r (@$regions)
	{
		print "$i...\n" if $verbose && $i % 5000 == 0;
		print $fout vulgar2line ($r), "\n";
	}
	close ($fout);
}


sub vulgar2line
{
	Carp::croak "obsolete function, call vulgarToLine instead\n";
}


sub vulgarToLine
{
	my $vulgar = $_[0];

	my @cols = qw(queryName queryStart queryEnd queryStrand targetName targetStart targetEnd targetStrand score blocks);
	my @values = map {$vulgar->{$_}} @cols;

	my $blocks = pop @values;

	my $line = join (" ", "vulgar:", @values);

	foreach my $b (@$blocks)
	{
		$line .= " " . join (" ", $b->{'type'}, $b->{'q'}, $b->{'t'});
	}
	return $line;
}

sub vulgarToBed
{
	my $vulgar = $_[0];

	my $queryName = $vulgar->{'queryName'};
	my $queryStart = $vulgar->{'queryStart'};
	my $queryEnd = $vulgar->{'queryEnd'};
	my $queryStrand = $vulgar->{'queryStrand'};
	my $targetName = $vulgar->{'targetName'};
	my $targetStart = $vulgar->{'targetStart'};
	my $targetEnd = $vulgar->{'targetEnd'};
	my $targetStrand = $vulgar->{'targetStrand'};
	my $score = $vulgar->{'score'};
	my @blocks = @{$vulgar->{'blocks'}};

	my $ncol = ($#blocks+1)*3;

	my @blockStarts;
	my @blockEnds;

	my $lastBlock = $blocks[$#blocks];

	if ($lastBlock->{'t'} != $lastBlock->{'q'})
	{
		#protein
		#$lastBlock->{'t'} += 3;
		#$lastBlock->{'q'} += 1;
	}

	my $chrom = $targetName;
	my $chromStart = $targetStart;
	my $chromEnd = $targetEnd;

	my $strand = $targetStrand;

	if ($strand eq '+')
	{
		$chromEnd -= 1;
	}
	else
	{
		$chromStart = $targetEnd;
		$chromEnd = $targetStart - 1;
		@blocks = reverse @blocks;
	}
	
	#my $blockIter = 0;
	my $blockStart = $chromStart;
	my $blockEnd = $chromStart - 1;

	for (my $i = 0; $i < @blocks; $i++)
	{
		my $b = $blocks[$i];

		if ($b->{"type"} eq 'M')
		{
			$blockEnd = $blockStart + $b->{"t"} - 1;

			my $blockStartPolish = $blockStart;
			my $blockEndPolish = $blockEnd;
			
			#consider split codon
			if ($i > 0 && $blocks[$i-1]->{'type'} eq 'S')
			{
				$blockStartPolish -= $blocks[$i-1]->{'t'};
			}

			if ($i < @blocks - 1 && $blocks[$i+1]->{'type'} eq 'S')
			{
				$blockEndPolish += $blocks[$i+1]->{'t'};
			}
			push @blockStarts, $blockStartPolish;
			push @blockEnds, $blockEndPolish;
		}
		$blockStart += $b->{"t"};
	}


	#add the last block
	#push @blockStarts, $blockStart;
	#push @blockEnds, $blockEnd;

	my @blockSizes;

	my $n = @blockStarts;
	for (my $i = 0; $i < @blockStarts; $i++)
	{
		$blockSizes[$i] = $blockEnds[$i] - $blockStarts[$i] + 1;
		$blockStarts[$i] -= $chromStart;
	}

	#$score /= (5 * ($queryEnd - $queryStart));

	return {chrom=>$chrom, 
			chromStart=>$chromStart, 
			chromEnd=>$chromEnd, 
			name=>$queryName, 
			score=>$score, 
			strand=>$strand,
			thickStart=>$chromStart,
			thickEnd=>$chromEnd,
			itemRgb=>"0,0,0",
			blockCount=>$n,
			blockSizes=>\@blockSizes,
			blockStarts=>\@blockStarts
	};
}



1;



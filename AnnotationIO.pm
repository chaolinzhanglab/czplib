package AnnotationIO;

=head1 NAME

AnnotationIO - read and write bio annotation files
subroutines starting with a hyphen should not be called outside

=cut

use strict;
use FileHandle;
use Data::Dumper;
use Carp ();

use Common;

=head2 readBedFile
BED file
=cut

#zero-based coordinates
#

my $debug = 0;

sub readBedFile
{
	my $in = $_[0];

	my $verbose = 0;

	if (@_ > 1)
	{
		$verbose = $_[1];
	}
	
	my $out = [];
	my $fd = new FileHandle;
	
	open ($fd, "<$in")||Carp::croak "can not open file $in to read\n";
	my $i = 0;
	
	while (my $line = <$fd>)
	{
		next if $line =~/^\s*$/;
		next if $line =~/^\#/;
		next if $line =~/^track name\=/;

		print "$i ...\n" if $verbose && int ($i / 10000) * 10000 - $i == 0;
		$i++;

		my $entry = line2bed ($line);
		#print join ("\t", $entry->{"chromStart"}, $entry->{"chromEnd"}), "\n";
		push @$out, $entry;
	}
	close ($fd);
	return $out;
}

sub line2bed
{
	my $line = $_[0];

	my @cols = split (/\s+/, $line);
	Carp::croak "less than three columns in line: $line\n" if @cols < 3;
	
	#	print join ("\t", @cols), "\n";
	my @colNames = qw (chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts);
	my $i;
	my $entry = {};
	for ($i = 0; $i < @colNames; $i++)
	{
		last if ($#cols < $i);
		$entry->{$colNames[$i]} = $cols[$i];
	}

	if (exists $entry->{"blockSizes"})
	{
		my $blockSizes = $entry->{"blockSizes"};
		my @bs = split (/\,/, $blockSizes);
		Carp::croak "in correct number of blocks at line $line\n" if $#bs != $entry->{"blockCount"} - 1;
		$entry->{"blockSizes"} = \@bs;
	}

	if (exists $entry->{"blockStarts"})
	{
		my $blockStarts = $entry->{"blockStarts"};
		my @bs = split (/\,/, $blockStarts);
		Carp::croak "in correct number of blocks at line $line\n" if $#bs != $entry->{"blockCount"} - 1;
		$entry->{"blockStarts"} = \@bs;
	}
		
	$entry->{"chromEnd"} -= 1;
	$entry->{"thickEnd"} -= 1 if (exists $entry->{"thickEnd"});
		
	Carp::croak "chromStart (" . $entry->{"chromStart"} . ")  > chromEnd (" . $entry->{"chromEnd"} . ")\n" 
	if ($entry->{"chromStart"} > $entry->{"chromEnd"});
		
		#print join ("\t", $entry->{"chromStart"}, $entry->{"chromEnd"}), "\n";
	return $entry;
}


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

		my $ret = line2vulgar ($line);
		push @regions, $ret;
	}
	close ($fin);
	return \@regions;

}


sub line2vulgar
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


sub gene2exon
{
    my $g = $_[0];
    my $nexon = $g->{"blockCount"};
    my @exons;
    for (my $i = 0; $i < $nexon; $i++)
    {
        my $exonStart = $g->{"chromStart"} + $g->{"blockStarts"}->[$i];
        my $exonEnd = $exonStart + $g->{"blockSizes"}->[$i] - 1;
        my $e = {
            chrom => $g->{"chrom"},
            chromStart=> $exonStart,
            chromEnd => $exonEnd,
            name => join (":", $g->{"name"}, $i),
            score => $g->{"score"},
            strand => $g->{"strand"}
        };
        push @exons, $e;
    }
    return \@exons;
}



sub readPslFile
{
	my ($in, $verbose) = @_;
	
	my $fin = new FileHandle;
	my @results;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	
	my $line = <$fin>;
	if ($line=~/^psLayout/)	# there is a header
	{
		while ($line = <$fin>)
		{
			chomp $line;
			last if ($line =~/^\-\-\-/);
		}
	}
	else
	{
		seek ($fin, 0, 0); # no header line, go to the very beginnning
	}

	#now ready to go
	
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
		print "$i ...\n" if $verbose && $i % 10000 == 0;		
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


#note Psl used absolute ccornidates
sub pslToBed
{
	my ($align, $gene) = @_;
	my $useQuery = 0;
	
	if (@_> 2)
	{
		$useQuery = $_[2];
	}
	
	my $chromStart = $align->{"tStart"};
	my $chromEnd = $align->{"tEnd"};
	my $name = $align->{"qName"};
	my @blockStarts;
    foreach my $s (@{$align->{"tStarts"}})
	{
		push @blockStarts, $s - $chromStart;
	}
			
	if ($useQuery)
	{
		$chromStart = $align->{"qStart"};
		$chromEnd = $align->{"qEnd"};
		$name = $align->{"tName"};
		@blockStarts = ();
		foreach my $s (@{$align->{"qStarts"}})
		{
			push @blockStarts, $s - $chromStart;
		}
	}
	
	my $chrom = $gene;
	my $score = Common::calcPslPercentIdentity ($align, 1);
	my $strand = '+';#always set to '+', the gene strand, although a few transcripts were aligned to the reverse strand
	my $blockCount = $align->{"blockCount"};
	my @blockSizes = @{$align->{"blockSizes"}};
	
	my $region = {chrom=>$gene,			#gene id
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
		blockStarts=>\@blockStarts
	};
	return $region;
}

sub printPslAlign
{
	my $align = $_[0];
	print join ("\t",
		$align->{"matches"},
		$align->{"misMatches"},
		$align->{"repMatches"},
		$align->{"nCount"},
		$align->{"qNumInsert"},
		$align->{"qBaseInsert"},
		$align->{"tNumInsert"},
		$align->{"tBaseInsert"},
		$align->{"strand"},
		$align->{"qName"},
		$align->{"qSize"},
		$align->{"qStart"},
		$align->{"qEnd"} + 1,
		$align->{"tName"}, #tName
		$align->{"tSize"},
		$align->{"tStart"},#tStart
		$align->{"tEnd"} + 1, #tEnd
		$align->{"blockCount"},
		join (",", @{$align->{"blockSizes"}}) . ",", #blockSizes
		join (",", @{$align->{"qStarts"}}) . ",", #qStarts
		join (",", @{$align->{"tStarts"}}) . ","),  #tStarts
		"\n";
}

sub writePslFile
{
	my ($aligns, $out) = @_;
	my $printHeader = 0;

	$printHeader = 1 if (@_ == 3 && $_[2] == 1);
	
	my $header = "psLayout version 3\n";
	$header .= "match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts\n";
	$header .= "     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count\n";
	$header .= "---------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
	my $fout = new FileHandle;

	open ($fout, ">$out") || Carp::croak "can not open file $out to write\n";
	print $fout $header if $printHeader;
	
	my $stdout = select ($fout);
	foreach my $align (@$aligns)
	{	
		printPslAlign ($align);
	}

	close ($fout);
	select ($stdout);
}


sub writeBedFile
{
	
	my ($regions, $header, $out, $append) = @_;
	my $fout = new FileHandle;
	if ($append && $append eq 'a')
	{
		open ($fout, ">>$out") || Carp::croak "can not open file $out to append\n";
	}
	else
	{
		open ($fout, ">$out") || Carp::croak "can not open file $out to write\n";
	}
	#my $stdout = select ($fout);
	#my @colNames = qw (chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts);
	
	if (@$regions <1)
	{
		#select ($stdout);
		close ($fout);
		return;
	}

	#my $colId;
	#for ($colId = @colNames - 1; $colId > 0; $colId--)
	#{
	#	last if (exists $regions->[0]->{$colNames[$colId]});
	#}
	
	#my $colNum = $colId + 1;	#total column number
	
	print $fout $header, "\n" if ($header =~/^track/);
	
	foreach my $r (@$regions)
	{
		#this is not a real copy
	
		print $fout printBedRegionToString ($r), "\n";
=annt
		my %rCopy = %$r;

		$rCopy{"chromEnd"} += 1;
		if (exists $rCopy{'thickEnd'})
		{
			$rCopy{'thickEnd'} += 1;
		}
		if (exists $rCopy{'blockCount'})
		{
			Carp::croak "no blockSizes\n" unless exists $rCopy {'blockSizes'};
			Carp::croak "no blockStarts\n" unless exists $rCopy {'blockStarts'};
			
			$rCopy{'blockSizes'} = join (",", @{$rCopy{'blockSizes'}});
			$rCopy{'blockStarts'} = join (",", @{$rCopy{'blockStarts'}});
		}
		
		#print STDOUT Dumper ($r), "\n";
		print join ("\t", $rCopy{"chrom"}, $rCopy{"chromStart"}, $rCopy{"chromEnd"});
		for (my $i = 3; $i < $colNum; $i++)
		{
			my $col = $colNames[$i];
			if (exists $rCopy{$col} && $rCopy{$col} ne '')
			{
				print "\t", $rCopy{$col};
			}
			else
			{
				print "\t", '.';
			}
		}
		print "\n";
=cut
	}

	close ($fout);
	#select ($stdout);
}


sub splitBedFileByChrom
{
	my ($inBedFile, $outDir, $verbose) = @_;

	Carp::croak "$inBedFile does not exist\n" unless -f $inBedFile;
	Carp::croak "$outDir does not exist\n" unless -d $outDir;

	my %tagCount;

	my $fin = new FileHandle;
	#print "reading tags from $inBedFile ...\n" if $verbose;
	open ($fin, "<$inBedFile") || Carp::croak "can not open file $inBedFile to read\n";

	#my %tagCount;
	#write tags according to chromosomes
	my $i = 0;
	my %fhHash;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		next if $line=~/^track name/;
	
		$i++;

		print "$i ...\n" if $i - int ($i / 50000) * 50000 == 0 && $verbose;
	
		$line =~/(\S+)\t/;

		my $chrom = $1;
		#print "chrom = $chrom\n";
		my $tmpFile = "$outDir/$chrom.bed";
	
		if (not exists $fhHash{$chrom})
		{
			my $fout = new FileHandle;
			open ($fout, ">>$tmpFile") || Carp::croak "can not open file $tmpFile to write\n";
			$fhHash{$chrom} = $fout;
		}
	
		#open ($fout, ">>$tmpFile") || Carp::croak "can not open file $tmpFile to write\n";
	
		my $fout = $fhHash{$chrom};
		print $fout $line, "\n";
		#close ($fout);

		if (exists $tagCount{$chrom})
		{
			$tagCount{$chrom}->{"n"} ++;
		}
		else
		{
			$tagCount{$chrom} = {f=>"$outDir/$chrom.bed", n=> 1};
		}
	}
	close ($fin);
	#close all file handles
	foreach my $chrom (keys %fhHash)
	{
		close ($fhHash{$chrom});
	}

	return \%tagCount;
}

#expend to 12 column format
sub bed2Full
{
	my $region = $_[0];
	Carp::croak "at least three columns should exist\n" 
	unless exists $region->{"chrom"} && exists $region->{"chromStart"} && exists $region->{"chromEnd"};
	
	$region->{"name"} = $region->{"chrom"} . ":" . $region->{"chromStart"} . "-" . ($region->{"chromEnd"} + 1)
	unless exists $region->{"name"};

	$region->{"score"} = 0 unless exists $region->{"score"};
	$region->{"strand"} = "+" unless exists $region->{"strand"};

	$region->{"thickStart"} = $region->{"chromStart"} unless exists $region->{"thickStart"};
	$region->{"thickEnd"} = $region->{"chromEnd"} unless exists $region->{"thickEnd"};
	$region->{"itemRgb"} = "0,0,0" unless exists $region->{"itemRgb"};

	$region->{"blockCount"} = 1 unless exists $region->{"blockCount"};
	my @blockSizes = ($region->{"chromEnd"} - $region->{"chromStart"} + 1);
	$region->{"blockSizes"} = \@blockSizes unless exists $region->{"blockSizes"};
	
	my @blockStarts = (0);
	$region->{"blockStarts"} = \@blockStarts unless exists $region->{"blockStarts"};

	return $region;
}

sub printBedRegion
{
	my $region = $_[0];
	my $str = printBedRegionToString ($region);
	print $str, "\n";
}

sub bed2line
{
	my $region = $_[0];
	return printBedRegionToString ($region);
}

#generate bed format into string
sub printBedRegionToString
{
	my $region = $_[0];
	
	my @colNames = qw (chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts);
	
	my $colNum = 12; #keys %$region;
	
	my %rCopy = %$region;
	$rCopy{"chromEnd"} += 1;
	if (exists $rCopy{'thickEnd'})
	{
		$rCopy{'thickEnd'} += 1;
	}
	
	if (exists $rCopy{'blockCount'})
	{
		Carp::croak "no blockSizes\n" unless exists $rCopy {'blockSizes'};
		Carp::croak "no blockStarts\n" unless exists $rCopy {'blockStarts'};
			
		$rCopy{'blockSizes'} = join (",", @{$rCopy{'blockSizes'}});
		$rCopy{'blockStarts'} = join (",", @{$rCopy{'blockStarts'}});
	}
	
	my $ret = join ("\t", $rCopy{"chrom"}, $rCopy{"chromStart"}, $rCopy{"chromEnd"});
	for (my $i = 3; $i < $colNum; $i++)
	{
		my $col = $colNames[$i];
		if (exists $rCopy{$col})
		{
			$ret .= "\t" . $rCopy{$col};
		}
		else
		{
			last;
			#Carp::croak "col=$col is not defined\n"; 
		}
	}
	return $ret;
}

sub readWigFile
{
	my $in = $_[0];
	my $fin = new FileHandle;
	my @ret;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		next if $line =~/^\#/;
		my ($chrom, $chromStart, $chromEnd, $score) = split (/\s/, $line);
		push @ret, {chrom=>$chrom, chromStart=>$chromStart, chromEnd=>$chromEnd-1, score=>$score};
	}
	close ($fin);
	return \@ret;
}

sub writeWigFile
{
	my ($regions, $header, $out) = @_;
	my $fout = new FileHandle;
	open ($fout, ">$out") || Carp::croak "can not open file $out to write\n";
	if ($header ne '')
	{
		print $fout $header, "\n";
	}

	foreach my $r (@$regions)
	{
		print $fout join ("\t", $r->{"chrom"},
				$r->{"chromStart"},
				$r->{"chromEnd"} + 1,
				$r->{"score"}), "\n";
	}
	close ($fout);
}

sub readSgrFile
{
	my $in = $_[0];
	my $fin = new FileHandle;
	my @ret;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	while (my $line =<$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		my ($chr, $pos, $score) = split ("\t", $line);
		push @ret, {chrom=>$chr, pos=>$pos, score=>$score};
	}
	close ($fin);
	return \@ret;
}

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

=head2 readMotifFile

Transfac motif file
references are ignored

readMotifFile (fileName, returnHash=0)

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

	my $fin = new FileHandle;
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

sub writeMotifFile
{
	my ($motifs, $out) = @_;

	my $fout = new FileHandle;
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


sub readMatchFile
{
	my $in = $_[0];

	my %siteHash;
	my $fin = new FileHandle;
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


#index PhastCons file
#the pointer indicate the position of the first line of each block

#->chromStart
#->chromEnd
#->chrom
#->pointer

sub indexBigPhastConsFile
{
	my $in = $_[0];
	my @ret;
	my $fin = new FileHandle;
	open ($fin, "<$in")|| Carp::croak "can not open file $fin to read\n";

	my $entry = 0;
	my $n = 0;
	while (my $line=<$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		if ($line=~/^fixedStep\schrom=(.*?)\sstart=(\d+)\sstep=(\d+)$/)
		{
			if ($entry)
			{
				$entry->{"chromEnd"} = $entry->{"chromStart"} + $n - 1;
				push @ret, $entry;
			}
			$entry = {chrom=>$1, chromStart=>$2, step=>$3, pointer=>tell($fin)};
			$n = 0;
		}
		else
		{
			next unless $entry;
			$n++;
		}
	}
	close ($fin);

	if ($entry)
	{
		$entry->{"chromEnd"} = $entry->{"chromStart"} + $n - 1;
		push @ret, $entry;
	}
	return {file=>Common::getFullPath($in), index=>\@ret};
}


sub indexBigFastaFile
{
	my $in = $_[0];
	my @ret;

	my $fin = new FileHandle;
	open ($fin, "<$in") || Carp::croak "can not open file $fin to read\n";

	#seek ($fin, 0, 0);
	my $currPos = tell ($fin);
	my $line = <$fin>;

	my $n = 0;
	while (1)
	{
		chomp $line;
		if ($line =~/^\s*$/)
		{
			$currPos = tell($fin);
			my $more = ($line =<$fin>);
			last unless $more;
			next;
		}
		
		if ($line =~/^\#/)
		{
			$currPos = tell($fin);
			my $more = ($line =<$fin>);
			last unless $more;
			next;
		}
		
		if ($line=~/^\>/)
		{
			$line = substr($line, 1);
			my @cols = split (/\s/, $line);
			my $id = shift @cols;
			my $entry = {id=>$id, pointer=>$currPos};
			push @ret, $entry;

			#if ($n - int($n / 500) * 500 == 0)
			#{
			#	print "$n...\n";
			#}
			#$n++;
		}
		$currPos = tell($fin);

		my $more = ($line =<$fin>);
		last unless $more;
	}
	close ($fin);
	return {file=>Common::getFullPath($in), index=>\@ret};
}

##---------------------------------------------------------
# readBigFastaFile
# in		: fasta file name
# seqInfo	: hash table {id=>id, pointer=>pointer}
#             when id is empty, we do not check consistency
#             sequence id
#
#----------------------------------------------------------
#
sub readBigFastaFile
{
	my ($in, $seqInfo) = @_;
	my $fin = new FileHandle;

	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	my $pointer = $seqInfo->{"pointer"};
	my $id = $seqInfo->{"id"};
	#print "readBigFastaFile: id=$id, pointer=$pointer\n";
	seek ($fin, $pointer, 0); #go to the point
	my $line = <$fin>;	#head line of the seq
	Carp::croak "can not find header line for seq $id at $pointer\n" unless $line=~/^\>/;
	my $id2 = substr($line, 1);
	my $desc = "";
	if ($id2 =~/^(\S+)\s+(\S.*?)$/)
	{
		$id2 = $1;
		$desc = $2;
	}
	
	if ($id ne '')
	{
		Carp::croak "the sequence id in fasta file ($id2) is not equal to the id ($id) in index\n" if $id ne $id2;
	}
	
	my $seq = "";
	while ($line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		next if $line =~/^\#/;
		
		last if ($line =~/^\>/);
		$seq .= $line;
	}	
	close ($fin);

	Carp::croak "the length of sequence $id is zero\n" if (length($seq) == 0);

	return {id=>$id, desc=>$desc, seq=>$seq};
}

##----------------------------------------------------------------
#
# Index phastcons file for efficient accessing
#
# Status: tested
# ----------------------------------------------------------------
# 
sub readPhastconsIndexFile
{
	my $in = $_[0];
	my $fin;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	my $line = <$fin>;
	chomp $line;
	my ($tmp, $file) = split ("=", $line);
	my @ret;
	while ($line =<$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		my ($chrom, $chromStart, $chromEnd, $step, $pointer) = split ("\t", $line);
		push @ret, {chrom=>$chrom, chromStart=>$chromStart, chromEnd=>$chromEnd, step=>$step, pointer=>$pointer};
	}
	close ($fin);
	return {file=>$file, index=>\@ret};
}

##----------------------------------------------------------------
# readBigPhastConsFile
#
# blockInfo{chrom=>, chromStart=>, chromEnd=>, step=>, pointer=>
# 	generated by readPhastconsIndexFile (above)
# Status: tested
# ----------------------------------------------------------------

sub readBigPhastConsFile
{
	my ($in, $blockInfo, $chromStart, $chromEnd) = @_;
	my $blockChromStart = $blockInfo->{"chromStart"};
	my $blockChromEnd = $blockInfo->{"chromEnd"};
	my $step = $blockInfo->{"step"};
	my @ret;
	
	return 0 if ($blockChromStart > $chromEnd || $chromStart > $blockChromEnd);	#no overlap
	
	my $start = ($blockChromStart >= $chromStart) ? $blockChromStart : $chromStart;
	my $end = ($blockChromEnd <= $chromEnd) ? $blockChromEnd : $chromEnd;

	my $fin = new FileHandle;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	seek ($fin, $blockInfo->{"pointer"}, 0);	#go to the first line of the block

	for (my $i = $blockChromStart; $i <= $blockChromEnd; $i+= $step)
	{
		my $line = <$fin>;
		chomp $line;
		if ($i >= $start && $i <= $end)
		{
			push @ret, $line;
		}
		elsif ($i > $end)
		{
			last;
		}
	}
	close ($fin);
	
	return {chrom=>$blockInfo->{"chrom"}, chromStart=>$start, chromEnd=>$end, step=>$step, scores=>\@ret};
}

#$seqs: reference to an array
#		each element should have {id=>id, desc=>desc, seq=>seq}
#
sub writeFastaFile
{
	my ($out, $seqs);
	
	my $fout = new FileHandle;
	open ($fout, ">$out") || Carp::croak "can not open file $out to write\n";
	foreach my $seq (@$seqs)
	{
		writeFastaSeq ($fout, $seq);
	}
	close ($fout);
}

##-------------------------------------------------------------
# write one fasta sequence
#
# Status: tested
# -------------------------------------------------------------
sub writeFastaSeq
{
	my ($fout, $seq) = @_;
	my $id = $seq->{"id"};
	$id .= "\t" . $seq->{"desc"} if (exists $seq->{"desc"} && length($seq->{"desc"}) > 0);
	print $fout ">$id\n";
	print $fout $seq->{"seq"}, "\n";
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

sub readsim4File
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
	
	my $fin = new FileHandle;
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


sub removeSmallGap
{
	my ($r, $gapSize) = @_;
	return $r unless exists $r->{"blockCount"} && $r->{"blockCount"} > 1;

	my @blockStarts;
	my @blockSizes;

	my $currBlockStart = 0; 
	my $currBlockSize = $r->{"blockSizes"}->[0];

	for (my $i = 1; $i < @{$r->{"blockStarts"}}; $i++)
	{
		my $currBlockEnd = $currBlockStart + $currBlockSize - 1;
		if ($r->{"blockStarts"}->[$i] > $currBlockEnd + $gapSize + 1) # a new block
		{
			push @blockStarts, $currBlockStart;
			push @blockSizes, $currBlockSize;
			$currBlockStart = $r->{"blockStarts"}->[$i];
			$currBlockSize = $r->{"blockSizes"}->[$i];
		}
		else
		{
			#extend the current block
			$currBlockSize = $r->{"blockStarts"}->[$i] + $r->{"blockSizes"}->[$i] - $currBlockStart;
		}
	}
	push @blockStarts, $currBlockStart;
	push @blockSizes, $currBlockSize;
	
	$r->{"blockCount"} = @blockStarts;
	$r->{"blockSizes"} = \@blockSizes;
	$r->{"blockStarts"} = \@blockStarts;

	return $r;
}



sub read1LQFile
{
	my $in = $_[0];
	my $fin = new FileHandle;

	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	while (my $line=<$fin>)
	{
		chomp $line;
		last if $line=~/^X\tY/;
	}

	my @ret;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		my @cols = split ("\t", $line);
		push @ret, {
			X=>$cols[0],
			Y=>$cols[1],
			SEQUENCE=>$cols[2],
			DESTYPE=>$cols[3],
			FEATURE=>$cols[4],
			QUALIFIER=>$cols[5],
			EXPOS=>$cols[6],
			PLEN=>$cols[7],
			POS=>$cols[8],
			CBASE=>$cols[9],
			PBASE=>$cols[10],
			TBASE=>$cols[11],
			IPBASE=>$cols[12],
			UNIT=>$cols[13],
			BLOCK=>$cols[14],
			ATOM=>$cols[15]
		};
	}
	close ($fin);
	return \@ret;
}

sub parseCelHeader
{
	my $fin = $_[0];
	my %ret;
	while (my $line =<$fin>)
	{
		chomp $line;
		last if $line=~/^\s*$/;
		my @cols = split ("=", $line);
		Carp::croak "incorrect name=value pair: $line\n" if @cols < 2;
		my $name = shift @cols;
		my $value = join ("=", @cols);
		$ret{$name} = $value;
	}
	return \%ret;
}

sub parseCelIntensity
{
	my $fin = $_[0];
	my $ret;
	my $line = <$fin>;
	my ($name, $numCells) = split ("=", $line);
	
	Carp::croak "NumberCells not found:$line\n" if $name ne 'NumberCells';
	$line = <$fin>;
	for (my $i = 0; $i < $numCells; $i++)
	{
		my $line = <$fin>;
		chomp $line;
		my @cols = split (/\s+/, $line);
		shift @cols if @cols > 5;
		Carp::croak "incorrect number of columns at line: $line\n" if @cols != 5;
		my ($x, $y, $mean, $stdev, $npixels) = @cols;
		$ret->[$x][$y] = {MEAN=>$mean, STDV=>$stdev};
	}
	return $ret;
}

sub parseCelMasksOrOutliers
{
	my $fin = $_[0];
	my @ret;
	my $line = <$fin>;
	my ($name, $numCells) = split ("=", $line);
	
	Carp::croak "NumberCells not found:$line\n" if $name ne 'NumberCells';
	$line = <$fin>;
	for (my $i = 0; $i < $numCells; $i++)
	{
		my $line = <$fin>;
		chomp $line;
		my @cols = split (/\s+/, $line);
		shift @cols if @cols > 2;
		Carp::croak "incorrect number of columns at line: $line\n" if @cols != 2;
		my ($x, $y) = @cols;
		push @ret, {X=>$x, Y=>$y};
	}
	return \@ret;
}

sub parseCelModified
{
	my $fin = $_[0];
	my @ret;
	my $line = <$fin>;
	my ($name, $numCells) = split ("=", $line);
	
	Carp::croak "NumberCells not found:$line\n" if $name ne 'NumberCells';
	$line = <$fin>;
	for (my $i = 0; $i < $numCells; $i++)
	{
		my $line = <$fin>;
		chomp $line;
		my @cols = split (/\s+/, $line);
		shift @cols if @cols > 3;
		Carp::croak "incorrect number of columns at line: $line\n" if @cols != 3;
		my ($x, $y, $origMean) = @cols;
		push @ret, {X=>$x, Y=>$y, ORIGMEAN=>$origMean};
	}
	return \@ret;
}


sub readCelFile
{
	my $in = $_[0];
	my %ret;
   	
	my $fin = new FileHandle;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	my $line = <$fin>;
	$line = <$fin>;
	chomp $line;
	my ($name, $value) = split ("=", $line);
	Carp::croak "only Cel file version 3 is recognized\n" if $value != 3;
	my $iter = 0;
	while ($line =<$fin>)
	{
		$iter++;
		#last if $iter > 10;
		chomp $line;
		next if $line=~/^\s*$/;
		#print $line, "\n";
		if ($line eq '[HEADER]')
		{
			#print "parsing header...\n";
			$ret{"header"} = parseCelHeader ($fin);
			Carp::croak "missing number of columns\n" unless exists $ret{"header"}->{"Cols"};
			Carp::croak "missing number of rows\n" unless exists $ret{"header"}->{"Rows"};
		}
		elsif ($line eq '[INTENSITY]')
		{
			#print "parsing intensity...\n";
			$ret{"intensity"} = parseCelIntensity ($fin);
			Carp::croak "incorrect number of columns\n" if @{$ret{"intensity"}} != $ret{"header"}->{"Cols"};
			Carp::croak "incorrect number of rows\n" if @{$ret{"intensity"}->[0]} != $ret{"header"}->{"Rows"};
		}
		elsif ($line eq '[MASKS]')
		{
			#print "parsing masks...\n";
			$ret{"masks"} = parseCelMasksOrOutliers ($fin);
		}
		elsif ($line eq '[OUTLIERS]')
		{
			$ret{"outliers"} = parseCelMasksOrOutliers ($fin);
		}
		elsif ($line eq '[MODIFIED]')
		{
			$ret{"modified"} = parseCelModified ($fin);
		}
	}
	close ($fin);
	return \%ret;	
}


sub readRNAplfoldFile
{
	my $inFile = $_[0];

	my $fin = new FileHandle;
	open ($fin, "<$inFile") || Carp::croak "can not open file $inFile to read\n";

	my @ret;
	my $seqLen = -1;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;

		if ($line =~/^\/sequence \{/)
		{
			my $seq = <$fin>;
			chomp $seq;
			$seqLen = length ($seq) - 1;

			#print "seqlen = $seqLen\n";
			for (my $i = 0; $i < $seqLen; $i++)
			{
				$ret[$i] = 0;
			}
		}
		elsif($line=~/^(\d+)\s+(\d+)\s+(\S+)\s+ubox$/)
		{
			my $i = $1 - 1;
			my $j = $2 - 1;
			my $p = $3 * $3;

			$ret[$i] += $p;
			$ret[$j] += $p;
		}
	}
	close ($fin);
	return \@ret;
}



sub readTreeFile
{
	my $in = $_[0];
	my $fin = new FileHandle;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";

	my $treeStr = "";
	
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		$line =~s/\s//g;

		$treeStr .= $line;
	}
	close ($fin);
	chop $treeStr if ($treeStr =~/\;$/);

	my $tokens = segmentToken ($treeStr);
	
	return parseTree ($tokens);
}

sub writeTreeFile
{
	my ($tree, $out) = @_;
	my $treeStr = Common::codeTree ($tree);
	my $fout = new FileHandle;
	open ($fout, ">$out") || Carp::croak "can not open file $out to  write\n";

	print $fout $treeStr, ";\n";

	close ($fout);
}

sub segmentToken
{
	my $treeStr = $_[0];
	
	##segment the tree text
	my @tokens; #= split (/\(|\,|\)/, $treeStr);
	
	my $tok = "";
	for (my $i = 0; $i < length ($treeStr); $i++)
	{
		my $c = substr($treeStr, $i, 1);
		if ($c eq '(')
		{
			push @tokens, $c;
		}
		elsif ($c eq ',')
		{
			push @tokens, $tok;
			push @tokens, $c;
			$tok = '';
		}
		elsif ($c eq ')')
		{
			push @tokens, $tok;
			push @tokens, $c;
			$tok = "";
		}
		else
		{
			$tok .= $c;
		}
	}
	if ($tok ne '')
	{
		push @tokens, $tok;
	}
	push @tokens, ";";
	return \@tokens;
}


sub parseTree
{
	my $tokens = $_[0];
	my @tokens = @$tokens;

	if (@tokens == 2)
	{
		return {iter => 0, id=>$tokens[0], blen=> 0};
	}
	
	my @nodes;
	
	my $root = {iter=> 0, id=>"", blen=> 0};
	push @nodes, $root;
	
	my $nodeIter = 1;
	my $currIter = $root->{"iter"};
	
	foreach my $tok (@tokens)
	{
	
		if ($tok eq '(')
		{
			my $leftChild = {iter=>$nodeIter++, parent=>$nodes[$currIter]};
			my $rightChild = {iter=>$nodeIter++, parent=>$nodes[$currIter]};
			push @nodes, $leftChild;
			push @nodes, $rightChild;
			
			$nodes[$currIter]->{"left"} = $leftChild;
			$nodes[$currIter]->{"right"} = $rightChild;

			print "tok = $tok, curr node = ", nodeInfo ($nodes[$currIter]), ", branch and go left child ", nodeInfo ($leftChild), "\n" if $debug;
			$currIter = $nodes[$currIter]->{"left"}->{"iter"};
		
		}
		elsif ($tok eq ',')
		{
			Carp::croak "no right sibling for node ", $nodes[$currIter]->{"iter"}, "\n" unless exists $nodes[$currIter]->{"parent"} && exists $nodes[$currIter]->{"parent"}->{"right"};
		
			print "tok = $tok, curr node = ", nodeInfo ($nodes[$currIter]), ", jump to right sibling ", nodeInfo ($nodes[$currIter]->{"parent"}->{"right"}), "\n" if $debug;
			$currIter = $nodes[$currIter]{"parent"}->{"right"}->{"iter"};
		}
		elsif ($tok eq ')')
		{
			#if ($currIter == 0) #done
			#{
			#	print "tok = $tok, curr node = ", nodeInfo ($nodes[$currIter]), ", done\n";
			#	return $root;
			#}
			Carp::croak "no parent for node ", $nodes[$currIter]->{"iter"}, "\n" unless exists $nodes[$currIter]->{"parent"};

			print "tok = $tok, curr node = ", nodeInfo ($nodes[$currIter]), ", jump to parent ", nodeInfo ($nodes[$currIter]->{"parent"}), "\n" if $debug;
			$currIter = $nodes[$currIter]->{"parent"}->{"iter"};
		}
		elsif ($tok eq ';')
		{
			print "tok = $tok, curr node = ", nodeInfo ($nodes[$currIter]), ", done\n" if $debug;
			#$root->{"nodes"} = \@nodes;
			return $root;
		}
		else
		{
			my ($id, $blen) = split (/\:/, $tok);
			$nodes[$currIter]->{"id"} = $id;
			$nodes[$currIter]->{"blen"} = $blen;

			print "tok = $tok, curr node = ", nodeInfo ($nodes[$currIter]), ", assign label $tok\n" if $debug;
		}
	}
}

#need more test
sub parseTreeRecursive
{
	my ($tokens, $currNode, $nodeIter) = @_;

	my @tokens = @$tokens;

	my $tok = shift @tokens;
	
	if ($tok eq '(')
	{
		my $leftChild = {iter=>$nodeIter++, parent=>$currNode};
		my $rightChild = {iter=>$nodeIter++, parent=>$currNode};
			
		$currNode->{"left"} = $leftChild;
		$currNode->{"right"} = $rightChild;

		print "curr node = ", nodeInfo ($currNode), ", branch and go left child ", nodeInfo ($leftChild), "\n" if $debug;
		
		parseTree (\@tokens, $leftChild, $nodeIter);
	}
	elsif ($tok eq ',')
	{
		Carp::croak "no right sibling for node ", $currNode->{"iter"}, "\n" unless exists $currNode->{"parent"} && exists $currNode->{"parent"}->{"right"};
		
		print "curr node = ", nodeInfo ($currNode), ", jump to right sibling ", nodeInfo ($currNode->{"parent"}->{"right"}), "\n" if $debug;
		parseTree (\@tokens, $currNode->{"parent"}->{"right"}, $nodeIter);
	}
	elsif ($tok eq ')')
	{
		if (@tokens == 0) #done
		{
		print "curr node = ", nodeInfo ($currNode), ", done\n";
		return;
		}
		Carp::croak "no parent for node ", $currNode->{"iter"}, "\n" unless exists $currNode->{"parent"};

		print "curr node = ", nodeInfo ($currNode), ", jump to parent ", nodeInfo ($currNode->{"parent"}), "\n" if $debug;
	    parseTree (\@tokens, $currNode->{"parent"}, $nodeIter);
	}
	else
	{
		my ($id, $blen) = split (/\:/, $tok);
		$currNode->{"id"} = $id;
		$currNode->{"blen"} = $blen;

		print "curr node = ", nodeInfo ($currNode), ", assign label $tok\n" if $debug;
		parseTree (\@tokens, $currNode, $nodeIter);
	}
}



sub nodeInfo
{
    my $node = $_[0];
    return $node->{"iter"} . "($node)";
}
		


1;



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
sub readBedFile
{
	my $in = $_[0];
	my $out = [];
	my $fd = new FileHandle;
	
	open ($fd, "<$in")||Carp::croak "can not open file $in to read\n";
	my $line;
	while ($line = <$fd>)
	{
		next if $line =~/^\s*$/;
		next if $line =~/^\#/;
		next if $line =~/^track name\=/;
		my @cols = split (/\s+/, $line);
		Carp::croak "less than three columns in $in\n" if @cols < 3;
		
		#	print join ("\t", @cols), "\n";
		my @colNames = qw (chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts);
		my $i;
		my $entry;
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
		push @$out, $entry;
	}
	close ($fd);
	return $out;
}


sub readPslFile
{
	my $in = $_[0];
	
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
	
	while ($line =<$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		
		my @cols = split (/\t/, $line);
		
		shift @cols if ($#cols == 21); #the first column is the bin id
		
		#print join ("|\n", @cols), "|\n";
		Carp::croak ("the number of columns is incorrect\n") if $#cols != 20;
		
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
	
	my ($regions, $header, $out) = @_;
	my $fout = new FileHandle;
	open ($fout, ">$out") || Carp::croak "can not open file $out to write\n";
	my $stdout = select ($fout);
	my @colNames = qw (chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts);
	
	if (@$regions <1)
	{
		select ($stdout);
		close ($fout);
		return;
	}

	my $colId;
	for ($colId = @colNames - 1; $colId > 0; $colId--)
	{
		last if (exists $regions->[0]->{$colNames[$colId]});
	}
	
	my $colNum = $colId + 1;	#total column number
	
	print $header, "\n" if ($header =~/^track/);
	
	foreach my $r (@$regions)
	{
		#this is not a real copy
		
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
	}
	close ($fout);
	select ($stdout);
}


sub printBedRegion
{
	my $region = $_[0];
	my $str = printBedRegionToString ($region);
	print $str, "\n";
}

#generate bed format into string
sub printBedRegionToString
{
	my $region = $_[0];
	
	my @colNames = qw (chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts);
	
	my $colNum = keys %$region;
	
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
			Carp::croak "col=$col is not defined\n"; 
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
sub readsim4File
{
	my $in = $_[0];
	my @ret;
	
	my ($seq1, $seq2, $seq2id, $strand) = ("", "", "", '+');
	my ($seq1len, $seq2len) = (0, 0);
	my $seq1blocks = [];
	my $seq2blocks = [];
	my $blockscore = [];
	
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
				push @ret, {seqfile1=>$seq1, seqfile2=>$seq2, seq2id=>$seq2id, 
						seq1len=>$seq1len, seq2len=>$seq2len,
						match=> $match, identity=>$identity, strand=>$strand, 
						seq1blocks=>$seq1blocks,
						seq2blocks=>$seq2blocks,
						blockscore=>$blockscore};
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
		else
		{
			Carp::croak "$in has been disrupted\n$line\n";
		}
	}
	if (@$seq2blocks > 0)
	{
		$identity = int ($identity/100 + 0.5);
		push @ret, {seqfile1=>$seq1, seqfile2=>$seq2, seq2id=>$seq2id, 
			seq1len=>$seq1len, seq2len=>$seq2len,
			match=> $match, identity=>$identity, strand=>$strand,
			seq1blocks=>$seq1blocks,
			seq2blocks=>$seq2blocks,
			blockscore=>$blockscore};
		
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
		print      ($seq1blocks->[$i]->{"start"} + 1) . "-" . ($seq1blocks->[$i]->{"end"} + 1) . "\t";
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


1;



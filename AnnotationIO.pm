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
		$entry->{"chromEnd"} -= 1;
		Carp::croak "chromStart (" . $entry->{"chromStart"} . ")  > chromEnd (" . $entry->{"chromEnd"} . ")\n" 
		if ($entry->{"chromStart"} > $entry->{"chromEnd"});
		
		#print join ("\t", $entry->{"chromStart"}, $entry->{"chromEnd"}), "\n";
		push @$out, $entry;
	}
	close ($fd);
	return $out;
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
		#print STDOUT Dumper ($r), "\n";
		print join ("\t", $r->{"chrom"}, $r->{"chromStart"}, $r->{"chromEnd"} + 1);
		for (my $i = 3; $i < $colNum; $i++)
		{
			my $col = $colNames[$i];
			if (exists $r->{$col} && $r->{$col} ne '')
			{
				print "\t", $r->{$col};
			}
			else
			{
				print "\t", '.';
			}
		}
		print "\n";
	}
	close ($fout);
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
		print join ("\t", 
				sprintf ("%4d", $pos->{"A"}), 
				sprintf ("%4d", $pos->{"C"}), 
				sprintf ("%4d", $pos->{"G"}), 
				sprintf ("%4d", $pos->{"T"})), "\n";
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
	return {file=>Common::getFullPath($in), index=>\@ret};
}

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

1;

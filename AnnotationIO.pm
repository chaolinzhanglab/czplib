package AnnotationIO;


use FileHandle;
use Data::Dumper;
use Carp ();
use strict;

sub readBedFile
{
	my $in = $_[0];
	my $out;
	my $fd = new FileHandle;
	open ($fd, "<$in")||Carp::croak "can not open file $in to read\n";
	my $line;
	while ($line = <$fd>)
	{
		next if $line =~/^\s*$/;
		next if $line =~/^track name\=/;
		my @cols = split (/\s+/, $line);
		my @colNames = qw (chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts);
		my $i;
		my $entry;
		for ($i = 0; $i < @colNames; $i++)
		{
			last if ($#cols < $i);
			$entry->{$colNames[$i]} = $cols[$i];
		}
		
		push @$out, $entry;
	}
	close ($fd);
	return $out;
}

1;

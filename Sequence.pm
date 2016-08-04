#
#===============================================================================
#
#         FILE:  Sequence.pm
#
#  DESCRIPTION:  Package to handle sequence data
#         BUGS:  ---
#        NOTES:  The package is spun off from the AnnotationIO.pm and Common.pm
#       AUTHOR:  Chaolin Zhang (cz), czhang@rockefeller.edu
#      COMPANY:  Rockefeller University
#      VERSION:  1.0
#      CREATED:  12/17/2010
#     REVISION:  ---
#===============================================================================


package Sequence;

use Bio::SeqIO;

use Common;
use Bed;

require Exporter;


our $VERSION = 1.02;


@ISA = qw (Exporter);

@EXPORT = qw (
	baseComp
	baseCount
	complement
	enumerateSeq
	generateRandomSeq
	indexBigFastaFile
	mutateSeq
	nibFrag
	readBigFastaFile
	revcom
	segmentStr
	wordcount
	seqEntropy
	writeFastaFile
	writeFastaSeq
	bedToSeq
	locateStopCodon
	isNMDTranscript
);


=head1 NAME

Sequence - subroutines to handle nucleotide sequences

subroutines starting with a hyphen should not be called outside

=cut

use strict;
use Data::Dumper;
use Carp;

use File::Temp qw(:mktemp);


=head2 indexBigFastaFile

=cut

sub indexBigFastaFile
{
	my $in = $_[0];
	my @ret;

	my $fin;
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

=head2 readBigFastaFile

my $fasta = readBigFastaFile ($inFile, $seqInfo);


inFile	: fasta file name
seqInfo	: hash table {id=>id, pointer=>pointer}
             when id is empty, we do not check consistency
             sequence id

=cut
sub readBigFastaFile
{
	my ($in, $seqInfo) = @_;
	my $fin;
	
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";
	my $pointer = $seqInfo->{"pointer"};
	my $id = $seqInfo->{"id"};
	#print "readBigFastaFile: id=$id, pointer=$pointer\n";
	seek ($fin, $pointer, 0); #go to the point
	my $line = <$fin>;	#head line of the seq
	chomp $line;
	
	if ($line!~/^\>/)
	{
		return 0;
		#Carp::croak "can not find header line for seq $id at $pointer\n";
	}

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


=head2 writeFastaFile

#$seqs: reference to an array
#		each element should have {id=>id, desc=>desc, seq=>seq}
#
=cut

sub writeFastaFile
{
	my ($out, $seqs);
	
	my $fout;
	open ($fout, ">$out") || Carp::croak "can not open file $out to write\n";
	foreach my $seq (@$seqs)
	{
		writeFastaSeq ($fout, $seq);
	}
	close ($fout);
}

=head2 writeFastaSeq

write one fasta sequence
 Status: tested

=cut
sub writeFastaSeq
{
	my ($fout, $seq) = @_;
	my $id = $seq->{"id"};
	$id .= "\t" . $seq->{"desc"} if (exists $seq->{"desc"} && length($seq->{"desc"}) > 0);
	print $fout ">$id\n";
	print $fout $seq->{"seq"}, "\n";
}


=head2 wordcount

count the number of all n-mer in one or more DNA/RNA sequences (case insensitive, i.e. sequeuence will be converted into upper case)

my $wordHash = wordcount ($seqStr, $wordSize)
my $wordHash = wordcount ($seqStr, $wordSize, {wordHash=>$oldWordHash)
my $wordHash = wordcount ($inFastaFile, $wordSize, {fasta=>1})

=cut

sub wordcount
{
	my ($seqStr, $wordSize, $opt) = @_;
	
	my $isFasta = exists $opt->{'fasta'} ? $opt->{'fasta'} : 0;
	my $wordHash = exists $opt->{'wordHash'} ? $opt->{'wordHash'} : {};
	my $noMask = exists $opt->{'noMask'} ? $opt->{'noMask'} : 0;
	
	if ($isFasta)
	{
		my $inFile = $seqStr;
		my $seqIO = Bio::SeqIO->new (-file =>$inFile, -format => 'Fasta');
		while (my $seq = $seqIO->next_seq())
		{
			$wordHash = wordcount ($seq->seq(), $wordSize, {wordHash=>$wordHash, noMask=>$noMask});
		}
    }
	else
	{
		my $seqLen = length ($seqStr);
		$seqStr = uc ($seqStr);
		for (my $i = 0; $i < length ($seqStr) - $wordSize + 1; $i++)
		{
			my $w = substr ($seqStr, $i, $wordSize);
			if ($noMask == 0)
			{
				next if $w=~/[^ACGTU]/;
			}
			$wordHash->{$w}++;
		}
	}
	return $wordHash;
}



#///////////////////////Sequence manipulation//////////////////////////
=head2 revcom

reverse complementary nucleotide sequence

my $rc = revcom ($seqStr); 

IUB codes, and different cases will be handled properly

=cut


sub revcom
{
	my $str = $_[0];
	return scalar reverse (complement($str));
	#return CORE::reverse (complement($str));
}


=head2 complement

=cut
sub complement
{
	my $str = $_[0];
	$str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	return $str;
}


sub enumerateSeq
{
	my $size = $_[0];
	my $nseq = 4 ** $size - 1;

	my @alphabets = ('A', 'C', 'G', 'T');
	my @sequences;

	for (my $i = 0; $i <= $nseq; $i++)
	{
		my @idx;
		my $x = $i;
		while (@idx < $size)
		{
			my $a = $x % 4;
			$x = int ($x / 4);
			push @idx, $a;
		}

		my $seq = join ("", reverse @alphabets[@idx]);
		push @sequences, $seq;
	}
	return \@sequences;
}

sub segmentStr
{
	my $str = $_[0];
	my $colNum = 60;
	$colNum = $_[1] if @_ > 1;
	

	#my ($str, $colNum) = @_;
	my $strNew = "";
	my $start = 0;

	while ($start < length ($str))
	{
		my $length = Common::min ($colNum, length($str) - $start);
		$strNew .= substr ($str, $start, $length) . "\n";
		#$str = substr ($str, $colNum);
		$start += $colNum;
	}
	chomp $strNew; #remove the last "\n"
	#$strNew .= $str;
	return $strNew;
}


=obsolete
#inefficient, should not be used
sub segmentStr2
{
	my $str = $_[0];
	my $colNum = 60;
	$colNum = $_[1] if @_ > 1;
	
	#my ($str, $colNum) = @_;
	my $strNew = "";
	while (length ($str) > $colNum)
	{
		$strNew .= substr ($str, 0, $colNum) . "\n";
		$str = substr ($str, $colNum);
	}
	$strNew .= $str;
	return $strNew;
}
=cut



=head2 mutationOffsetTable

=cut

my %mutationOffsetTable = (
		'A'=>['C', 'G', 'T'],
		'C'=>['G', 'T', 'A'],
		'G'=>['T', 'A', 'C'],
		'T'=>['A', 'C', 'G'],
		'N'=>['N', 'N', 'N'],
		'a'=>['c', 'g', 't'],
		'c'=>['g', 't', 'a'],
		'g'=>['t', 'a', 'c'],
		't'=>['a', 'c', 'g'],
		'n'=>['n', 'n', 'n']);


=head2 mutateSeq

my $mutant = mutateSeq ($seqStr, $numMutation);

mutate $numMutation nucleotides in the sequence $seqStr

=cut

sub mutateSeq
{
	my ($seqStr, $numMutation) = @_;
	
	$numMutation = max (0, $numMutation); #can not be negative
	$numMutation = min (length($seqStr), $numMutation); #can not be more than the length of the sequence

	if ($numMutation > 0)
	{
		my $mutPositions = sampleSeq (0, length($seqStr), $numMutation);
		foreach my $pos (@$mutPositions)
		{
			my $offset = int (rand(3));
			my $base = substr ($seqStr, $pos, 1);
			my $mutant = exists $mutationOffsetTable{$base} ? $mutationOffsetTable {$base}->[$offset] : 'N';
			
				#print "seqStr=$seqStr, len=", $seq->length(), ", pos=$pos, base=$base, offset=$offset, mutant=$mutant\n";
			substr($seqStr, $pos, 1) = $mutant;
		}
	}
	return $seqStr;
}


=head2 generateRandomSeq

generate a random nucleotide sequence

my $seqStr = generateRandomSeq ($seqLen)

=cut
sub generateRandomSeq
{
	my $seqLen = $_[0];
	my @ret;
	my @alphabet = qw(A C G T);

	foreach my $i (0 .. ($seqLen -1))
	{
		my $idx = int (rand (4));
		push @ret , $alphabet[$idx];
	}
	return join("", @ret);
}



sub nibFrag
{
	my ($nibFrag, $chromNib, $chromStart, $chromEnd, $strand, $name, $cacheDir) = @_;

	Carp::croak "$cacheDir does not exists\n" unless -d $cacheDir;

	$strand = 'm' if $strand eq '-';

	my $tmpFile = mktemp ("XXXXXX");

	my $nibFragTmp = "$cacheDir/$tmpFile.fa";
	my $cmd = "$nibFrag -masked -name=\"$name\" $chromNib $chromStart $chromEnd $strand $nibFragTmp";
	system ($cmd);
	my $seqIO = Bio::SeqIO->new (-file=>$nibFragTmp, format=>'fasta');
	my $seq=$seqIO->next_seq();
	unlink $nibFragTmp;
	
	return $seq;
}


#base composition
sub baseComp
{
	my $seqs = $_[0];
	my ($a, $c, $g, $t) = (0, 0, 0, 0);
	foreach my $seqStr (@$seqs)
	{
		my $a2 = ($seqStr=~tr/aA//);
		my $c2 = ($seqStr=~tr/cC//);
		my $g2 = ($seqStr=~tr/gG//);
		my $t2 = ($seqStr=~tr/tTuU//);
		$a += $a2;
		$c += $c2;
		$g += $g2;
		$t += $t2;
	}
	my $n = $a + $c +$g + $t;
	$a /= $n;
	$c /= $n;
	$g /= $n;
	$t /= $n;
	return {A=>$a, C=>$c, G=>$g, T=>$t, N=>$n};
}

#base composition
sub baseCount
{
	my $seqs = $_[0];
	my ($a, $c, $g, $t) = (0, 0, 0, 0);
	foreach my $seqStr (@$seqs)
	{
		my $a2 = ($seqStr=~tr/aA//);
		my $c2 = ($seqStr=~tr/cC//);
		my $g2 = ($seqStr=~tr/gG//);
		my $t2 = ($seqStr=~tr/tTuU//);
		$a += $a2;
		$c += $c2;
		$g += $g2;
		$t += $t2;
	}
	my $n = $a + $c +$g + $t;
	return {A=>$a, C=>$c, G=>$g, T=>$t, N=>$n};
}


sub maskSeqInRegions
{
	my ($seqStr, $regions) = @_;
	#print "seq = $seqStr\n";
	print "length of seq = ", length ($seqStr), "\n";
	foreach my $r (@$regions)
	{
		my $start = $r->{"chromStart"};
		my $end = $r->{"chromEnd"};

		#print "start = $start, end = $end\n";

		$start = 0 if $start < 0;
		$end = length ($seqStr) - 1 if $end >= length ($seqStr);

		my $seg = substr ($seqStr, $start, $end - $start + 1);
		$seg =~tr/a-z/n/;
		$seg =~tr/A-Z/N/;

		substr ($seqStr, $start, $end - $start + 1) = $seg;
	}
	return $seqStr;	
}


sub seqEntropy
{
	my ($seqStr, $wordSize) = @_;

	my $wordCount = wordcount ($seqStr, $wordSize);
	my @freq = values %$wordCount;
	my $s = sum (\@freq);
	@freq = map {$_/$s} @freq;
	
	my $entropy = entropy (\@freq);
	return $entropy;
}



=head2 bedToSeq

my $seqStr = bedToSeq ($chromSeq, $bed)

$chromSeq = {id=>, seq=>};

=cut
sub bedToSeq
{
	my ($chromSeq, $bed) = @_;
	bedToFull ($bed) unless exists $bed->{'blockStarts'};
	my $name = $bed->{'name'};

	my $seqStr = "";
	foreach (my $i = 0; $i < $bed->{'blockCount'}; $i++)
	{
		my $blockStart = $bed->{'chromStart'} + $bed->{'blockStarts'}->[$i];
		my $blockSeqStr = substr ($chromSeq->{"seq"}, $blockStart, $bed->{'blockSizes'}->[$i]);
		if (length ($blockSeqStr) != $bed->{'blockSizes'}->[$i])
		{
			warn "incorrect size for transcript $name, block $i\n";
			return "";
		}
		$seqStr .= $blockSeqStr;
	}

	return $bed->{'strand'} eq '+' ? $seqStr : revcom ($seqStr);
}


=head2 locateStopCodon

find the stop codon
return:
 the last position of the stop codon, -1 if not found

=cut
sub locateStopCodon
{
	my ($seqStr, $start, $checkStartCodon) = @_;
	$checkStartCodon = 0 unless $checkStartCodon;

	$seqStr = uc($seqStr);
	
	if ($checkStartCodon)
	{
		return -1 unless substr ($seqStr, $start, 3) eq 'ATG'; 
	}

	my $pos = $start;
	my $stop = -1;
	while (my $codon = substr ($seqStr, $pos, 3))
	{
		last if length($codon) < 3;
		if ($codon eq 'TAA' || $codon eq 'TAG' || $codon eq 'TGA')
		{
			$stop = $pos + 2;
			last;
		}
		$pos += 3;
	}
	return $stop;
}

=head2 isNMDTranscript

determine if a transcript is an NMD transcript

my $status = isNMDTranscript ($tsSeqStr, $tsBed);

$tsSeqStr: mRNA sequences
$tsBed: bed of transcript, start codon has to be determined in advance

return: 
	1: nmd
	0: not nmd
	-1: cannot be determined, because,e.g., the transcript is in complete and no stop codon is found
	
	it will also change score, thickStart, thickEnd of $tsBed

=cut

sub isNMDTranscript
{
	my ($tsSeqStr, $tsBed) = @_;
	bedToFull ($tsBed) unless exists $tsBed->{'blockStarts'};
	
	#get ORF and NMD of the transcript
	my $startCodon = $tsBed->{'strand'} eq '+' ? $tsBed->{'thickStart'} : $tsBed->{'thickEnd'};
	my $startCodonOnTs = genomeToTranscript ($tsBed, $startCodon);
	Carp::croak "failed to map start codon on to transcript:", Dumper ($tsBed), "\n" if $startCodonOnTs < 0;

	#look for stop codon
	my $checkStartCodon = 1;
	my $stopCodonOnTs = locateStopCodon ($tsSeqStr, $startCodonOnTs, $checkStartCodon);

	#stop codon not found
	if ($stopCodonOnTs < 0)
	{
		$tsBed->{'score'} = -1;
		return -1;
	}	

	my $stopCodon = transcriptToGenome ($tsBed, $stopCodonOnTs);
	$tsBed->{'strand'} eq '+' ? $tsBed->{'thickEnd'} = $stopCodon : $tsBed->{'thickStart'} = $stopCodon;

	#ORF OK but no splice, so no NMD
	if ($tsBed->{'blockCount'} <= 2)
	{
		$tsBed->{'score'} = 0; #
		return 0;
	}

	#ORF found and there are multiple exons
	
	#determine the start of the last exon junction
	my $lastExonJunctionStart = $tsBed->{'strand'} eq '+' ? 
		$tsBed->{'chromStart'} + $tsBed->{'blockStarts'}->[$tsBed->{'blockCount'}-2] + $tsBed->{'blockSizes'}->[$tsBed->{'blockCount'}-2] - 1 : 
		$tsBed->{'chromStart'} + $tsBed->{'blockStarts'}->[1];
	
	my $lastExonJunctionStartOnTs = genomeToTranscript ($tsBed, $lastExonJunctionStart);
	
	$tsBed->{'score'} = $lastExonJunctionStartOnTs - $stopCodonOnTs >= 50 ? 1 : 0;
	return $tsBed->{'score'};
}




1;



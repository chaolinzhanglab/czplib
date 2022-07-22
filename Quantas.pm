#
#===============================================================================
#
#         FILE:  Quantas.pm
#
#  DESCRIPTION:  Package for Quantas
#         BUGS:  ---
#        NOTES:  
#       AUTHOR:  Chaolin Zhang (cz), cz2294@columbia.edu
#      COMPANY:  Columbia University
#      VERSION:  1.0
#      CREATED:  11/25/2014
#     REVISION:  ---
#===============================================================================

package Quantas;


require Exporter;

our $VERSION = 1.01;

@ISA = qw (Exporter);

@EXPORT = qw (
	readConfigFile
	readASConfigFile
	readExprConfigFile
	readExprDataFile
	readASDataFile
	readBed6DataFile
);



=head1 NAME

QuantAS - subroutines that deal with QuantAS files

subroutines starting with a hyphen should not be called outside

=cut

use strict;
use warnings;
use Data::Dumper;
use Carp;


use Common;


=head2
readConfigFile - obsolete now

=cut

sub readConfigFile
{
    my ($configFile, $base, $type, $suffix) = @_;
    my $fin;
    open ($fin, "<$configFile") || Carp::croak "cannot open file $configFile to read\n";
    my $i = 0;
    my %groups;

    while (my $line = <$fin>)
    {
        chomp $line;
        next if $line=~/^\s*$/;
        next if $line=~/^\#/;
        my ($sampleName, $groupName) = split (/\t/, $line);
		$sampleName .= $suffix if $suffix;

        $groups{$groupName}->{"id"} = $i++ unless exists $groups{$groupName};
        push @{$groups{$groupName}->{"samples"}}, $sampleName;

		#check whether input file exists
        my $inputFile = $base ne '' ? "$base/$sampleName" : $sampleName;
        if (-d $inputFile)
        {
			Carp::croak "undefined type?\n" unless $type;
            $inputFile = "$inputFile/$type.count.txt";
        }
        Carp::croak "Input file $inputFile does not exist\n" unless -f $inputFile;
    }
    close ($fin);
    return \%groups;
}

sub readASConfigFile
{
    my ($configFile, $base, $type) = @_;
    my $fin;
    open ($fin, "<$configFile") || Carp::croak "cannot open file $configFile to read\n";
    my $i = 0;
    my %groups;

    while (my $line = <$fin>)
    {
        chomp $line;
        next if $line=~/^\s*$/;
        next if $line=~/^\#/;
        my ($sampleName, $groupName) = split (/\t/, $line);
        $groups{$groupName}->{"id"} = $i++ unless exists $groups{$groupName};
        push @{$groups{$groupName}->{"samples"}}, $sampleName;

		#check whether input file exists
        my $inputFile = $base ne '' ? "$base/$sampleName" : $sampleName;
        if (-d $inputFile)
        {
			Carp::croak "undefined type?\n" unless $type;
            $inputFile = "$inputFile/$type.count.txt";
        }
        Carp::croak "Input file $inputFile does not exist\n" unless -f $inputFile;
    }
    close ($fin);
    return \%groups;
}

sub readExprConfigFile
{
    my ($configFile, $base, $suffix) = @_;
    my $fin;
    open ($fin, "<$configFile") || Carp::croak "cannot open file $configFile to read\n";
    my $i = 0;
    my %groups;

    while (my $line = <$fin>)
    {
        chomp $line;
        next if $line=~/^\s*$/;
        next if $line=~/^\#/;
        my ($sampleName, $groupName) = split (/\t/, $line);

		$sampleName .= $suffix if $suffix;

        $groups{$groupName}->{"id"} = $i++ unless exists $groups{$groupName};
        push @{$groups{$groupName}->{"samples"}}, $sampleName;

		#check whether input file exists
        my $inputFile = $base ne '' ? "$base/$sampleName" : $sampleName;
        Carp::croak "Input file $inputFile does not exist\n" unless -f $inputFile;
    }
    close ($fin);
    return \%groups;
}


=head2 readExprDataFile

=cut

sub readExprDataFile
{
    my ($inputFile, $pseudoCount) = @_;

	$pseudoCount = 1 unless defined $pseudoCount;

    my $fin;
    my @data;
    my @geneInfo;
    open ($fin, "<$inputFile") || Carp::croak "cannot open file $inputFile to read\n";

    my $totalTagNum = 0;
    while (my $line = <$fin>)
    {
        chomp $line;
        next if $line =~/^\s*$/;
        next if $line =~/^\#/;

        my ($geneId, $symbol, $tagNum, $exonLen, $RPKM) = split (/\t/, $line);

        push @geneInfo, [$geneId, $symbol, $exonLen];
        push @data, [$tagNum,$RPKM];

        $totalTagNum += $tagNum;
    }
    close ($fin);

    #recalculate RPKM by adding pseudo count
    for (my $i = 0; $i < @geneInfo; $i++)
    {
        my $exonLen = $geneInfo[$i][2];
        my $tagNum = $data[$i][0];

        $tagNum = $pseudoCount if $tagNum < $pseudoCount;
        $data[$i][1] = $tagNum * 1e9 / $exonLen / $totalTagNum;
    }
	
	#geneInfo should be removed in the future
    return {info=>\@geneInfo, geneInfo=>\@geneInfo, data=>\@data};
}


sub readASDataFile
{
    my ($inputFile, $type) = @_;

    my $fin;
    my @data;
    my @ASInfo;
    open ($fin, "<$inputFile") || Carp::croak "cannot open file $inputFile to read\n";
    while (my $line = <$fin>)
    {
        chomp $line;
        next if $line =~/^\s*$/;
        next if $line =~/^\#/;

        my @cols = split (/\t/, $line);
        my (@infoCols, @dataCols);

        if ($type eq 'cass' || $type eq 'iret' || $type eq 'mutx' || $type eq 'taca' | $type eq 'ss')
        {
            @infoCols = @cols[0..7];
            @dataCols = @cols[8..$#cols];
            pop @dataCols if $type eq 'taca' && @dataCols > 5;
            #to exclude the last column of taca in the new format
            #04/17/2014
        }
        elsif ($type eq 'alt3' || $type eq 'alt5')
        {
            @infoCols = @cols[0..7];
            push @infoCols, $cols[10];
            @dataCols = @cols[8..9];
            push @dataCols, @cols[11..$#cols];
        }
        elsif ($type eq 'apat')
        {
            @infoCols = @cols[0..7];
            @dataCols = @cols[8..9];

            push @infoCols, @cols[10..11];
        }
        elsif ($type eq 'apa')
        {
            #polyA seq data
            @infoCols = @cols[0..7];
            push @infoCols, @cols[10..$#cols];
            @dataCols = @cols[8..9];
        }
        elsif ($type eq 'snv')
        {
            @infoCols = @cols[0..7];
            push @infoCols, @cols[10..11];
            @dataCols = @cols[8..9];
            push @dataCols, @cols[12..$#cols];
        }
        else
        {
            Carp::croak "incorrect AS type: $type\n";
        }

        push @ASInfo, \@infoCols;
        push @data, \@dataCols;
    }
    close ($fin);

	#the ASInfo should be removed
    return {info=>\@ASInfo, ASInfo=>\@ASInfo, data=>\@data};
}

sub readBed6DataFile
{
    my ($inputFile) = @_;

    my $fin;
    my @data;
    my @info;
    open ($fin, "<$inputFile") || Carp::croak "cannot open file $inputFile to read\n";

    my $totalTagNum = 0;
    while (my $line = <$fin>)
    {
        chomp $line;
        next if $line =~/^\s*$/;
        next if $line =~/^\#/;

		my @cols = split (/\t/, $line);
		Carp::croak "require 6 columns in bed\n" unless @cols == 6;

		my $score = $cols[4];
        push @info, [@cols[0..3], $cols[5]];
        push @data, $score;
    }
    close ($fin);

    return {info=>\@info, data=>\@data};
}

1;



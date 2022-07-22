#!/usr/bin/perl -w;
#
package PXSXconfig;
use strict;
use Carp;

my $ESEhome = "/home/zhang/zhangc/src/ESE3/PXSX";
sub getDB
{
	#print "getDB called ...\n";
	my $dbName = $_[0];
	open (FD, "<$ESEhome/DB/$dbName.dat") || Carp::croak "can not open file $dbName.dat to read\n";
	Carp::croak ("unknown matrix database name $@\n") if ($@);
	my $db;
	$db->{"desc"} = <FD>;
	chomp $db->{"desc"};
	my $line;
	while ($line = <FD>)
	{
		chomp $line;
		next if $line=~/^$/;
		$line=~tr/a-z/A-Z/;
		$db->{"word"}->{$line} = 1;
		$db->{"size"} = length ($line);
	}
	close (FD);
	return $db;	
}
sub getDBName
{
	my @dbArray;
	my $dbName;
	opendir (DB, "$ESEhome/DB") || die "could not open DB dir\n";
	my @dbNames = readdir (DB); 
	close (DB);
	foreach $dbName(@dbNames)
	{
		next unless $dbName =~/^(.*?)\.dat$/;
		$dbName = $1;
		
		my $dbDesc = `head -n 1 $ESEhome/DB/$dbName.dat`;
		chomp $dbDesc;
		#print "name = $dbName, desc = $dbDesc\n";
		push @dbArray, {name=>$dbName, desc=>$dbDesc};
		
	}
	return \@dbArray;
}

sub getVersion {"1.00";}
sub getReleaseTime {"11/06/2005";}
sub getURL {"http://rulai.cshl.edu/tools/ESE"};
sub getRef {"Copyright (c) 2005 - 2006, Chaolin Zhang";};
sub getAuthor {"PXEX " . getVersion () . " implemented by Chaolin Zhang (zhangc\@cshl.edu)";}
		
1;



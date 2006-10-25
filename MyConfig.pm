# create on Nov 22, 2004
# by Chaolin Zhang
#
# Common subroutines

package MyConfig;
use strict;


my $projHome = "/mnt/raygun";

sub getDefaultCache
{
	my $base = "cache";
   	$base = $_[0] if (@_ > 0);
	
	my $cacheHome = ".";
	$cacheHome = $ENV{'CACHEHOME'} if exists $ENV{'CACHEHOME'};
	my $cache = "$cacheHome/$base" . "_" . time();
	return $cache;
}

sub getGenomeDir
{
	my $species =$_[0];
	return "$projHome/data/mm5"	if ($species eq 'mm5');
	return "$projHome/data/mm6" if ($species eq 'mm6');
	return "$projHome/data/hg17" if ($species eq 'hg17');
	return "$projHome/data/hg18" if ($species eq 'hg18');
}

1;



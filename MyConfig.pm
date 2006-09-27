# create on Nov 22, 2004
# by Chaolin Zhang
#
# Common subroutines

package MyConfig;
use strict;


sub getDefaultCache
{
	my $base = "cache";
   	$base = $_[0] if (@_ > 0);
	
	my $cacheHome = ".";
	$cacheHome = $ENV{'CACHEHOME'} if exists $ENV{'CACHEHOME'};
	my $cache = "$cacheHome/$base" . "_" . time();
	return $cache;
}

1;



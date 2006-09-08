# $Id: NoVersion.pm,v 1.1 2006/09/08 01:55:38 zhangc Exp $

=head1 NAME

WWW::Search::Null::NoVersion - class for testing WWW::Search

=head1 SYNOPSIS

  use WWW::Search;
  my $oSearch = new WWW::Search('Null::NoVersion');

=head1 DESCRIPTION

This class is a specialization of WWW::Search that has no $VERSION.

This module is for testing the WWW::Search module.

=head1 AUTHOR

Martin Thurn <mthurn@cpan.org>

=cut

package WWW::Search::Null::NoVersion;

use strict;

use vars qw( @ISA );
@ISA = qw( WWW::Search );

1;


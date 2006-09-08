

=head1 NAME

WWW::Search::Null::Count - class for testing WWW::Search clients

=head1 SYNOPSIS

=for example begin 

  require WWW::Search;
  my $iCount = 4;
  my $oSearch = new WWW::Search('Null::Count',
                                '_null_count' => $iCount,
                               );
  $oSearch->native_query('Makes no difference what you search for...');
  my @aoResults = $oSearch->results;
  # ...You get $iCount results.

=for example end 

=for example_testing
is(scalar(@aoResults), $iCount, 'got the right number of results');
is($oSearch->approximate_result_count, $iCount, 'got the right approx_results');
my $oResult = shift @aoResults;
is($oResult->url, "url1", 'url');
is(scalar(@{$oResult->related_urls}), $iCount, 'got N related_urls');
is(scalar(@{$oResult->related_titles}), $iCount, 'got N related_titles');
is(scalar(@{$oResult->urls}), $iCount+1, 'got N+1 urls');
my $raURL = $oResult->urls;
# diag("sURL =$sURL=");
is(scalar(@{$raURL}), $iCount+1, 'got N+1 urls in arrayref');
$oSearch = new WWW::Search('Null::Count');
@aoResults = $oSearch->results;

=head1 DESCRIPTION

This class is a specialization of WWW::Search that returns some hits,
but no error message.  The number of hits returned can be controlled
by adding a '_null_count' hash entry onto the call to
WWW::Search::new().  The default is 5.

This module might be useful for testing a client program without
actually being connected to any particular search engine.

=head1 AUTHOR

Martin Thurn <mthurn@cpan.org>

=cut

package WWW::Search::Null::Count;

use WWW::Search::Result;
use strict;

use vars qw( @ISA $VERSION );
@ISA = qw( WWW::Search );
$VERSION = do { my @r = (q$Revision: 1.1 $ =~ /\d+/g); sprintf "%d."."%03d" x $#r, @r };

sub _native_setup_search
  {
  my ($self, $native_query, $native_opt) = @_;
  # print STDERR " + ::Null::Count::_native_setup_search()\n";
  if (! defined $self->{_null_count})
    {
    # print STDERR " +   setting default _null_count to 5\n";
    $self->{_null_count} = 5;
    } # if
  } # _native_setup_search


sub _native_retrieve_some
  {
  my $self = shift;
  # print STDERR " + ::Null::Count::native_retrieve_some()\n";
  my $response = new HTTP::Response(200,
                                    "This is a test of WWW::Search");
  $self->{response} = $response;
  my $iCount = $self->{_null_count};
  # print STDERR " +   iCount is $iCount\n";
  $self->_elem('approx_count', $iCount);
  for my $i (1..$iCount)
    {
    my $oResult = new WWW::Search::Result;
    $oResult->url(qq{url$i});
    $oResult->title(qq{title$i});
    $oResult->description("description$i");
    $oResult->change_date("yesterday");
    $oResult->index_date("today");
    $oResult->raw(qq{<A HREF="url$i">});
    $oResult->score(100-$i*2);
    $oResult->normalized_score(1000-$i*20);
    $oResult->size($i*2*1024);
    $oResult->source('WWW::Search::Null::Count');
    $oResult->company('Dub Dub Dub Search, Inc.');
    $oResult->location('Ashburn, VA');
    if ($i % 2)
      {
      $oResult->urls("url$i", map { "url$i.$_" } (1..$iCount));
      $oResult->related_urls(map { "url-r$i.$_" } (1..$iCount));
      my @aoTitles = map { "title-r$i.$_" } (1..$iCount);
      $oResult->related_titles(\@aoTitles);
      }
    else
      {
      for my $j (1..$iCount)
        {
        $oResult->add_url(qq{url$i.$j});
        $oResult->add_related_url(qq{url-r$j});
        $oResult->add_related_title(qq{title-r$i});
        } # for $j
      } # else
    push(@{$self->{cache}}, $oResult);
    } # for $i
  return 0;
  } # native_retrieve_some


1;


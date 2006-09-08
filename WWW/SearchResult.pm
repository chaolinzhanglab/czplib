# SearchResult.pm
# by John Heidemann
# Copyright (C) 1996 by USC/ISI
# $Id: SearchResult.pm,v 1.1 2006/09/08 01:55:28 zhangc Exp $
#
# Copyright (c) 1996 University of Southern California.
# All rights reserved.
#
# Redistribution and use in source and binary forms are permitted
# provided that the above copyright notice and this paragraph are
# duplicated in all such forms and that any documentation, advertising
# materials, and other materials related to such distribution and use
# acknowledge that the software was developed by the University of
# Southern California, Information Sciences Institute.  The name of the
# University may not be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.


=head1 NAME

WWW::SearchResult - class for results returned from WWW::Search

=head1 SYNOPSIS

    require WWW::Search;
    require WWW::SearchResult;
    $search = new WWW::Search;
    $search->native_query(WWW::Search::escape_query($query));
    # Get first result:
    $result = $search->next_result();

=head1 DESCRIPTION

A framework for returning the results of C<WWW::Search>.

=head1 SEE ALSO

L<WWW::Search>

=head1 REQUIRED RESULTS

The particular fields returned in a result are backend- (search
engine-) dependent.  However, all search engines are required to
return a url and title.  (This list may grow in the future.)

=head1 METHODS AND FUNCTIONS

=cut

#####################################################################

package WWW::SearchResult;

require LWP::MemberMixin;
@ISA = qw(LWP::MemberMixin);
use Carp ();
$VERSION = do { my @r = (q$Revision: 1.1 $ =~ /\d+/g); sprintf "%d."."%03d" x $#r, @r };

=head2 new

To create a new WWW::SearchResult, call

    $result = new WWW::SearchResult();

=cut

sub new
  {
  my $class = shift;
  my $self = bless { }, $class;
  $self->{urls} = ();
  return $self;
  } # new

=head2 url

Returns the primary URL.  Note that there may be a list of urls, see
also methods C<urls> and C<add_url>.  Nothing special is guaranteed
about the primary URL other than that it is the first one returned by
the back end.

Every result is required to have at least one URL.

=cut

sub url
  {
  my $self = shift;
  if (@_ < 1)
    {
    # No arguments, just return the current value:
    return ${$self->{urls}}[0];
    } # if no args
  unshift @{$self->{urls}}, $_[0];
  return $self->{urls}->[0];
  } # url

sub _elem_array
  {
  my $self = shift;
  my $elem = shift;
  if (@_ < 1)
    {
    # No arguments
    return wantarray ? @{$self->{$elem}} : $self->{$elem};
    } # if
  if (ref($_[0]))
    {
    # Trying to assign an arrayref:
    $self->{$elem} = $_[0];
    }
  else
    {
    # Trying to set to a scalar (or list of scalars); make sure it's
    # an array even if they give one element:
    $self->{$elem} = undef;
    push @{$self->{$elem}}, @_;
    }
  # Always return array refrence
  return $self->{$elem};
  } # _elem_array

sub _add_elem_array
  {
  my $self = shift;
  my $elem = shift;
  # No matter how many they're adding:
  push(@{$self->{$elem}}, @_);
  } # _add_elem_array


=head2 add_url

Add a URL to the list.

=cut

sub add_url { return shift->_add_elem_array('urls', @_); }

=head2 urls

Return a reference to the list of urls.
There is also a primary URL (C<url>).

=cut

sub urls { return shift->_elem_array('urls', @_); }

=head2 add_related_url

Add a URL to the related_url list.

=cut

sub add_related_url { return shift->_add_elem_array('related_urls', @_); }

=head2 related_urls

Return a reference to the list of related urls.

=cut

sub related_urls { return shift->_elem_array('related_urls', @_); }

=head2 add_related_title

Add a title to the list or related titles.

=cut

sub add_related_title { return shift->_add_elem_array('related_titles', @_); }

=head2 related_titles

Return a reference to the list of related titles.

=cut

sub related_titles { return shift->_elem_array('related_titles', @_); }


=head2 title, description, score, change_date, index_date, size, raw

Set or get attributes of the result.

None of these attributes is guaranteed to be provided by 
a given backend.  If an attribute is not provided
its method will return C<undef>.

Typical contents of these attributes:

=over 8

=item title

The title of the hit result (typically that provided by the 'TITLE'
HTML tag).

=cut

sub title { return shift->_elem('title', @_); }

=item description

A brief description of the result, as provided (or not) by the search engine.
Often the first few sentences of the document.

=cut

sub description { return shift->_elem('description', @_); }

=item source

Source is either the base url for this result (as listed on the search
engine's results page) or another copy of the full url path of the
result.  It might also indicate the source site address where the
resource was found, for example, 'http://www.cnn.com' if the search
result page said "found at CNN.com".

This value is backend-specific; in fact very few backends set this
value.

=cut

sub source { return shift->_elem('source', @_); }

=item score

A backend specific, numeric score of the search result.
The exact range of scores is search-engine specific.
Usually larger scores are better, but this is no longer required.
See normalized_score for a backend independent score.

=cut

sub score { return shift->_elem('score', @_); }

=item normalized_score

This is intended to be a backend-independent score of the search
result.  The range of this score is between 0 and 1000.  Higher values
indicate better quality results.

This is not really implemented since no one has created an
backend-independent ranking algorithm.

=cut

sub normalized_score { return shift->_elem('normalized_score', @_); }

=item change_date

When the result was last changed.  Typically this is the modification
time of the destination web page.

=cut

sub change_date { return shift->_elem('change_date', @_); }

=item index_date

When the search engine indexed the result.

=cut

sub index_date { return shift->_elem('index_date', @_); }

=item size

The approximate size of the result, in bytes.  This is only an
approximation because search backends often report the size as
"18.4K"; the best we can do with that number is return it as the value
of 18.4 * 1024.

=cut

sub size { return shift->_elem('size', @_); }

=item raw

The raw HTML for the entire result.  Raw should be exactly the raw
HTML for one entry.  It should not include list or table setup
commands (like ul or table tags), but it may include list item or
table data commands (like li, tr, or td).  Whether raw contains a list
entry, table row, br-separated lines, or plain text is search-engine
dependent.  In fact, many backends do not even return it at all.

=back

=cut

sub raw { return shift->_elem('raw', @_); }

=head2 Others

More attributes of the result.  Backend-specific.

=cut

sub company { return shift->_elem('company', @_); }
sub location { return shift->_elem('location', @_); }
sub bid_count { return shift->_elem('bids', @_); }
sub bid_amount { return shift->_elem('bid', @_); }

=head1 AUTHOR

WWW::SearchResult was written by John Heidemann.
WWW::SearchResult is maintained by Martin Thurn.

=cut

1;

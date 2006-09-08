# PubMed.pm
# by Jim Smyser
# Copyright (C) 2000 by Jim Smyser 
# $Id: PubMed.pm,v 1.1 2006/09/08 01:55:31 zhangc Exp $

#redeveloped by Chaolin Zhang
#Copyright (c) 2002 by Chaolin Zhang
#Email: zcl98@mails.tsinghua.edu.cn
#May, 2002

package WWW::Search::PubMed;

=head1 NAME

WWW::Search::PubMed - class for searching National Library of Medicine 

=head1 SYNOPSIS

use WWW::Search;

$query = "lung cancer treatment"; 
$search = new WWW::Search('PubMed');
$search->native_query(WWW::Search::escape_query($query));
$search->maximum_to_retrieve(100);
while (my $result = $search->next_result()) {

$url = $result->url;
$title = $result->title;
$desc = $result->description;

print <a href=$url>$title<br>$desc<p>\n"; 
} 

=head1 DESCRIPTION

WWW::Search class for searching National Library of Medicine
(PubMed). If you never heard of PubMed, Medline or don't know
the difference between a Abstract and Citation -- you then
can live without this backend.

This class exports no public interface; all interaction should
be done through WWW::Search objects.

=head1 AUTHOR

C<WWW::Search::PubMed> is written and maintained by Jim Smyser
<jsmyser@bigfoot.com>.

=head1 COPYRIGHT

WWW::Search Copyright (c) 1996-1998 University of Southern California.
All rights reserved. PubMed.pm by Jim Smyser.                                           
                                                               
THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

=cut
#'

#####################################################################
require Exporter;
@EXPORT = qw();
@EXPORT_OK = qw();
@ISA = qw(WWW::Search Exporter);
$VERSION = '1.0';

use Carp ();
use WWW::Search(qw(generic_option strip_tags));

require WWW::SearchResult;

sub native_setup_search 
{
    my($self, $native_query, $native_options_ref) = @_;
    $self->{_debug} = $native_options_ref->{'search_debug'};
    $self->{_debug} = 2 if ($native_options_ref->{'search_parse_debug'});
    $self->{_debug} = 0 if (!defined($self->{_debug}));
    $self->{agent_e_mail} = 'jsmyser@bigfoot.com';
    $max =  $self->maximum_to_retrieve;
    $self->user_agent('user');
    $self->{_next_to_retrieve} = 1;
    $self->{'_num_hits'} = 0;
    if (!defined($self->{_options})) 
	{
		$self->{_options} = 
		{
 			'search_url' => 'http://www.ncbi.nlm.nih.gov/entrez/query.fcgi',
			'db' => 'PubMed',
			'orig_db' => 'PubMed',
			'term' => $native_query,
			'cmd' => 'search',
			'cmd_current' => '',
			'WebEnv' => '',
			'dispmax' => '10'
		};
		#modified by Chaolin Zhang here. allow the user to override options
    }
    my $options_ref = $self->{_options};    
	if (defined($native_options_ref))
    {
        # Copy in new options.
        foreach (keys %$native_options_ref)
        {
        	$options_ref->{$_} = $native_options_ref->{$_};
        } 
	} 
        # Process the options.
    my($options) = '';
    foreach (sort keys %$options_ref)
    {
        next if (generic_option($_));
        $options .= $_ . '=' . $options_ref->{$_} . '&';
    }
    chop $options;#the string $options is not added to search_url?
    $self->{_next_url} = $self->{_options}{'search_url'} . "?". $options;

	print STDERR "options= $options\n" if $self->{_debug};
	print STDERR "_next_url = ", $self->{_next_url}, "\n" if $self->{_debug};
} 

# private
sub native_retrieve_some 
{
	

    my ($self) = @_;
    
    print STDERR "Entering sub routine WWW::Search::PubMed::native_retrieve_some...\n"
	if $self->{_debug};
    # Fast exit if already done:
    return undef if (!defined($self->{_next_url}));
   	
	print STDERR "_next_url=", $self->{_next_url}, "\n" if $self->{_debug};
    
	# If this is not the first page of results, sleep so as to not
    # overload the server:
    
	print STDERR "_next_to_retrieve = ", $self->{'_next_to_retrieve'}, "\n"
	if $self->{_debug};
    $self->user_agent_delay if 1 < $self->{'_next_to_retrieve'};
            
    # Get some if were not already scoring somewhere else:
    my($response) = $self->http_request('GET', $self->{_next_url});
        
    $self->{response} = $response;
    if (!$response->is_success)
    {
		print STDERR "Http response is not successful, return undef\n" 
		if $self->{_debug};
        return undef;
    }
    $self->{'_next_url'} = undef;
	
	print STDERR "Http response is successful, go on\n" 
    if $self->{_debug};

	# parse the output
    my ($ITEMS, $PAGES, $CURPG, $HITS, $DESC) = qw(IT PG CP HI DE);
    my $hits_found = 0;
    my $state = $ITEMS;
    my $hit;
	my $total_page = 1;
	my $curr_page = 1;
    #print STDERR join ("\n", $self->split_lines($response->content()))
	#if $self->{_debug};    

    # parse the output
	foreach ($self->split_lines($response->content()))
    {
        next if m@^\s+$@; # short circuit for blank lines
        
		if ($state eq $ITEMS && m|>Items .*?of ([\d,]+)</|i) 
        {
    	    $self->approximate_result_count($1);
			print STDERR "total_items = $1\n"
			if $self->{_debug};
        	$state = $PAGES;
        }
		elsif ($state eq $PAGES && m@type="BUTTON" value="Page">.*?type="text" value="(\d+)"></td>$@i) 
        {
    	    #$self->approximate_result_count($1);
    	    $curr_page = $1;
			print STDERR "current page=$curr_page\n";
			
        	$state = $CURPG;
        }
		elsif (($state eq $PAGES || $state eq $CURPG) && m|^<td><div class="medium2">&nbsp;of&nbsp;(\d+)</div></td>$|i) 
  		{
			$total_page = $1;
			print STDERR "total_page = $total_page\n"
			if $self->{_debug};
			#if $self->{_debug};
			$state = $HITS;
		}
   		elsif (($state eq $PAGES || $state eq $HITS) && m@^<td width="100%"><.*?><a href="(.*?)">(.+)</a></font></td>$@i) 
        {
			$hit = ();
        	my ($url, $title) = ($1,$2);
        	$hit = new WWW::SearchResult;
        	$hits_found++;
        	$url =~ s/dopt=Abstract/dopt=XML/g;
        	$url =~ s/amp;//g;
        	$hit->add_url($url);
        	$hit->title($title);

			if ($self->{_debug})
			{
				print STDERR "[**PubMed Entry match**]\n";
				print STDERR "url=$url\n";
				print STDERR "title=$title\n";
			}
			print STDERR "$url\n";
        	$state = $DESC;
        } 
    	elsif ($state eq $DESC && m|^<td colspan="2">(.+)</font></td>$|i) 
        {
        	$desc = $1;
        	$desc =~ s/<font size=\"-1\">//g;
        	$desc =~ s/<br>/\n/g;
			$hit->description($desc);
			if ($self->{_debug})
			{
				print STDERR "description=$desc\n";
			}
			if (defined($hit)) 
			{
            	push(@{$self->{cache}}, $hit);
        	}
        	$state = $HITS;
        }
		else
		{
		}
	}

    if ($curr_page >= $total_page && $state eq $HITS)
	{
		#last page! this page has at least one hits!!
		print  STDERR "query end\n"
		if $self->{_debug};
       	$self->{_next_url} = undef;
	}
	elsif ($hits_found < $self->{_options}{'dispmax'})
	{
		#some error occur, retrieve this page again
		#most likely because of network problem such as time out
		my $options_ref = $self->{_options};
		$options_ref->{'showndispmax'} = (defined ($options_ref->{'dispmax'})) ? $options_ref->{'dispmax'} : 10;
    
    	my($options) = '';
    	foreach (sort keys %$options_ref)
    	{
       		next if (generic_option($_));
       		$options .= $_ . '=' . $options_ref->{$_} . '&';
    	}
    	chop $options;#the string $options is not added to search_url?
    	$self->{_next_url} = $self->{_options}{'search_url'} . "?". $options . "&page+$curr_page.x=5&page+$curr_page.y=3";
		print  STDERR "url of next page=$self->{_next_url}\n"
		if $self->{_debug};
	}
	else
	{
		$curr_page++;
		my $options_ref = $self->{_options};
    	$options_ref->{'showndispmax'} = (defined ($options_ref->{'dispmax'})) ? $options_ref->{'dispmax'} : 10;
    	my($options) = '';
    	foreach (sort keys %$options_ref)
    	{
       		next if (generic_option($_));
       		$options .= $_ . '=' . $options_ref->{$_} . '&';
    	}
    	chop $options;#the string $options is not added to search_url?
    	$self->{_next_url} = $self->{_options}{'search_url'} . "?". $options . "&page+$curr_page.x=5&page+$curr_page.y=3";
		print  STDERR "url of next page=$self->{_next_url}\n"
		if $self->{_debug};
	};
    return $hits_found;
}
1;

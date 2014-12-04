#
#===============================================================================
#
#         FILE:  PhyloTree.pm
#
#  DESCRIPTION:  Package to handle Phylogenetic trees
#         BUGS:  ---
#        NOTES:  The package is spun off from the AnnotationIO.pm
#       AUTHOR:  Chaolin Zhang (cz), czhang@rockefeller.edu
#      COMPANY:  Rockefeller University
#      VERSION:  1.0
#      CREATED:  12/16/10 20:58:53
#     REVISION:  ---
#===============================================================================


package PhyloTree;


require Exporter;


our $VERSION = 1.01;

@ISA = qw (Exporter);

@EXPORT = qw (
	codeTree
	copyTree
	getLeafNodes
	getNodes
	nodeInfo
	parseTree
	readTreeFile
	releaseNode
	releaseTree
	removeLeafNodes
	removeSpecies
	segmentToken
	subTree
	totalBranchLength
	writeTreeFile
);

=head1 NAME

PyloTree - handle phylogenetic trees


=cut

use strict;
use warnings;

use Data::Dumper;
use Carp;

my $debug = 0;

=head2 nodeInfo

return information of a node of the tree

my $str = nodeInfo ($node)

=cut

sub nodeInfo
{
    my $node = $_[0];
    return $node->{"iter"} . "($node)";
}
		
		
=head2 readTreeFile

my $tree = readTreeFile ($inFile)

=cut


sub readTreeFile
{
	my $in = $_[0];
	my $fin;
	open ($fin, "<$in") || Carp::croak "can not open file $in to read\n";

	my $treeStr = "";
	
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		$line =~s/\s//g;

		$treeStr .= $line;
	}
	close ($fin);
	chop $treeStr if ($treeStr =~/\;$/);

	my $tokens = segmentToken ($treeStr);
	
	return parseTree ($tokens);
}


=head2 writeTreeFile

writeTreeFile ($tree, $outFile)

=cut

sub writeTreeFile
{
	my ($tree, $out) = @_;
	my $treeStr = codeTree ($tree);
	my $fout;
	open ($fout, ">$out") || Carp::croak "can not open file $out to  write\n";

	print $fout $treeStr, ";\n";

	close ($fout);
}

=head2 segmentToken

my $tokens = segmentToken ($treeStr)

=cut
sub segmentToken
{
	my $treeStr = $_[0];
	
	##segment the tree text
	my @tokens; #= split (/\(|\,|\)/, $treeStr);
	
	my $tok = "";
	for (my $i = 0; $i < length ($treeStr); $i++)
	{
		my $c = substr($treeStr, $i, 1);
		if ($c eq '(')
		{
			push @tokens, $c;
		}
		elsif ($c eq ',')
		{
			push @tokens, $tok;
			push @tokens, $c;
			$tok = '';
		}
		elsif ($c eq ')')
		{
			push @tokens, $tok;
			push @tokens, $c;
			$tok = "";
		}
		else
		{
			$tok .= $c;
		}
	}
	if ($tok ne '')
	{
		push @tokens, $tok;
	}
	push @tokens, ";";
	return \@tokens;
}


=head2 parseTree

my $tree = parseTree ($tokens)

=cut

sub parseTree
{
	my $tokens = $_[0];
	my @tokens = @$tokens;

	if (@tokens == 2)
	{
		return {iter => 0, id=>$tokens[0], blen=> 0};
	}
	
	my @nodes;
	
	my $root = {iter=> 0, id=>"", blen=> 0};
	push @nodes, $root;
	
	my $nodeIter = 1;
	my $currIter = $root->{"iter"};
	
	foreach my $tok (@tokens)
	{
	
		if ($tok eq '(')
		{
			my $leftChild = {iter=>$nodeIter++, parent=>$nodes[$currIter]};
			my $rightChild = {iter=>$nodeIter++, parent=>$nodes[$currIter]};
			push @nodes, $leftChild;
			push @nodes, $rightChild;
			
			$nodes[$currIter]->{"left"} = $leftChild;
			$nodes[$currIter]->{"right"} = $rightChild;

			print "tok = $tok, curr node = ", nodeInfo ($nodes[$currIter]), ", branch and go left child ", nodeInfo ($leftChild), "\n" if $debug;
			$currIter = $nodes[$currIter]->{"left"}->{"iter"};
		
		}
		elsif ($tok eq ',')
		{
			Carp::croak "no right sibling for node ", $nodes[$currIter]->{"iter"}, "\n" unless exists $nodes[$currIter]->{"parent"} && exists $nodes[$currIter]->{"parent"}->{"right"};
		
			print "tok = $tok, curr node = ", nodeInfo ($nodes[$currIter]), ", jump to right sibling ", nodeInfo ($nodes[$currIter]->{"parent"}->{"right"}), "\n" if $debug;
			$currIter = $nodes[$currIter]{"parent"}->{"right"}->{"iter"};
		}
		elsif ($tok eq ')')
		{
			#if ($currIter == 0) #done
			#{
			#	print "tok = $tok, curr node = ", nodeInfo ($nodes[$currIter]), ", done\n";
			#	return $root;
			#}
			Carp::croak "no parent for node ", $nodes[$currIter]->{"iter"}, "\n" unless exists $nodes[$currIter]->{"parent"};

			print "tok = $tok, curr node = ", nodeInfo ($nodes[$currIter]), ", jump to parent ", nodeInfo ($nodes[$currIter]->{"parent"}), "\n" if $debug;
			$currIter = $nodes[$currIter]->{"parent"}->{"iter"};
		}
		elsif ($tok eq ';')
		{
			print "tok = $tok, curr node = ", nodeInfo ($nodes[$currIter]), ", done\n" if $debug;
			#$root->{"nodes"} = \@nodes;
			return $root;
		}
		else
		{
			my ($id, $blen) = split (/\:/, $tok);
			$nodes[$currIter]->{"id"} = $id;
			$nodes[$currIter]->{"blen"} = $blen;

			print "tok = $tok, curr node = ", nodeInfo ($nodes[$currIter]), ", assign label $tok\n" if $debug;
		}
	}
}

=head2 parseTreeRecursive
need more test

=cut
sub parseTreeRecursive
{
	my ($tokens, $currNode, $nodeIter) = @_;

	my @tokens = @$tokens;

	my $tok = shift @tokens;
	
	if ($tok eq '(')
	{
		my $leftChild = {iter=>$nodeIter++, parent=>$currNode};
		my $rightChild = {iter=>$nodeIter++, parent=>$currNode};
			
		$currNode->{"left"} = $leftChild;
		$currNode->{"right"} = $rightChild;

		print "curr node = ", nodeInfo ($currNode), ", branch and go left child ", nodeInfo ($leftChild), "\n" if $debug;
		
		parseTree (\@tokens, $leftChild, $nodeIter);
	}
	elsif ($tok eq ',')
	{
		Carp::croak "no right sibling for node ", $currNode->{"iter"}, "\n" unless exists $currNode->{"parent"} && exists $currNode->{"parent"}->{"right"};
		
		print "curr node = ", nodeInfo ($currNode), ", jump to right sibling ", nodeInfo ($currNode->{"parent"}->{"right"}), "\n" if $debug;
		parseTree (\@tokens, $currNode->{"parent"}->{"right"}, $nodeIter);
	}
	elsif ($tok eq ')')
	{
		if (@tokens == 0) #done
		{
		print "curr node = ", nodeInfo ($currNode), ", done\n";
		return;
		}
		Carp::croak "no parent for node ", $currNode->{"iter"}, "\n" unless exists $currNode->{"parent"};

		print "curr node = ", nodeInfo ($currNode), ", jump to parent ", nodeInfo ($currNode->{"parent"}), "\n" if $debug;
	    parseTree (\@tokens, $currNode->{"parent"}, $nodeIter);
	}
	else
	{
		my ($id, $blen) = split (/\:/, $tok);
		$currNode->{"id"} = $id;
		$currNode->{"blen"} = $blen;

		print "curr node = ", nodeInfo ($currNode), ", assign label $tok\n" if $debug;
		parseTree (\@tokens, $currNode, $nodeIter);
	}
}



=head2 copyTree

my $to = copyTree ($from)

=cut
sub copyTree
{
	my $rootFrom = $_[0];
	my $rootTo = {};

	$rootTo = $_[1] if (@_ > 1);

	$rootTo->{"iter"} = $rootFrom->{"iter"};
	$rootTo->{"id"} = $rootFrom->{"id"};
	$rootTo->{"blen"} = $rootFrom->{"blen"};

	if (exists $rootFrom->{"left"})
	{
		
		my $leftChildTo = {};
		my $rightChildTo = {};
		$rootTo->{"left"} = $leftChildTo;
		$leftChildTo->{"parent"} = $rootTo;
		copyTree ($rootFrom->{"left"}, $rootTo->{"left"});
				
		$rootTo->{"right"} = $rightChildTo;
		$rightChildTo->{"parent"} = $rootTo;
		copyTree ($rootFrom->{"right"}, $rootTo->{"right"});
	}
	
	return $rootTo;
}


=head2 codeTree

return the Newick format text of the tree

my $str = codeTree ($tree)

=cut

sub codeTree
{
	my $root = $_[0];
	my $str = "";
	$str = $_[1] if (@_ > 1);

	if ($root == 0)
	{
		return ""; #empty
	}
	if (exists $root->{"left"})
	{
		$str .= "(";
		$str = codeTree ($root->{"left"}, $str);
		$str .= ",";
		$str = codeTree ($root->{"right"}, $str);
		$str .= ")";
		$str .= ":" . $root->{"blen"} if exists $root->{"parent"};
	}
	else #leaf
	{
		$str .= $root->{"id"};
	    $str .=	":" . $root->{"blen"} if exists $root->{"parent"};
	}
	return $str;
}

=head2 getNodes

my $nodes = getNodes ($tree)
will return all nodes of the tree

or

getNodes ($tree, \@nodes) 
will add nodes of $tree to @nodes

=cut
sub getNodes
{
	my $root = $_[0];
	my $nodes = [];
	$nodes = $_[1] if @_ > 1;

	push @$nodes, $root;
	if ($root->{"left"})
	{
		getNodes ($root->{"left"}, $nodes);
		getNodes ($root->{"right"}, $nodes);
	}
	return $nodes;
}

=head2 subTree

my $tree2 = subTree($tree, \@species);

get the minimal subtree that include all the given species 

=cut
sub subTree
{
	my ($tree, $species) = @_;
	#my $treeCpy = copyTree ($tree);
	my $leaves = getLeafNodes ($tree);

	my %speciesToKeep;

	foreach my $s (@$species)
	{
		Carp::croak "species $s does not exist in the tree\n" unless exists $leaves->{$s};
		$speciesToKeep{$s} = 1;
	}
	
	my @leavesToRm;
	foreach my $s (keys %$leaves)
	{
		push @leavesToRm, $leaves->{$s} unless exists $speciesToKeep{$s};
	}
	return removeLeafNodes ($tree, \@leavesToRm);
}


=head2 removeSpecies

my $tree2 = removeSpecies ($tree, \@species_to_remove)

=cut

sub removeSpecies
{
	my ($tree, $species) = @_;
	#my $treeCpy = copyTree ($tree);
	my $leaves = getLeafNodes ($tree);

	my @leavesToRm;
	foreach my $s (@$species)
	{
		Carp::croak "species $s does not exist in the three\n" unless exists $leaves->{$s};
		push @leavesToRm, $leaves->{$s}; # if exists $leaves->{$s};
	}
	return removeLeafNodes ($tree, \@leavesToRm);
}

=head2 removeLeafNodes

remove a set of leaf nodes

my $tree2 = removeLeafNodes ($tree, \@leaf_nodes_to_remove)

=cut

sub removeLeafNodes
{
	my ($root, $nodes) = @_;
	
	my $i = 0;
	foreach my $n (@$nodes)
	{
		Carp::croak "The node is not a leaf node: ", Dumper ($n), "\n" if exists $n->{"left"};
		
		print $i++,",  node to remove: ", $n->{"id"}, "\n" if $debug;
		
		if (not exists $n->{"parent"})
		{#single node tree
			return 0;
		}
		elsif (not exists $n->{"parent"}->{"parent"})
		{
			#parent is the root
			my $branch = ($n->{"parent"}->{"left"}->{"iter"} == $n->{"iter"})? 'right' : 'left';
			
			$root = $root->{$branch};
			$root->{"blen"} = 0;
			delete $root->{"parent"};
		}
		elsif ($n->{"parent"}->{"left"}->{"iter"} == $n->{"iter"})
		{
			#delete left leaf
			my $rightSibling = $n->{"parent"}->{"right"};
			$rightSibling->{"blen"} += $n->{"parent"}->{"blen"};
			my $grandParent = $n->{"parent"}->{"parent"};
			
			my $branch = ($grandParent->{"left"}->{"iter"} == $n->{"parent"}->{"iter"})? 'left' : 'right';
			$grandParent->{$branch} = $rightSibling;
			$rightSibling->{"parent"} = $grandParent;
			
		}
		else
		{
			#delete right leaf
			print "delete right leaf\n" if $debug;
			print "left sibling: ", $n->{"parent"}->{"left"}->{"id"}, "\n" if $debug;
			my $leftSibling = $n->{"parent"}->{"left"};
			$leftSibling->{"blen"} += $n->{"parent"}->{"blen"};
			my $grandParent = $n->{"parent"}->{"parent"};
			my $branch = ($grandParent->{"left"}->{"iter"} == $n->{"parent"}->{"iter"})? 'left' : 'right';
			$grandParent->{$branch} = $leftSibling;
			$leftSibling->{"parent"} = $grandParent;
		}

		releaseNode ($n->{"parent"}) if exists $n->{"parent"};
		releaseNode ($n);

		print codeTree ($root), "\n\n" if $debug;
	}
	return $root;
}

=head2 releaseNode

release the memory of a list of nodes

releaseNode (\%nodes)

=cut

sub releaseNode
{
	my $node = $_[0];
	foreach my $k (keys %$node)
	{
		delete $node->{$k};
	}
	return;
}

=head2 releaseTree

releaseTree ($tree)
=cut
sub releaseTree
{
	my $root = $_[0];
	if (exists $root->{"left"})
	{
		releaseTree ($root->{"left"});
		releaseTree ($root->{"right"});
	}
	
	releaseNode ($root);
}

=head2 totalBranchLengh

calculate the total branch length of a tree

my $tbl = totalBranchLength ($tree)

=cut

sub totalBranchLength
{
	my $root = $_[0];
	print "current node: ", nodeInfo ($root), "id=", $root->{"id"}, ", blen=", $root->{"blen"}, "\n" if $debug;
	my $tbn = 0;
	#$root->{"blen"} if exists $root->{"parent"};
	
	if (exists $root->{"left"})
	{
		$tbn += $root->{"left"}->{"blen"} + totalBranchLength ($root->{"left"});
		$tbn += $root->{"right"}->{"blen"} + totalBranchLength ($root->{"right"});
	}
	
	#$tbn = $root->{"blen"} if exists $root->{"parent"};
	return $tbn;
}

=head2 getLeafNodes

my $leafNodes = getLeafNodes ($tree)

=cut
sub getLeafNodes
{
	my $root = $_[0];
	my $nodes = getNodes ($root); 
	#print Dumper ($nodes), "\n";
	my %leaves;
	foreach my $n (@$nodes)
	{
		if (not exists $n->{"left"})
		{
			Carp::croak "no id for leaf node ", Dumper ($n), "\n" unless exists $n->{"id"} && $n->{"id"} ne '';
			$leaves{$n->{"id"}} = $n;
		}
	}
	return \%leaves;
}



1;



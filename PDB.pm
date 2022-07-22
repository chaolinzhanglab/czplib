#
#===============================================================================
#
#         FILE:  PDB.pm
#
#  DESCRIPTION:  Package to handle PDB files
#         BUGS:  ---
#        NOTES:  
#       AUTHOR:  Chaolin Zhang (cz), cz2294@columbia.edu
#      COMPANY:  Columbia University
#      VERSION:  1.0
#      CREATED:  02/10/2017
#     REVISION:  ---
#===============================================================================

package PDB;


require Exporter;

our $VERSION = 1.00;

@ISA = qw (Exporter);

@EXPORT = qw (
	three2one
	getPrimarySeq
	extractComplexInfo
	res2res
);



=head1 NAME

PDB - subroutines that deal with PDB files

subroutines starting with a hyphen should not be called outside

=cut

use strict;
use warnings;
use Data::Dumper;
use Carp;

use Bio::Structure::IO;


my %_t2o = (
    'ALA' => 'A',
    'ARG' => 'R',
    'ASN' => 'N',
    'ASP' => 'D',
    'CYS' => 'C',
    'GLN' => 'Q',
    'GLU' => 'E',
    'GLY' => 'G',
    'HIS' => 'H',
    'ILE' => 'I',
    'LEU' => 'L',
    'LYS' => 'K',
    'MET' => 'M',
    'PHE' => 'F',
    'PRO' => 'P',
    'SEC' => 'U',
    'SER' => 'S',
    'THR' => 'T',
    'TRP' => 'W',
    'TYR' => 'Y',
    'VAL' => 'V',
    'ASX' => 'B',
    'GLX' => 'Z'
);




=head2 three2one

convert amino acid in three letter format to one letter abbreviation

=cut

sub three2one
{
    my $aa = $_[0];
    return exists $_t2o{$aa} ? $_t2o{$aa} : "";
}



=head2 getPrimarySeq

#return amino acid or nucleotide sequence and type
{seq=>$seqStr, type=>$type}

#0 otherwise

=cut

sub getPrimarySeq
{
    my ($struct, $chain) = @_;

    my @seq;
    for my $res ($struct->get_residues ($chain))
    {
        my $resid = $res->id;
		#print $resid, "\n";
        my ($r, $n) = split("-", $resid);
		next if $r eq 'HOH';
		next if $r eq '6U0' || $r eq 'UVP' || $r eq 'GTP' || $r eq 'G5J'; #just workaround
        push @seq, $r; # assume HETATOM has been removed
    }

    return 0 if @seq < 1;

    my $type = "";
    my $seqStr = "";

	#print "chain=", $chain->id, "\n";

	#print Dumper (\@seq), "\n";

    if (length ($seq[0]) == 3)
    {
        if (three2one($seq[0]))
        {
            $type = "protein";
            foreach my $res (@seq)
            {
                if (length ($res) != 3)
				{
					#can be small molecule, eg. 1CVJ
					warn "unknown amino acid: $res\n";
					next;

					#if ($res eq 'DA' || $res eq 'DC' || $res eq 'DG' || $res eq 'DT')
                	#{
					#	#DNA adduct?
					#	return 0;
					#}
					#else
					#{
                    #	Carp::croak "unknown amino acid: $res in ", Dumper (\@seq), "\n";
					#}
                }

                if (my $aa = three2one ($res))
                {
                    $seqStr .= $aa;
                }
                elsif ($res eq 'UNK' || $res eq 'MSE' || $res eq 'AMP' || $res eq 'SO4' || $res eq 'DIO' || $res eq '1PE' 
						|| $res eq 'ADP' || $res eq 'ALF' || $res eq 'GOL' || $res eq 'ANP' || $res eq 'M2M' || $res eq 'EDO' 
						|| $res eq 'IPH' || $res eq 'IPA' || $res eq 'ACT' || $res eq 'UNX' || $res eq 'SAH' || $res eq 'TRS'
						|| $res eq 'CME' || $res eq 'PO4' || $res eq 'BEF')
                {
					warn "unknown amino acid: $res\n";
					next;
                    #return 0;
                }
                else
                {
                    Carp::croak "unknown amino acid: $res in ", Dumper (\@seq), "\n";
                }
            }
            #$seqStr = join ("", map {$three2one{$_}} @seq);
        }
        elsif ($seq[0] eq 'UNK')
        {
            return 0;
        }
        else
        {
            Carp::croak "unknown amino acid: $seq[0]\n";
        }
    }
    elsif (length ($seq[0]) == 1)
    {
        if ($seq[0] =~/[ACGUI]/i)
        {
            $type = "rna";

            foreach my $b (@seq)
            {
                if (length($b) != 1)
                {
                    if ($b eq 'DA' || $b eq 'DC' || $b eq 'DG' || $b eq 'DT' || $b eq 'DU')
                    {
						$type="dna_rna_hybrid";
						$seqStr .= substr($b, 1, 1);
						next;
                    }
					elsif ($b eq '5BU' || three2one($b) ne '')
					{
						#some weired stuff
						#return 0;
						warn "some weird stuff: $b\n";
						next;
					}
					else
					{
						warn "unknown nucleotide: $b\n";
						next;
                    	#Carp::croak "unknown nucleotide: $b in ", Dumper (\@seq), "\n";
					}
                }
	
                if ($b=~/[ACGUI]/i)
                {
                    $seqStr .= $b;
                }
                elsif ($b eq 'N' || $b eq 'K')
                {
					warn "unknown nucleotide: $b\n";
                    next;
                    #return 0;
                }
                else
                {
                    Carp::croak "unknown nucleotide: $b in ", Dumper (\@seq), "\n";
                }
            }
            #$seqStr = join ("", @seq);
        }
        elsif ($seq[0] eq 'N')
        {
			Carp::croak "unknown nucleotide: $seq[0]\n", Dumper (\@seq), "\n";
            return 0;
        }
        else
        {
            Carp::croak "unknown nucleotide: $seq[0]\n";
        }
    }
    elsif (length ($seq[0]) == 2)
    {
        my $b = substr ($seq[0], 1, 1);
        if ($b =~/[ACGTIU]/i)
        {
            $type = "dna";

            foreach my $bb (@seq)
            {
                if (length($bb) != 2)
                {
                    if ($bb eq 'A' || $bb eq 'C' || $bb eq 'G' || $bb eq 'U')
                    {
						$type="dna_rna_hybrid";
						$seqStr .= $bb;
						next;
                    }
					else
					{
                    	Carp::croak "unknown nucleotide: $bb in ", Dumper (\@seq), "\n";
					}
                }

                my $b = substr ($bb, 1, 1);
                if ($b=~/[ACGTIU]/i)
                {
                    $seqStr .= $b;
                }
                else
                {
                    Carp::croak "unknown nucleotide: $b in ", Dumper (\@seq), "\n";
                }
            }
            #$seqStr = join ("", @seq);
        }
    }

    if ($type eq 'protein' || $type eq 'rna' || $type eq 'dna' || $type eq 'dna_rna_hybrid')
    {
        return {seq=>$seqStr, type=>$type};
    }
    else
    {
        return 0;
    }
}


=head2 extractComplexInfo

#return reference of complex structure if legitimate
{
id=>
complex=>{chainid}->{db, acc, name, type, seq}
}

#return 0 otherwise

=cut

sub extractComplexInfo
{
	my $struct = $_[0];
    
	#my $inFile = $_[0];
    #my $io = Bio::Structure::IO->new(-file => $inFile, -format => 'PDB');
    #my $struct = $io->next_structure;

    my $structid = $struct->id;
    #Carp::croak Dumper ($struct), "\n";

    my @chains = $struct->get_chains;
    return 0 if @chains < 2;

    my %complex;

    for my $chain ($struct->get_chains)
    {
        my $chainid = $chain->id;

        my $seqInfo = getPrimarySeq ($struct, $chain);

		if ($seqInfo == 0)
		{
			warn ("unable to parse sequence in structure $structid, chain $chainid\n");
        	return 0;
        }

		my $type = $seqInfo->{'type'};
        my $seq = $seqInfo->{'seq'};

        $complex{$chainid} = {type=>$type, seq=>$seq};
    }

    my $annot = $struct->annotation;
    my ($dbref) = $annot->get_Annotations('dbref');
	
	if (not defined $dbref)
	{
		warn ("dbref not found in struct $structid\n");
    	return 0;
	}
    my @cols = split (/\s+/, $dbref->value);
	my $ncols = @cols;
	if ($ncols % 9 != 0)
	{
		warn "error in parsing dbref in struct $structid\n";
		return 0;
	}

    for (my $i = 0; $i < @cols; $i += 9)
    {
        my $chainid = $cols[$i+1];

		if (not exists $complex{$chainid})
		{
			warn "error in parsing dbref in struct $structid: chain $chainid does not exist\n";
			return 0;
		}

        my $db = $cols[$i+4];
        my $acc = $cols[$i+5];
        my $name = $cols[$i+6];

        $complex{$chainid}{'db'} = $db;
		$complex{$chainid}{'acc'} = $acc;
		$complex{$chainid}{'name'} = $name;
    }

	my ($compnd) = $annot->get_Annotations('compnd');
	
	my @molecules;
	
	my $compndStr = $compnd->value;
	$compndStr=~s/\s+$//g;
	$compndStr .= ";" unless $compndStr=~/;$/;
	#in some weired cases, there is no ; at the end, so we have to add this for correct parsing below

	#print "compnd=$compndStr\n";

	my $nchains = 0;

	#while ($compndStr=~m/MOL_ID: (\d+);\s+MOLECULE: (\S.*?);*\s*CHAIN: (\S.*?);/g)
	#while ($compndStr=~m/MOL_ID: (\d+);\s+MOLECULE: (\S.*?);*\s*CHAIN: (\w{1,2}[;|\,\s.*?;])/g)
	#print $compndStr, "\n";

	while ($compndStr=~m/MOL_ID: (\d+);\s+MOLECULE:\s+(\S.*?);*\s*CHAIN: (\w{1,2}(?:\s*\,\s*\w{1,2}){0,})\,?;/g)
	{
		my ($id, $desc, $chainId) = ($1, $2, $3);
		chop $chainId if $chainId=~/[;\,]$/;
		#print "id=$id, desc=$desc, chainId=$chainId\n";
	
		my %mol = (id=>$id, desc=>$desc, chain=>$chainId);
		if ($mol{'chain'}=~/\,$/)
		{
			#due to inconsistancy in some record
			chop $mol{'chain'};
		}

		my @chains = split(/\s*\,\s*/, $mol{'chain'});
		$nchains += @chains;
	
		if (not exists $complex{$chains[0]})
		{
			warn "chain $chains[0] does not exist in structure $structid. Detected chains: \n", join("\t", sort keys %complex), "\n";
			return 0;
		}

		$mol{'type'} = $complex{$chains[0]}{'type'};
		$mol{'seq'} = $complex{$chains[0]}{'seq'};

		$mol{'db'} = $complex{$chains[0]}{'db'};
		$mol{'acc'} = $complex{$chains[0]}{'acc'};
		$mol{'name'} = $complex{$chains[0]}{'name'};

		push @molecules, \%mol;
	}
	#Carp::croak Dumper (\@molecules), "\n";

	
	if ($nchains != keys %complex) #this can happy when there are many chains and the meta data was truncated by bioperl. rediculous! we ignore those stucture for now
	{
		warn "inconsistent number of chains in structure $structid. Detected chains: \n", join("\t", sort keys %complex), "\n", "Chains in compnd: \n", join(", ", map{$_->{"chain"}} @molecules), "\n";
		return 0;
	}

	#Carp::croak Dumper (\@molecules), "\n";

    my %ret = (id=>$structid, molecules=>\@molecules, chains=>\%complex);
    return \%ret;
}



=head2 res2res

calculate residue to residue distance by the pair of atoms closest to each other

=cut


sub res2res
{
	my ($struct, $res1, $res2) = @_;

	my $d = -1;
	my ($atomId1, $atomId2);

	foreach my $atom1 ($struct->get_atoms ($res1))
	{
		my ($x1, $y1, $z1) = $atom1->xyz;
		foreach my $atom2 ($struct->get_atoms ($res2))
		{
			my ($x2, $y2, $z2) = $atom2->xyz;
			my $d2 = sqrt(($x1-$x2)**2 + ($y1-$y2)**2 + ($z1-$z2)**2);

			next if $d >= 0 && $d2 > $d;
					
			$d = $d2;
			$atomId1 = $atom1->id;
			$atomId2 = $atom2->id;
		}
	}
	
	return {d=>$d, atom1=>$atomId1, atom2=>$atomId2};
}

1;



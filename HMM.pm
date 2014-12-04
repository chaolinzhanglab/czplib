#!/usr/bin/perl -w
#

package HMM;

use strict;
use warnings;
use Data::Dumper;
use Common;


our $VERSION = 1.01;


=head1 NAME

HMM - Hidden Markov Model

implemented according to Rabiner 1989

=cut


my $verbose = 1;
my $eps = 1e-40;
my $dblMax = 1e40;



=head2 setVerbose

setVerbose (1); #verbose mode

=cut


sub setVerbose
{
	$verbose = $_[0];
}


=head2 viterbi

my $states = viterbi ($obs, $A, $e, $pi, $logTransformed = 0);

$obs: a single or multiple observation sequences
$A: transition matrix
$e: emission matrix or emission function that can be called as &$e ($i, $o)
$pi: initial state distribution
$logTransformed: whether $A, $B, $pi are log are already log transformed

return:
 a single or multiple state sequences depending on the input observation sequences

=cut

sub viterbi
{
	my ($obs, $A, $e, $pi, $logTransformed) = @_;

	#Carp::croak "obs=", Dumper ($obs), "\n";
	#Carp::croak "logTransformed = $logTransformed\n";
	
	$logTransformed = 0 unless $logTransformed;
	
	if (ref($obs->[0]) eq 'ARRAY')
	{
		my @multi_q;
		foreach my $o (@$obs)
		{
			my $q = viterbi ($o, $A, $e, $pi, $logTransformed);
			push @multi_q, $q;
		}
		return \@multi_q;
	}
	
	my $N = @$pi; # the number of state
	my @delta;
	my @phi;

	my @q; #state

	my $T = @$obs;
	
	#initialization
	for (my $i = 0; $i < $N; $i++)
	{
		my $o = $obs->[0];
		my $B = 0;

		if (ref ($e) eq 'ARRAY')
		{
			$B = $e->[$i][$o];
		}
		else
		{
			$B = &$e($i, $o);
		}
		
		$delta[$i][0] = -$dblMax;
		
		if ($logTransformed)
		{
			$delta[$i][0] = $pi->[$i] + $B;
		}
		else
		{
			$delta[$i][0] = log ($pi->[$i]) + log ($B) unless ($pi->[$i] == 0 || $B == 0);
		}
		$phi[$i][0] = 0;
	}
	
	my $t = 0;
	#print "d[$t]=", join ("\t", @{$delta[$t]}), ", phi[$t]=", join ("\t", @{$phi[$t]})  , "\n";
	
	#recursion
	
	for (my $t = 1; $t < $T; $t++)
	{
		#print "$t of $T ...\n" if $verbose && $t - int ($t / 5000) * 5000 == 0;

		my $o = $obs->[$t];
		for (my $j = 0; $j < $N; $j++)
		{
			my %v;
			for (my $i = 0; $i < $N; $i++)
			{
				$v{$i} = -$dblMax;
				if ($logTransformed)
				{
					$v{$i} = $delta[$i][$t-1] + $A->[$i][$j];
				}
				else
				{
					$v{$i} = $delta[$i][$t-1] + log ($A->[$i][$j]) unless $A->[$i][$j] == 0;
				}
			}

			my @w = sort {$v{$b} <=> $v{$a}} keys %v;

			$phi[$j][$t] = $w[0];
			my $B = 0;
			if (ref ($e) eq 'ARRAY')
			{
				$B = $e->[$j][$o];
			}
			else
			{
				$B = &$e($j, $o);
			}

			$delta[$j][$t] = -$dblMax;
			if ($logTransformed)
			{
				$delta[$j][$t] = $v{$phi[$j][$t]} + $B;
			}
			else
			{
				$delta[$j][$t] = $v{$phi[$j][$t]} + log ($B) unless $B == 0;
			}
		}

		#print "d[$t]=", join ("\t", @{$delta[$t]}), ", phi[$t]=", join ("\t", @{$phi[$t]})  , "\n";
		
	}

	#termination
	
	my %v;

	for (my $i = 0; $i < $N; $i++)
	{
		$v{$i} = $delta[$i][$T - 1];
	}

	my @w = sort {$v{$b} <=> $v{$a}} keys %v;
	
	$q[$T-1] = $w[0];
	my $P = $v{$q[$T-1]};

	for (my $t = $T - 2; $t >= 0; $t--)
	{
		$q[$t] = $phi[$q[$t+1]][$t+1];
	}
	
	#release memory
	for (my $i = 0; $i < $N; $i++)
	{
		$delta[$t] = [];
		$phi[$t] = [];
	}
	@delta = ();
	@phi = ();
	
	return \@q;
}


=head2 forward

my $alpha = forward ($obs, $A, $e, $pi, $retDetail = 0);

$obs: a single or multiple  observation sequence
$A: reference to the transition matrix
$e: reference to the emission routine, &$e ($state, $obs), or tabularized emission matrix
$pi: initial state distribution 

return: $alpha->[$t][$j] or {alpha=>$alpha, scaling=>$scaling, logP=>$logP}

=cut

sub forward
{
	my ($obs, $A, $e, $pi) = @_;

	my $retDetail = 0;
	if (@_ > 4)
	{
		$retDetail = $_[4];
	}
	
	if (ref ($obs->[0]) eq 'ARRAY')
	{
		my @multi_ret;
		
		foreach my $o (@$obs)
		{
			my $ret = forward ($o, $A, $e, $pi, $retDetail);
			push @multi_ret, $ret;
		}
		return \@multi_ret;
	}
	
	my $N = @$pi;
	my $T = @$obs;

	my @alpha;

	my @scaling;

	for (my $i = 0; $i < $N; $i++)
	{
		$alpha[$i][$T - 1] = 0;
	}
	$scaling[$T - 1] = 0;
	
	#initialization
	my $c = 0;
	for (my $i = 0; $i < $N; $i++)
	{
		my $o = $obs->[0];
		my $B = 0;

		#print "i= $i, o=$o\n";
		if (ref ($e) eq 'ARRAY')
		{
			$B = $e->[$i][$o];
		}
		else
		{
			$B = &$e ($i, $o);
		}
	
		$alpha[$i][0] = $pi->[$i] * $B;

		$c+= $alpha[$i][0];
	}

	$scaling[0] = 1 /$c;

	for (my $i = 0; $i < $N; $i++)
	{
		$alpha[$i][0] /= $c;
	}

	#induction
	for (my $t = 0; $t < $T - 1; $t++)
	{
		#print "$t of $T ...\n" if $verbose && $t - int ($t / 50000) * 50000 == 0;

		my $o = $obs->[$t+1];
	
		my $c = 0;
		for (my $j = 0; $j < $N; $j++)
		{
			$alpha[$j][$t+1] = 0;
			
			for (my $i = 0; $i < $N; $i++)
			{
				$alpha[$j][$t+1] += $alpha[$i][$t] * $A->[$i][$j];
			}

			#print "j=$j, o=$o\n";

			my $B = 0;
			if (ref ($e) eq 'ARRAY')
			{
				$B = $e->[$j][$o];
			}
			else
			{
				$B = &$e ($j, $o);
			}
			$alpha[$j][$t+1] *= $B;

			$c += $alpha[$j][$t+1];
		}

		$scaling[$t+1] = 1 / $c;

		for (my $j = 0; $j < $N; $j++)
		{
			$alpha[$j][$t+1] /= $c;
		}
	}
	
	#print "alpha = ", Dumper (\@alpha), "\n";
	#print "scaling = ", Dumper (\@scaling), "\n";

	if ($retDetail)
	{
		my $logP = 0;
		for (my $j = 0; $j < $N; $j++)
		{
			$logP += $alpha[$j][$T -1];
		}
		$logP = log($logP);

		for (my $t = 0; $t < $T; $t++)
		{
			if ($scaling[$t] == 0)
			{
				$logP = -$dblMax;
				last;
			}
			$logP -= log ($scaling[$t]);
		}

		return {alpha=>\@alpha, scaling=>\@scaling, logP=>$logP};
	}
	else
	{
		@scaling = ();
		return \@alpha;
	}
}

=head2 backward

my $beta = backward ($obs, $A, $e, $pi, $scaling);

$obs: the observation sequence
$A: reference to the transition matrix
$e: reference to the emission routine, &$e ($state, $obs), or tabularized emission matrix
$pi: initial state distribution 
$scaling: the scaling coefficient obtained by the forward subroutine, for numerical stability

return: $beta->[$t][$j]

=cut

sub backward
{
	my ($obs, $A, $e, $pi, $scaling) = @_;


	if (ref ($obs->[0]) eq 'ARRAY')
	{
		my @multi_ret;
		foreach my $o (@$obs)
		{
			my $ret = backward ($o, $A, $e, $pi);
			push @multi_ret, $ret;
		}
		return \@multi_ret;
	}
	
	my $N = @$pi;
	my $T = @$obs;

	my @beta;

	#initialization
	for (my $i = 0; $i < $N; $i++)
	{
		$beta[$i][$T-1] = 1;
	}

	#induction
	
	for (my $t = $T - 2; $t >= 0; $t--)
	{
		#print "$t of $T ...\n" if $verbose && $t - int ($t / 5000) * 5000 == 0;
		
		my $o = $obs->[$t+1];
		for (my $i = 0; $i < $N; $i++)
		{
			$beta[$i][$t] = 0;

			for (my $j = 0; $j < $N; $j++)
			{
				my $B = 0;
				if (ref ($e) eq 'ARRAY')
				{
					$B = $e->[$j][$o];
				}
				else
				{
					$B = &$e ($j, $o);
				}
	
				$beta[$i][$t] += $A->[$i][$j] * $B * $beta[$j][$t+1];
			}

			#scaling using the coefficient obtained by forward algorithm
			$beta[$i][$t] *= $scaling->[$t];
		}
	}

	#print "beta=", Dumper (\@beta), "\n";

	return \@beta;
	
}

=head2 posteriorDecode

my $ret = posteriorDecode ($obs, $A, $e, $pi, $retDetail);

$obs: the observation sequence
$A: reference to the transition matrix
$e: reference to the emission routine, &$e ($state, $obs), or tabularized emission matrix
$pi: initial state distribution 
$retDetail: whether to return more information (optional)

return: $gamma->[$t][$j]
 or {alpha=>$alpha, beta=>$beta, gamma=>\@gamma}
=cut

sub posteriorDecode
{
	my ($obs, $A, $e, $pi) = @_;

	my $retDetail = 0;

	if (@_ > 4)
	{
		$retDetail = $_[4];
	}
	
	if (ref ($obs->[0]) eq 'ARRAY')
	{
		my @multi_ret;
		foreach my $o (@$obs)
		{
			my $ret = posteriorDecode ($o, $A, $e, $pi, $retDetail);
			push @multi_ret, $ret;
		}
		return \@multi_ret;
	}


	my $N = @$pi;
	my $T = @$obs;

	my $ret = forward ($obs, $A, $e, $pi, 1);
	
	my $alpha = $ret->{"alpha"};
	my $scaling = $ret->{"scaling"};

	my $beta = backward ($obs, $A, $e, $pi, $scaling);

	my @gamma;
	for (my $t = 0; $t < $T; $t++)
	{
		#print "$t of $T ...\n" if $verbose && $t - int ($t / 5000) * 5000 == 0;
		
		my $P = 0;
		for (my $i = 0; $i < $N; $i++)
		{
			$P += $alpha->[$i][$t] * $beta->[$i][$t];
		}

		for (my $i = 0; $i < $N; $i++)
		{
			$gamma[$i][$t] = $alpha->[$i][$t] * $beta->[$i][$t] / $P;
		}
	}

	if ($retDetail)
	{
		return {alpha=>$alpha, beta=>$beta, gamma=>\@gamma, scaling=>$scaling};
	}
	else
	{
		for (my $i = 0; $i < $N; $i++)
		{
			$alpha->[$i] = [];
			$beta->[$i] = [];
		}
		$alpha = [];
		$beta = [];
		$scaling = [];
		return \@gamma;
	}
}

=obsolete
# do one round of baumwelch update
# 
# obs: the observation sequence
# A: reference to the transition matrix
# e: reference to the emission routine, &$e ($state, $obs), or tabularized emission matrix
# pi: initial state distribution 
# return: $gamma->[$t][$j]
#

sub _baumwelch
{
	my ($obs, $A, $e, $pi) = @_;

	my $multi_obs;

	if (ref ($obs->[0]) eq 'ARRAY')
	{
		$multi_obs = $obs;
	}
	else
	{
		$multi_obs = [$obs];
	}
	

	my $alphabetSize = 0;

	if (@_ > 4)
	{
		$alphabetSize = $_[4];
	}
	else
	{
		$alphabetSize = Common::max (@$obs) + 1;
	}
	
	my $N = @$pi;

	
	
	my $T = @$obs;


	#copy parameters
	my @A2;
	my @B2;
	
	my @pi2 = @$pi;
	for (my $i = 0; $i < $N; $i++)
	{
		@{$A2[$i]} = @{$A->[$i]};
		for (my $j = 0; $j < $alphabetSize; $j++)
		{
			if (ref ($e) eq 'ARRAY')
			{
				$B2[$i][$j] = $e->[$i][$j];
			}
			else
			{
				$B2[$i][$j] = &$e ($i, $j);
			}
		}
	}
	
	print "do posterior decoding ...\n" if $verbose;

	my $ret = posteriorDecode ($obs, \@A2, \@B2, \@pi2, 1);

	my $alpha = $ret->{"alpha"};
	my $beta = $ret->{"beta"};
	my $gamma = $ret->{"gamma"};
	
	
	my @theta;

	print "calculate theta ...\n" if $verbose;
	for (my $t = 0; $t < $T - 1; $t++)
	{
		#print "$t of $T ...\n" if $verbose && $t - int ($t / 5000) * 5000 == 0;
		
		my $o = $obs->[$t+1];
		my $P = 0;
		for (my $i = 0; $i < $N; $i++)
		{
			for (my $j = 0; $j < $N; $j++)
			{
				$theta[$i][$j][$t] = $alpha->[$i][$t] * $A2[$i][$j] * $B2[$j][$o] * $beta->[$j][$t+1];
				$P += $theta[$i][$j][$t];
			}
		}
		for (my $i = 0; $i < $N; $i++)
		{
			for (my $j = 0; $j < $N; $j++)
			{
				$theta[$i][$j][$t] /= $P; 
			}
		}
	}
	

	#update parameters
	
	print "update pi ...\n" if $verbose;
	
	for (my $i = 0; $i < $N; $i++)
	{
		$pi2[$i] = $gamma->[$i][0];
	}
	
	my @gammaS;

	for (my $j = 0; $j < $N; $j++)
	{
		$gammaS[$j] = 0;
		for (my $t = 0; $t < $T - 1; $t++)
		{
			$gammaS[$j] += $gamma->[$j][$t];
		}
	}
	
	print "update transition matrix ...\n" if $verbose;

	for (my $i = 0; $i < $N; $i++)
	{
		for (my $j = 0; $j < $N; $j++)
		{
			$A2[$i][$j] = 0;
			for (my $t = 0; $t < $T - 1; $t++)
			{
				$A2[$i][$j] += $theta[$i][$j][$t];
			}
			$A2[$i][$j] /= $gammaS[$i];
		}
	}

	print "update emission matrix ...\n" if $verbose;
	
	for (my $j = 0; $j < $N; $j++)
	{
		$gammaS[$j] += $gamma->[$j][$T-1];
	}
	
	my @gammaS2;

	for (my $j = 0; $j < $N; $j++)
	{
		for (my $k = 0; $k < $alphabetSize; $k++)
		{
			$gammaS2[$j][$k] = 0;
		}
	}
	
	for (my $j = 0; $j < $N; $j++)
	{
		for (my $t = 0; $t < $T; $t++)
		{
			my $o = $obs->[$t];
			$gammaS2[$j][$o] += $gamma->[$j][$t];
		}
	}
	
	for (my $j = 0; $j < $N; $j++)
	{
		for (my $k = 0; $k < $alphabetSize; $k++)
		{
			$B2[$j][$k] = $gammaS2[$j][$k] / $gammaS[$j];
		}
	}

	#release memory
	for (my $i = 0; $i < $N; $i++)
	{
		$alpha->[$i] = [];
		$beta->[$i] = [];
		$gamma->[$i] = [];
		for (my $j = 0; $j < $N; $j++)
		{
			$theta[$i][$j] = [];
		}
		$theta[$i] = [];
	}
	$alpha = [];
	$beta = [];
	$gamma = [];
	@theta = ();

	return {A=>\@A2, B=>\@B2, pi=>\@pi2};
}
=cut


=head2 baumwelch
 do one round of baumwelch update

my $gamma = baumwelch ($mobs, $A, $e, $pi);

$mobs: a single or multiple observation sequences
$A: reference to the transition matrix
$e: reference to the emission routine, &$e ($state, $obs), or tabularized emission matrix
$pi: initial state distribution, averaged from all sequences

return: $gamma->[$t][$j]

=cut

sub baumwelch
{
	my ($mobs, $A, $e, $pi) = @_;

	my $multi_obs;

	if (ref ($mobs->[0]) eq 'ARRAY')
	{
		$multi_obs = $mobs;
	}
	else
	{
		$multi_obs = [$mobs];
	}
	

	my $alphabetSize = 0;

	if (@_ > 4)
	{
		$alphabetSize = $_[4];
	}
	else
	{
		$alphabetSize = 0;
		foreach my $obs (@$multi_obs)
		{
			$alphabetSize = Common::max (@$obs, $alphabetSize) + 1;
		}
	}
	
	my $N = @$pi;

	#copy parameters
	my @A2;
	my @B2;
	
	my @pi2; # = @$pi;
	

	my @Anom; # the nominator of A
	my @Bnom; # the nominator of B
	my @Aden; # the denominator of A and B
	my @Bden;
	
	for (my $i = 0; $i < $N; $i++)
	{
		$pi2[$i] = 0;
		for (my $j = 0; $j < $N; $j++)
		{
			$Anom[$i][$j] = 0;
		}

		for (my $k = 0; $k < $alphabetSize; $k++)
		{
			$Bnom[$i][$k] = 0;
			if (ref ($e) eq 'ARRAY')
			{
				$B2[$i][$k] = $e->[$i][$k];
			}
			else
			{
				$B2[$i][$k] = &$e ($i, $k);
			}
		}
	}
	
	for (my $s = 0; $s < @$multi_obs; $s++)
	{
		print "##processing observation sequence #$s...\n" if $verbose;

		my $obs = $multi_obs->[$s];
		my $T = @$obs;
		
		print "do posterior decoding ...\n" if $verbose;

		my $piTmp = $pi;
		$piTmp = $pi->[0] if (ref ($pi->[0]) eq 'ARRAY');

		my $ret = posteriorDecode ($obs, $A, $e, $piTmp, 1);

		my $alpha = $ret->{"alpha"};
		my $beta = $ret->{"beta"};
		my $gamma = $ret->{"gamma"};
		my $scaling = $ret->{"scaling"};
		#update pi

		##my @piNew;
		for (my $j = 0; $j < $N; $j++)
		{
			##$piNew[$j] = $gamma->[$j][0];
			$pi2[$j] += $gamma->[$j][0];
		}
	
		##$pi2[$s] = \@piNew;
		

		#my $P = 0; # the posterior probability of observing a sequence

		#for (my $j = 0; $j < $N; $j++)
		#{
		#	$P += $alpha->[$j][$T-1];
		#}
		
		
		#print "P=$P\n";
		#print "calculate theta ...\n" if $verbose;
		
		for (my $t = 0; $t < $T - 1; $t++)
		{
			#print "$t of $T ...\n" if $verbose && $t - int ($t / 5000) * 5000 == 0;
		
			my $o = $obs->[$t+1];
			my $k = $obs->[$t];

			for (my $i = 0; $i < $N; $i++)
			{
				for (my $j = 0; $j < $N; $j++)
				{
					$Anom[$i][$j] += $alpha->[$i][$t] * $A->[$i][$j] * $B2[$j][$o] * $beta->[$j][$t+1];
				}

				$Aden[$i] += $alpha->[$i][$t] * $beta->[$i][$t] / $scaling->[$t];
				$Bden[$i] = $Aden[$i];
				$Bnom[$i][$k] += $alpha->[$i][$t] * $beta->[$i][$t] / $scaling->[$t];
			}
			
		}

		my $t = $T - 1;
		my $k = $obs->[$t];

		for (my $i = 0; $i < $N; $i++)
		{
			$Bden[$i] += $alpha->[$i][$t] * $beta->[$i][$t] / $scaling->[$t];
			$Bnom[$i][$k] += $alpha->[$i][$t] * $beta->[$i][$t] / $scaling->[$t];
		}

		#release memory
		for (my $i = 0; $i < $N; $i++)
		{
			$alpha->[$i] = [];
			$beta->[$i] = [];
			$gamma->[$i] = [];
		
		}
		$alpha = [];
		$beta = [];
		$gamma = [];
	}
	
	

	#update parameters
	
	print "update pi ...\n" if $verbose;
	
	for (my $j = 0; $j < $N; $j++)
	{
		$pi2[$j] /= @$multi_obs;
	}
	#@pi2 = @{$pi2[0]} unless ref ($mobs->[0]) eq 'ARRAY';
	
	#my @pi3 = @pi2;
	
	print "update transition matrix ...\n" if $verbose;

	for (my $i = 0; $i < $N; $i++)
	{
		#my $s = 0;
		for (my $j = 0; $j < $N; $j++)
		{
			$A2[$i][$j] = $Anom[$i][$j] / $Aden[$i];
			#$s+= $A2[$i][$j];
		}
		
		#for (my $j = 0; $j < $N; $j++)
		#{
		#	$A2[$i][$j] /= $s;
		#}
	}

	print "update emission matrix ...\n" if $verbose;
	
	for (my $j = 0; $j < $N; $j++)
	{
		#my $s = 0;
		for (my $k = 0; $k < $alphabetSize; $k++)
		{
			$B2[$j][$k] = $Bnom[$j][$k] / $Bden[$j];
			#$s+= $B2[$j][$k];
		}
		
		#for (my $k = 0; $k < $alphabetSize; $k++)
		#{
		#	$B2[$j][$k] /= $s;
		#}
	}

	#release memory
	for (my $i = 0; $i < $N; $i++)
	{
		$Anom[$i] = [];
		$Bnom[$i] = [];
		
	}
	@Anom = ();
	@Bnom = ();
	@Aden = ();
	@Bden = ();
	return {A=>\@A2, B=>\@B2, pi=>\@pi2};
}



1;







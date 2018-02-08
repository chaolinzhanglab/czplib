#
#===============================================================================
#
#         FILE:  SGE.pm
#
#  DESCRIPTION:  interface of SGE
#         BUGS:  ---
#        NOTES:  The package is spun off from the Common.pm
#       AUTHOR:  Chaolin Zhang (cz), czhang@rockefeller.edu
#      COMPANY:  Rockefeller University
#      VERSION:  1.0
#      CREATED:  12/17/10
#     REVISION:  ---
#===============================================================================



package SGE;

require Exporter;


our $VERSION = 1.01;

@ISA = qw (Exporter);

@EXPORT = qw (
	checkSGEJobStatus
	waitUntilSGEJobsDone
);


=head1 NAME

SGE - subroutines to handle Sun Grid Engine (SGE)

=cut

use strict;
use warnings;

use Data::Dumper;
use Carp;


=head2 waitUntilSGEJobsDone
#return the status of unfinished jobs among a specified list

waitUntilSGEJobsDone ($jobIds, $verbose=0, $user="");

$jobIds: job ids to monitor
$verbose:
$user: 


=cut
sub waitUntilSGEJobsDone
{
	my ($jobIds, $verbose, $user) = @_;
	my $total = @$jobIds;
	my $secondSlept = 0;
	
	while (1)
	{
		my $status = checkSGEJobStatus ($jobIds, $user);

		return 1 unless keys %$status > 0;

		my $summary = $status->{'summary'};
		
		my $njobs = 0;
		foreach my $stat (keys %$summary)
		{
			Carp::croak "detect failed jobs: ", Dumper ($status), "\n" unless $stat eq 'r' || $stat eq 't' || $stat eq 'qw';
			$njobs += $summary->{$stat};
		}

		return 1 if $njobs == 0;
		#bug fix. different programs/runs might interfere with each other
		#04/06/2014 Chaolin Zhang

		#TODO: check whether each job finished correctly by running qacct


		#my $n = keys %$status;
		#$n--;
		my $date = `date`;
		chomp $date;
			
		print "$njobs tasks of $total jobs are not finished yet at $date ...\n" if $verbose && $secondSlept % 60 == 0;
		sleep (10); #10 seconds
		$secondSlept += 10;
	}
}


=head2 checkSGEJobStatus

my $status = checkSGEJobStatus ($jobIds, $user);

$status->{$jobId}{$taskId} : is the status of a specific task of job id
$taskId is assigned for array jobs, and otherwise 1

$status->{'summary'}->{$status}: is the number of jobs/tasks of with $status
	
	
=cut

sub checkSGEJobStatus
{
	my ($jobIds, $user) = @_;
	Carp::croak "no job id specified in:", Dumper ($jobIds), "\n" unless @$jobIds > 0;
	
	my %jobHash = map {$_=> 1} @$jobIds;
	my %jobStatus;
	my $cmd = "qstat";
	$cmd .= " -u $user" if $user;

	my @qstat = `$cmd`;
	return {} unless @qstat > 0;
	
	#remove title rows
	shift @qstat; shift @qstat;
	
	my %summary;
	
	foreach my $line (@qstat)
	{
		chomp $line;

		$line=~s/^\s*//;
		my @cols = split (/\s+/, $line);
		
		my $id = $cols[0];
		my $u = $cols[3];
		my $status = $cols[4];
		my $taskId = $cols[$#cols];
	
		#print $taskId, "\n";
		if ($user)
		{
			next unless $u eq $user;
		}
		next unless exists $jobHash {$id};

		my $ntask = 1;
		if ($taskId && $taskId=~/^(\d+)-(\d+):/)
		{
			$ntask = $2-$1 + 1;
		}
		$taskId = 1 unless $taskId;

		$summary{$status} += $ntask;
		$jobStatus{$id}{$taskId} = $status;
	}
		
	$jobStatus{'summary'} = \%summary;
	return \%jobStatus;
}





1;



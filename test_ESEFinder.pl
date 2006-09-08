  use ESEfinder;
  use strict;

  my $seq; # a Bio::PrimarySeqI or Bio::SeqI object

  $seq = Bio::PrimarySeq->new
      (-primery_id => 'test',
       -seq=>'atgcatgctaggtgtgtgttttgtgggttgtactagctagtgat'
       );

  my $ese_finder = Bio::Tools::Analysis::DNA::ESEfinder->
      new(-seq => $seq);

  # run ESEfinder prediction on a DNA sequence
  $ese_finder->run();

  die "Could not get a result"
      unless $ese_finder->status =~ /^COMPLETED/;

	  print $ese_finder->result;      # print raw prediction to STDOUT

  foreach my $feat ( $ese_finder->result('Bio::SeqFeatureI') ) {

      # do something to SeqFeature
      # e.g. print as GFF
      print $feat->gff_string, "\n";
      # or store within the sequence - if it is a Bio::SeqI
	  #$seq->add_SeqFeature($feat)

  }

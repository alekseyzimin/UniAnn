#!/usr/bin/env perl
my @narray=("A","C","G","T");
#initialize code hashes
$n=0;
$n2=0;
$n3=0;
for($i=0;$i<4;$i++){
  $code{$narray[$i]}=$n;
  $n++;
  for($j=0;$j<4;$j++){
    $code2{"$narray[$i]$narray[$j]"}=$n2;
    $n2++;
    for($k=0;$k<4;$k++){
      $code3{"$narray[$i]$narray[$j]$narray[$k]"}=$n3;
      $n3++;
    }
  }
}
my $donor_length;
my $acceptor_length;
my $use_negatives=1;

#we load the genome sequences
open(FILE,$ARGV[0]);
while(my $line=<FILE>){
  chomp($line);
  if($line=~ /^>/){
    if(not($scf eq "")){
      $genome_seqs{$scf}=$seq;
      $seq="";
    }
    my @f=split(/\s+/,$line);
    $scf=substr($f[0],1);
  }else{
    $seq.=$line;
  } 
}   
$genome_seqs{$scf}=$seq if(not($scf eq ""));

#we load SNAP HMMs
if(-e $ARGV[1]){
  open(FILE,$ARGV[1]);
  $line=<FILE>;
  if($line =~ /^zoeHMM/){#check format
    while($line=<FILE>){
      chomp($line);
      if($line=~/^Donor 0HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<4;$j++){
            $donor_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
        $donor_length=$i;
      }elsif($line=~/^Donor 1HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<16;$j++){
            $donor_hmm_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^Donor 2HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line); 
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<64;$j++){
            $donor_hmm2_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^Acceptor 0HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<4;$j++){
            $acceptor_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
        $acceptor_length=$i;
      }elsif($line=~/^Acceptor 1HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<16;$j++){
            $acceptor_hmm_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^Acceptor 2HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<64;$j++){
            $acceptor_hmm2_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^Start/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NNN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<4;$j++){
            $coding_start_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^SDonor/){
        while($line=<FILE>){
          chomp($line);
          my @f=split(/\s+/,$line);
          last if($line eq "NNNNNN");
          $sdonor{$f[0]}=$f[1];
        }
      }elsif($line=~/^SAcceptor/){
        while($line=<FILE>){
          chomp($line);
          my @f=split(/\s+/,$line);
          last if($line eq "NNNNNN");
          $sacceptor{$f[0]}=$f[1];
        }
      }
    }
  }
}

#we load scores for the negative model
if( -e $ARGV[2]){
  $use_negatives=1;
  open(FILE,$ARGV[2]);
  $line=<FILE>;
  if($line =~ /^zoeHMM/){#check format
    while($line=<FILE>){
      chomp($line);
      if($line=~/^Donor 0HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<4;$j++){
            $donor_nfreq[$i][$j]=$f[$j];
          }
          $i++;
        }
        $donor_length=$i;
      }elsif($line=~/^Donor 1HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<16;$j++){
            $donor_hmm_nfreq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^Donor 2HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<64;$j++){
            $donor_hmm2_nfreq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^Acceptor 0HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<4;$j++){
            $acceptor_nfreq[$i][$j]=$f[$j];
          }
          $i++;
        }
        $acceptor_length=$i;
      }elsif($line=~/^Acceptor 1HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<16;$j++){
            $acceptor_hmm_nfreq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^Acceptor 2HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<64;$j++){
            $acceptor_hmm2_nfreq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }
    }
  }
}

#load psauron scores
open(FILE,$ARGV[3]);
$line=<FILE>;
$line=<FILE>;
$line=<FILE>;
$line=<FILE>;
while($line=<FILE>){
  chomp($line);
  my @f=split(/,/,$line);
  $scores{$f[0]}=join(",",@f[9..14]);
}

open(FILEGT,">out.gt.txt");
open(FILEAG,">out.ag.txt");
open(FILEPS,">out.ps.txt");
for my $g(keys %genome_seqs){
  #only doing forward for now!!!
  print "DEBUG scaffold $g $donor_length $acceptor_length\n";
  my $seq_fwd=uc($genome_seqs{$g});
  my $seq_rev=$seq_fwd;
  my @don_fwd_pos=();
  my @acc_fwd_pos=();
  my %acceptor_fwd_hmm2_score=();
  my %don_fwd_hmm2_score=();
  $seq_rev=~tr/ACGTNacgtn/TGCANtgcan/;
  $seq_rev=reverse($seq_rev);
  @psauron_scores=split(/,/,$scores{$g});
  @psauron_frame0=split(/;/,$psauron_scores[0]);
  @psauron_frame1=split(/;/,$psauron_scores[1]);
  @psauron_frame2=split(/;/,$psauron_scores[2]);
  my $mult=30;
  my $off=0.2;
  for(my $i=0;$i<$#psauron_frame0;$i++){
    $psauron_frame0[$i]=log($psauron_frame0[$i]*$mult+1e-6)/log($mult);
    #$psauron_frame0[$i]=($psauron_frame0[$i]-$off)/(1-$off);
  }
  for(my $i=0;$i<$#psauron_frame1;$i++){
    $psauron_frame1[$i]=log($psauron_frame1[$i]*$mult+1e-6)/log($mult);
    #$psauron_frame1[$i]=($psauron_frame1[$i]-$off)/(1-$off);
  }
  for(my $i=0;$i<$#psauron_frame2;$i++){
    $psauron_frame2[$i]=log($psauron_frame2[$i]*$mult+1e-6)/log($mult);
    #$psauron_frame2[$i]=($psauron_frame2[$i]-$off)/(1-$off);
  }
  my $j=0;
  #insert large negative score for an in frame stop 
  for(my $i=0;$i<length($genome_seqs{$g})-3;$i+=3){
    splice @psauron_frame0,$j,0,(-1e6) if(is_stop(substr($genome_seqs{$g},$i,3)));
    splice @psauron_frame1,$j,0,(-1e6) if(is_stop(substr($genome_seqs{$g},$i+1,3)));
    splice @psauron_frame2,$j,0,(-1e6) if(is_stop(substr($genome_seqs{$g},$i+2,3)));
    $j++;
  }
      
  #find donors fwd
  while ($seq_fwd =~ /GT/g) {
    push @don_fwd_pos, pos($seq_fwd) - 2 if(pos($seq_fwd)>$donor_length-5);  # subtract length of "GT" (2) to get start index
  }
  #find acceptors fwd  
  while ($seq_fwd =~ /AG/g) {
    push @acc_fwd_pos, pos($seq_fwd) - 2 if(pos($seq_fwd)>$acceptor_length-3);  # subtract length of "AG" (2) to get start index
  }
  
  for $pos(@don_fwd_pos){
    my $donor_seq=substr($seq_fwd,$pos-3,$donor_length);
    my $donor_hmm2_score=0;
    my $donor_hmm2_nscore=0;
    
    for(my $i=0;$i<($donor_length-2);$i++){
      $donor_hmm2_score+=$donor_hmm2_freq[$i][$code3{substr($donor_seq,$i,3)}] if(defined($code3{substr($donor_seq,$i,3)}));
    }
    $donor_hmm2_score+=$donor_hmm_freq[0][$code2{substr($donor_seq,0,2)}] if(defined($code2{substr($donor_seq,0,2)}));
    if($use_negatives == 0){
      $donor_hmm2_score=$donor_hmm2_score;
      #print "DEBUG $pos $donor_seq $donor_score $donor_hmm_score $donor_hmm2_score\n";
    }else{
      for(my $i=0;$i<($donor_length-2);$i++){
        $donor_hmm2_nscore+=$donor_hmm2_nfreq[$i][$code3{substr($donor_seq,$i,3)}] if(defined($code3{substr($donor_seq,$i,3)}));
      }
      $donor_hmm2_nscore+=$donor_hmm_nfreq[0][$code2{substr($donor_seq,0,2)}] if(defined($code2{substr($donor_seq,0,2)}));

      $donor_hmm2_score=($donor_hmm2_score-$donor_hmm2_nscore);
      #print "$pos $donor_seq $donor_score $donor_hmm_score $donor_hmm2_score\n";
      $don_fwd_hmm2_score{$pos}=$donor_hmm2_score;
    }
  }
  for $pos(@acc_fwd_pos){
    my $acceptor_seq=substr($seq_fwd,$pos-($acceptor_length-5),$acceptor_length);
    my $acceptor_hmm2_score=0;
    my $acceptor_hmm2_nscore=0;
    
    for(my $i=0;$i<($acceptor_length-2);$i++){
      $acceptor_hmm2_score+=$acceptor_hmm2_freq[$i][$code3{substr($acceptor_seq,$i,3)}] if(defined($code3{substr($acceptor_seq,$i,3)}));
    }   
    $acceptor_hmm2_score+=$acceptor_hmm_freq[0][$code2{substr($acceptor_seq,0,2)}] if(defined($code2{substr($acceptor_seq,0,2)}));
    if($use_negatives == 0){
      $acceptor_hmm2_score=$acceptor_hmm2_score/2;
    }else{
      for(my $i=0;$i<($acceptor_length-2);$i++){
        $acceptor_hmm2_nscore+=$acceptor_hmm2_nfreq[$i][$code3{substr($acceptor_seq,$i,3)}] if(defined($code3{substr($acceptor_seq,$i,3)}));
      }
      $acceptor_hmm2_nscore+=$acceptor_hmm_nfreq[0][$code2{substr($acceptor_seq,0,2)}] if(defined($code2{substr($acceptor_seq,0,2)}));

      $acceptor_hmm2_score=($acceptor_hmm2_score-$acceptor_hmm2_nscore)/2;
      #print "$pos $acceptor_seq $acceptor_score $acceptor_hmm_score $acceptor_hmm2_score\n";
      $acceptor_fwd_hmm2_score{$pos}=$acceptor_hmm2_score;
    }
  }
  
  my ($p0,$p1,$p2)=(0,0,0);
  my ($pp0,$pp1,$pp2)=(0,0,0);
  for(my $i=0;$i<length($seq_fwd);$i++){
    ($p0,$p1,$p2)=(0,0,0);
    $p0=$psauron_frame0[int($i/3)] if defined($psauron_frame0[int($i/3)]);
    $p1=$psauron_frame1[int(($i-1)/3)] if ($i>0);
    $p2=$psauron_frame2[int(($i-2)/3)] if ($i>1);
    if($i%3==0 || $i%3==1){
      $p0=$pp0 if($p0==-1e6);
    }
    if($i%3==1 || $i%3==2){
      $p1=$pp1 if($p1==-1e6);
    }
    if($i%3==2 || $i%3==0){
      $p2=$pp2 if($p2==-1e6);
    }

    my $max_score=$p0;
    $max_score=$p1 if($p1>$max_score);
    $max_score=$p2 if($p2>$max_score);
    my $scoreN=0.1-$max_score;
    my $scoreI0=0.1-$max_score;
    my $scoreI1=0.1-$max_score;
    my $scoreI2=0.1-$max_score;

    if($p0==-1e6){#stop in frame 0,1 or 2
      $scoreN=10;
      $scoreI0=-1;
    }
    if($p1==-1e6){#stop in frame 0,1 or 2
      $scoreN=10;
      $scoreI1=-1;
    }
    if($p2==-1e6){#stop in frame 0,1 or 2
      $scoreN=10;
      $scoreI2=-1;
    }
    printf FILEPS "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n",$i,$scoreN,$p0,$p1,$p2,$scoreI0,$scoreI1,$scoreI2,substr($seq_fwd,$i,1);
    $pp0=$p0 if($p0 > -1e6);
    $pp1=$p1 if($p1 > -1e6);
    $pp2=$p2 if($p2 > -1e6);
  }
  my $splice_t=-1.5;
  for(my $i=0;$i<length($seq_fwd);$i++){
    if(defined($don_fwd_hmm2_score{$i})){
      $don_fwd_hmm2_score{$i}=($don_fwd_hmm2_score{$i}>$splice_t)?$don_fwd_hmm2_score{$i}*70:-1e3;
      #print FILEGT "$g $i\t",$don_fwd_hmm2_score{$i},"\n";
      print FILEGT "$i\t",$don_fwd_hmm2_score{$i},"\n";
    }
  }
  for(my $i=0;$i<length($seq_fwd);$i++){
    if(defined($acceptor_fwd_hmm2_score{$i})){
      $acceptor_fwd_hmm2_score{$i}=($acceptor_fwd_hmm2_score{$i}>($splice_t*2))?$acceptor_fwd_hmm2_score{$i}*60:-1e3;
      #print FILEAG "$g $i\t",$acceptor_fwd_hmm2_score{$i},"\n";
      print FILEAG "$i\t",$acceptor_fwd_hmm2_score{$i},"\n";
    }
  }
}

sub is_stop{
  if(uc($_[0]) eq "TAG" || uc($_[0]) eq "TAA" || uc($_[0]) eq "TGA"){
    return 1;
  }else{
    return 0;
  }
}

sub is_start{
  if(uc($_[0]) eq "ATG"){
    return 1;
  }else{
    return 0;
  }
}

sub check_start_stop{
  my $candidate_exon=@_[0];
  my $start_coord=-1;
  my $stop_coord=-1;
  for(my $k=0;$k<length($seq);$k+=3){
    $start_coord=$k if(is_start(substr($candidate_exon,$k,3)) && $start_coord==-1);
    if(is_stop(substr($candidate_exon,$k,3))){
      $stop_coord=$k;
      last;
    }
  }
  return ($start_coord,$stop_coord);
}

sub find_intervals {
    my ($vals, $T, $L, $M) = @_;
    my @intervals;
    my ($start, $end);

    # -------------------------
    # Step 1: find raw intervals above threshold T
    # -------------------------
    for my $i (0 .. $#$vals) {
        if ($vals->[$i] > $T) {
            $start = $i if !defined $start;
            $end   = $i;
        } else {
            if (defined $start) {
                push @intervals, [$start, $end];
                undef $start;
                undef $end;
            }
        }
    }
    push @intervals, [$start, $end] if defined $start;

    # -------------------------
    # Step 2: merge intervals with gaps < L
    # -------------------------
    my @merged;
    for my $iv (@intervals) {
        if (!@merged) {
            push @merged, $iv;
            next;
        }
        my ($ps, $pe) = @{$merged[-1]};
        my ($cs, $ce) = @$iv;

        if ($cs - $pe - 1 < $L) {
            # merge
            $merged[-1]->[1] = $ce;
        } else {
            push @merged, $iv;
        }
    }

    # -------------------------
    # Step 3: filter out intervals shorter than M
    # -------------------------
    my @filtered = grep {
        my ($s, $e) = @$_;
        ($e - $s + 1) >= $M
    } @merged;

    return @filtered;
}



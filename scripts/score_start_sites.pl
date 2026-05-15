#!/usr/bin/env perl
#
my $model_file=$ARGV[0];
my $start_length=0;
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

if(-e $model_file){
  open(FILE,$model_file);
  $line=<FILE>;
  if($line =~ /^zoeHMM/){#check format
    while($line=<FILE>){
      chomp($line);
      if($line=~/^ATG 0HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NNN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<4;$j++){
            $start_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
        $start_length=$i;
      }elsif($line=~/^ATG 1HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NNN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<16;$j++){
            $start_hmm_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^ATG 2HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NNN TRM/);
          chomp($line); 
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<64;$j++){
            $start_hmm2_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }
    }
  }
}

#we load the genome sequences
my $seq="";
while(my $line=<STDIN>){
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
print "DEBUG start pattern length = $start_length\n";
open(FILEATG,">out.atg.txt");
for my $g(keys %genome_seqs){
  #only doing forward for now!!!
  my $seq_fwd=uc($genome_seqs{$g});
  my @start_fwd_pos=();
  #find starts fwd
  while ($seq_fwd =~ /ATG/g) {
    push @start_fwd_pos, pos($seq_fwd) - 3 if(pos($seq_fwd)>$start_length-2);  # subtract length of "ATG" (2) to get start index
  }
  for $pos(@start_fwd_pos){
    #ATG position fixed on the start seq
    my $start_seq=substr($seq_fwd,$pos-19,$start_length);
    my $start_hmm2_score=0;
    #print "DEBUG $start_seq\n";
    for(my $i=0;$i<($start_length-2);$i++){
      $start_hmm2_score+=$start_hmm2_freq[$i][$code3{substr($start_seq,$i,3)}] if(defined($code3{substr($start_seq,$i,3)}));
    }
    $start_hmm2_score+=$start_hmm_freq[0][$code2{substr($start_seq,0,2)}] if(defined($code2{substr($start_seq,0,2)}));
    $start_hmm2_score-=500 if($start_hmm2_score<13);
    print FILEATG "$pos\t",$start_hmm2_score,"\n";
  }
}







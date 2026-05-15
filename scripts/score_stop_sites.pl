#!/usr/bin/env perl
#
my $model_file=$ARGV[0];
my $stop_length=0;
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
            $stop_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
        $stop_length=$i;
      }elsif($line=~/^ATG 1HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NNN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<16;$j++){
            $stop_hmm_freq[$i][$j]=$f[$j];
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
            $stop_hmm2_freq[$i][$j]=$f[$j];
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
print "DEBUG stop pattern length = $stop_length\n";
open(FILESTOP,">out.stop.txt");
for my $g(keys %genome_seqs){
  #only doing forward for now!!!
  my $seq_fwd=uc($genome_seqs{$g});
  my @stop_fwd_pos=();
  #find stops fwd
  while ($seq_fwd =~ /TAA|TAG|TGA/g) {
    push @stop_fwd_pos, pos($seq_fwd) - 3 if(pos($seq_fwd)>$stop_length-2);  # subtract length of "ATG" (2) to get stop index
  }
  for $pos(@stop_fwd_pos){
    #ATG position fixed on the stop seq
    my $stop_seq=substr($seq_fwd,$pos-19,$stop_length);
    my $stop_hmm2_score=0;
    #print "DEBUG $stop_seq\n";
    for(my $i=0;$i<($stop_length-2);$i++){
      $stop_hmm2_score+=$stop_hmm2_freq[$i][$code3{substr($stop_seq,$i,3)}] if(defined($code3{substr($stop_seq,$i,3)}));
    }
    $stop_hmm2_score+=$stop_hmm_freq[0][$code2{substr($stop_seq,0,2)}] if(defined($code2{substr($stop_seq,0,2)}));
    $stop_hmm2_score=-1000 if($stop_hmm2_score<7);
    print FILESTOP "$pos\t",$stop_hmm2_score,"\n";
  }
}







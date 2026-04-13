#!/usr/bin/env perl
#
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

my $fasta=$ARGV[0];
my $gff=$ARGV[1];
my $score_seq=$ARGV[2];
my $seq="";
my $start=-1;
open (FILE,"gffread -W -w /dev/stdout -g $fasta $gff |");
while($line=<FILE>){
  chomp($line);
  if($line =~/^>/){
    if($start>-1){
      push(@seqs,uc(substr($seq,$start-20,25))) if(length(substr($seq,$start-20,25))==25);
    }
    $seq="";
    if($line=~/\sCDS=(\d+)-(\d+)\s/){
      $start=$1;
    }else{
      $start=-1;
    }
  }else{
    $seq.=$line;
  }
}

my $start_length=25;
foreach $start_seq (@seqs){
 #print "Training: $start_seq\n";
 for(my $i=0;$i<$start_length;$i++) {$start_pwm[$i][$code{substr($start_seq,$i,1)}]++ if(defined($code{substr($start_seq,$i,1)}));}
 for(my $i=0;$i<($start_length-1);$i++) {$start2_pwm[$i][$code2{substr($start_seq,$i,2)}]++ if(defined($code2{substr($start_seq,$i,2)}));}
 for(my $i=0;$i<($start_length-2);$i++) {$start3_pwm[$i][$code3{substr($start_seq,$i,3)}]++ if(defined($code3{substr($start_seq,$i,3)}));}
 $w++;
}

my $score_floor_value=1e-10;
for(my $i=0;$i<$start_length;$i++){
  for(my $j=0;$j<4;$j++){
    $start_pwm[$i][$j]=log($start_pwm[$i][$j]/$w*4+$score_floor_value);
  }
}

for(my $i=0;$i<$start_length-1;$i++){
  for(my $j=0;$j<16;$j++){
    $start2_pwm[$i][$j]=log($start2_pwm[$i][$j]/$w*16+$score_floor_value);
  }
} 

for(my $i=0;$i<$start_length-2;$i++){
  for(my $j=0;$j<64;$j++){
    $start3_pwm[$i][$j]=log($start3_pwm[$i][$j]/$w*64+$score_floor_value);
  }
} 

#we load the genome sequences
open(FILE,$score_seq);
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


open(FILEST,">out.atg.txt");
for my $g(keys %genome_seqs){
  #only doing forward for now!!!
  print "DEBUG scaffold $g $start_length $acceptor_length\n";
  my $seq_fwd=uc($genome_seqs{$g});
  my $seq_rev=$seq_fwd;
  my @start_fwd_pos=();
  #find starts fwd
  while ($seq_fwd =~ /ATG/g) {
    push @start_fwd_pos, pos($seq_fwd) - 3 if(pos($seq_fwd)>23);  # subtract length of "ATG" (2) to get start index
  }
  for $pos(@start_fwd_pos){
    my $start_seq=substr($seq_fwd,$pos-19,$start_length);
    my $start_hmm2_score=0;
    my $start_hmm2_nscore=0;
    
    for(my $i=0;$i<($start_length-2);$i++){
      $start_hmm2_score+=$start3_pwm[$i][$code3{substr($start_seq,$i,3)}] if(defined($code3{substr($start_seq,$i,3)}));
    }
    $start_hmm2_score+=$start2_pwm[0][$code2{substr($start_seq,0,2)}] if(defined($code2{substr($start_seq,0,2)}));
    print "$pos $start_seq $start_score $start_hmm_score $start_hmm2_score\n";
    $startfwd_hmm2_score{$pos}=$start_hmm2_score;
  }
}







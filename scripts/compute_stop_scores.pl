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

my $stop_length=30;
my $fasta=$ARGV[0];
my $gff=$ARGV[1];
my $seq="";
my $stop=-1;
open (FILE,"gffread -W -w /dev/stdout -g $fasta $gff |");
while($line=<FILE>){
  chomp($line);
  if($line =~/^>/){
    if($stop>-1){
      push(@seqs,uc(substr($seq,$stop-20,$stop_length))) if(length($seq)>=$stop-20+$stop_length);
    }
    $seq="";
    if($line=~/\sCDS=(\d+)-(\d+)\s/){
      $stop=$2-2;
    }else{
      $stop=-1;
    }
  }else{
    $seq.=$line;
  }
}

foreach $stop_seq (@seqs){
 #print "Training: $stop_seq\n";
 for(my $i=0;$i<$stop_length;$i++) {$stop_pwm[$i][$code{substr($stop_seq,$i,1)}]++ if(defined($code{substr($stop_seq,$i,1)}));}
 for(my $i=0;$i<($stop_length-1);$i++) {$stop2_pwm[$i][$code2{substr($stop_seq,$i,2)}]++ if(defined($code2{substr($stop_seq,$i,2)}));}
 for(my $i=0;$i<($stop_length-2);$i++) {$stop3_pwm[$i][$code3{substr($stop_seq,$i,3)}]++ if(defined($code3{substr($stop_seq,$i,3)}));}
 $w++;
}

my $score_floor_value=1e-10;
print "zoeHMM\n";
print "ATG 0HMM\n";
for(my $i=0;$i<$stop_length;$i++){
  for(my $j=0;$j<4;$j++){
    printf("%.3f ", log($stop_pwm[$i][$j]/$w*4+$score_floor_value));
  }
  print "\n";
}
print "NNN TRM\n";

print "ATG 1HMM\n";
for(my $i=0;$i<$stop_length-1;$i++){
  for(my $j=0;$j<16;$j++){
    printf("%.3f ", log($stop2_pwm[$i][$j]/$w*16+$score_floor_value));
  }
  print "\n";
} 
print "NNN TRM\n";

print "ATG 2HMM\n";
for(my $i=0;$i<$stop_length-2;$i++){
  for(my $j=0;$j<64;$j++){
    printf("%.3f ", log($stop3_pwm[$i][$j]/$w*64+$score_floor_value));
  }
  print "\n";
} 
print "NNN TRM\n";


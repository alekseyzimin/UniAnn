#!/usr/bin/env perl
$intron_stop=-2;
$n=0;
open(FILE,$ARGV[0]);
while($line=<FILE>){
  chomp($line);
  if(not(substr($line,0,1) eq ">")){
    $seq.=uc($line);
  }
}
print STDERR "Loaded sequences\n";
while($seq =~ /TAA|TAG|TGA/g){
  $stops{$-[0]}=1;
}
print STDERR "Found stops\n";
while($line=<STDIN>){
  chomp($line);
  @F=split(/\t/,$line);
  for($i=3;$i<7;$i++){
    $F[$i]=log($F[$i]/(1-$F[$i])+1e-10);
    $F[$i]=1/(1+exp(-$F[$i]/7));
  }
  for($i=0;$i<3;$i++){
    $coding_frame[($n-$i-1)%3]=log(exp(1)*$F[$i+3]*(1-$F[6]+0.15)+1e-10);
  }
  $nc_score=log(exp(1)*$F[6]+1e-10);
  if($stops{($n-2)}){
    $coding_frame[($n-2)%3]=-1000000;
    $intron_score=$intron_stop;
  }else{
    $intron_score=$nc_score;
  }
  printf "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n",$n,$nc_score,$coding_frame[0],$coding_frame[1],$coding_frame[2],$intron_score,substr($seq,$n,1);
  $n++;
}
print STDERR "Done processing scores\n";

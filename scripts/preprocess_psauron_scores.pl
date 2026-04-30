#!/usr/bin/env perl
my $stop_value=-1e6;

#usage on empty
unless(defined($ARGV[0]) && defined($ARGV[1])){
  die "Usage:\npreprocess_psauron_scores.pl genome.fa psauron_score.csv\nMake sure that psauron was run with -a option!\n";
}
  

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

#load psauron scores
open(FILE,$ARGV[1]);
$line=<FILE>;
$line=<FILE>;
$line=<FILE>;
$line=<FILE>;
while($line=<FILE>){
  chomp($line);
  my @f=split(/,/,$line);
  my $g=$f[0];
  $psauron_scores_0f{$g}=$f[9];
  $psauron_scores_1f{$g}=$f[10];
  $psauron_scores_2f{$g}=$f[11];
  $psauron_scores_0r{$g}=$f[12];
  $psauron_scores_1r{$g}=$f[13];
  $psauron_scores_2r{$g}=$f[14];
}
$now=localtime(); 
print "DEBUG $now: psauron scores loaded\n";

open(FILEPS,">out.ps.txt");
for my $g(keys %genome_seqs){
  #only doing forward for now!!!
  $now=localtime();
  print "DEBUG $now: starting scaffold $g\n";
  my $seq_fwd=uc($genome_seqs{$g});
  @psauron_frame0=split(/;/,$psauron_scores_0f{$g});
  @psauron_frame1=split(/;/,$psauron_scores_1f{$g});
  @psauron_frame2=split(/;/,$psauron_scores_2f{$g});
   
  my $mult=30;
  my $lmult=log($mult);
  #my $off=0.2;
  $_ = log($_*$mult+1e-6)/$lmult for @psauron_frame0;
  $_ = log($_*$mult+1e-6)/$lmult for @psauron_frame1;
  $_ = log($_*$mult+1e-6)/$lmult for @psauron_frame2;

  $now=localtime();
  print "DEBUG $now: log transform done for $g\n";
  my $j=0;
  
  #find all stops and insert large negative score for an in frame stop 
  my %stops;
  while($seq_fwd =~ /TAA|TAG|TGA/g){
    $stops{$-[0]}=1;
  }
  $now=localtime();
  print "DEBUG $now: stops located for $g\n";

  my %stops_f0;
  my %stops_f1;
  my %stops_f2;
  my $j=0;
  my $n_stops_f0=0;
  my $n_stops_f1=0;
  my $n_stops_f2=0;
  for(my $i=0;$i<length($seq_fwd)-3;$i+=3){
    if($stops{$i}){
      $stops_f0{$j-$n_stops_f0}++;
      $n_stops_f0++;
    }elsif($stops{$i+1}){
      $stops_f1{$j-$n_stops_f1}++;
      $n_stops_f1++;
    }elsif($stops{$i+2}){
      $stops_f2{$j-$n_stops_f2}++;
      $n_stops_f2++;
    }
    $j++;
  }
  $now=localtime();
  print "DEBUG $now: stops in frames $g\n";
  @psauron_frame0_wstops=insert_before_positions(\@psauron_frame0,\%stops_f0);
  @psauron_frame1_wstops=insert_before_positions(\@psauron_frame1,\%stops_f1);
  @psauron_frame2_wstops=insert_before_positions(\@psauron_frame2,\%stops_f2);
  $now=localtime();
  print "DEBUG $now: stops inserted for $g\n";
  #replace scores by averages between the stops
  #@psauron_frame0_ave=average_between_stops(\@psauron_frame0);
  #@psauron_frame1_ave=average_between_stops(\@psauron_frame1);
  #@psauron_frame2_ave=average_between_stops(\@psauron_frame2);
      
  my ($p0,$p1,$p2)=(0,0,0);
  my ($pp0,$pp1,$pp2)=(0,0,0);
  for(my $i=0;$i<length($seq_fwd);$i++){
    ($p0,$p1,$p2)=(0,0,0);
    $p0=$psauron_frame0_wstops[int($i/3)] if defined($psauron_frame0_wstops[int($i/3)]);
    $p1=$psauron_frame1_wstops[int(($i-1)/3)] if ($i>0);
    $p2=$psauron_frame2_wstops[int(($i-2)/3)] if ($i>1);
    
    if($i%3==0 || $i%3==1){
      $p0=$pp0 if($p0==-1e6);
    }
    if($i%3==1 || $i%3==2){
      $p1=$pp1 if($p1==-1e6);
    }
    if($i%3==2 || $i%3==0){
      $p2=$pp2 if($p2==-1e6);
    }
    
    #my $min_diff=0.5;
    #my @scores_sorted= sort {$a<=>$b} ($p0,$p1,$p2);
    #if($scores_sorted[2]-$scores_sorted[1]>$min_diff && $scores_sorted[2]>0){
    #  $p0+=1-$scores_sorted[2] if($p0>-1e6 && $scores_sorted[2]>0);
    #  $p1+=1-$scores_sorted[2] if($p1>-1e6 && $scores_sorted[2]>0);
    #  $p2+=1-$scores_sorted[2] if($p2>-1e6 && $scores_sorted[2]>0);
    #  $scores_sorted[2]=1;
    #}else{
    #  $p0-=$scores_sorted[1] if($p0>-1e6 && $scores_sorted[1]>0);
    #  $p1-=$scores_sorted[1] if($p1>-1e6 && $scores_sorted[1]>0);
    #  $p2-=$scores_sorted[1] if($p2>-1e6 && $scores_sorted[1]>0);
    #}

    my $max_score=$p0;
    $max_score=$p1 if($p1 > $max_score);
    $max_score=$p2 if($p2 > $max_score);

    my $scoreN=0.1-$max_score;
    my $scoreI0=0.1-$max_score;
    my $scoreI1=0.1-$max_score;
    my $scoreI2=0.1-$max_score;


    if($p0==-1e6){#stop in frame 0,1 or 2
      $scoreN=20;
      $scoreI0=-1;
    }
    if($p1==-1e6){#stop in frame 0,1 or 2
      $scoreN=20;
      $scoreI1=-1;
    }
    if($p2==-1e6){#stop in frame 0,1 or 2
      $scoreN=20;
      $scoreI2=-1;
    }
    printf FILEPS "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n",$i,$scoreN,$p0,$p1,$p2,$scoreI0,$scoreI1,$scoreI2,substr($seq_fwd,$i,1);
    $pp0=$p0 if($p0 > -1e6);
    $pp1=$p1 if($p1 > -1e6);
    $pp2=$p2 if($p2 > -1e6);
  }
  $now=localtime();
  print "DEBUG $now: output complete for $g\n"; 
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

sub average_between_stops {
    my ($arr_ref, $marker) = @_;
    $marker //= $stop_value;   # default marker

    my @arr = @$arr_ref;
    my $n = @arr;

    my $start = undef;   # start index of a block
    my $sum   = 0;
    my $count = 0;

    for (my $i = 0; $i <= $n; $i++) {

        # Case 1: inside a block and we hit a marker or end of array
        if (defined $start && ($i == $n || $arr[$i] == $marker)) {
            my $avg = $count ? $sum / $count : $marker;

            # replace values in the block
            for my $j ($start .. $i-1) {
                $arr[$j] = $avg;
            }

            # reset block
            undef $start;
            $sum = 0;
            $count = 0;
        }

        # Case 2: start of a new block
        if ($i < $n && $arr[$i] != $marker) {
            if (!defined $start) {
                $start = $i;
            }
            $sum += $arr[$i];
            $count++;
        }
    }

    return @arr;
}

sub insert_before_positions {
    my ($arr_ref, $insert_ref) = @_;
    # $insert_ref: position => value OR position => [values...]

    my @arr = @$arr_ref;
    my %ins = %$insert_ref;

    my @out;
    for my $i (0 .. $#arr) {

        if(exists $ins{$i}) {
          for(my $j=0;$j<$ins{$i};$j++){
            push @out, $stop_value;
          }
        }

        push @out, $arr[$i];
    }

    return @out;
}


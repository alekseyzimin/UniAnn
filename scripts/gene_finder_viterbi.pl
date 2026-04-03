#!/usr/bin/env perl
use strict;
use warnings;
use List::Util qw(max);

#------------------------------------------------------------
# States:
# 0: N   (noncoding)
# 1: E0  2: E1  3: E2  (exon frames)
# 4: I0  5: I1  6: I2  (intron frames)
#------------------------------------------------------------

my @states = (0..6);
my %state_name = (
    0 => 'N',
    1 => 'E0', 2 => 'E1', 3 => 'E2',
    4 => 'I0', 5 => 'I1', 6 => 'I2',
);

my $MIN_INTRON = 30;   # minimum intron length in bases
my $MIN_EXON = 40;   # minimum exon length in bases

#------------------------------------------------------------
# Command line:
#   script.pl seq.fasta emissions.txt gt.txt ag.txt
#
# emissions.txt:
#   pos \t N \t E0 \t E1 \t E2 \t I0 \t I1 \t I2
#
# gt.txt, ag.txt:
#   pos \t log_score
#   (pos = coordinate of G in GT, or A in AG; 0-based)
#------------------------------------------------------------

my ($f_fasta, $f_emit, $f_gt, $f_ag) = @ARGV;
die "Usage: $0 seq.fasta emissions.txt gt.txt ag.txt\n" unless $f_ag;

my $seq = read_fasta($f_fasta);
my @seq = split //, $seq;
my $L   = @seq;

my @emit     = load_emissions($f_emit, $L);
my @gt_score = load_sparse_scores($f_gt, $L);
my @ag_score = load_sparse_scores($f_ag, $L);
#------------------------------------------------------------
# Safety check: verify GT/AG coordinates match the sequence
#------------------------------------------------------------

for my $pos (0 .. $L-2) {

    # Check GT
    if ($gt_score[$pos] > -1e8) {   # means a real score exists
        my $dinuc = $seq[$pos] . $seq[$pos+1];
        if ($dinuc ne 'GT') {
            warn "WARNING: GT score at position $pos but sequence has $dinuc\n";
        }
    }

    # Check AG
    if ($ag_score[$pos] > -1e8) {
        my $dinuc = $seq[$pos] . $seq[$pos+1];
        if ($dinuc ne 'AG') {
            warn "WARNING: AG score at position $pos but sequence has $dinuc\n";
        }
    }
}


#------------------------------------------------------------
# Baseline transition matrix (log-space), excluding GT/AG overrides
#------------------------------------------------------------

my %trans;
foreach my $s (@states) {
    foreach my $t (@states) {
        $trans{$s}{$t} = -1e9; # impossible by default
    }
}

my $max_log_prob=1;

# Noncoding self
$trans{0}{0} = log(1);

# Exon frame cycling: none
$trans{1}{1} = $max_log_prob;
$trans{2}{2} = $max_log_prob;
$trans{3}{3} = $max_log_prob;

# Intron states: no cycling, only self-transition
$trans{4}{4} = $max_log_prob;   # I0 -> I0
$trans{5}{5} = $max_log_prob;   # I1 -> I1
$trans{6}{6} = $max_log_prob;   # I2 -> I2


#------------------------------------------------------------
# Viterbi with GT/AG splice overrides, stop codons, min intron length,
# and coordinate-dependent emissions
#------------------------------------------------------------

my (@dp, @bt, @intron_len, @exon_len);

# Initialization at position 0 — force start in E0
for my $s (@states) {
    my $start_prob = ($s == 1) ? 0.0 : -1e9;  # start in E0
    my $e = $emit[0][$s];
    $dp[0][$s] = $start_prob + $e;
    $bt[0][$s] = -1;
    $intron_len[0][$s] = is_intron($s) ? 1 : 0;
    $exon_len[0][$s] = is_exon($s) ? 1 : 0;
}

for (my $i = 1; $i < $L; $i++) {
    my $b_prev = $seq[$i-1];
    my $b      = $seq[$i];

    for my $to (@states) {
        my $emit_log = $emit[$i][$to];
        #fix for TAG stop shere AG is acceptor
        $emit_log=$emit[$i+1][$to] if($emit[$i][$to]<=-1e6 && $b_prev eq 'A' && $b eq 'G');
        my $best = -1e18;
        my $best_from = -1;

        for my $from (@states) {
            my $log_t = $trans{$from}{$to};

            #--------------------------------------------------------
            # Exon → Intron only at GT AND only if exon length ≥ MIN_EXON
            #--------------------------------------------------------
            if (is_exon($from) && is_intron($to) &&
                $b_prev eq 'G' && $b eq 'T') {
              my $len = $exon_len[$i-1][$from]-2;
              print STDERR "DEBUG at $i trying transition $state_name{$from} $state_name{$to} score $dp[$i-1][$from] emit $emit_log length $len\n";
              $log_t = -1e9;
              #my $mod=$i%3;
              #$mod=2-$mod if($mod>0);
              $log_t = $gt_score[$i-1] if ($to-4 == $from-1 && $len >= $MIN_EXON);
              
              print STDERR "DEBUG probability $log_t\n";
            }

            #--------------------------------------------------------
            # Intron -> exon only at AG AND only if intron length ≥ MIN_INTRON
            #--------------------------------------------------------
            if (is_intron($from) && is_exon($to) &&
                $b_prev eq 'A' && $b eq 'G') {

              my $len = $intron_len[$i-1][$from]+2;
# Only allow transition if intron length is compatible with frame transition
              print STDERR "DEBUG at $i trying transition $state_name{$from} $state_name{$to} score $dp[$i-1][$from] emit $emit_log length $len\n";
              $log_t = -1e9;
              if ($len >= $MIN_INTRON) {
                if(($state_name{$from} eq 'I0' && $state_name{$to} eq 'E0') || ($state_name{$from} eq 'I1' && $state_name{$to} eq 'E1') || ($state_name{$from} eq 'I2' && $state_name{$to} eq 'E2')){
#intron length must be divisible by 3
                  $log_t = $ag_score[$i-1] if ( $len%3 == 0 );
                } elsif(($state_name{$from} eq 'I0' && $state_name{$to} eq 'E1') || ($state_name{$from} eq 'I1' && $state_name{$to} eq 'E2') || ($state_name{$from} eq 'I2' && $state_name{$to} eq 'E0')){
                  $log_t = $ag_score[$i-1] if ( $len%3 == 1 );
                } elsif(($state_name{$from} eq 'I0' && $state_name{$to} eq 'E2') || ($state_name{$from} eq 'I1' && $state_name{$to} eq 'E0') || ($state_name{$from} eq 'I2' && $state_name{$to} eq 'E1')){
                  $log_t = $ag_score[$i-1] if ( $len%3 == 2 );
                }
              }
              print STDERR "DEBUG probability $log_t\n";
            }

            #--------------------------------------------------------
            # Exon -> Noncoding after STOP codon (TAA, TAG, TGA), frame-correct
            #--------------------------------------------------------
            if (is_exon($from) && $to == 0 && $i >= 2) {
                my $codon = $seq[$i-2] . $seq[$i-1] . $seq[$i];
                if ($codon eq 'TAA' || $codon eq 'TAG' || $codon eq 'TGA') {
                  print STDERR "DEBUG at $i trying transition $state_name{$from} $state_name{$to} score $dp[$i-1][$from] emission $emit_log\n";
                    my $frame = $from - 1;   # E0=1→0, E1=2→1, E2=3→2
                    if ( (($i-2) % 3) == $frame ) {
                        # force exon -> N at stop codon
                        $log_t = $max_log_prob;
                    }
                  print STDERR "DEBUG probability $log_t\n";
                }
            }

            #--------------------------------------------------------
            # Noncoding ->exon after START codon (ATG)
            #--------------------------------------------------------
            if (is_exon($to) && $from == 0 && $i >= 2) {
                my $codon = $seq[$i-2] . $seq[$i-1] . $seq[$i];
                if ($codon eq 'ATG') {
                    print STDERR "DEBUG at $i trying transition $state_name{$from} $state_name{$to} score $dp[$i-1][$from] emission $emit_log\n";
                    my $frame = $to-1;   # E0=1→0, E1=2→1, E2=3→2
                    if ( (($i-2) % 3) == $frame ) {
                        # allow N -> exon at start codon
                        $log_t = $max_log_prob; 
                    }
                print STDERR "DEBUG probability $log_t\n";
                }
            }

            my $cand = $dp[$i-1][$from] + $log_t + $emit_log;
            if ($cand > $best) {
                $best = $cand;
                $best_from = $from;
            }
        }

        $dp[$i][$to] = $best;
        $bt[$i][$to] = $best_from;

        # Track intron length
        if (is_intron($to)) {
            if ($best_from >= 0 && is_intron($best_from)) {
                $intron_len[$i][$to] = $intron_len[$i-1][$best_from] + 1;
            } else {
                $intron_len[$i][$to] = 1;
            }
        } else {
            $intron_len[$i][$to] = 0;
        }
        # Track exon length
        if (is_exon($to)) {
          if ($best_from >= 0 && is_exon($best_from)) {
            $exon_len[$i][$to] = $exon_len[$i-1][$best_from] + 1;
          } else {
            $exon_len[$i][$to] = 1;
          }
        } else {
          $exon_len[$i][$to] = 0;
        }
    }
}

# Termination
my $best_final = -1e18;
my $best_state = -1;
for my $s (@states) {
    if ($dp[$L-1][$s] > $best_final) {
        $best_final = $dp[$L-1][$s];
        $best_state = $s;
    }
}

# Backtrace
my @path_states;
my $cur = $best_state;
for (my $i = $L-1; $i >= 0; $i--) {
    $path_states[$i] = $cur;
    $cur = $bt[$i][$cur];
    last if $cur < 0;
}

my @path_labels = map { $state_name{$_} } @path_states;

#------------------------------------------------------------
# GFF3 writer
#------------------------------------------------------------

my $seqid = get_fasta_header($f_fasta);

print "##gff-version 3\n";

my $current_state = $path_labels[0];
my $start = 0;
my $index=0;

for my $i (1..$#path_labels) {
    if ($path_labels[$i] ne $current_state) {
        write_gff_feature($seqid, $current_state, $start, $i-1);
        $current_state = $path_labels[$i];
        $start = $i;
    }
}
write_gff_feature($seqid, $current_state, $start, $#path_labels);

#print "##FASTA\n";
#print ">$seqid\n$seq\n";

#------------------------------------------------------------
# Helpers
#------------------------------------------------------------

sub is_exon {
    my ($s) = @_;
    return ($s == 1 || $s == 2 || $s == 3);
}

sub is_intron {
    my ($s) = @_;
    return ($s == 4 || $s == 5 || $s == 6);
}

sub read_fasta {
    my ($file) = @_;
    open my $fh, "<", $file or die "Cannot open FASTA $file: $!";
    my $seq = '';
    while (<$fh>) {
        chomp;
        next if /^>/;
        $seq .= $_;
    }
    close $fh;
    $seq =~ s/\s+//g;
    return uc($seq);
}

sub load_emissions {
    my ($file, $L) = @_;
    my @emit;
    open my $fh, "<", $file or die "Cannot open emissions $file: $!";
    while (<$fh>) {
        next if /^#/;
        chomp;
        next unless length;
        my ($pos, @vals) = split /\t/;
        die "Need 7 emission values per line in $file\n"
            unless @vals >= 7;
        for my $s (0..6) {
            $emit[$pos][$s] = $vals[$s];
        }
    }
    close $fh;

    for my $i (0..$L-1) {
        for my $s (0..6) {
            $emit[$i][$s] = -1e9 unless defined $emit[$i][$s];
        }
    }
    return @emit;
}

sub load_sparse_scores {
    my ($file, $L) = @_;
    my @scores = ((-1e9) x $L);  # default: impossible

    open my $fh, "<", $file or die "Cannot open scores $file: $!";
    while (<$fh>) {
        next if /^#/;
        chomp;
        next unless length;
        my ($pos, $score) = split /\t/;
        $scores[$pos] = $score + 0.0;
    }
    close $fh;

    return @scores;
}

sub get_fasta_header {
    my ($file) = @_;
    open my $fh, "<", $file or die $!;
    while (<$fh>) {
        if (/^>(\S+)/) {
            close $fh;
            return $1;
        }
    }
    close $fh;
    return "sequence";
}

sub write_gff_feature {
    my ($seqid, $state, $start0, $end0) = @_;
    my @f=split(/\//,$f_fasta);
    my $tname=$f[-1];
    my $start = $start0;  # 0-based
    my $end   = $end0 ;
    $start0-- if($start0 == 0);

    my $type;
    if ($state =~ /^E/) {
        $type = "CDS";
        $start = $start0+2;
    } elsif ($state =~ /^I/) {
        $end = $end0 + 2;
        $type = "intron";
    } else {
        $type = "region";   # noncoding
        $end = $end0 + 2;
        $index++; #global
    }
    $best_final=int($best_final*100)/100;
    print join("\t",
        $seqid,
        "HMM",
        $type,
        $start,
        $end,
        $best_final,
        ".",
        ".",
        "Parent=$tname.$index;state=$state"
    ), "\n";
}


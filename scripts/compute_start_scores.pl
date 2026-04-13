#!/usr/bin/perl
use strict;
use warnings;

my ($fasta, $gff) = @ARGV;
die "Usage: $0 genome.fasta annotation.gff\n" unless $gff;

# -----------------------------
# Load genome FASTA
# -----------------------------
my %seq;
my $id;

open FA, "<", $fasta or die $!;
while (<FA>) {
    chomp;
    if (/^>(\S+)/) {
        $id = $1;
        $seq{$id} = "";
    } else {
        $seq{$id} .= uc($_);
    }
}
close FA;

# -----------------------------
# Reverse complement
# -----------------------------
sub revcomp {
    my $s = reverse($_[0]);
    $s =~ tr/ACGTN/TGCAN/;
    return $s;
}

# -----------------------------
# Parse GFF: collect exons and CDS per transcript
# -----------------------------
my %tx;   # $tx{id}{chr,strand,exons=>[],cds=>[]}

open GFF, "<", $gff or die $!;
while (<GFF>) {
    next if /^#/;
    chomp;
    my @f = split /\t/;
    next unless @f >= 9;

    my ($chr,$src,$type,$start,$end,$score,$strand,$phase,$attr) = @f;

    my ($idtag) = $attr =~ /ID=([^;]+)/;
    my ($parent) = $attr =~ /Parent=([^;]+)/;

    # mRNA/transcript defines the container
    if ($type eq "mRNA" || $type eq "transcript") {
        $tx{$idtag}{chr}    = $chr;
        $tx{$idtag}{strand} = $strand;
    }

    # Exons
    if ($type eq "exon" && $parent) {
        push @{$tx{$parent}{exons}}, [$start, $end];
    }

    # CDS
    if ($type eq "CDS" && $parent) {
        push @{$tx{$parent}{cds}}, [$start, $end];
    }
}
close GFF;

# -----------------------------
# Build transcript sequences
# -----------------------------
foreach my $tid (keys %tx) {

    my $chr    = $tx{$tid}{chr};
    my $strand = $tx{$tid}{strand};
    my @exons  = sort { $a->[0] <=> $b->[0] } @{$tx{$tid}{exons}};

    next unless @exons;

    # Concatenate exon sequences
    my $tseq = "";
    for my $e (@exons) {
        my ($s,$e2) = @$e;
        $tseq .= substr($seq{$chr}, $s-1, $e2-$s+1);
    }

    # Reverse complement if needed
    $tseq = revcomp($tseq) if $strand eq "-";

    # -----------------------------
    # Find translation start site
    # -----------------------------
    my @cds = sort { $a->[0] <=> $b->[0] } @{$tx{$tid}{cds}};
    next unless @cds;

    # First CDS genomic coordinate
    my ($cds_start, $cds_end) =
        ($strand eq "+") ? @{$cds[0]} : @{$cds[-1]};

    # Convert genomic CDS start to transcript coordinate
    my $tpos = 0;
    my $cds_tpos;

    for my $e (@exons) {
        my ($s,$e2) = @$e;
        my $elen = $e2 - $s + 1;

        if ($strand eq "+") {
            if ($cds_start >= $s && $cds_start <= $e2) {
                $cds_tpos = $tpos + ($cds_start - $s);
                last;
            }
        } else {
            if ($cds_end >= $s && $cds_end <= $e2) {
                $cds_tpos = $tpos + ($e2 - $cds_end);
                last;
            }
        }

        $tpos += $elen;
    }

    next unless defined $cds_tpos;

    # -----------------------------
    # Extract ±20 bp around start
    # -----------------------------
    my $win_start = $cds_tpos - 20;
    my $win_end   = $cds_tpos + 20;

    $win_start = 0 if $win_start < 0;
    $win_end   = length($tseq)-1 if $win_end >= length($tseq);

    my $subseq = substr($tseq, $win_start, $win_end - $win_start + 1);

    print ">$tid|start=$cds_tpos|len=", length($subseq), "\n";
    print "$subseq\n";
}


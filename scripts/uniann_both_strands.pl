#!/usr/bin/env perl
use strict;
use warnings;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);
use File::Spec;
use File::Temp qw(tempdir tempfile);

# uniann.sh calls this helper when -a / --all-prob is enabled.
# Prepare both orientations for the existing forward pipeline and combine their GFFs.
# The decoder and its scoring rules remain unchanged.
my ($wrapper, $fasta, $psauron, $sites, $mult, @flags) = @ARGV;
# Resolve paths before child processes change directory.
# Derive the output name from the FASTA pathname supplied by the caller.
$fasta = File::Spec->rel2abs($fasta);
($wrapper, $psauron, $sites) = map { abs_path($_) // die "Cannot resolve $_\n" }
    ($wrapper, $psauron, $sites);
die "Multiplier must be between 1 and 100\n"
    unless defined($mult) && $mult =~ /\A\d+(?:\.\d*)?(?:[eE][+-]?\d+)?\z/
        && $mult >= 1 && $mult <= 100;

# Read the FASTA and create its reverse complement.
# Both-strand mode accepts one FASTA sequence.
# Keep its header so both orientations use the same PSAURON and site-score ID.
open my $fh, '<', $fasta or die "Cannot read $fasta: $!\n";
my ($header, $seq) = ('', '');
while (<$fh>) {
    if (/^>/) {
        die "Both-strand mode requires a single FASTA sequence\n" if length $header;
        $header = $_;
    } else {
        s/\s+//g;
        die "FASTA sequence before header\n" if length($_) && !length($header);
        $seq .= $_;
    }
}
close $fh;
die "FASTA sequence is empty\n" unless length $seq;
my ($seqid) = $header =~ /^>(\S+)/;
die "FASTA header is empty\n" unless defined $seqid;
my $length = length $seq;
# Reverse base order and complement all IUPAC symbols, including lowercase.
my $reverse = reverse $seq;
$reverse =~ tr/ACGTRYMKBDHVSWNacgtrymkbdhvswn/TGCAYRKMVHDBSWNtgcayrkmvhdbswn/;

# Prepare the PSAURON frame columns for each orientation.
# Use the original PSAURON CSV for the plus run.
# For the minus run, adapt the matching sequence's frame columns for the existing emission preprocessor.
open $fh, '<', $psauron or die "Cannot read $psauron: $!\n";
my @csv = <$fh>;
close $fh;
my @reverse_csv = @csv;
my $matched = 0;
# PSAURON's first four lines contain metadata and the column header.
for my $i (4 .. $#csv) {
    my $line = $csv[$i];
    $line =~ s/[\r\n]+\z//;
    my @cells = split /,/, $line, -1;
    next unless $cells[0] eq $seqid;
    die "Duplicate PSAURON scores for $seqid\n" if $matched++;
    for my $column (12 .. 14) {
        die "Missing reverse_all_prob scores for $seqid (run PSAURON with -a)\n"
            unless defined($cells[$column]) && length($cells[$column]);
        for my $prob (split /;/, $cells[$column], -1) {
            die "Invalid reverse_all_prob probability for $seqid\n"
                unless $prob =~ /\A(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\z/
                    && $prob >= 0 && $prob <= 1;
        }
    }
    # Zero-based columns 9..11 are forward frames; 12..14 are reverse frames.
    # Reverse arrays already describe rc[0:], rc[1:], rc[2:], in that order.
    # Copy the arrays without reversing values or permuting frames.
    # The existing preprocessor aligns them to the reverse-complemented FASTA.
    @cells[9 .. 11] = @cells[12 .. 14];
    $reverse_csv[$i] = join(',', @cells) . "\n";
}
die "Missing PSAURON scores for $seqid\n" unless $matched;

# Split site scores by strand and map minus positions to the reverse complement.
# Site rows contain: sequence ID, 1-based position, strand, type, motif, probability.
# A seventh probability column is also accepted.
# Positions identify the first motif base in the indicated direction, using original FASTA coordinates.
open $fh, '<', $sites or die "Cannot read $sites: $!\n";
my (%rows, %has_donor);
while (<$fh>) {
    next if /^\s*(?:#|$)/;
    my @f = split;
    next unless @f >= 6 && ($f[2] eq '+' || $f[2] eq '-');
    next unless $f[0] eq $seqid;
    die "Site position outside FASTA sequence\n"
        unless $f[1] =~ /\A\d+\z/ && $f[1] >= 1 && $f[1] <= $length;
    my $strand = $f[2];
    # Mirror uniann.sh's probability-column selection for its scaling factor.
    # A donor must have a positive log(probability * multiplier + 1e-10) for that factor.
    # Track this per strand so strands without a qualifying donor can be skipped.
    my $factor_prob = @f > 6 ? $f[6] : $f[5];
    $has_donor{$strand} = 1 if $f[3] eq 'donor'
        && $factor_prob =~ /\A(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\z/
        && $factor_prob <= 1 && $factor_prob * $mult + 1e-10 > 1;
    # The forward pipeline only accepts '+' rows.
    # Reflect minus positions into reverse-complement coordinates and relabel their strand.
    # Keep the site types and scores unchanged.
    $f[1] = $length - $f[1] + 1 if $strand eq '-';
    $f[2] = '+';
    push @{$rows{$strand}}, join("\t", @f) . "\n";
}
close $fh;

# Run the existing forward pipeline independently for each strand.
# uniann.sh writes fixed out.* filenames and reuses an existing out.ps.txt.
# Separate directories prevent either strand from reusing the other's scores.
my $work = tempdir('uniann-both-XXXXXX', TMPDIR => 1, CLEANUP => 1);
my @combined;
for my $strand ('+', '-') {
    unless ($has_donor{$strand}) {
        warn "Skipping $strand strand: no donor with a positive rescaled score\n";
        next;
    }
    my $dir = "$work/" . ($strand eq '+' ? 'plus' : 'minus');
    mkdir $dir or die "Cannot create $dir: $!\n";
    my $input = "$dir/" . basename($fasta);
    write_file($input, $header, $strand eq '+' ? $seq : $reverse, "\n");
    write_file("$dir/psauron.csv", $strand eq '+' ? @csv : @reverse_csv);
    write_file("$dir/sites.txt", @{$rows{$strand}});
    # Change directory in a child so the parent retains its original context.
    # Forward -n and -v, but not -a, which would invoke this helper recursively.
    my $pid = fork();
    die "Cannot fork: $!\n" unless defined $pid;
    if (!$pid) {
        chdir $dir or die "Cannot enter $dir: $!\n";
        exec 'bash', '-o', 'pipefail', $wrapper, '-f', $input,
            '-p', "$dir/psauron.csv", '-s', "$dir/sites.txt", '-m', $mult, @flags;
        die "Cannot run UniAnn: $!\n";
    }
    # Run strands sequentially and require a completed GFF before merging.
    waitpid($pid, 0);
    die "UniAnn failed on $strand strand\n" if $? || !-f "$input.uniann.gff";
    # Map reverse GFF coordinates back and update IDs and Parent references.
    my $suffix = $strand eq '+' ? 'f' : 'r';
    open my $gff, '<', "$input.uniann.gff" or die "Cannot read child GFF: $!\n";
    while (<$gff>) {
        next if /^#/ || /^\s*$/;
        chomp;
        my @f = split /\t/, $_, -1;
        die "Invalid UniAnn GFF row\n" unless @f == 9;
        if ($strand eq '-') {
            # Reflect the 1-based inclusive interval, swapping its endpoints.
            # Other feature fields, including CDS phase, stay with that feature.
            @f[3, 4] = ($length - $f[4] + 1, $length - $f[3] + 1);
            $f[6] = '-';
        }
        # Both runs share an ID namespace.
        # Append f or r to IDs and every member of Parent lists to keep links consistent.
        $f[8] =~ s/(\A|;)(ID|Parent)=([^;]+)/$1 . $2 . '=' . join(',', map { $_ . $suffix } split \/,\/, $3)/eg;
        push @combined, join("\t", @f) . "\n";
    }
    close $gff;
}
# Write the combined GFF after all eligible strands succeed.
# Stage beside the final output so rename is atomic and a failed run leaves the previous GFF intact.
my ($out, $pending) = tempfile('.uniann-gff-XXXXXX', DIR => dirname($fasta), UNLINK => 1);
print {$out} "##gff-version 3\n", @combined or die "Cannot write combined GFF: $!\n";
close $out or die "Cannot close combined GFF: $!\n";
# File::Temp starts with private permissions; apply the usual output permissions before publishing.
chmod(0666 & ~umask, $pending) or die "Cannot set combined GFF permissions: $!\n";
rename $pending, "$fasta.uniann.gff" or die "Cannot publish combined GFF: $!\n";
print "Output gff file is $fasta.uniann.gff\n";

sub write_file {
    my ($path, @content) = @_;
    open my $out, '>', $path or die "Cannot write $path: $!\n";
    print {$out} @content or die "Cannot write $path: $!\n";
    close $out or die "Cannot close $path: $!\n";
}

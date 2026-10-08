"""Small deterministic inputs covering UniAnn's transition boundaries."""
from pathlib import Path
import random


def cases():
    single = 'C' * 30 + 'ATG' + 'AAA' * 40 + 'TAA' + 'C' * 35
    yield 'single_exon', single, False
    yield 'early_start', single[30:], False
    yield 'short_single_exon', 'C' * 30 + 'ATG' + 'AAA' * 8 + 'TAA' + 'C' * 40, False
    yield 'unknown_bases', 'N' * 40 + single + 'N' * 30, False
    yield 'all_unknown', 'N' * 100, False
    for length in (37, 38, 39, 40, 41, 42):
        intron = 'GT' + 'C' * (length - 4) + 'AG'
        seq = 'C' * 30 + 'ATG' + 'AAA' * 20 + intron + 'AAA' * 40 + 'TAA' + 'C' * 35
        yield f'intron_{length}', seq, False
    yield 'lowercase', seq.lower(), False
    yield 'tag_acceptor_overlap', 'C' * 30 + 'ATG' + 'AAA' * 20 + 'GT' + 'C' * 36 + 'TAG' + 'AAA' * 40 + 'TAA', False
    for seed in range(6):
        rng = random.Random(seed)
        yield f'random_{seed}', ''.join(rng.choices('ACGTN', k=200)), True


def write_inputs(directory, sequence, random_scores=False):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    sequence_upper = sequence.upper()
    (directory / 'sequence.fa').write_text('>fixture\n' + sequence + '\n')
    rng = random.Random(712)
    with (directory / 'emissions.txt').open('w') as out:
        for position, base in enumerate(sequence_upper):
            if base == 'N':
                scores = [0, -1e6, -1e6, -1e6, 0]
            elif random_scores:
                scores = [rng.choice((-3, -1, 0, 0.25, 1, 4)) for _ in range(5)]
            else:
                scores = [0.1, 1.2, 0.7, 0.3, -0.2]
            out.write(str(position) + '\t' + '\t'.join(map(str, scores)) + '\n')
    motifs = {'gt': ('GT',), 'ag': ('AG',), 'atg': ('ATG',), 'stop': ('TAA', 'TAG', 'TGA')}
    for kind, alternatives in motifs.items():
        with (directory / (kind + '.txt')).open('w') as out:
            for position in range(len(sequence)):
                if any(sequence_upper.startswith(motif, position) for motif in alternatives):
                    score = rng.choice((-30, 0, 20, 300)) if random_scores else 300
                    out.write(f'{position}\t{score}\n')
    return [directory / name for name in ('sequence.fa', 'emissions.txt', 'gt.txt', 'ag.txt', 'atg.txt', 'stop.txt')]

"""Check wrapper flag forwarding without invoking external gene-model tools."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


INPUT_ARGS = ['sequence.fa', 'out.ps.txt', 'out.gt.txt', 'out.ag.txt', 'out.atg.txt', 'out.stop.txt']


def decoder_args(wrapper_flags):
    """Run the wrapper with stub tools and return the arguments the decoder saw."""
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        shutil.copyfile(ROOT / 'scripts' / 'uniann.sh', work / 'uniann.sh')
        (work / 'uniann').write_text('#!/bin/sh\nprintf "%s\\n" "$@" > decoder.args\n')
        (work / 'gffread').write_text('#!/bin/sh\ncat\n')
        for name in ('uniann', 'gffread'):
            (work / name).chmod(0o755)
        (work / 'sequence.fa').write_text('>fixture\nGT\n')
        (work / 'psauron.csv').write_text('unused: emission input is already prepared\n')
        (work / 'out.ps.txt').write_text('0 0 0 0 0 0\n')
        (work / 'sites.tsv').write_text('fixture\t1\t+\tdonor\tGT\t0.99\n')
        env = dict(os.environ, PATH=str(work) + os.pathsep + os.environ['PATH'])
        subprocess.run(['bash', 'uniann.sh', *wrapper_flags,
                        '-f', 'sequence.fa', '-p', 'psauron.csv', '-s', 'sites.tsv'],
                       cwd=work, env=env, capture_output=True, check=True, start_new_session=True)
        return (work / 'decoder.args').read_text().splitlines()


def both_strand_results(wrapper_flags):
    """Exercise the real both-strand helper and child wrappers with stub tools."""
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        for name in ('uniann.sh', 'uniann_both_strands.pl'):
            shutil.copyfile(ROOT / 'scripts' / name, work / name)
        (work / 'preprocess_psauron_scores.pl').write_text(
            '#!/bin/sh\nprintf "0 0 0 0 0 0\\n" > out.ps.txt\n')
        (work / 'uniann').write_text(
            '#!/usr/bin/env python3\n'
            'import json, os, pathlib, sys\n'
            'args = sys.argv[1:]\n'
            'sequence = pathlib.Path(args[0]).read_text().splitlines()[1]\n'
            'with open(os.environ["DECODER_LOG"], "a") as out:\n'
            '    out.write(json.dumps({"args": args, "sequence": sequence}) + "\\n")\n'
            'print("fixture\\tstub\\tgene\\t1\\t3\\t.\\t+\\t.\\tID=gene1")\n')
        (work / 'gffread').write_text('#!/bin/sh\ncat\n')
        for name in ('uniann', 'gffread', 'preprocess_psauron_scores.pl'):
            (work / name).chmod(0o755)
        (work / 'sequence.fa').write_text('>fixture\nACGTAA\n')
        # The helper validates reverse probabilities and prepares each orientation.
        cells = ['fixture'] + ['0'] * 8 + ['0.9'] * 6
        (work / 'psauron.csv').write_text('metadata\n' * 4 + ','.join(cells) + '\n')
        (work / 'sites.tsv').write_text(
            'fixture\t1\t+\tdonor\tGT\t0.99\n'
            'fixture\t6\t-\tdonor\tGT\t0.99\n')
        env = dict(os.environ, PATH=str(work) + os.pathsep + os.environ['PATH'],
                   DECODER_LOG=str(work / 'decoder.jsonl'))
        subprocess.run(['bash', 'uniann.sh', '-a', *wrapper_flags,
                        '-f', 'sequence.fa', '-p', 'psauron.csv', '-s', 'sites.tsv'],
                       cwd=work, env=env, capture_output=True, check=True, start_new_session=True)
        calls = [json.loads(line) for line in (work / 'decoder.jsonl').read_text().splitlines()]
        return calls, (work / 'sequence.fa.uniann.gff').read_text().splitlines()


class WrapperTests(unittest.TestCase):
    def test_noviterbi_skips_dump_and_keeps_following_flag(self):
        self.assertEqual(decoder_args(['-n', '-m', '2']), INPUT_ARGS + ['--no-dp-dump'])

    def test_default_keeps_dump(self):
        self.assertEqual(decoder_args([]), INPUT_ARGS)

    def test_both_strands_forward_dump_flag_to_both_decoders(self):
        for flags, expected in (([], []), (['-n', '-m', '2'], ['--no-dp-dump'])):
            with self.subTest(flags=flags):
                calls, gff = both_strand_results(flags)
                self.assertEqual([call['sequence'] for call in calls], ['ACGTAA', 'TTACGT'])
                for call in calls:
                    self.assertEqual(call['args'][1:6], INPUT_ARGS[1:])
                    self.assertEqual(call['args'][6:], expected)
                self.assertEqual(gff, [
                    '##gff-version 3',
                    'fixture\tstub\tgene\t1\t3\t.\t+\t.\tID=gene1f',
                    'fixture\tstub\tgene\t4\t6\t.\t-\t.\tID=gene1r',
                ])


if __name__ == '__main__':
    unittest.main()

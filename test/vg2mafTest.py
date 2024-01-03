import sys
import os
import unittest
import subprocess

class Vg2mafTest(unittest.TestCase):
    def setUp(self):
        if os.system('vg 2> /dev/null') != 256:
            raise RuntimeError('vg not found in PATH')
        if os.system('../vg2maf 2> /dev/null') != 256:
            raise RuntimeError('../vg2maf not found')

        self.tiny_vg = 'tiny.vg'

    def reverse_comp(s):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return "".join([complement[c] for c in s.reversed()])
    
    def maf2path(self, maf_name):
        """ return a map of pathname->sequence. """
        out_map = {}
        s_lines = []
        with open(maf_name, 'r') as maf_file:
            for line in maf_file:
                if line.startswith('s'):
                    s_lines.append(line.rstrip().split())

        # nothing below makes sense unless rows are sorted by path position
        s_lines = sorted(s_lines, key=lambda x : (x[1], int(x[2])))
        
        for toks in s_lines:
            name = toks[1]
            strand = toks[4]
            seq_len = toks[5]
            length = int(toks[3])
            seq = toks[6].replace('-', '')
            self.assertEqual(len(seq), length)
            self.assertTrue(strand in ['+', '-'])
            if strand == '-':
                seq = Vg2mafTest.reverse_comp(seq)
            if name not in out_map:
                out_map[name] = ''
            out_map[name] += seq
        return out_map

    def vg2path(self, vg_name):
        """ return a map of pathname->sequence"""
        paths_name = 'vg-paths.txt'
        subprocess.check_call('vg paths -x {} -F > {}'.format(vg_name, paths_name), shell=True)
        out_map = {}
        name = None
        with open(paths_name, 'r') as paths_file:
            for line in paths_file:
                if line.startswith('>'):
                    name = line.rstrip()[1:]
                    assert name not in out_map
                    out_map[name] = ''
                else:
                    assert name
                    out_map[name] += line.rstrip()
        return out_map
                                        
    def test_simple_forward(self):
        """ do a real easy conversion """
        subprocess.check_call(['vg', 'index', self.tiny_vg, '-j', 'tiny.dist'])
        out_maf_name = 'tiny.maf'
        with open(out_maf_name, 'w') as out_maf_file:
            subprocess.check_call(['vg2maf', self.tiny_vg, '-d', 'tiny.dist', '-r', 'x'], stdout=out_maf_file)

        vg_paths = self.vg2path(self.tiny_vg)
        maf_paths = self.maf2path(out_maf_name)

        self.assertEqual(len(vg_paths), len(maf_paths))
        for path_name, path_seq in vg_paths.items():
            self.assertEqual(maf_paths[path_name], path_seq)

        

if __name__ == '__main__':
    unittest.main()

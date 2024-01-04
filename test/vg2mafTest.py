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
        return "".join([complement[c] for c in reversed(s)])
    
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

    def flip_path(self, vg_name, path_name, out_vg_name, out_path_name):
        """ make a graph where a backwards version of path_name is added as out_path_name """
        path_gaf_name = path_name + '.gaf'
        with open(path_gaf_name, 'w') as path_gaf_file:
            subprocess.check_call(['vg', 'paths', '-x', self.tiny_vg, '-A', '-Q', path_name], stdout=path_gaf_file)
        out_path_gaf_name = out_path_name + '.gaf'
        with open(path_gaf_name, 'r') as path_gaf_file, open(out_path_gaf_name, 'w') as out_path_gaf_file:
            for line in path_gaf_file:
                toks = line.rstrip().split()
                # rename the name (not checking for degenerate cases)
                toks[0] = toks[0].replace(path_name, out_path_name)
                path = toks[5]
                steps = []
                step = ""
                for c in path:
                    if c in ['>', '<']:
                        if step:
                            steps.append(step)
                        # flip the sign
                        step = '>' if c == '<' else '<'
                    else:
                        step += c
                if step:
                    steps.append(step)
                # flip the steps
                toks[5] = "".join(reversed(steps))
                out_path_gaf_file.write('\t'.join(toks) + '\n')

        # add the flipped path back to the graph
        with open(out_vg_name, 'w') as out_vg_file:
            subprocess.check_call(['vg', 'augment', vg_name, out_path_gaf_name, '-F', '-B'], stdout=out_vg_file)

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

    def test_simple_reverse(self):
        """ do an easy conversion on reverse strand """
        rev_vg_name = 'tiny_rev.vg'
        self.flip_path(self.tiny_vg, 'x', rev_vg_name, 'x_rev')

        subprocess.check_call(['vg', 'index', self.tiny_vg, '-j', 'tiny.dist'])
        out_maf_name = 'tiny_rev.maf'
        with open(out_maf_name, 'w') as out_maf_file:
            subprocess.check_call(['vg2maf', rev_vg_name, '-d', 'tiny.dist', '-r', 'x'], stdout=out_maf_file)

        vg_paths = self.vg2path(rev_vg_name)
        maf_paths = self.maf2path(out_maf_name)

        self.assertEqual(len(vg_paths), len(maf_paths))
        for path_name, path_seq in vg_paths.items():
            self.assertEqual(maf_paths[path_name], path_seq)

    def test_snp_forward(self):
        """ test a single forward strand snp """
        # Manually make this alignment
        # x CAAATAAG
        # y --AAGA--
        snp_gam_json_name = 'snp.gam.json'
        with open(snp_gam_json_name, 'w') as snp_gam_json_file:
            snp_gam_json_file.write('{"name": "snp", "path": {"mapping": [{"position": {"node_id": "1", "offset": "2"}, "edit": [{"from_length": 2, "to_length": 2}, {"from_length": 1, "to_length": 1, "sequence": "G"}, {"from_length": 1, "to_length": 1}]}]}, "sequence": "AAGA"}\n')
        # convert to gam and index it
        snp_gam_name = 'snp.gam'
        with open(snp_gam_name, 'w') as snp_gam_file:
            subprocess.check_call(['vg', 'view', '-JaG', snp_gam_json_name], stdout=snp_gam_file)
        subprocess.check_call(['vg', 'gamsort', snp_gam_name, '-i', snp_gam_name + '.gai'], stdout=subprocess.DEVNULL)

        # convert vg+gam to maf
        subprocess.check_call(['vg', 'index', self.tiny_vg, '-j', 'tiny.dist'])
        out_maf_name = 'snp.maf'
        with open(out_maf_name, 'w') as out_maf_file:
            subprocess.check_call(['vg2maf', self.tiny_vg, '-d', 'tiny.dist', '-r', 'x', '-g', snp_gam_name], stdout=out_maf_file)

        lines_by_offset = {}
        with open(out_maf_name, 'r') as out_maf_file:
            for line in out_maf_file:
                if line.startswith('s'):
                    toks = line.rstrip().split()
                    offset = int(toks[2])
                    if offset not in lines_by_offset:
                        lines_by_offset[offset] = []
                    lines_by_offset[offset].append(toks)

        # there are 10 nodes in the graph
        self.assertEqual(len(lines_by_offset), 10)

        # there should be 2 lines for node 1
        self.assertEqual(len(lines_by_offset[0]), 2)

        self.assertEqual(['s', 'x', '0', '8', '+', '50', 'CAAATAAG'], lines_by_offset[0][0])
        self.assertEqual(['s', 'snp', '0', '4', '+', '4', '--AAGA--'],lines_by_offset[0][1])

    def test_deletion_forward(self):
        """ test a single forward strand deletion """
        # Manually make this alignment
        # x CAAATAAG
        # y -AA--A--
        snp_gam_json_name = 'del.gam.json'
        with open(snp_gam_json_name, 'w') as snp_gam_json_file:
            snp_gam_json_file.write('{"name": "del", "path": {"mapping": [{"position": {"node_id": "1", "offset": "1"}, "edit": [{"from_length": 2, "to_length": 2}, {"from_length": 2, "to_length": 0}, {"from_length": 1, "to_length": 1}]}]}, "sequence": "AAA"}\n')
        # convert to gam and index it
        snp_gam_name = 'del.gam'
        with open(snp_gam_name, 'w') as snp_gam_file:
            subprocess.check_call(['vg', 'view', '-JaG', snp_gam_json_name], stdout=snp_gam_file)
        subprocess.check_call(['vg', 'gamsort', snp_gam_name, '-i', snp_gam_name + '.gai'], stdout=subprocess.DEVNULL)

        # convert vg+gam to maf
        subprocess.check_call(['vg', 'index', self.tiny_vg, '-j', 'tiny.dist'])
        out_maf_name = 'del.maf'
        with open(out_maf_name, 'w') as out_maf_file:
            subprocess.check_call(['vg2maf', self.tiny_vg, '-d', 'tiny.dist', '-r', 'x', '-g', snp_gam_name], stdout=out_maf_file)

        lines_by_offset = {}
        with open(out_maf_name, 'r') as out_maf_file:
            for line in out_maf_file:
                if line.startswith('s'):
                    toks = line.rstrip().split()
                    offset = int(toks[2])
                    if offset not in lines_by_offset:
                        lines_by_offset[offset] = []
                    lines_by_offset[offset].append(toks)

        # there are 10 nodes in the graph
        self.assertEqual(len(lines_by_offset), 10)

        # there should be 2 lines for node 1
        self.assertEqual(len(lines_by_offset[0]), 2)

        self.assertEqual(['s', 'x', '0', '8', '+', '50', 'CAAATAAG'], lines_by_offset[0][0])
        self.assertEqual(['s', 'del', '0', '3', '+', '3', '-AA--A--'],lines_by_offset[0][1])

    def test_insertion_forward(self):
        """ test a single forward strand insertion """
        # Manually make this alignment
        # x CAA---ATAAG
        # y CAAGGGATAAG
        snp_gam_json_name = 'ins.gam.json'
        with open(snp_gam_json_name, 'w') as snp_gam_json_file:
            snp_gam_json_file.write('{"name": "ins", "path": {"mapping": [{"position": {"node_id": "1", "offset": "0"}, "edit": [{"from_length": 3, "to_length": 3}, {"from_length": 0, "to_length": 3, "sequence": "GGG"}, {"from_length": 5, "to_length": 5}]}]}, "sequence": "CAAGGGATAAG"}\n')
        # convert to gam and index it
        snp_gam_name = 'ins.gam'
        with open(snp_gam_name, 'w') as snp_gam_file:
            subprocess.check_call(['vg', 'view', '-JaG', snp_gam_json_name], stdout=snp_gam_file)
        subprocess.check_call(['vg', 'gamsort', snp_gam_name, '-i', snp_gam_name + '.gai'], stdout=subprocess.DEVNULL)

        # convert vg+gam to maf
        subprocess.check_call(['vg', 'index', self.tiny_vg, '-j', 'tiny.dist'])
        out_maf_name = 'ins.maf'
        with open(out_maf_name, 'w') as out_maf_file:
            subprocess.check_call(['vg2maf', self.tiny_vg, '-d', 'tiny.dist', '-r', 'x', '-g', snp_gam_name], stdout=out_maf_file)

        lines_by_offset = {}
        with open(out_maf_name, 'r') as out_maf_file:
            for line in out_maf_file:
                if line.startswith('s'):
                    toks = line.rstrip().split()
                    offset = int(toks[2])
                    if offset not in lines_by_offset:
                        lines_by_offset[offset] = []
                    lines_by_offset[offset].append(toks)

        # there are 10 nodes in the graph
        self.assertEqual(len(lines_by_offset), 10)

        # there should be 2 lines for node 1
        self.assertEqual(len(lines_by_offset[0]), 2)

        self.assertEqual(['s', 'x', '0', '8', '+', '50', 'CAA---ATAAG'], lines_by_offset[0][0])
        self.assertEqual(['s', 'ins', '0', '11', '+', '11', 'CAAGGGATAAG'],lines_by_offset[0][1])


    def test_multi_insertion_forward(self):
        """ test a triple forward strand insertion """
        # Manually make this alignment (ie all insertions unaligned)
        # x CAA--------ATAAG
        # y CAAGGG-----ATAAG
        # y CAA---CC---ATAAG
        # y CAA-----GGGATAAG        
        snp_gam_json_name = 'ins3.gam.json'
        with open(snp_gam_json_name, 'w') as snp_gam_json_file:
            snp_gam_json_file.write('{"name": "ins1", "path": {"mapping": [{"position": {"node_id": "1", "offset": "0"}, "edit": [{"from_length": 3, "to_length": 3}, {"from_length": 0, "to_length": 3, "sequence": "GGG"}, {"from_length": 5, "to_length": 5}]}]}, "sequence": "CAAAGGGTAAG"}\n')
            snp_gam_json_file.write('{"name": "ins2", "path": {"mapping": [{"position": {"node_id": "1", "offset": "0"}, "edit": [{"from_length": 3, "to_length": 3}, {"from_length": 0, "to_length": 2, "sequence": "CC"}, {"from_length": 5, "to_length": 5}]}]}, "sequence": "CAAACCTAAG"}\n')
            snp_gam_json_file.write('{"name": "ins3", "path": {"mapping": [{"position": {"node_id": "1", "offset": "0"}, "edit": [{"from_length": 3, "to_length": 3}, {"from_length": 0, "to_length": 3, "sequence": "GGG"}, {"from_length": 5, "to_length": 5}]}]}, "sequence": "CAAAGGGTAAG"}\n')
        # convert to gam and index it
        snp_gam_name = 'ins3.gam'
        with open(snp_gam_name, 'w') as snp_gam_file:
            subprocess.check_call(['vg', 'view', '-JaG', snp_gam_json_name], stdout=snp_gam_file)
        subprocess.check_call(['vg', 'gamsort', snp_gam_name, '-i', snp_gam_name + '.gai'], stdout=subprocess.DEVNULL)

        # convert vg+gam to maf
        subprocess.check_call(['vg', 'index', self.tiny_vg, '-j', 'tiny.dist'])
        out_maf_name = 'ins3.maf'
        with open(out_maf_name, 'w') as out_maf_file:
            subprocess.check_call(['vg2maf', self.tiny_vg, '-d', 'tiny.dist', '-r', 'x', '-g', snp_gam_name], stdout=out_maf_file)

        lines_by_offset = {}
        with open(out_maf_name, 'r') as out_maf_file:
            for line in out_maf_file:
                if line.startswith('s'):
                    toks = line.rstrip().split()
                    offset = int(toks[2])
                    if offset not in lines_by_offset:
                        lines_by_offset[offset] = []
                    lines_by_offset[offset].append(toks)

        # there are 10 nodes in the graph
        self.assertEqual(len(lines_by_offset), 10)

        # there should be 2 lines for node 1
        self.assertEqual(len(lines_by_offset[0]), 4)

        self.assertEqual(['s', 'x', '0', '8', '+', '50', 'CAA--------ATAAG'], lines_by_offset[0][0])
        self.assertEqual(['s', 'ins1', '0', '11', '+', '11', 'CAAGGG-----ATAAG'], lines_by_offset[0][1])
        self.assertEqual(['s', 'ins2', '0', '10', '+', '10', 'CAA---CC---ATAAG'], lines_by_offset[0][2])
        self.assertEqual(['s', 'ins3', '0', '11', '+', '11', 'CAA-----GGGATAAG'], lines_by_offset[0][3])

        
        

if __name__ == '__main__':
    unittest.main()

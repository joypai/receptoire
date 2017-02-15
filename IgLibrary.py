import re
from operator import itemgetter
from collections import Counter
from itertools import groupby
from Bio.Seq import Seq
from Bio import SeqIO, BiopythonWarning, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SeqRecord import SeqRecord
import warnings
warnings.simplefilter('ignore', BiopythonWarning)

REGIONS = ['FR1-IMGT', 'CDR1-IMGT', 'FR2-IMGT', 'CDR2-IMGT', 'FR3-IMGT', 'CDR3-IMGT']

class AntibodyLibrary:
    def __init__(self, name, igblast_file, changeo_file, chain):
        self.name = name
        self.chain = chain
        self.igblast_file = igblast_file
        self.changeo_file = changeo_file
        self.entries = self._create_sequences()       # list of AntibodySequence objects
        self.heavy_seqs, self.light_seqs = self._split_heavy_light()

    def _create_sequences(self):
        """ create AntibodySequence object for each sequence """
        igblast_output = open(self.igblast_file).readlines()
        chunks = [ list(y) for x,y in groupby(igblast_output, lambda x: not re.match('# IGBLASTN', x)) ]
        chunks = list(filter(lambda x: not x[0].startswith('# IGBLASTN') and len(x) > 0,  chunks))
        entries = [ AntibodySequence(info, self.chain) for info in chunks if any("# Alignment" in s for s in info) ]

        return entries

    def _split_heavy_light(self):
        heavy, light = [], []

        for seq in self.entries:
            if seq.chain_type == 'VH': heavy.append(seq)
            else: light.append(seq)

        return heavy, light

    def report_entries(self):
        print('heavy chain sequences: ')
        for entry in self.heavy_seqs:
            print(entry.get_id())
        print('light chain sequences: ')
        for entry in self.light_seqs:
            print(entry.get_id())

    def access_entries(self):
        for entry in self.entries:
            print(entry.get_id())
            entry.display()



class ClonalAntibodyLibrary(AntibodyLibrary):
    def __init__(self, name, igblast_file, changeo_file, chain):
        AntibodyLibrary.__init__(self, name, igblast_file, changeo_file, chain)
        self.clones = self._make_clones()
        self.functional_clones = self._filter_functional()

    def _make_clones(self):
        clone_info = [ itemgetter(*[0,45])(entry.strip().split('\t')) for entry in open(self.changeo_file).readlines() ]

        groups = {}
        for x,y in clone_info:
            groups.setdefault(y,[]).append(x)

        clone_db = { x:y for x,y in groups.items() if x != "CLONE" }        # skip header entry
        clones, sizes = [], []

        for clone_id, seq_ids in clone_db.items():
            # make clone
            clone = Clone(clone_id, seq_ids)
            clone.add_seqs([ entry for entry in self.entries if entry.entry_id in seq_ids ])
            clones.append(clone)
            sizes.append(len(seq_ids))

        sorted_clones = [ x for (y,x) in sorted(zip(sizes, clones), reverse=True) ]

        return sorted_clones

    def _filter_functional(self):
        clones_filtered, sizes = [], []
        for c in self.clones:
            seq_ids = [ s.entry_id for s in c.sequences if s.prod == "Yes" ]
            if len(seq_ids) == 0: continue
            clone = Clone(c.clone_id, seq_ids)
            clone.add_seqs([ entry for entry in self.entries if entry.entry_id in seq_ids ])
            clones_filtered.append(clone)
            sizes.append(len(seq_ids))

        sorted_clones_filtered = [ x for (y,x) in sorted(zip(sizes, clones_filtered), reverse=True) ]
        return sorted_clones_filtered





"""
class ReferenceAntibodyLibrary(AntibodyLibrary):
    def __init__(self, name, igblast_file, changeo_file):
        AntibodyLibrary.__init__(self, name, igblast_file, changeo_file)
        self.reference;
        self.sec_ref;
"""


class Clone:
    def __init__(self, clone_id, sequence_ids):
        self.clone_id = clone_id
        self.sequence_ids = sequence_ids          # list of AntibodySequence objects

    def add_seqs(self, seq_objs):
        self.sequences = seq_objs

        # modify AntibodySequence object to include clone information
        for s in seq_objs:
            s.assign_clone(self.clone_id)

    def output_clone(self):
        print("Clone %s:" % self.clone_id)
        print("\n".join([s for s in self.sequence_ids]))

    def __lt__(self, other):
        return len(self.sequence_ids) < len(other.sequence_ids)

    def __gt__(self, other):
        return other.__lt__(self)

    def __eq__(self, other):
        return self.clone_id == other.clone_id

    def __ne__(self, other):
        return not self.__eq__(other)

    def print_vdjs(self):
        print("Clone %s:" % self.clone_id)
        print("seq_id\tV\tD\tJ")

        for seq in self.sequences:
            print("%s\t%s" % (seq.entry_id, "\t".join(seq.extract_vdj())))

    def get_vdjs(self):
        vdjs = {}
        vdj_count = { 'V': [] , 'D': [], 'J': [] }

        for seq in self.sequences:
            vdj = seq.extract_vdj()
            vdjs[seq.entry_id] = vdj

            vdj_count['V'].append(vdj[0])
            vdj_count['D'].append(vdj[1])
            vdj_count['J'].append(vdj[2])

        # get consensus germline vdj
        germline_vdj = []
        for r in ['V','D','J']:
            count = Counter(vdj_count[r])
            consensus = count.most_common()[0][0]
            germline_vdj.append(consensus)

        vdjs['germline'] = tuple(germline_vdj)

        return vdjs

    def output_region_fastas(self):
        filenames = []

        # get germline sequence
        germline = self.get_germline()
        for region in REGIONS:
            nt_records = [ SeqRecord(Seq(germline[region][0]), id='germline', description=('| %s' % region)) ]
            aa_records = [ SeqRecord(Seq(germline[region][1].replace('*','Z')), id='germline', description=('| %s' % region)) ]
            nt_records += [ SeqRecord(Seq(s.regions[region][0]), id=s.get_id(), description='| '+region) for s in self.sequences ]
            aa_records += [ SeqRecord(Seq(s.regions[region][1].replace('*','Z')), id=s.get_id(),
                                  description='| '+region) for s in self.sequences ]

            file_basename = self.clone_id+"_"+region
            with open(file_basename+"_nt.fa",'w') as out_nt:
                SeqIO.write(nt_records, out_nt, "fasta")
            with open(file_basename+"_aa.fa",'w') as out_aa:
                SeqIO.write(aa_records, out_aa, "fasta")

            filenames += [file_basename+"_nt.fa", file_basename+"_aa.fa"]

        return filenames

    def get_germline(self):
        germlines = { x: {} for x in REGIONS }
        for s in self.sequences:
            gl = s.get_germline(True)
            for r in gl:
                if gl[r] in germlines[r]: germlines[r][gl[r]] += 1
                else: germlines[r][gl[r]] = 1

        consensus_gl = {}
        for r in germlines:
            consensus_gl[r] = max(germlines[r].items(), key=itemgetter(1))[0]

        return consensus_gl


class AntibodySequence:

    def __init__(self, info, chain):
        self.entry_id = info[0].strip().split('# Query: ')[1]
        self._fill_features(info[4:], chain)

        self.regions, self.germline_regions = self._extract_regions()
        self.nt_mismatches, self.aa_mismatches = self._calculate_mismatches()

    def _fill_features(self, info, chain):
        """ Parse IgBlast output and extract information about sequence and regions """
        sections = [ list(y) for x,y in groupby(info, lambda x: not re.match('#', x)) ]

        for heading, content in zip(*[iter(sections)] * 2):
            if heading[0].startswith('# V-(D)-J rearrangement'):
                #if chain == 'IgH':
                try:
                    self.top_v, self.top_d, self.top_j, self.chain_type, self.stop_codon, self.inframe, self.prod, _ = content[0].split('\t')
                except ValueError:
                    self.top_v, self.top_j, self.chain_type, self.stop_codon, self.inframe, self.prod, _ = content[0].split('\t')
                    self.top_d = 'None'
            elif heading[0].startswith('# V-(D)-J junction'):
                junction_parts = [ x for x in content[0].split('\t')[1:len(content[0].split('\t'))-2]
                                if x != 'N/A' and x.find('(') == -1 ]
                self.junction = "".join(junction_parts)
            elif heading[0].startswith('# Sub-region sequence details'):
                self.cdr3_nt, self.cdr3_aa = content[0].strip().split('\t')[1:]
            elif heading[0].startswith('# Alignment summary'):
                self.v_insertions = 0
                self.v_deletions = 0
                self.details = {}
                self.first_region = None

                for region_content in content[:-2]:
                    region_content = region_content.strip()
                    region, start_idx, end, length, matches, mismatches, num_gaps, pc = [ int(x) if i not in [0,7] else x
                                                    for i,x in enumerate(region_content.split('\t')) ]
                    region = region.split(' ')[0]       # simplify CDR3-IMGT (germline) to CDR3-IMGT

                    if self.first_region is None and length >= 3:        # check beginning frame
                        self.start = start_idx
                        self.first_region = region
                        self.overflow = length % 3

                    dels = length - (end-start_idx+1)
                    self.v_insertions += num_gaps - dels
                    self.v_deletions += dels
                    if length < 3 and region != "CDR3-IMGT":
                        self.details[region] = ['NA'] * 6
                    else:
                        self.details[region] = [start_idx-self.start, end-self.start, length, matches, mismatches, num_gaps]

            elif heading[0].startswith('# Hit table'):
                self.seq = ''
                self.germseq = ''
                seen = set()

                for hit in content:
                    if hit[0] not in seen and hit[0] in ['V','J']:
                        seq, germseq = hit.split('\t')[14:16]
                        self.seq += seq
                        self.germseq += germseq

                        if hit[0] == 'V': # add junction / (D) sequence
                            self.seq += self.junction
                            self.germseq += 'N' * len(self.junction)
                        seen.add(hit[0])

        if not hasattr(self,'details'):
            print(self.entry_id, "has no details")
            self.details = { r:['NA']*6 for r in REGIONS }
        # fill empty regions
        for r in set(REGIONS)-set(self.details.keys()):
            self.details[r] = ['NA','NA','NA','NA','NA','NA']

        """
        # update overflow
        start, end, length = self.details[self.first_region][0:3]
        num_deletions = self.seq[start:start+length].count('-')
        num_insertions = self.germseq[start:start+length].count('-')
        #self.overflow = self.overflow - abs(num_deletions - num_insertions) % 3
        if self.entry_id == 'ApoLZ_Exp138_Plate4_IgH_21-GMr':
            print('overflow before:', self.overflow)
            print("indels:",num_deletions, num_insertions)
        """

        return self

    def _extract_regions(self):
        region_seqs = {}
        region_germseqs = {}
        insertions, deletions = 0, 0
        for region, details in sorted(self.details.items(), key=lambda i:REGIONS.index(i[0])):
            start, end, length = details[:3]

            if length == "NA":
                region_seqs[region] = ('N'*6, 'X'*2)
                region_germseqs[region] = ('N'*6, 'X'*2)
                continue

            if region == self.first_region:
                full_aa_seq = str(Seq(self.seq[self.overflow:].replace('-','')).translate())
                full_aa_germseq = str(Seq(self.germseq[self.overflow:].replace('-','')).translate())
                self.aa_seq = full_aa_seq
                self.aa_germseq = full_aa_germseq
                seq_aln = self.seq[self.overflow:length]
                deletions += seq_aln.count('-')
                seq = seq_aln.replace('-','')
                germseq = self.germseq[self.overflow:length]
                insertions += germseq.count('-')
                germseq = germseq.replace('-','')

            elif re.match('CDR3', region):
                if hasattr(self, 'cdr3_nt'):
                    seq = self.cdr3_nt
                    germseq = self.germseq[start+deletions:start+deletions+len(self.cdr3_nt)].replace('-','')
                else:
                    seq_aln = self.seq[start+deletions:]
                    seq = seq_aln.replace('-','')
                    germseq = self.germseq[start+deletions:].replace('-','')
            else:
                seq_aln = self.seq[start+deletions:start+deletions+length]
                germseq = self.germseq[start+deletions:start+deletions+length]
                deletions += seq_aln.count('-')
                insertions += germseq.count('-')
                seq = seq_aln.replace('-','')
                germseq = germseq.replace('-','')

            aa_seq, full_aa_seq = full_aa_seq[:round(len(seq)/3)], full_aa_seq[round(len(seq)/3):]
            aa_germseq, full_aa_germseq = full_aa_germseq[:round(len(germseq)/3)], full_aa_germseq[round(len(germseq)/3):]

            region_seqs[region] = (seq,aa_seq)
            region_germseqs[region] = (germseq,aa_germseq)

        return region_seqs, region_germseqs

    def get_germline(self, regions=False):
        if regions:
            return self.germline_regions
        else:
            return (self.germseq, self.aa_germseq)

    def assign_clone(self, clone_id):
        self.clone = clone_id

    def extract_vdj(self):
        return self.top_v, self.top_d, self.top_j


    def _calculate_mismatches(self):
        """ determine number of mismatches in regions up to CDR3 """

        if hasattr(self,'cdr3_nt'):         # short sequence, no CDR3, no need to truncate
            seq = self.seq[self.overflow:].split(self.cdr3_nt)[0]
            germseq = self.germseq[self.overflow:self.overflow+len(seq)]
            aa_seq = self.aa_seq.split(self.cdr3_aa)[0]
            aa_germseq = self.aa_germseq[:len(aa_seq)]
        elif self.regions['CDR3-IMGT'][1] != 'XX':
            cdr3_nt_seq, cdr3_aa_seq = self.regions['CDR3-IMGT']
            cdr3_nt_germseq, cdr3_aa_germseq = self.germline_regions['CDR3-IMGT']
            if self.entry_id == 'plate3_-01-49-2mRK':
                print(self.entry_id)
                print('cdr3s:',cdr3_nt_seq, cdr3_aa_seq)
                print('cdr3s germ:',cdr3_nt_germseq, cdr3_aa_germseq)
                print('germseq', self.germseq[self.overflow:])
            seq = self.seq[self.overflow:].split(cdr3_nt_seq)[0]
            germseq = self.germseq[self.overflow:].split(cdr3_nt_germseq)[0]
            aa_seq = self.aa_seq.split(cdr3_aa_seq)[0]
            aa_germseq = self.aa_germseq.split(cdr3_aa_germseq)[0]
        else:       # remove CDR3
            seq = self.seq[self.overflow:]
            germseq = self.germseq[self.overflow:]
            aa_seq = self.aa_seq
            aa_germseq = self.aa_germseq

        assert len(seq) == len(germseq), 'nt lengths not equal\n'+self.entry_id+"\n"+seq+"\n"+germseq+"\n"+str(len(seq))+" "+str(len(germseq))
        nt_mismatches = sum(1 for x,y in zip(seq, germseq) if x != y and y != "N")

        # aa
        if len(aa_seq) == len(aa_germseq):
            aa_mismatches = sum(1 for x,y in zip(aa_seq, aa_germseq) if x != y and y != "X")
        else:
            aln = pairwise2.align.globalxs(aa_seq, aa_germseq, -10, -.1)
            a1, a2, _ , _, _ = aln[0]
            aa_mismatches = sum(1 for x,y in zip(a1, a2) if x != y and y != "X")

          #  print(self.entry_id,"\n",a1,"\n",a2,sep="")

        return nt_mismatches, aa_mismatches

    def get_id(self):
        return self.entry_id

    def display(self):
        print("seq %s" % self.seq)
        print("germline seq %s" % self.germseq)
        print("overflow %d" % self.overflow)
        if hasattr(self, 'cdr3_aa'): print("CDR3 %s" % self.cdr3_aa)



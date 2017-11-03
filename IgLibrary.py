import sys
import re
from operator import itemgetter
from collections import Counter
from itertools import groupby
from Bio.Seq import Seq
from Bio import SeqIO, BiopythonWarning, pairwise2
from Bio.SeqRecord import SeqRecord
import warnings
warnings.simplefilter('ignore', BiopythonWarning)

REGIONS = ['FR1-IMGT', 'CDR1-IMGT', 'FR2-IMGT', 'CDR2-IMGT', 'FR3-IMGT', 'CDR3-IMGT']
V_fasta = { 'human': '/rugpfs/fs0/nuss_lab/scratch/jpai/software/ncbi-igblast-1.6.1/IMGT_Human_IGV.fasta',
           'mouse': '/rugpfs/fs0/nuss_lab/scratch/jpai/software/ncbi-igblast-1.6.1/IMGT_Mouse_IGV.fasta' }
V_tcr_fasta = { 'human': '/rugpfs/fs0/nuss_lab/scratch/jpai/software/ncbi-igblast-1.6.1/IMGT_Human_TRBV.fasta',
           'mouse': '/rugpfs/fs0/nuss_lab/scratch/jpai/software/ncbi-igblast-1.6.1/IMGT_Mouse_TRBV.fasta' }

class AntibodyLibrary:
    def __init__(self, name, igblast_file, chain, organism, seq_type):
        self.name = name
        self.chain = chain
        self.seq_type = seq_type
        self.igblast_file = igblast_file
        self.entries = self._create_sequences()       # list of AntibodySequence objects
        self.heavy_seqs, self.light_seqs = self._split_heavy_light()
        self._fill_fwr1_beginning(organism, seq_type)

    def _create_sequences(self):
        """ create AntibodySequence object for each sequence """
        igblast_output = open(self.igblast_file).readlines()
        igblast_records = [ list(y) for x,y in groupby(igblast_output, lambda x: not re.match('# IGBLASTN', x)) ]
        igblast_records = list(filter(lambda x: not x[0].startswith('# IGBLASTN') and len(x) > 0,  igblast_records))
        entries = [ AntibodySequence(info, self.chain) for info in igblast_records if any("# Alignment" in s for s in info) ]

        return entries

    def _fill_fwr1_beginning(self, organism, seq_type):
        if seq_type == "Ig":
            V_db = SeqIO.to_dict(SeqIO.parse(V_fasta[organism],"fasta"), key_function = lambda rec: rec.description.split('|')[1])
        else:
            V_db = SeqIO.to_dict(SeqIO.parse(V_tcr_fasta[organism],"fasta"), key_function = lambda rec: rec.description.split('|')[1])

        for entry in self.entries:
            if not entry.full_fwr1() and (entry.entry_id.find('PGT') == -1 and entry.entry_id.find('10-1074') == -1):
                entry.fill_fwr1_beginning(str(V_db[entry.extract_vdj()[0].split(',')[0]].seq))
            else:
                entry.fill_fwr1_beginning()

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
    def __init__(self, name, igblast_file, changeo_file, chain, organism, seq_type):
        AntibodyLibrary.__init__(self, name, igblast_file, chain, organism, seq_type)
        self.changeo_file = changeo_file
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

    def _output_summary(self):
        with open('results.txt','w') as outf:
            for s in self.entries:
                outf.write("\t".join([s.entry_id,s.seq[s.overflow:],s.aa_seq,s.top_v.split(',')[0],s.top_d.split(',')[0],s.top_j.split(',')[0],s.prod]) + "\n")


class ReferenceAntibodyLibrary(AntibodyLibrary):
    def __init__(self, name, igblast_file, chain, organism, seq_type, ref_id, sec_ref_ids):
        AntibodyLibrary.__init__(self, name, igblast_file, chain, organism, seq_type)
        self.ref_id = ref_id
        self.sec_ref_ids = sec_ref_ids
        self.ref_obj = [ s for s in self.entries if s.entry_id == self.ref_id ][0]

    def get_vdjs(self):
        vdjs = {}
        vdj_count = { 'V': [] , 'D': [], 'J': [] }

        for seq in self.entries:
            vdj = seq.extract_vdj()
            vdjs[seq.entry_id] = vdj

            vdj_count['V'].append(vdj[0])
            vdj_count['D'].append(vdj[1])
            vdj_count['J'].append(vdj[2])

        return vdjs

    def _filter_sequences(self, sameV):
        if sameV:   # include only functional antibody sequences with same V-call as reference in alignment
            ref_V = self.ref_obj.top_v
            filtered_seqs = [ s for s in self.entries if s.entry_id in self.sec_ref_ids or s.prod == 'Yes' \
                             and s.top_v.split('*')[0] == ref_V.split('*')[0] ]
        else:   # include only functional antibody sequences in alignment
            filtered_seqs = [ s for s in self.entries if s.entry_id in self.sec_ref_ids or s.prod == 'Yes']

        fs = sorted(filtered_seqs, key=lambda x:x.entry_id)

        return fs

    def _write_id_order(self, seqs):
        with open(self.name+"_ids.tmp", 'w') as f:
            f.write("\n".join([ s.entry_id for s in sorted(seqs, key=lambda x: x.entry_id) ]))

        return self.name+"_ids.tmp"

    def output_region_fastas(self, filter_functional=False, sameV=False):
        filenames = []

        if filter_functional or sameV:
            seqs = self._filter_sequences(sameV)
        else:
            seqs = self.entries

        id_filename = self._write_id_order(seqs)

        for region in REGIONS:
            nt_records = [ SeqRecord(Seq(s.regions[region][0]), id=s.get_id(), description='| '+region) for s in seqs ]
            aa_records = [ SeqRecord(Seq(s.regions[region][1].replace('*','Z')), id=s.get_id(),
                                description='| '+region) for s in seqs ]

            file_basename = self.name+"_"+region
            with open(file_basename+"_nt.fa",'w') as out_nt:
                SeqIO.write(nt_records, out_nt, "fasta")
            with open(file_basename+"_aa.fa",'w') as out_aa:
                SeqIO.write(aa_records, out_aa, "fasta")

            filenames += [file_basename+"_nt.fa", file_basename+"_aa.fa"]

        self.filtered_sequences = seqs

        return id_filename, filenames



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
        self.outfile = open('seqs.txt','a')
        self.regions, self.germline_regions = self._extract_regions()
        self.nt_mismatches, self.aa_mismatches = self._calculate_mismatches()

    def fill_fwr1_beginning(self, germline_v=None):
        if germline_v is None:
            self.fwr1_beginning_nt, self.fwr1_beginning_aa = "", ""
            self.full_fwr1_input = self.seq
            return

        num_missing_bases = self.overflow + self.vdj_start
        germline_v = germline_v.replace('.','').upper()
        self.fwr1_beginning_nt = germline_v[:num_missing_bases-1]
        self.fwr1_beginning_aa = str(Seq(self.fwr1_beginning_nt).translate())
        self.full_fwr1_input = germline_v[:self.vdj_start-1] + self.seq

    def full_fwr1(self):
        return True if self.vdj_start == 1 else False

    def add_vdj_start(self, position):
        self.vdj_start = 1 if isinstance(position, bool) and position else int(position)

    def _fill_features(self, info, chain):
        """ Parse IgBlast output and extract information about sequence and regions """
        sections = [ list(y) for x,y in groupby(info, lambda x: not re.match('#', x)) ]
        for heading, content in zip(*[iter(sections)] * 2):
            if heading[0].startswith('# V-(D)-J rearrangement'):
                try:
                    self.top_v, self.top_d, self.top_j, self.chain_type, self.stop_codon, self.inframe, self.prod, _ = content[0].split('\t')
                except ValueError:
                    self.top_v, self.top_j, self.chain_type, self.stop_codon, self.inframe, self.prod, _ = content[0].split('\t')
                    self.top_d = 'None'
            elif heading[0].startswith('# V-(D)-J junction'):
                junction_parts = [ x for x in content[0].strip().split('\t')[1:-1] if x != 'N/A' and x.find('(') == -1 ]
                try:
                    junction_overlaps = [ m for m in re.findall('.*\((.*)\).*',content[0]) ]
                except AttributeError:
                    junction_overlaps = []
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

                    if self.first_region is None:
                        if not hasattr(self, 'overflow'):
                            self.overflow = length % 3
                        if length >= 3:        # check beginning frame
                            self.start = start_idx
                            self.first_region = region

                    dels = length - (end-start_idx+1)
                    self.v_insertions += num_gaps - dels
                    self.v_deletions += dels
                    if length < 3 and region != "CDR3-IMGT":
                        self.details[region] = ['NA'] * 6
                    elif region == "CDR3-IMGT" and start_idx-self.start < self.details['FR3-IMGT'][1]:      # check for cdr3 overlapping with fr3
                        cdr3_length = end-self.start - self.details['FR3-IMGT'][1]
                        if cdr3_length < 3:
                            self.details[region] = ['NA'] * 6
                        else:
                            self.details[region] = [self.details['FR3-IMGT'][1]+1, end-self.start, cdr3_length, matches, mismatches, num_gaps]
                    else:
                        self.details[region] = [start_idx-self.start, end-self.start, length, matches, mismatches, num_gaps]


            elif heading[0].startswith('# Hit table'):
                self.seq = ''
                self.germseq = ''
                seen = set()

                for hit in content:
                    if hit[0] not in seen and hit[0] in ['V','J']:
                        s_start, _, _, _, seq, germseq = hit.split('\t')[10:16]

                        if self.chain_type != "VH" and hit[0] == "J" and len(junction_overlaps) != 0 and seq.startswith(junction_overlaps[0]):
                            seq = seq[len(junction_overlaps[0]):]
                            germseq = germseq[len(junction_overlaps[0]):]
                        self.seq += seq
                        self.germseq += germseq

                        if hit[0] == 'V': # add junction / (D) sequence
                            self.vseq = seq
                            self.seq += self.junction
                            self.germseq += 'N' * len(self.junction)
                            self.vdj_start = int(s_start)
                        seen.add(hit[0])

        if not hasattr(self,'details'):
            print(self.entry_id, "has no details")
            self.details = { r:['NA']*6 for r in REGIONS }
        # fill empty regions
        for r in set(REGIONS)-set(self.details.keys()):
            self.details[r] = ['NA','NA','NA','NA','NA','NA']

        # update overflow
        start, end, length = self.details[self.first_region][0:3]
        num_deletions = self.seq[start:start+length].count('-')
        num_insertions = self.germseq[start:start+length].count('-')
        if num_deletions != 0 and num_insertions != 0:              # deletion and insertion in first region
            self.overflow = (self.overflow + num_deletions) % 3
        elif num_insertions != 0:
            self.overflow = (self.overflow - num_insertions) % 3

        #print(self.entry_id,"\t",self.overflow,sep="")

        return self


    def _extract_regions(self):
        region_seqs = {}
        region_germseqs = {}
        insertions, deletions = 0, 0
        for region, details in sorted(self.details.items(), key=lambda i:REGIONS.index(i[0])):
            start, end, length = details[:3]
            aa_seq = None

            if length == 'NA':
                if region != 'CDR3-IMGT' or not hasattr(self, 'cdr3_nt'):
                    region_seqs[region] = ('N'*6, 'X'*2)
                    region_germseqs[region] = ('N'*6, 'X'*2)
                    continue

            if region == self.first_region:
                full_aa_seq = str(Seq(self.seq[self.overflow:].replace('-','')).translate())
                full_aa_germseq = str(Seq(self.germseq[self.overflow:].replace('-','')).translate())

                # rare codon
                self.outfile.write(self.entry_id+"\t"+full_aa_seq+"\t"+self.seq[self.overflow:].replace('-','')+"\n")

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
                    aa_seq = self.cdr3_aa
                    if start == 'NA':
                        a = self.seq.find(self.cdr3_nt)
                        germseq = self.germseq[a:a+len(self.cdr3_nt)].replace('-','')
                    else:
                        germseq = self.germseq[start+deletions:start+deletions+len(self.cdr3_nt)].replace('-','')
                else:
                    seq_aln = self.seq[start+deletions:]
                    self.cdr3_start = start+deletions
                    self.cdr3_aa_start = start
                    self.cdr3_aligned = seq_aln
                    self.cdr3_aligned_germ = self.germseq[start+deletions:]
                    seq = seq_aln.replace('-','')
                    germseq = self.germseq[start+deletions:].replace('-','')
            else:
                if re.match('FR3', region) and hasattr(self, 'cdr3_nt'):
                    # correct for shorter FR3 using CDR3 sequence info
                    cdr_start = self.seq.find(self.cdr3_nt)
                    length += len(self.seq[start+deletions+length:cdr_start])

                seq_aln = self.seq[start+deletions:start+deletions+length]
                germseq = self.germseq[start+deletions:start+deletions+length]
                deletions += seq_aln.count('-')
                insertions += germseq.count('-')
                seq = seq_aln.replace('-','')
                germseq = germseq.replace('-','')

            if not aa_seq:      # CDR3 aa sequence not defined by igblast
                aa_seq, full_aa_seq = full_aa_seq[:round(len(seq)/3)], full_aa_seq[round(len(seq)/3):]
            if aa_seq == '': aa_seq = 'X'           # not enough amino acids due to frameshift
            aa_germseq, full_aa_germseq = full_aa_germseq[:round(len(germseq)/3)], full_aa_germseq[round(len(germseq)/3):]
            #print(self.entry_id, germseq, aa_germseq)

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

        if hasattr(self,'cdr3_nt'):       # remove CDR3 as defined by IgBlast
            seq = self.seq[self.overflow:].split(self.cdr3_nt)[0]
            germseq = self.germseq[self.overflow:self.overflow+len(seq)]
            aa_seq = self.aa_seq.split(self.cdr3_aa)[0]
            aa_germseq = self.aa_germseq[:len(aa_seq)]
        elif self.regions['CDR3-IMGT'][1] != 'XX':              # remove sequence after FR3
            seq = self.seq[:self.cdr3_start]
            germseq = self.germseq[:self.cdr3_start]
            aa_seq = self.aa_seq[:round(self.cdr3_aa_start/3)]
            aa_germseq = self.aa_germseq[:round(self.cdr3_aa_start/3)]
            if aa_seq == '': aa_seq = 'X'
        else:         # short sequence, no CDR3, no need to truncate
            seq = self.seq[self.overflow:]
            germseq = self.germseq[self.overflow:]
            aa_seq = self.aa_seq
            aa_germseq = self.aa_germseq

        assert len(seq) == len(germseq), 'nt lengths not equal\n'+self.entry_id+"\n"+seq+"\n"+germseq+"\n"+str(len(seq))+" "+str(len(germseq))
        nt_mismatches = sum(1 for x,y in zip(seq, germseq) if x != y and y not in ("N","-") and x not in ("N","-"))

        # aa
        aln = pairwise2.align.globalxs(aa_seq, aa_germseq, -10, -.1)
        a1, a2, _ , _, _ = aln[0]
        aa_mismatches = sum(1 for x,y in zip(a1, a2) if x != y and y not in ("X","-") and x not in ("X","-"))

        return nt_mismatches, aa_mismatches

    def _make_chart_seq(self):
        self.inframe

        new_seq = ""
        codon = ""
        self.v_seq
        for x,y in zip(self.seq, self.germseq):
            if x == "-" and y in ["A","T","C","G"]: # deletion
                new_seq += "D"
            elif y == "-" and x in ["A","T","C","G"]: # insertion
                new_seq += "I"
                codon += x
            else:
                new_seq += "N"

    def get_id(self):
        return self.entry_id

    def display(self):
        print("seq %s" % self.seq)
        print("germline seq %s" % self.germseq)
        print("overflow %d" % self.overflow)
        if hasattr(self, 'cdr3_aa'): print("CDR3 %s" % self.cdr3_aa)



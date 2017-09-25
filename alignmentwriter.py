import xlsxwriter
from xlsxwriter.utility import xl_range

REGIONS = ['FR1-IMGT', 'CDR1-IMGT', 'FR2-IMGT', 'CDR2-IMGT', 'FR3-IMGT', 'CDR3-IMGT']

# amino acid equivalency
EQUIVALENCY = {
    'G':['G','A','V','L','I'], 'A':['G','A','V','L','I'], 'V':['G','A','V','L','I'], 'L':['G','A','V','L','I'], 'I':['G','A','V','L','I'],
    'S':['S','T'], 'T':['S','T'],
    'C':['C','M'], 'M':['C','M'],
    'D':['D','N','E','Q'], 'N':['D','N','E','Q'], 'E':['D','N','E','Q'], 'Q':['D','N','E','Q'],
    'R':['R','K','H'], 'K':['R','K','H'], 'H':['R','K','H'],
    'F':['F','Y','W'], 'Y':['F','Y','W'], 'W':['Y','F','W']
}

# positions important for neutralization by PGT
NEUT_POS_HC = [ 55, 70, 102, 104, 111 ]
NEUT_POS_LC = [ 24, 32, 49, 92, 95 ]


def set_up_workbook(outfile):
    workbook = xlsxwriter.Workbook(outfile)
    global CENTER, CENTER_BOLD, BOLD, RED, HIGHLIGHT, RED_AND_HIGHLIGHT
    CENTER = workbook.add_format({ 'align':'center' })
    CENTER_BOLD = workbook.add_format({ 'align':'center', 'bold':1 })
    BOLD = workbook.add_format({ 'bold':1 })
    RED = workbook.add_format({ 'color':'red', 'align':'center' })
    HIGHLIGHT = workbook.add_format({ 'bg_color':'#FFFCC', 'align':'center' })
    RED_AND_HIGHLIGHT = workbook.add_format({ 'color':'red', 'bg_color':'#FFFCC', 'align':'center' })

    return workbook


def initialize_dicts(stats, vdjs, ref_id=None, sec_ref_ids=[]):
    SEQ_IDS = [ref_id] + list(sec_ref_ids)

    for s in sorted(vdjs):
        if s not in SEQ_IDS: SEQ_IDS.append(s)
        for d in stats:
            d[s] = 0

    return SEQ_IDS


def output_ref_excel(library, ref_id, sec_ref_ids, chain, filtered=False, align_cdr3=False):
    #outfile = library.name + "_alignment_filtered.xlsx" if filtered else library.name + "_alignment.xlsx"
    outfile = library.name + "_alignment.xlsx"
    nt_files = [ library.name+"_"+x+"_nt_aligned.fa" for x in REGIONS ]
    aa_files = [ library.name+"_"+x+"_aa_aligned.fa" for x in REGIONS ]

    important_pos = NEUT_POS_HC if chain == 'IgH' else NEUT_POS_LC
    important_pos = adjust_important_columns(aa_files, important_pos)

    # make workbook with 2 worksheets: one for nucleotide alignment, one for amino acid alignment
    workbook = set_up_workbook(outfile)
    nt_worksheet = workbook.add_worksheet('nucleotide')
    aa_worksheet = workbook.add_worksheet('amino acid')

    seqs = library.filtered_sequences if filtered else library.entries
    vdjs = { s.entry_id: s.extract_vdj() for s in seqs}

    write_worksheet(nt_worksheet, nt_files, 'nt', ref_id, vdjs, align_cdr3, sec_ref_ids)
    write_worksheet(aa_worksheet, aa_files, 'aa', ref_id, vdjs, align_cdr3, sec_ref_ids, important_pos)


def output_clone_excel(clone_id, vdjs, library_name):
    """
    vdjs = { sequence_id: (v_call, d_call, j_call) }
    """
    outfile = library_name + "_clone" + clone_id + "_alignment.xlsx"
    nt_files = [ clone_id+"_"+x+"_nt_aligned.fa" for x in REGIONS ]
    aa_files = [ clone_id+"_"+x+"_aa_aligned.fa" for x in REGIONS ]

    ref_id = 'germline'

    # make workbook with 2 worksheets: one for nucleotide alignment, one for amino acid alignment
    workbook = set_up_workbook(outfile)
    nt_worksheet = workbook.add_worksheet('nucleotide')
    aa_worksheet = workbook.add_worksheet('amino acid')

    write_worksheet(nt_worksheet, nt_files, 'nt', ref_id, vdjs, True)
    write_worksheet(aa_worksheet, aa_files, 'aa', ref_id, vdjs, True)


def write_worksheet(worksheet, aln_files, seq_type, ref_id, vdjs, align_cdr3, sec_ref_ids=[], important_pos=[]):
    headers = ['sequence', 'V', 'D', 'J', 'matches (with germline)', 'mismatches (with germline)', 'gaps', 'insertions']
    if seq_type == 'aa':
        headers += ['equiv mismatches', 'nonequiv mismatches',
                                            'neutralization position mismatches',
                                            'neutralization equiv mismatches',
                                            'neutralization nonequiv mismatches']
        col_offset = last_header_col = 13
    else:
        col_offset = last_header_col = 8

    worksheet.write_row('A1', headers, BOLD)
    worksheet.freeze_panes(2,1)         # freeze first column and row
    worksheet.set_column(0,0,20)        # widen seq name column

    match, mismatch, gap, insertion, equiv, nonequiv = {}, {}, {}, {}, {}, {}
    mismatches_specific, equiv_specific, nonequiv_specific = {}, {}, {}
    equiv_count_by_pos = {}
    stats = (match, mismatch, gap, insertion, equiv, nonequiv, mismatches_specific, equiv_specific, nonequiv_specific)

    SEQ_IDS = initialize_dicts(stats, vdjs, ref_id, sec_ref_ids)

    # output region alignments
    last_col = last_header_col
    for i in range(len(REGIONS)):
        region = REGIONS[i]
        region_file = aln_files[i]
        last_col = output_region(region_file, last_col, region, worksheet, seq_type, col_offset,
                                 stats, equiv_count_by_pos, important_pos, sec_ref_ids, align_cdr3)

    num_seqs = len(match)+1
    if seq_type == 'nt':
        worksheet.merge_range(num_seqs,0,num_seqs,last_header_col-1,'number of mismatches',CENTER_BOLD)
        worksheet.merge_range(num_seqs+1,0,num_seqs+1,last_header_col-1,'number of sequences',CENTER_BOLD)
    else:
        worksheet.merge_range(num_seqs,0,num_seqs,last_header_col-1,'number of equivalent mismatches',CENTER_BOLD)
        worksheet.merge_range(num_seqs+1,0,num_seqs+1,last_header_col-1,'total number of mismatches',CENTER_BOLD)
        worksheet.merge_range(num_seqs+2,0,num_seqs+2,last_header_col-1,'number of sequences',CENTER_BOLD)

    # report number of match/mismatch/etc. (row-wise)
    for i, sid in enumerate(SEQ_IDS):
        seq_details = [sid] + list(vdjs[sid]) + [ match[sid], mismatch[sid], gap[sid], insertion[sid] ]
        if seq_type == 'aa':
            if sid == ref_id:
                seq_details += [0,0,0,0,0]
            elif sid in sec_ref_ids:
                seq_details += [ mismatch[sid], 0, mismatches_specific[sid], mismatches_specific[sid], 0 ]
            else:
                seq_details += [ equiv[sid], nonequiv[sid], mismatches_specific[sid], equiv_specific[sid], nonequiv_specific[sid] ]

        worksheet.write_row(i+1, 0, seq_details)

    length_full_seq = last_col-1

    # report number of mutations at each position (col-wise)
    for i in range(last_header_col,length_full_seq+1):
        col_range = xl_range(2+len(sec_ref_ids),i,num_seqs-1,i)

        if seq_type == 'nt':
            worksheet.write_formula(num_seqs,i,'=COUNTIF(%s,"<>.") - COUNTIF(%s,"N") - COUNTIF(%s,"-")' % (col_range, col_range, col_range))
            worksheet.write_formula((num_seqs+1),i,'=COUNTA(%s) - COUNTIF(%s,"N") - COUNTIF(%s,"-")' % (col_range, col_range, col_range))
        else:
            worksheet.write((num_seqs),i,equiv_count_by_pos[i])
            worksheet.write_formula(num_seqs+1,i,'=COUNTIF(%s,"<>.") - COUNTIF(%s,"X") - COUNTIF(%s,"-")' % (col_range, col_range, col_range))
            worksheet.write_formula((num_seqs+2),i,'=COUNTA(%s) - COUNTIF(%s,"X") - COUNTIF(%s,"-")' % (col_range, col_range, col_range))

    worksheet.set_column(last_header_col, length_full_seq, 2)
    worksheet.activate()


def output_region(region_file,last_col,region_name,worksheet,seq_type,col_offset,stats,equiv_count_by_pos,important_pos,sec_ref_ids, align_cdr3):
    match, mismatch, gap, insertion, equiv, nonequiv, mismatches_specific, equiv_specific, nonequiv_specific = stats
    ref_seqs_by_pos = {}
    row = 1
    col = last_col

    for line in open(region_file):
        id, seq = line.strip().split('\t')
        seq = seq.replace('Z','*')
        if row == 1:
            ref_seq = seq

        if not align_cdr3:
            if region_name == 'CDR3-IMGT':       # don't output alignment, just CDR3 sequence
                worksheet.write_string(row,col,seq.replace('-','').replace('.',''))
                row += 1
                continue

        for pos in range(0,len(seq)):
            write_equiv, write_equiv_specific = False, False
            if row == 1:            # reference sequences, all matches
                ref_seqs_by_pos[pos] = set()
                equiv_count_by_pos[col+pos] = 0
            elif row <= len(sec_ref_ids)+1 and seq_type == 'aa':
                # secondary reference sequences: for equivalency calculations
                if seq[pos] != '-': ref_seqs_by_pos[pos].add(seq[pos])       # add current aa to list of accepted ref aa

            if seq[pos] == '.' or seq[pos] == ref_seq[pos]:
                match[id] += 1
            else:
                if ref_seq[pos] == "-" and seq[pos] != "-":
                    insertion[id] += 1
                elif ref_seq[pos] != "-" and seq[pos] == "-":
                    gap[id] += 1
                elif (seq_type == 'aa' and seq[pos] != "X" and ref_seq[pos] != "X") or (seq_type == 'nt' and seq[pos] != "N" and ref_seq[pos] != "N"):
                    mismatch[id] += 1
                    if (col+pos-col_offset) in important_pos:
                        mismatches_specific[id] += 1

            # find equivalent matches with references for input sequences
            if seq_type == 'aa':# and row > len(sec_ref_ids)+1:
                accepted = ref_seqs_by_pos[pos].copy()
                for r in ref_seqs_by_pos[pos]:
                    if r in EQUIVALENCY:
                        accepted |= set(EQUIVALENCY[r])

                if len(accepted) > 0:
                    if seq[pos] != '.' and seq[pos] != '-' and seq[pos] in accepted:
                        equiv[id] += 1

                        if (col+pos-col_offset) in important_pos:
                            equiv_specific[id] += 1
                            write_equiv_specific = True
                        else:
                            write_equiv = True
                    elif seq[pos] != '.' and seq[pos] != '-' and seq[pos] != "X":
                        nonequiv[id] += 1
                        if (col+pos-col_offset) in important_pos:
                            nonequiv_specific[id] += 1

            # write equivalent mismatches in red
            if write_equiv_specific and row > len(sec_ref_ids)+1:
                worksheet.write_string(row,col+pos,seq[pos],RED_AND_HIGHLIGHT)
                equiv_count_by_pos[col+pos] += 1
            elif write_equiv and row > len(sec_ref_ids)+1:
                worksheet.write_string(row,col+pos,seq[pos],RED)
                equiv_count_by_pos[col+pos] += 1
            elif row > 1 and row <= len(sec_ref_ids)+1 and (col+pos-col_offset) in important_pos and seq_type == 'aa':
                worksheet.write_string(row,col+pos,seq[pos],HIGHLIGHT)
            else:
                worksheet.write_string(row,col+pos,seq[pos],CENTER)

            end = col+pos
        row +=1

    if region_name != "CDR3-IMGT" or align_cdr3:
        worksheet.merge_range(0,col,0,end,region_name,CENTER_BOLD)
    else:
        worksheet.merge_range(0,col,0,col+1,region_name,CENTER_BOLD)
        end = col
        equiv_count_by_pos[end] = 0

    return end+1


def adjust_important_columns(aln_files, important_pos):
    """ adjust neutralization columns based on any gaps/insertions """
    seq = ""
    for aln_file in aln_files:
        f = open(aln_file)

        beginning = f.tell()
        for line in f:
            seq += line.split('\t')[1].strip()
            f.seek(beginning)
            break       # read first line (germline reference) only
        f.close()

    for pos in range(0,len(seq)):
        if seq[pos] == '-':
            tmp = important_pos
            important_pos = [ i+1 if i >= pos else i for i in tmp ]

    return important_pos


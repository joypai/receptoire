import operator
from collections import OrderedDict
import xlsxwriter
from xlsxwriter.utility import xl_range
import sys

SORTED_REGIONS = { r:i for i,r in enumerate(['CDR3-IMGT', 'FR1-IMGT','CDR1-IMGT','FR2-IMGT','CDR2-IMGT','FR3-IMGT']) }
CLONE_COLORS = ['#FF0000','#FF6666','#FF6600','#FFCC99','#FFFF00','#FFFFCC',            # red, orange, yellow
          '#008000','#00FF00', '#CCFFCC', '#00FFFF','#99CCFF','#0000FF',                # green, blue
          '#800080','#660099','#FF00FF','#CC99FF','#FF99CC','#E0E0E0']                  # purple, pink, grey

def write_headers(worksheet, clone, bold):
    field_names = ['seq_id','productive','in-frame','stop-codon','V','D','J','full_input','VDJ_start','V_insertions','V_deletions',
                   'nt_mismatches (up to CDR3)', 'aa_mismatches (up to CDR3)', 'CDR3_nt','CDR3_nt_length (V-region)','CDR3_matches (V-region)',
                   'CDR3_mismatches (V-region)','CDR3_gaps (V-region)','CDR3_aa','CDR3_aa_length',
                   'FR1_nt','FR1_nt_length','FR1_matches','FR1_mismatches','FR1_gaps','FR1_aa','FR1_aa_length',
                   'CDR1_nt','CDR1_nt_length','CDR1_matches','CDR1_mismatches','CDR1_gaps','CDR1_aa','CDR1_aa_length',
                   'FR2_nt','FR2_nt_length','FR2_matches','FR2_mismatches','FR2_gaps','FR2_aa','FR2_aa_length',
                   'CDR2_nt','CDR2_nt_length','CDR2_matches','CDR2_mismatches','CDR2_gaps','CDR2_aa','CDR2_aa_length',
                   'FR3_nt','FR3_nt_length','FR3_matches','FR3_mismatches','FR3_gaps','FR3_aa','FR3_aa_length',
                   'junction','junction_length','germline']
    worksheet.set_column(0,0,33)
    if clone:
        field_names = ['clone_id', 'num_seqs'] + field_names + ['alignment_file']
        worksheet.set_column(0,1,8)
        worksheet.set_column(2,2,32)

    worksheet.write_row(0,0,field_names,bold)
    worksheet.freeze_panes(1,1)
    worksheet.set_zoom(65)


def write_summary(library, chain, filter_prod=False, ref_id=None, sec_ref_ids=[], clone=True):
    if clone:
        outfile = library.name+"_clone_filtered_summary.xlsx" if filter_prod else library.name+"_clone_summary.xlsx"
        colors = CLONE_COLORS
    else:
        outfile = library.name+"_ref_filtered_summary.xlsx" if filter_prod else library.name+"_ref_summary.xlsx"
        colors = ['#FFFFFF']

    # make workbook and formats
    workbook = xlsxwriter.Workbook(outfile)
    seq_sheet = workbook.add_worksheet(chain)
    chart_sheet = workbook.add_worksheet(chain+' chart')
    bold = workbook.add_format({ 'bold':1 })                                            # for headers

    write_headers(seq_sheet, clone, bold)
    seq_sheet.activate()

    if clone:
        total_num_seqs = output_by_clone(workbook, seq_sheet, library.clones, library.name)
        write_chart(workbook, chart_sheet, library.clones, total_num_seqs)

        if filter_prod:
            seq_filtered_sheet = workbook.add_worksheet(chain+' filtered')
            chart_filtered_sheet = workbook.add_worksheet(chain+' filtered chart')
            write_headers(seq_filtered_sheet, clone, bold)
            total_num_seqs = output_by_clone(workbook, seq_filtered_sheet, library.functional_clones, library.name)
            write_chart(workbook, chart_filtered_sheet, library.functional_clones, total_num_seqs)
    else:
        total_num_seqs = output_all(workbook, seq_sheet, library, library.name, False)

        if filter_prod:
            seq_filtered_sheet = workbook.add_worksheet(chain+' filtered')
            write_headers(seq_filtered_sheet, clone, bold)
            total_num_seqs = output_all(workbook, seq_filtered_sheet, library, library.name, True)


def output_by_clone(workbook, worksheet, clones, library_name):
    center_format = workbook.add_format({ 'align':'center', 'valign':'center' })        # for clone_id and num_seqs columns
    link_format = workbook.add_format({ 'align':'center', 'valign':'vcenter', 'color':'blue', 'underline':1 })
    fill_form = workbook.add_format({ 'align':'fill' })

    total_num_seqs, color_count = 0, 0
    row = 1
    for clone in clones:
        color = CLONE_COLORS[color_count]
        clone_color_format = workbook.add_format({ 'bg_color':color })
        clone_info_format = workbook.add_format({ 'bg_color':color, 'valign':'vcenter', 'align':'center' })

        first_row = row
        seq_objs = clone.sequences
        num_seqs = len(clone.sequence_ids)
        total_num_seqs += num_seqs
        if num_seqs > 1:     # use colors to output clone
            worksheet.merge_range(row,0,row+num_seqs-1,0,clone.clone_id,clone_info_format)
            worksheet.merge_range(row,1,row+num_seqs-1,1,num_seqs,clone_info_format)
            form = clone_color_format
        else:
            color = '#FFFFFF'
            worksheet.write(row,0,clone.clone_id,center_format)
            worksheet.write(row,1,num_seqs,center_format)
            form = None

        for s in seq_objs:
            seq_entry = [ s.entry_id, s.prod, s.inframe, s.stop_codon, s.top_v, s.top_d, s.top_j,
                         s.seq, s.vdj_start, s.v_insertions, s.v_deletions, s.nt_mismatches, s.aa_mismatches ]
            for r in OrderedDict(sorted(s.details.items(), key=lambda x:SORTED_REGIONS.get(x[0]))):
                seq_entry += [s.regions[r][0]] + s.details[r][2:] + [s.regions[r][1], len(s.regions[r][1])]
            seq_entry += [ s.junction, len(s.junction), s.germseq ]

            worksheet.write_row(row, 2, seq_entry, form)
            row += 1
        col = len(seq_entry) + 2

        # write link to clone alignment file
        if num_seqs > 1:
            #col = len(seq_entry) + 2
            clone_file = "clone_alignments/" + library_name + "_clone" + clone.clone_id + "_alignment.xlsx"
            worksheet.merge_range(first_row,col,row-1,col,'alignment',link_format)
            worksheet.write_url(row-1,col,clone_file,link_format,'alignment')
        else:
            worksheet.write(row-1, col, ' ')            # write empty cell to prevent overflow from last column

        color_count = 0 if color_count == len(CLONE_COLORS)-1 else color_count+1
        row += 1    # skip line between clones

    return total_num_seqs


def output_all(workbook, worksheet, library, library_name, filter_prod):
    row = 1
    for s in library.entries:
        if filter_prod and s.prod != "Yes": continue
        seq_entry = [ s.entry_id, s.prod, s.inframe, s.stop_codon, s.top_v, s.top_d, s.top_j,
                        s.seq, s.vdj_start, s.v_insertions, s.v_deletions, s.nt_mismatches, s.aa_mismatches ]
        for r in OrderedDict(sorted(s.details.items(), key=lambda x:SORTED_REGIONS.get(x[0]))):
            seq_entry += [s.regions[r][0]] + s.details[r][2:] + [s.regions[r][1], len(s.regions[r][1])]
        seq_entry += [ s.junction, len(s.junction), s.germseq ]

        worksheet.write_row(row, 0, seq_entry)
        row += 1


def write_chart(workbook, worksheet, sorted_clones, total_num_seqs):
    worksheet_name = worksheet.get_name()
    chart_title = "heavy chain clones" if worksheet_name.find("IgH") != -1 else "light chain clones"
    worksheet.set_column(0,0,10)

    # assign colors to clones in chart
    clones, counts, color_points = [], [], []
    single_count, color_count = 0, 0

    for c in sorted_clones:
        clone_id = c.clone_id
        num_seqs = len(c.sequence_ids)
        if num_seqs > 1:
            clones.append(clone_id)
            counts.append(num_seqs)
            color_points.append({'fill':{'color':CLONE_COLORS[color_count]},'border':{'color':'#606060'}})
        else:
            single_count += 1

        # set color count for next clone
        color_count = 0 if color_count == len(CLONE_COLORS)-1 else color_count+1

    clones.append('singles')
    counts.append(single_count)
    color_points.append({'fill': {'color':'white'}, 'border': {'color':'#606060'}})

    # write clone count information
    worksheet.write_row(0,0,('clone_id', 'num_seqs'))
    worksheet.write_column('A2',clones)
    worksheet.write_column('B2',counts)

    sum_range = xl_range(1,1,len(clones),1)
    worksheet.write(len(clones)+2,0,'total samples')
    worksheet.write(len(clones)+2,1,'=SUM(%s)' % sum_range)
    worksheet.insert_textbox('G11',str(total_num_seqs),
                             {'font': {'size':12}, 'width': 60, 'height': 30, 'border': {'none':True}})

    # create clone breakdown chart
    clone_chart = workbook.add_chart({'type':'doughnut'})
    clone_chart.set_title({'name':chart_title})
    clone_chart.add_series({
        'categories': [worksheet_name,1,0,len(clones),0],
        'values': [worksheet_name,1,1,len(clones),1],
        'points': color_points
    })
    clone_chart.set_legend({'layout': {'x':0.8, 'y': 0.2, 'width': 0.3, 'height':0.7}})

    worksheet.insert_chart('D3',clone_chart)

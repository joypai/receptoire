#!/usr/bin/env python3
import sys, click, os, subprocess
from IgLibrary import ClonalAntibodyLibrary, ReferenceAntibodyLibrary
from alignmentwriter import output_clone_excel, output_ref_excel
from summarywriter import write_summary

IGBLAST='/rugpfs/fs0/nuss_lab/scratch/jpai/software/ncbi-igblast-1.6.1'
DIST_MODEL = { 'mouse': 'm1n_compat', 'human': 'hs1f_compat'}


# helper functions:
def run_igblastn(input_fasta, model_organism, seq_type, align5end):
    curdir = os.getcwd()
    os.chdir(IGBLAST)
    output_file = os.path.splitext(input_fasta)[0]

    cmd = "{igblastn}/bin/igblastn -query {infile} -out {outfile} \
            -germline_db_V {igblastn}/database/{organism}_{db_type}_v \
            -germline_db_D {igblastn}/database/{organism}_{db_type}_d \
            -germline_db_J {igblastn}/database/{organism}_{db_type}_j \
            -auxiliary_data {igblastn}/optional_file/{organism}_gl.aux \
            -domain_system imgt -ig_seqtype {seq_type} -outfmt '7 std qseq sseq btop' \
            -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 \
            -show_translation {extend} \
            -organism {organism}".format(igblastn=IGBLAST, organism=model_organism,
                                         infile=input_fasta, outfile=output_file,
                                         seq_type=seq_type, db_type="ig" if seq_type == "Ig" else "trb",
                                         extend="-extend_align5end" if align5end else "")

    click.secho("Running Igblast ... ", fg="blue", bold=True, nl=False)
    return_code = subprocess.call(cmd, shell=True)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert return_code == 0, 'Error in call to Igblast. Exiting.'
    click.secho("Done", fg="blue", bold=True)
    os.chdir(curdir)

    return output_file


@click.group()
@click.option('--input-file', '-i', type=click.Path(exists=True), required=True, help='input fasta')
@click.option('--model_organism', '-m', type=click.Choice(['mouse', 'human']),
              prompt="Please enter model organism used",
              help='organism from which sequences were obtained')
@click.option('--seq_type', '-s', type=click.Choice(['Ig', 'TCR']), default='Ig',
              help='sequence type')
@click.option('--extend_5end', '-e', is_flag=True,
              help='invoke IgBLAST parameter to extend alignment for 5\' end (-extend_align5end)')
@click.pass_context
def run_analysis(ctx, input_file, model_organism, seq_type, extend_5end):
    subprocess.call(['date'])
    ctx.obj['INPUT'] = input_file
    ctx.obj['ORGANISM'] = model_organism
    ctx.obj['SEQ_TYPE'] = seq_type
    ctx.obj['EXTEND'] = extend_5end

    # determine if heavy or light chain based on input fasta file name
    ctx.obj['CHAIN'] = 'IgH' if input_file.find('IgH') != -1 else 'IgK' if input_file.find('IgK') != -1 else 'IgL'



@run_analysis.command(help='perform clonal analysis on antibody sequences')
@click.option('--outdir', '-o', default='results', help='output directory name')
@click.option('--cut-off', '-c', type=float, help='clone cutoff distance')
@click.option('--filter_functional', '-f', is_flag=True, help='output clones for functional sequences')
@click.pass_context
def clone(ctx, outdir, cut_off, filter_functional):
    input_fasta = os.path.abspath(ctx.obj['INPUT'])
    chain = ctx.obj['CHAIN']
    model_organism = ctx.obj['ORGANISM']
    seq_type = ctx.obj['SEQ_TYPE']
    library_name = os.path.basename(input_fasta).split('.fasta')[0]

    cur_dir = os.path.dirname(input_fasta)
    os.chdir(cur_dir)

    click.secho('--> Running clonal analysis', fg='green')
    click.echo('input file:\t\t%s' % input_fasta)
    click.echo('sequence type:\t\t%s' % seq_type)
    click.echo('organism:\t\t%s' % model_organism)
    if cut_off: click.echo('specified cut-off:\t%f' % cut_off)
    click.echo('filter nonfunctional:\t%s' % filter_functional)
    click.echo('extend 5\' end:\t\t%s' % ctx.obj['EXTEND'])

    igblast_out = run_igblastn(input_fasta, model_organism, seq_type, ctx.obj['EXTEND'])
    library = ClonalAntibodyLibrary(library_name, igblast_out, chain, model_organism, seq_type)

    click.secho('Creating clone alignment files: ', fg='blue', bold=True, nl=False)
    if not os.path.exists(cur_dir+"/"+outdir):
        os.makedirs(cur_dir+"/"+outdir)
    if not os.path.exists(os.path.join(cur_dir,outdir,'clone_alignments')):
        os.makedirs(os.path.join(cur_dir,outdir,'clone_alignments'))
    for c in library.clones:
        if (len(c.sequence_ids) < 2): continue      # no need to align singlets

        region_filenames = c.output_region_fastas()
        id_order_file = c.clone_id+"_ids.tmp"
        with open(id_order_file, 'w') as f:
            f.write('germline\n')
            f.write("\n".join([ x for x in sorted(c.sequence_ids) ]))

        # align with clustal
        for r in region_filenames:
            cmd = "perl /data04-scratch/toliveira/jpai/antibody_pipeline/align_region_clustal.pl "+r+" "+id_order_file +" True -c "
            return_code = subprocess.call(cmd, shell=True, stdout=open(c.clone_id+"_"+r+"_align.txt","w"))#, stderr=subprocess.PIPE)

        # create clone alignment file
        output_clone_excel(c.clone_id, c.get_vdjs(), library.name)
        click.secho(' %s' % c.clone_id, nl=False)

    subprocess.call("mv *.xlsx "+os.path.join(cur_dir,outdir,'clone_alignments'), shell=True, stderr=subprocess.PIPE)
    subprocess.call("rm *.tmp; rm *.fa; rm *align.txt", shell=True, stderr=subprocess.PIPE)

    # create summary excel file
    click.secho('\nCreating summary file: ', fg='blue', bold=True, nl=False)
    write_summary(library, chain, filter_functional)
    click.secho(library_name+"_clone_summary.xlsx")

    # organize output files
    if not os.path.exists(cur_dir+'/intermediate_files'): os.makedirs(cur_dir+'/intermediate_files')
    subprocess.call("mv *.tab "+cur_dir+"/intermediate_files", shell=True, stderr=subprocess.PIPE)
    subprocess.call("mv *.pdf "+cur_dir+"/intermediate_files", shell=True, stderr=subprocess.PIPE)
    subprocess.call(["mv", igblast_out, cur_dir+"/intermediate_files"])
    subprocess.call("mv *.xlsx "+cur_dir+"/"+outdir, shell=True, stderr=subprocess.PIPE)

    click.secho("Done", fg="green")
    subprocess.call(['date'])


@run_analysis.command(help='align antibody sequences against reference')
@click.option('--reference', '-r', type=click.Path(exists=True), required=True, help="reference sequence to align sequences to")
@click.option('--sec-ref', '-s', type=click.Path(exists=True), default=None, help="secondary references for AA equivalency")
@click.option('--filter_functional', '-f', is_flag=True, help="only align functional sequences")
@click.option('--filter_same_v', '-v', is_flag=True, help="only align sequences with same V gene as reference")
@click.option('--align_cdr3', '-a', is_flag=True, help="align CDR3")
@click.pass_context
def reference(ctx, reference, sec_ref, filter_functional, filter_same_v, align_cdr3):
    input_fasta = os.path.abspath(ctx.obj['INPUT'])
    chain = ctx.obj['CHAIN']
    model_organism = ctx.obj['ORGANISM']
    seq_type = ctx.obj['SEQ_TYPE']
    library_name = os.path.basename(input_fasta).split('.fasta')[0]

    cur_dir = os.path.dirname(input_fasta)
    os.chdir(cur_dir)

    ref_id = open(reference).readline().strip().strip('>')

    click.secho('--> Running analysis against reference', fg='green')
    click.echo('input file:\t\t%s' % input_fasta)
    click.echo('sequence type:\t\t%s' % seq_type)
    click.echo('organism:\t\t%s' % model_organism)
    click.echo('reference id:\t\t%s' % ref_id)

    if sec_ref is not None:
        sec_ref_ids = [ x.strip().strip('>') for x in open(sec_ref).readlines() if x.find('>') != -1 ] if sec_ref else []
        click.echo('secondary reference id:\t%s' % ",".join(sec_ref_ids))
    else:
        sec_ref_ids = []

    click.echo('filter productive:\t%s' % filter_functional)
    click.echo('filter same V-call:\t%s' % filter_same_v)
    click.echo('align CDR3:\t\t%s' % align_cdr3)
    click.echo('extend 5\' end:\t\t%s' % ctx.obj['EXTEND'])

    with open('input_ref_combined.fasta','w') as outfile:
        if sec_ref is not None:
            subprocess.call(['cat', reference, sec_ref, ctx.obj['INPUT']], stdout=outfile)
        else:
            subprocess.call(['cat', reference, ctx.obj['INPUT']], stdout=outfile)

    igblast_out = run_igblastn(cur_dir+'/input_ref_combined.fasta', model_organism, seq_type, ctx.obj['EXTEND'])
    library = ReferenceAntibodyLibrary(library_name, igblast_out, chain, model_organism, seq_type, ref_id, sec_ref_ids)

    click.secho('Creating alignment file: ', fg='blue', bold=True, nl=False)
    id_filename, region_filenames = library.output_region_fastas(filter_functional, filter_same_v)

    # align with clustal
    for r in region_filenames:
        cmd = "perl /data04-scratch/toliveira/jpai/antibody_pipeline/align_region_clustal.pl {rfile} {order} {cdr3} \
            -a {ref} {srefs}".format(rfile=r,order=id_filename,cdr3=align_cdr3,ref=ref_id,srefs=" ".join(sec_ref_ids))
        return_code = subprocess.call(cmd, shell=True, stdout=open("ref_align.txt","w"))#, stderr=subprocess.PIPE))

    click.secho(library_name+"_alignment.xlsx")

    # create clone alignment file
    output_ref_excel(library, ref_id, sec_ref_ids, chain, filter_functional, align_cdr3)
    subprocess.call("rm *.tmp; rm *.fa; rm ref_align.txt", shell=True)

    # create summary excel file
    click.secho('Creating summary file: ', fg='blue', bold=True, nl=False)
    write_summary(library, chain, filter_functional, ref_id, clone=False)
    click.secho(library_name+"_ref_summary.xlsx")

    click.secho("Done", fg="green")
    subprocess.call(['date'])

if __name__ == '__main__':
    run_analysis(obj={})

#!/usr/bin/env python3
import sys, click, os, subprocess
from IgLibrary import ClonalAntibodyLibrary
from alignmentwriter import output_clone_excel
from summarywriter import write_summary

# tools:
#CLUSTALW2='~/software/clustalw2'
IGBLAST='/data_alpha/home/jpai/software/ncbi-igblast-1.6.1'

DIST_MODEL = { 'mouse': 'm1n', 'human': 'hs1f'}


# helper functions:
def run_igblastn(input_fasta, model_organism):
    curdir = os.getcwd()
    os.chdir(IGBLAST)
    output_file = os.path.splitext(input_fasta)[0]

    cmd = "{igblastn}/bin/igblastn -query {infile} -out {outfile} \
            -germline_db_V {igblastn}/database/{organism}_ig_v \
            -germline_db_D {igblastn}/database/{organism}_ig_d \
            -germline_db_J {igblastn}/database/{organism}_ig_j \
            -auxiliary_data {igblastn}/optional_file/{organism}_gl.aux \
            -domain_system imgt -ig_seqtype Ig -outfmt '7 std qseq sseq btop' \
            -show_translation -organism {organism}".format(igblastn=IGBLAST, organism=model_organism, infile=input_fasta, outfile=output_file)

    click.secho("Running Igblast ... ", fg="blue", bold=True, nl=False)
    return_code = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert return_code == 0, 'Error in call to Igblast. Exiting.'
    click.secho("Done", fg="blue", bold=True)
    os.chdir(curdir)

    return output_file


def changeo(input_fasta, igblast_out, organism, *args):
    # parse Igblast output
    click.secho("Parsing Igblast output", fg="blue", bold=True)
    cmd = "MakeDb.py igblast -s {infile} -i {outfile} -r {igblastn} \
        --regions --scores --failed".format(infile=input_fasta, outfile=igblast_out, igblastn=IGBLAST)

    return_code = subprocess.call(cmd, shell=True, stderr=subprocess.PIPE)
    click.echo("\n")
    #assert return_code == 0, 'Error in call to MakeDb.py. Exiting.'

    # check for CDR3 region; remove bad sequences with missing CDR3 that escaped changeo filtering
    filtered_output = [ line for line in open(igblast_out+"_db-pass.tab") if len(line.strip().split('\t')) == 45 and line.strip().split('\t')[44] != "" ]

    with open(igblast_out+".tab","w") as f:
        f.write("".join(filtered_output))

    # determine clone cut-off
    click.secho("Determining clone cut-off distance", fg="blue", bold=True)
    if args[0] is None:      # no clone cut-off distance specified; use dynamic clone clustering
        cmd = "Rscript ~/scripts/cluster_cutoff.R {outfile}.tab {outfile}_dist.pdf {model} | cut -d ' ' -f2".format(outfile=igblast_out, model=DIST_MODEL[organism])
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        cutoff = float(p.communicate()[0])
        click.secho("chosen cut-off distance: ", fg="yellow", nl=False)
    else:                   # use specified clone distance cut-off
        cutoff=float(args[0])
        click.secho("specified cut-off distance: ", fg="yellow", nl=False)
    click.echo(cutoff)

    # identify clones
    cmd = "DefineClones.py bygroup -d {outfile}.tab --model {clone_model} \
        --act set --sym min --norm len --nproc 1 --failed \
        --dist {cutoff}".format(outfile=igblast_out, cutoff=cutoff, clone_model=DIST_MODEL[organism])
    return_code = subprocess.call(cmd, shell=True, stderr=subprocess.PIPE)
    assert return_code == 0, 'Error in call to DefineClones.py. Exiting.'

    """
    # reconstruct germline
    cmd = "CreateGermlines.py -d {outfile}_clone-pass.tab -r {igblastn} \
        -g dmask --cloned --failed".format(outfile=igblast_out, igblastn=IGBLAST)
    return_code = subprocess.call(cmd, shell=True, stderr=subprocess.PIPE)
    assert return_code == 0, 'Error in call to CreateGermlines.py. Exiting.'

    #return igblast_out+'_clone-pass_germ-pass.tab'
    """
    return igblast_out+'_clone-pass.tab'



@click.group()
@click.option('--input-file', '-i', type=click.Path(exists=True), required=True, help='input fasta')
@click.pass_context
def cli(ctx, input_file):
    subprocess.call(['date'])
    ctx.obj['INPUT'] = input_file

    # determine if heavy or light chain based on input fasta file name
    ctx.obj['CHAIN'] = 'IgH' if input_file.find('IgH') != -1 else 'IgK' if input_file.find('IgK') != -1 else 'IgL'


@cli.command(help='perform clonal analysis on antibody sequences')
@click.option('--model-organism', '-m', type=click.Choice(['mouse', 'human']),
              prompt="Please enter model organism used",
              help='organism from which sequences were obtained')
@click.option('--outdir', '-o', default='results')
@click.option('--cut-off', '-c', type=float,
              help='clone cutoff distance')
@click.pass_context
def clone(ctx, model_organism, outdir, cut_off):
    input_fasta = ctx.obj['INPUT']
    chain = ctx.obj['CHAIN']
    library_name = os.path.basename(input_fasta).split('.fasta')[0]

    cur_dir = os.path.dirname(input_fasta)
    os.chdir(cur_dir)

    click.secho('--> Running clonal analysis', fg='green')
    click.echo('input file: %s' % ctx.obj['INPUT'])
    click.echo('organism: %s' % model_organism)
    if cut_off: click.echo('specified cut-off: %f' % cut_off)

    igblast_out = run_igblastn(ctx.obj['INPUT'], model_organism)
    changeo_out = changeo(input_fasta, igblast_out, model_organism, cut_off)
    library = ClonalAntibodyLibrary(library_name, igblast_out, changeo_out, chain)
    library.report_entries()

    click.secho('Creating clone alignment files: ', fg='blue', bold=True, nl=False)
    if not os.path.exists(cur_dir+"/"+outdir): os.makedirs(cur_dir+"/"+outdir)
    if not os.path.exists(os.path.join(cur_dir,outdir,'clone_alignments')): os.makedirs(os.path.join(cur_dir,outdir,'clone_alignments'))
    for c in library.clones:
        if len(c.sequence_ids) < 2: continue

        region_filenames = c.output_region_fastas()
        id_order_file = c.clone_id+"_ids.tmp"
        with open(id_order_file, 'w') as f:
            f.write('germline\n')
            f.write("\n".join([ x[:30] for x in sorted(c.sequence_ids) ]))  # truncate to 30 chars to match clustal

        # align with clustal
        for r in region_filenames:
            cmd = "perl ~/scripts/align_region_clustal.pl "+r+" "+id_order_file +" -c "
            return_code = subprocess.call(cmd, shell=True)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # create clone alignment file
        output_clone_excel(c.clone_id, c.get_vdjs(), library.name)
        click.secho(' %s' % c.clone_id, nl=False)

    subprocess.call("mv *.xlsx "+os.path.join(cur_dir,outdir,'clone_alignments'), shell=True)
    subprocess.call("rm *.tmp; rm *.fa", shell=True)

    # create summary excel file
    click.secho('\nCreating summary file', fg='blue', bold=True)
    write_summary(library, chain, True)

    # organize output files
    if not os.path.exists(cur_dir+'/intermediate_files'): os.makedirs(cur_dir+'/intermediate_files')
    subprocess.call("mv *.tab "+cur_dir+"/intermediate_files", shell=True)
    subprocess.call("mv *.pdf "+cur_dir+"/intermediate_files", shell=True)
    subprocess.call(["mv", igblast_out, cur_dir+"/intermediate_files"])
    subprocess.call("mv *.xlsx "+cur_dir+"/"+outdir, shell=True)

    click.secho("Done", fg="green")
    subprocess.call(['date'])


@cli.command(help='align antibody sequences against reference')
@click.option('--model_organism', '-m', type=click.Choice(['mouse', 'human']),
            help='organism from which sequences were obtained')
@click.option('--reference', '-r', type=click.Path(exists=True))
@click.option('--sec-ref', '-s', multiple=True)
@click.pass_context
def reference(ctx, model_organism, reference, sec_ref):
    input_fasta = ctx.obj['INPUT']
    chain = ctx.obj['CHAIN']

    cur_dir = os.path.dirname(input_fasta)
    os.chdir(cur_dir)

    click.secho('--> Running analysis against reference', fg='green')
    click.echo('input file: %s' % ctx.obj['INPUT'])

    for r in sec_ref: print(r)

if __name__ == '__main__':
    cli(obj={})

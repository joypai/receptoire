## Set-up

#### Requirements
* python packages: biopython, xlsxwriter, click, changeo
* igblast 1.6.1 (and blast database)
* clustalw 2.1
* antibody pipeline python and perl scripts


#### Install required modules and scripts
`pip install biopython xlsxwriter click changeo
echo 'export PATH=/data04-scratch/toliveira/jpai/antibody_pipeline:$PATH' >> ~/.bashrc;
source ~/.bashrc`

## Running the pipeline

Convert trace files to fasta
`batch_ab1tofasta.py {ab1_data_dir} {output_dir} {output_fasta_name}`

example: 
`cd /data04-scratch/toliveira/jpai/amelia_escolano/antibody/alivaMAb/LC;
batch_ab1tofasta.py IP1_P1/data IP1_P1 IP1_P1_LC.fasta`


### Clone analysis
`antibody_analysis.py -i {input_fasta} -m {organism} [-e] clone -c {cut_off} -o {output_dir}`

required arguments:
* input_fasta
* organism: mouse or human
optional parameters:
* cut_off: clone distance
    * default=0.05
* output_dir
    * default=results
flags:
* -e: extend 5’ end of alignment (igblastn parameter extend_align5end)

example:
`cd /data04-scratch/toliveira/jpai/amelia_escolano/antibody/alivaMAb/LC/IP1_P1;
antibody_analysis.py -i IP1_P1_LC.fasta -m human clone -c 0.05`


### Alignment to reference
`antibody_analysis.py -i {input_fasta} -m {organism} [-e] reference -r {primary_reference_fasta} -s {secondary_reference_fasta} [-a -f -v]`

required arguments:
* input_fasta
* organism: mouse or human
* primary_reference_fasta: reference to align sequences to
optional parameters:
* secondary_reference_fasta: additional references to use for AA equivalency
flags:
* -a: align CDR3 region
* -e: extend 5’ end of alignment (igblastn parameter extend_align5end)
* -f: only align functional sequences
* -v: only align sequences with same V gene as reference

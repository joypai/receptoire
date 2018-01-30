# receptoire: antibody repertoire analysis 

## Set-up

#### Requirements
* python packages: biopython, xlsxwriter, click
* igblast 1.6.1 (and blast database)
* clustalw 2.1

#### Install required modules and scripts
```shell
pip install biopython xlsxwriter click
echo 'export PATH=/path/to/receptoire/:$PATH' >> ~/.bashrc;
source ~/.bashrc
```

## Running the program

Convert trace files to fasta
```shell
batch_ab1tofasta.py {ab1_data_dir} {output_dir} {output_fasta_name}
```

### Clone analysis
```shell
receptoire.py -i {input_fasta} -m {organism} [-e] clone -c {cut_off=0.05} -s {seq_type=Ig} -o {output_dir=results}
```

required arguments:
* input_fasta
* organism: sequence origin [ mouse | human ]

optional parameters:
* cut_off: clone distance
    * default = 0.05
* seq_type: receptor type [ Ig | TCR ]
   * default = Ig
* output_dir
    * default = results
flags:
* -e: extend 5’ end of alignment (igblastn parameter extend_align5end)


### Alignment to reference
```shell
receptoire.py -i {input_fasta} -m {organism} [-e] reference -r {primary_reference_fasta} -s {secondary_reference_fasta} [-a -f -v]
```

required arguments:
* input_fasta
* organism: sequence origin [ mouse | human ]
* primary_reference_fasta: reference to align sequences to

optional parameters:
* secondary_reference_fasta: additional references to use for AA equivalency
flags:
* -a: align CDR3 region
* -e: extend 5’ end of alignment (igblastn parameter extend_align5end)
* -f: only align functional sequences
* -v: only align sequences with same V gene as reference

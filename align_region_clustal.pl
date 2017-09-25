# USAGE
# align all:    perl parse_clustal.pl -a ref_id [ sec_ref_id1 sec_ref_id2 ... ]
# align clone:  perl parse_clustal.pl -c

BEGIN { $ENV{CLUSTALDIR} = '/rugpfs/fs0/nuss_lab/scratch/jpai/software/clustalw-2.1-linux-x86_64-libcppstatic/' }

use strict;
use warnings;
use IO::Handle;
use Bio::AlignIO;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::SeqIO;
use Bio::SimpleAlign;
use Data::Printer;
use feature qw(say);

open ERROR, '>', "error.log" or die $!;
STDERR->fdopen( \*ERROR,'w') or die $!;

my %seq_ids;

sub align_region {
    die "ERROR: Could not open alignment input fasta ('$ARGV[0]')\n$1" if ! -e -s $ARGV[0];
    my $in_fasta = Bio::SeqIO->new( -file => $ARGV[0],
                                    -format => "fasta" );
    my @sequences;
    while ( my $seq = $in_fasta->next_seq() ) {
        push(@sequences, $seq);
        say $seq->id;
        #$seq_ids{( substr $seq->id, 0, 30 )} = $seq->id;       # preserve sequence names to fix clustal truncation
        $seq_ids{$seq->id} = $seq->id;       # preserve sequence names to fix clustal truncation
    }

    # set up clustalw environment and align
    my $aligner = Bio::Tools::Run::Alignment::Clustalw->new();
    my $aln = $aligner->align(\@sequences);

    return $aln;
}

# extract primary reference type
my $template;
if ($ARGV[3] eq '-a') {
    $template = $ARGV[4];       # ref_seq name passed in as argument
}
elsif ($ARGV[3] eq '-c') {
    $template = 'germline';
}
else {
    say "not a valid reference option";
    exit;
}
#say $template;

# get sequence ids to preserve alignment order across regions
my $id_order_file = $ARGV[1];

my $aln = align_region();

# slice alignment to match reference residues
my ($first_col, $last_col);
if ( $ARGV[0] =~ /FR1/ ) {
    $first_col = $aln->column_from_residue_number( $template, 1 );         # find first residue
}
else {
    $first_col = 1;
}

# find last residue for cdr3 splice
if ( $ARGV[0] =~ /CDR3/ ) {
    my $ref_seq = $aln->select_noncont_by_name($template);
    my $no_gaps = $ref_seq->remove_gaps;
    my $num_res = $no_gaps->length;
    $ref_seq->gap_char('-');
    my $res = $ref_seq->num_residues;
    my $s = ($ref_seq->get_seq_by_id($template))->seq;
    my $gaps = ($s =~ s/\./\./g);
    my $real_numres = $num_res - $gaps;
    $last_col = $aln->column_from_residue_number($template,$real_numres);
}
else {
    $last_col = $aln->length;
}

$aln = $aln->slice($first_col,$last_col,1);


# output alignment
my $outfname = ( split /.fa/, $ARGV[0] )[0] . "_aligned.fa";
say $outfname;
open( my $outfile, '>', $outfname ) or die "Could not open file '$outfname' $!";

my @refs;
if (scalar(@ARGV) > 4) {
    @refs = ($template,@ARGV[5..scalar(@ARGV)-1]);         # secondary reference sequences for equivalency analysis
    #say "refs: @refs[1..scalar(@refs)]";

    # make new alignment object with only reference sequence
    my $ref_aln = $aln->select_noncont_by_name(@refs);
    foreach my $seq ($ref_aln->each_seq) {
        my $ref_seq_out = $seq->seq;
        $ref_seq_out =~ tr/\./-/;
        say $outfile $seq->id . "\t" . $ref_seq_out;
    }
}

# replace query sequence matches with dots
$aln = $aln->sort_by_list($id_order_file);
$aln = $aln->set_new_reference($template);
if ( ($ARGV[0] !~ /CDR3/) || $ARGV[2] eq "True" ) {
    $aln->match();
}


foreach my $seq ($aln->each_seq) {
    if ( grep $_ eq $seq->id, @refs ) {
        next
    }
    else {
        say $outfile $seq_ids{$seq->id} . "\t" . $seq->seq;
    }
}

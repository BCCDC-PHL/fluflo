#!/usr/bin/env nextflow

//enable domain-specific-language 2
nextflow.enable.dsl=2

/**
---------------------------------------------------------------------------------
function definition
---------------------------------------------------------------------------------
*/

def helpMe() {
  log.info """

Overview:
Nextflow pipeline for generation of phylogenetic tree of influenza.
Replaces snakefile.

Usage:
nextflow run main.nf -profile conda [OPTIONS]

Mandatory arguments:
 --work_dir                     User's directory that contains input 'config' & 'data' folders.
                                Ensure trailing / is present in absolute path. If not specified,
                                will run on a test dataset located at /home/jess.cal/flu/TestSet_FluTree20_HA/

Optional arguments:
 --seqs                         Multi-fasta file containing consensus sequences of interest [./data/sequences.fasta]
 --ref                          Reference genome used to align reads to during guided assembly [./config/Ref.gb]
 --meta                         File containing metadata for sequences under analysis [./data/metadata.csv]
 --drop_strains                 Excluded strains/ samples [./config/dropped_strains.txt]
 --colors                       Colors used in final auspice visualization [./config/colors.csv]
 --lat_long                     Sample latitudes and longitudes [./config/lat_longs.csv]
 --auspice                      Specifications for visualization in auspice (ex. title) [./config/auspice_config.json]
 --divergence_units             Units used to measure divergence in phylogeny refining step [mutations]
 --version                      Current fluflo version number
 --help                         This usage statement
        """
}

def version() {
  log.info """
  fluflo version: ${workflow.manifest.version}
  """
}

//displays help upon request
if (params.help) {
  helpMe()
  exit 0 //stop running
}

//version upon request
if (params.version) {
  version()
  exit 0
}

/**
---------------------------------------------------------------------------------
program introduction
---------------------------------------------------------------------------------
*/

// this prints program header with mandatory input
log.info """

🅕 🅛 🅤 Ⓕ Ⓛ Ⓞ
\n========================================
data directory: ${params.work_dir}

"""

process align {

  tag "Aligning sequences to ${params.ref} & filling gaps with N"
  publishDir "${params.work_dir}/results/", mode: 'copy'

  input:
  tuple file(sequences), file(reference)

  output:
  path("aligned.fasta")

  """
  augur align \
  --sequences ${sequences} \
  --reference-sequence ${reference} \
  --output aligned.fasta \
  --fill-gaps \
  --nthreads 32
  """
}

process tree {
  tag "Building IQTREE - GTR+I+R method - 10 iterations - FAST mode"
  publishDir "${params.work_dir}/results/", mode: 'copy'

  input:
  file(aln)

  output:
  path("aligned.fasta.treefile")

  """
  iqtree -alninfo -ninit 10 -n 10 -me 0.05 -nt auto -s ${aln} -m GTR+I+R
  """
}

process refine {
  tag "Refining phylogeny with Augur"
  publishDir "${params.work_dir}/results/", mode: 'copy'

  input:
  tuple file(tree), file(msa)

  output:
  tuple path("tree.nwk"), path("branch_lengths.json")

  """
  augur refine \
      --tree ${tree} \
      --alignment ${msa} \
      --metadata ${params.meta} \
      --timetree \
      --keep-root \
      --divergence-units ${params.divergence_units} \
      --output-tree tree.nwk \
      --output-node-data branch_lengths.json
  """
}

process ancestral {
    tag "Reconstructing ancestral sequences and mutations"
    publishDir "${params.work_dir}/results/", mode: 'copy'

    input:
    tuple file(refine_tree), file(branch_len), file(msa)

    output:
    path("nt_muts.json")

        """
        augur ancestral \
            --tree ${refine_tree} \
            --alignment ${msa} \
            --output-node-data nt_muts.json \
            --keep-overhangs \
            --keep-ambiguous
        """
}

process translate {
  tag "Translating amino acid sequences"
  publishDir "${params.work_dir}/results", mode: 'copy'

  input:
  tuple file(node_data), file(refine_tree), file(branch_len), file(ref)

  output:
  file("aa_muts.json")

  """
  augur translate \
      --tree ${refine_tree} \
      --ancestral-sequences ${node_data} \
      --reference-sequence ${ref} \
      --output-node-data aa_muts.json \
      """
}

process export {
  tag "Exporting data files for auspice"
  publishDir "${params.work_dir}/auspice", mode: 'copy'

  input:
  tuple file(refine_tree), file(branch_len), file(nt_muts),\
  file(aa_muts)

  output:
  file("flu_na.json")

      """
      augur export v2 \
          --tree ${refine_tree} \
          --metadata ${params.meta} \
          --node-data ${branch_len} ${nt_muts} ${aa_muts} \
          --colors ${params.colors} \
          --lat-longs ${params.lat_long} \
          --minify-json \
          --auspice-config ${params.auspice} \
          --output flu_na.json
      """
}

workflow {
  seq_ch = Channel.fromPath(params.seqs, checkIfExists:true)
  ref_ch = Channel.fromPath(params.ref, checkIfExists:true)

  align(seq_ch.combine(ref_ch)) | tree
  msa_ch = align.out
  refine(tree.out.combine(align.out))
  ancestral(refine.out.combine(msa_ch))
  translate(ancestral.out.combine(refine.out.combine(ref_ch)))
  export(refine.out.combine(ancestral.out.combine(translate.out)))
}

workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: 'jessica.caleta@bccdc.ca', subject: 'Flu Tree Analysis Complete on Sabin', body: msg)
}

/**

process clean {
  tag "Removing Nextflow work directory?"
    shell:
        "rm -rfv work"
}
*/

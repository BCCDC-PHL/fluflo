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
 --ref                          Reference genome used to align reads to during guided assembly [./config/Ref.gb (default); *.fasta file also accepted here]
 --ref_anno                     Reference genome annotation file; only required when using a FASTA under --ref [*.gb / *.gff file accepted here]
 --meta                         File containing metadata for sequences under analysis [./data/metadata.csv]
 --drop_strains                 Excluded strains/ samples [./config/dropped_strains.txt]
 --colors                       Colors used in final auspice visualization [./config/colors.csv]
 --lat_long                     Sample latitudes and longitudes [./config/lat_longs.csv]
 --auspice                      Specifications for visualization in auspice (ex. title) [./config/auspice_config.json]
 --divergence_units             Units used to measure divergence in phylogeny refining step [mutations]
 --recursion_limit              Augur recursion threshold used to set env variable [10000]
 --conda_cache                  File system path where Conda env is to be stored [fluflo-main/work/]
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

/**
---------------------------------------------------------------------------------
process definition
---------------------------------------------------------------------------------
*/

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
  --nthreads ${task.cpus}
  """
}

process tree {
  tag "Building IQTREE - model: ${params.sub_model} - iterations: ${params.bootstrap}"
  publishDir "${params.work_dir}/results/", mode: 'copy'

  input:
  file(aln)

  output:
  path("aligned.fasta.treefile"), emit: treefile
  path("aligned.fasta.contree"), optional: true, emit: contree
  
  """
  iqtree -alninfo \
  -ninit 10 \
  ${params.bootstrap} \
  -me 0.05 \
  -nt ${task.cpus} \
  -s ${aln} \
  -m GTR+I+R 
  """
}

process refine {
  tag "Refining phylogeny with Augur"
  publishDir "${params.work_dir}/results/", mode: 'copy'

  input:
  tuple file(tree), file(msa), file(metadata)

  output:
  tuple path("tree.nwk"), path("branch_lengths.json")

  script:
  root_value = params.root_name ? "--root ${params.root_name}" : "--keep-root" 
  """
  augur refine \
      --tree ${tree} \
      --alignment ${msa} \
      --metadata ${metadata} \
      --timetree \
      ${root_value} \
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

process fix_aa_json {
  tag "Fixing aa_muts.json when using a GFF3 file for augur translate."
  publishDir "${params.work_dir}/results/", mode: 'copy'

  input:
  tuple file(aa_muts_json), file(nt_muts_json)

  output:
  file("aa_muts_fix.json")

  """
  fix_aa_muts.py --aa_json ${aa_muts_json} --nt_json ${nt_muts_json} --outpath aa_muts_fix.json
  """
}

process export {
  tag "Exporting data files for auspice"
  publishDir "${params.work_dir}/auspice", mode: 'copy'

  input:
  tuple file(refine_tree), file(branch_len), file(nt_muts), file(aa_muts) 
  file(metadata)
  tuple file(auspice_config), file(colors), file(lat_long)

  output:
  file("flu.json")

  """
  export AUGUR_RECURSION_LIMIT=${params.recursion_limit}

  augur export v2 \
      --tree ${refine_tree} \
      --metadata ${metadata} \
      --node-data ${branch_len} ${nt_muts} ${aa_muts} \
      --colors ${colors} \
      --lat-longs ${lat_long} \
      --minify-json \
      --auspice-config ${auspice_config} \
      --output flu.json
  """
}

/**
---------------------------------------------------------------------------------
workflow
---------------------------------------------------------------------------------
*/

workflow {
  seq_ch = Channel.fromPath(params.seqs, checkIfExists:true)
  ref_ch = Channel.fromPath(params.ref, checkIfExists:true)
  meta_ch = Channel.fromPath(params.meta, checkIfExists:true)
  config_ch = Channel.fromPath([params.auspice, params.colors, params.lat_long], checkIfExists:true).collect()


  // Catch invalid reference input combinations 
  ref_gb_format = (params.ref =~ /.+\.[Gg]b$/)

  if (!ref_gb_format && params.ref_anno == 'NO_FILE' ){                         // Cannot have an empty --ref_anno parameter if reference is in non-GenBank format
    error "ERROR: Parameter --ref_anno (.gff3 or .gb format) must be specified if non-GenBank formatted reference is provided under --ref."
  }
  if (params.ref_anno != 'NO_FILE' && !(params.ref_anno =~ /.+\.gff.?|.+\.[Gg]b/) ){     // Can only have .gff or .gb formats in the --ref_anno parameter
    error "ERROR: Parameter --ref_anno must be in either .gff or .gb (GenBank) format."
  }
  
  // Load the ref_anno_ch channel appropriately 
  if (ref_gb_format){                                                         // Copy the ref_ch channel if in GenBank format (ref_ch can be reused as ref_anno_ch)
    ref_anno_ch = ref_ch
  }else{                                                                      // Load new channel from scratch if different reference annotation format specified
    ref_anno_ch = Channel.fromPath(params.ref_anno, checkIfExists:true)
  }

  align(seq_ch.combine(ref_ch)) | tree
  msa_ch = align.out
  tree_ch = tree.out.contree.ifEmpty(tree.out.treefile)
  refine(tree_ch.combine(msa_ch).combine(meta_ch))
  ancestral(refine.out.combine(msa_ch))
  translate(ancestral.out.combine(refine.out).combine(ref_anno_ch))
  
  ch_aa_muts = translate.out

  if (params.ref_anno != 'NO_FILE' && params.ref_anno =~ /.+\.gff.?/ ) {        // If gff annotation format used, augur translate outputs need to be fixed (causes downstream schema error)
    ch_aa_muts = fix_aa_json(ch_aa_muts.combine(ancestral.out))
  }
  export(refine.out.combine(ancestral.out).combine(ch_aa_muts), meta_ch, config_ch)
}

/**
---------------------------------------------------------------------------------
optional notification of completion
---------------------------------------------------------------------------------
*/
/**
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

    sendMail(to: '', subject: 'Flu Tree Analysis Complete on Sabin', body: msg)
}



process clean {
  tag "Removing Nextflow work directory?"
    shell:
        "rm -rfv work"
}
*/

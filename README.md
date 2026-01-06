![image](/pics/fluflo_logo.png)

## Introduction

Nextflow pipeline for generation of phylogenetic trees to be visualized with Auspice. 
FLUFLO is written by JMC and adapted from a snakefile by Kimia Kamelian which generates 
phylogenies with the Augur bioinformatic toolkit (MAFFT, IQ-TREE, TimeTree..) that can be visualized in Auspice from [Nextstrain](https://docs.nextstrain.org/projects/auspice/en/stable/index.html)

The original intent of the pipeline was for Influenza A sequences, with the flexibility
afforded by Nextflow to be applied to other pathogens, used with various workload managers (SGE vs. SLURM),
and easily adjusted.

## Table of Contents

- [Introduction](#introduction)
- [Quick-Start Guide](#quick-start-guide)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Input](#input)
- [Output](#output)
- [Workflow](#workflow)
- [References](#references)

## Quick-Start Guide

Change into project directory:
```
cd /home/user/flu/fluflo/
```
Run FLUFLO pipeline:
```
nextflow run main.nf -profile conda --work_dir /home/user/flu/input_data/
```
For details on available arguments, enter:
```
nextflow run main.nf --help
```

## Dependencies

[Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) is required to build an environment with required workflow dependencies.

This bioinformatic pipeline requires Nextflow:
```
conda install -c bioconda nextflow
```
or download and add the nextflow executable to a location in your user $PATH variable:
```
curl -fsSL get.nextflow.io | bash
mv nextflow ~/bin/
```
Nextflow requires Java v8.0+, so check that it is installed:
```
java -version
```
The OS-independent conda environment activated upon running fluflo is specified in the
```environment.yml``` file of the project directory and is built when 
```-profile conda``` is included in the command line. Nextflow will save
the environment to the project directory by default. Alternatively, the 
necessary conda environment can be saved to a different shared location 
accesible to compute nodes by adding ```--conda_cache /path/to/new/location/```.

## Installation

To copy the program into a directory of your choice, from desired directory run:
```
git clone https://github.com/BCCDC-PHL/fluflo.git
cd fluflo
nextflow run main.nf \
-profile conda \
--work_dir /home/user/flu/input_data/
```
or run directly using:
```
nextflow run BCCDC-PHL/fluflo \
-profile conda \
--work_dir /home/user/flu/input_data/
```

## Input

The pipeline requires the following files which should be present in the config
and data folders of the directory containing sequences to be analyzed. These
are named the same within different directories - the only thing that needs to be changed
each run is the input directory, which can be specified with the --work_dir flag on the
command line.

- Multi-fasta file containing consensus sequences of interest [./data/sequences.fasta]
- One of the two following configurations is required to describe the reference genome:
  1. A single GenBank file passed to `--ref` describing the reference genome used for both alignment and amino acid annotation [./config/Ref.gb] OR 
  2. A reference sequence in FASTA format under `--ref` AND a reference annotation in GFF format under `--ref_anno` [⚠️](#warning).

- File containing metadata for sequences under analysis [./data/metadata.csv]
- Excluded strains/ samples [./config/dropped_strains.txt]
- Colors used in final auspice visualization [./config/colors.csv]
- Sample latitudes and longitudes [./config/lat_longs.csv]
- Specifications for visualization in auspice (ex. title) [./config/auspice_config.json]

#### Warning: 
When generating phylogenies from concatenated sequences with `fluflo >=v1.0.0`, the *.gff3 file passed with `--ref_anno` must not contain lines specifying regions that conflict with the assumption of a single contiguous sequence as specified by the sequence-region line in the header.

**causes error:**
![error](/pics/concat_gff_example_error.png) 

**runs successfully:**
![success](/pics/concat_gff_example_correct.png) 

### Bootstrapping:

Visualization of bootstrap support is implemented in fluflo v1.0.1. Bootstrap parameters are passed to IQ-Tree via the `--bootstrap` parameter. 

![big-o](/pics/big-o.png)
Trees generated with `-b 100` that exceeded ~ 1000 samples timed out after 7 days running with 64GB memory. Trees with 20 nonparametric bootstrap replicates (`-b 20`) also became intractable beyond ~ 5000 tips. The number of trees was too low to accurately determine best fit and is only included to demonstrate the cutoff of samples before requiring a reduction in number of bootstrap replicates. Ten replicates worked across all trees in the dataset but is insufficient to conclude clade support. Note the [difference in interpretation](https://iqtree.github.io/doc/Frequently-Asked-Questions#how-do-i-interpret-ultrafast-bootstrap-ufboot-support-values) of UltraFast bootstrap support values compared to standard support values. 

#### Warning: 
Note that using UFBoot (`-B`) and SH-aLRT (`--alrt`) in conjuction as recommended by IQ-Tree results in incorrect parsing by `augur refine` because node names are replaced by SH-aLRT/UFBoot support values (see [Issue #856](https://github.com/nextstrain/augur/issues/856) of Nextstrain/augur). 

This was confirmed by analyzing 5 trees from the fluflo_v0.4.0:fluflo_v1.0.1 verification dataset using `-B 1000 --alrt 1000`. The resulting trees either failed (3/5) or were produced and could be visualized but without branch support values (2/5) (ie. the nodes were labelled incorrectly but matched with JSON node-data files like nt_muts.json, aa_mut.json, branch_lengths.json). As a result of the lack of reliability displayed on branches and lack of reliability in successful generation of trees, using `-B` & `--alrt` together is not recommended until parsing is fixed. Using `--alrt` independently was not tested, but would likely function properly, as it is a single value like nonparametric or UFBoot support values as opposed to two (SH-aLRT/UFBoot), which are parsed as node names instead of confidence values. 

## Output

The output directories are `results`, `auspice`, and `reports`. The default location of these folders is the same as the input directory set by `--work_dir` but can be written to an alternative directory using `--out_dir`.

`results`:
- aligned.fasta
- aligned.fasta.treefile
- aligned.fasta.contree (only when bootstrapping)
- branch_lengths.json
- tree.nwk
- nt_muts.json
- aa_muts.json

*NOTE: results folder generated by fluflo has less files than the results folder
produced with the original snakefile because they are extraneous to the analysis,
thus, designed not to be captured in Nextflow channels and copied to the results
directory.

`auspice`:
- flu_na.json

`reports`:
- fluflo_usage.html
- fluflo_timeline.html
- fluflo_dag.html


## Workflow

```mermaid
flowchart TD
    p0((sequences.fasta))
    p1((Ref.gb/
    reference.fasta))
    p2((metadata.csv))
    p3((auspice_config.json,
    colors.csv,
    lat_longs.csv))
    p4([collect])
    p5((ref_anno.gb/
    ref_anno.gff3))
    p6([combine])
    p7[align]
    p8[tree]
    p11([combine])
    p12([combine])
    p13[refine]
    p14([combine])
    p15[ancestral]
    p16([combine])
    p17([combine])
    p18[translate]
    p19([combine])
    p20[fix_aa_json]
    p21([combine])
    p22([combine])
    p23[export]
    p24(( ))
    p0 -->|seq_ch| p6
    p1 -->|ref_ch| p6
    p2 -->|meta_ch| p12
    p3 --> p4
    p4 -->|config_ch| p23
    p5 -->|ref_anno_ch| p17
    p6 --> p7
    p7 --> p8
    p8 -->|tree_ch| p11
    p7 --> p11
    p11 --> p12
    p12 --> p13
    p13 --> p14
    p7 --> p14
    p14 --> p15
    p15 --> p16
    p13 --> p16
    p16 --> p17
    p17 --> p18
    p18 --> p19
    p15 --> p19
    p19 --> p20
    p20 --> p22
    p13 --> p21
    p15 --> p21
    p21 --> p22
    p22 --> p23
    p2 -->|meta_ch| p23
    p23 --> p24

```


## References

1. Hadfield, J. et al. NextStrain: Real-time tracking of pathogen evolution. Bioinformatics 34, 4121–3 (2018).

2. Huddleston J, Hadfield J, Sibley TR, Lee J, Fay K, Ilcisin M, Harkins E, Bedford T, Neher RA, Hodcroft EB, (2021). Augur: a bioinformatics toolkit for phylogenetic analyses of human pathogens. Journal of Open Source Software, 6(57), 2906, https://doi.org/10.21105/joss.02906

3. Katoh, K., Misawa, K., Kuma, K., & Miyata, T. (2002). MAFFT: a novel method for rapid
multiple sequence alignment based on fast Fourier transform. Nucleic Acids Research,
30(14), 3059–3066. https://doi.org/10.1093/nar/gkf436

4. Nguyen, L.-T., Schmidt, H. A., Haeseler, A. von, & Minh, B. Q. (2014). IQ-TREE: A Fast and
Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. Molecular Biology and Evolution, 32(1), 268–274. https://doi.org/10.1093/molbev/msu300

5. Sagulenko, P., Puller, V., & Neher, R. A. (2018). TreeTime: Maximum-likelihood phylodynamic analysis. Virus Evolution, 4(1). https://doi.org/10.1093/ve/vex042

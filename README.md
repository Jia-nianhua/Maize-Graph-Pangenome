# Maize-Graph-Pangenome

# MaizeGP downstream analysis tutorial

## 0. Introduction

MaizeGP 1.0 is a graph-based pangenome resource constructed from 48 diverse maize inbred lines, integrating a graph-based genome structure, a comprehensive structural variant (SV) catalog, a curated transposable element (TE) library, and multiple integrated ATAC-seq datasets. This tutorial provides step-by-step example scripts for core downstream analyses using MaizeGP datasets, including graph-based genotyping in resequencing populations, genotyping of new genome assemblies using the graph, population-level SV genotyping, integration of new assemblies into the graph, TE annotation, and SV-based GWAS.

These examples are designed as a starting point for applying MaizeGP to user-specific research questions. All code blocks are ready to run with minimal modification.

![Main Figure](https://github.com/Jnhcau/Maize-Graph-Pangenome/blob/main/image/mainFig.png)

**Key resources**

- **MaizeGP Portal:** An interactive web platform for locus-level graph exploration and data download — visit https://www.maizepan.cn. The portal enables users to browse bubble-resolved SVs, inspect genome-specific paths through the graph, and download versioned datasets without command-line setup.

- **48pan.gfa.gz:** Graph-based pangenome structure constructed from 48 maize assemblies in GFA format.

- **48pansv.vcf.gz:** Structural variant (SV) catalog derived from 48 maize genomes, including insertions, deletions, and complex variants.

- **Bubble Structural Variant Data:** Bubble-resolved structural variants extracted from the graph pangenome, representing alternative haplotypes across 48 genomes.

- **ATAC-seq Accessible Chromatin Peaks:** Integrated ATAC-seq peaks from multiple public datasets, mapped to the B73 v5 reference genome.

- **Sequence Presence/Absence Variation (PAV):** Genome-specific sequences identified by aligning 47 maize assemblies to the B73 reference genome.

- **PAV Candidate Genes:** Genes associated with presence/absence variation, identified by projecting PAV sequences onto annotated coding regions.

- **Orthologous Groups (OG):** Orthologous gene groups identified across 48 maize genomes.

- **GitHub repository:** All scripts used in this tutorial, along with the full pangenome construction pipeline, are available at https://github.com/Jia-nianhua/Maize-Graph-Pangenome.

## 1. Software requirements and environment setup

The analyses in this tutorial span several tool categories. Before you begin, make sure the following are installed and available on your `PATH`.

| Category             | Tools                              | Purpose in this tutorial                                     |
| -------------------- | ---------------------------------- | ------------------------------------------------------------ |
| Graph-based analysis | `vg`, `minigraph`                  | Pangenome graph indexing, read mapping, variant calling, and graph construction |
| Variant processing   | `bcftools`, `samtools`, `bedtools` | VCF/BCF manipulation, BAM handling                           |
| SV genotyping        | `Paragraph`                        | Genotyping known SVs from short-read BAM files               |
| TE annotation        | `RepeatMasker`                     | Annotating sequences against the 48-pan-TE library           |
| GWAS                 | `plink`, `GEMMA`                   | PCA, kinship matrix estimation, and mixed-model association  |

**Recommended environment**

We recommend creating a dedicated conda environment to avoid version conflicts:

```
conda create -n maizegp -c bioconda vg minigraph bcftools samtools bedtools repeatmasker plink gemma
conda activate maizegp
```

---

## 2. SV genotyping in resequencing populations using the graph pangenome

This section describes how to use the MaizeGP graph (`48pan.gfa.gz`) to genotype structural variants in resequencing populations. The workflow consists of four stages: (1) index the pangenome graph, (2) map short reads to the graph with `vg giraffe`, (3) call variants from graph-based read support, and (4) merge, filter, and export genotypes across samples.

**Input requirements:**

- The graph file `48pan.gfa.gz` from MaizeGP
- Paired-end short reads for each sample (FASTQ format, gzipped)
- Sufficient memory (≥ 20 GB recommended) and threads for large cohorts

**Stage 1: Build graph indices**

```
mkdir -p index gam pack vcf
vg autoindex --workflow giraffe --gfa graph/48pan.gfa --prefix index/maizegp_48pan --threads 32
vg snarls index/maizegp_48pan.xg > index/maizegp_48pan.snarls
```

The `vg autoindex` command generates all necessary indices (`.gbz`, `.min`, `.dist`, `.xg`) for the giraffe mapper. The `.snarls` file defines the sites (bubbles) in the graph where variation will be called.

**Stage 2: Map reads and call variants per sample**

Repeat the following for each sample, replacing `sample1` with the actual sample name:

```
vg giraffe -Z index/maizegp_48pan.gbz -m index/maizegp_48pan.min -d index/maizegp_48pan.dist -f reads/sample1_R1.fq.gz -f reads/sample1_R2.fq.gz -t 32 > gam/sample1.gam
vg stats -a gam/sample1.gam > gam/sample1.gam.stat

vg pack -x index/maizegp_48pan.xg -g gam/sample1.gam -Q 5 -o pack/sample1.pack -t 32

vg call -s index/maizegp_48pan.xg -r index/maizegp_48pan.snarls -k pack/sample1.pack -s sample1 -a -A -t 32 > vcf/sample1.vcf
```

A few details on the parameters: `-Q 5` in `vg pack` sets the minimum mapping quality for reads contributing to the coverage pack; `-a` enables calls at all graph sites (not just those with alternate alleles); `-A` produces genotype likelihoods for all alleles.

**Stage 3: Merge across samples and filter**

```
bgzip -f vcf/sample1.vcf
bcftools index -c vcf/sample1.vcf.gz

ls vcf/*.vcf.gz > vcf.list
bcftools merge -l vcf.list -Oz -o maizegp.population.sv.vcf.gz
bcftools index -c maizegp.population.sv.vcf.gz

bcftools query -f '%CHROM\t%POS\t%ID[\t%GT]\n' maizegp.population.sv.vcf.gz > maizegp.population.sv.genotype.txt

bcftools +fill-tags maizegp.population.sv.vcf.gz -Oz -o maizegp.population.sv.tagged.vcf.gz -- -t MAF,F_MISSING
bcftools index -c maizegp.population.sv.tagged.vcf.gz

bcftools view -i 'F_MISSING < 0.2 && MAF > 0.05' maizegp.population.sv.tagged.vcf.gz -Oz -o maizegp.population.sv.filtered.vcf.gz
bcftools index -c maizegp.population.sv.filtered.vcf.gz
```

The filter in the last step retains sites with < 20% missing genotypes and minor allele frequency > 5%. Adjust these thresholds according to your population size and analysis goals.

**Output:** A population-level filtered VCF (`maizegp.population.sv.filtered.vcf.gz`) ready for downstream association or population-genetic analyses, and a genotype matrix (`maizegp.population.sv.genotype.txt`) for custom processing.

![Main Figure](https://github.com/Jnhcau/Maize-Graph-Pangenome/blob/main/image/1.png)
---

## 3. SV genotyping of new assemblies using the graph pangenome

If you have a newly assembled genome (in FASTA format) and want to genotype SVs relative to the MaizeGP pangenome, this section shows how to map the assembly to the graph with `minigraph` and extract variant calls.

This approach is faster than read-based mapping for assemblies, because it aligns the contig sequences directly to the graph without going through short-read alignment. Typical runtime for a maize-sized assembly is under one hour with 30 threads.

```
minigraph -cxasm --call -t 30 pangenome.sv.gfa.gz sample.fa > sample.bed

paste *.bed | ./k8 mgutils.js merge -s <(./agc listset maizepan) - | gzip > maizepan.sv.bed.gz
./k8 mgutils-es6.js merge2vcf -r0 maizepan.sv.bed.gz > maizepan.sv.vcf
```

**Output:** A BED file describing SVs for each assembly (first command), which can be merged across assemblies into a multi-sample VCF (subsequent commands using k8/mgutils helpers).

![Main Figure](https://github.com/Jnhcau/Maize-Graph-Pangenome/blob/main/image/2.png)
---

## 4. Adding a new assembly to the MaizeGP graph

To expand the pangenome graph, a new genome assembly can be incorporated using minigraph in construction mode. The new assembly (NEW.fa) is aligned against the existing graph (48pan.gfa), producing an augmented graph (NEW.gfa) that preserves the original graph structure while integrating novel sequences and structural variation from the added genome.

```
minigraph -cxggs -t 16 48pan.gfa NEW.fa > NEW.gfa
```

**Input requirements:**

- The existing graph in GFA format (`48pan.gfa`)
- A new genome assembly in FASTA format (`NEW.fa`)
- Thread count (`-t`) should be adjusted to your available cores

**Output:** An augmented GFA graph (`NEW.gfa`) that includes the new assembly. This graph can then be used for all downstream analyses described in Section 2&3.

![Main Figure](https://github.com/Jnhcau/Maize-Graph-Pangenome/blob/main/image/3.png)
---

## 5. SV genotyping using the 48-pan SV catalog

When you have a predefined set of SVs (e.g., the 48-pan SV catalog in VCF format) and short-read sequencing data, Paragraph provides an alternative genotyping strategy. Unlike graph-based genotyping, Paragraph uses local reassembly around each SV breakpoint.

**Input requirements:**

- A VCF of known SVs (`48pansv.vcf.gz`)
- Sorted, indexed BAM files for each sample
- A reference genome FASTA (`B73NAMV5.fa`)
- Paragraph installed (see its documentation for build instructions)

```
idxdepth -b sample.bam -r B73NAMV5.fa -o sample.txt

/opt/paragraph/bin/multigrmpy.py -i 48pansv.vcf.gz -m sample.txt -r B73NAMV5.fa -o sample --threads 128

ls *.vcf.gz > vcf.list
bcftools merge -l vcf.list -Oz -o maizegp.population.sv.vcf.gz
bcftools index -c maizegp.population.sv.vcf.gz
```

`idxdepth` computes per-locus depth statistics from the BAM file; `multigrmpy.py` then performs localized read reassembly and genotype likelihood calculation at each SV site. The final merge step combines all individual VCFs into a population-level genotype matrix.
**Output:** A multi-sample SV genotype VCF (`maizegp.population.sv.vcf.gz`) ready for downstream association (Section 8)

![Main Figure](https://github.com/Jnhcau/Maize-Graph-Pangenome/blob/main/image/4.png)
---

## 6. Merging and extending the SV catalog with user-provided SV calls

If you have your own set of SV calls (from any caller or pipeline) and want to combine them with the MaizeGP 48-pan SV catalog, SURVIVOR provides a straightforward merging workflow. The merged SV set can then be used as input for Paragraph-based genotyping (Section 5), enabling you to genotype both the MaizeGP-discovered SVs and your own SVs in a single pass.

**Workflow overview:**

1. Prepare your SV calls in VCF format, ensuring chromosome names are consistent with the MaizeGP catalog (B73 v5 reference).
2. Merge your VCF with `48pansv.vcf.gz` using SURVIVOR. The input file list (`merged_sv.list`) should contain the paths to both VCFs, one per line.
3. The merged, nonredundant VCF becomes the new input panel for Section 5.

```
SURVIVOR merge merged_sv.list 1000 2 1 1 1 50 merged_sv_catalog.vcf
```

**Parameter notes:**

| Value  | Description                                                  |
| ------ | ------------------------------------------------------------ |
| `1000` | Maximum distance (bp) between breakpoints for two SVs to be considered the same event |
| `2`    | Minimum number of supporting callers — only output SVs identified in at least 2 samples |
| `1`    | Take the type into account — only merge SVs of the same type (deletion / insertion / inversion, etc.) |
| `1`    | Take the strands of SVs into account — only merge SVs with the same strand orientation |
| `1`    | Estimate distance based on the size of SV (1 = yes, 0 = no)  |
| `50`   | Minimum SV length (bp); SVs shorter than this are discarded  |

**Output:** A nonredundant merged SV catalog (`merged_sv_catalog.vcf`) that combines the MaizeGP panel with your own SV calls. Use this as the input VCF (`-i`) in Section 5 to genotype the expanded panel.

![Main Figure](https://github.com/Jnhcau/Maize-Graph-Pangenome/blob/main/image/5.png)
---

## 7. TE annotation of genome assemblies using the 48-pan-TE library

The MaizeGP project includes a high-quality, manually curated transposable element library (`TElib.clean.fa`) built from the 48-genome pangenome. This library can be used with RepeatMasker to annotate TE content in newly assembled or user-provided genome assemblies.

**Input requirements:**

- The 48-pan-TE library (`TElib.clean.fa`)
- A genome assembly in FASTA format (`assembly.fa`)

```
RepeatMasker -e rmblast -pa 60 -qq -lib TElib.clean.fa assembly.fa
```

Key parameters: `-e rmblast` selects the RMBlast search engine, `-pa 60` uses 60 parallel threads, and `-qq` runs in quick mode (faster, suitable for large query sets against a target library).

**Output:** RepeatMasker produces a `.out` annotation file, a `.tbl` summary table, and a `.masked` FASTA with TE bases replaced by Ns. The `.out` file is the primary result for downstream interpretation.

<img src="https://github.com/Jnhcau/Maize-Graph-Pangenome/blob/main/image/6.png" width="50%">
---

## 8. SV-based GWAS

Genome-wide association studies using the SV genotypes generated in Sections 2 or 5 can identify SVs associated with phenotypic variation. This section provides a standard workflow using PLINK for data preparation and GEMMA for mixed-linear-model association.

**Input requirements:**

- PLINK-format genotype files (`.bed`, `.bim`, `.fam`) — these can be converted from the filtered SV VCF using PLINK's `--vcf` option
- A phenotype file in GEMMA format (`multi_trait.gemma.pheno`)

```
plink --bfile SNP --pca 3 --out pca

awk 'NR==FNR { pc[$2]=$3" "$4" "$5; next } { print 1, pc[$2] }' pca.eigenvec PLINK.fam > covariates.txt

gemma -bfile SNP -gk 2 -p multi_trait.gemma.pheno
gemma -bfile SNP -k output/result.sXX.txt -lmm 2 -p multi_trait.gemma.pheno -c covariates.txt -n 1
```

The first `gemma` command (`-gk 2`) computes a standardized relatedness matrix (kinship). The second runs a univariate linear mixed model (`-lmm 2`) with the kinship matrix and population-structure covariates (top 3 PCs). The `-n 1` flag specifies column 1 of the phenotype file as the trait to analyze.

**Output:** GEMMA writes association statistics (effect sizes, standard errors, p-values) to a results directory under `output/`. Manhattan and Q-Q plots can be generated from these outputs using R packages such as `qqman`.

<img src="https://github.com/Jnhcau/Maize-Graph-Pangenome/blob/main/image/7.png" width="50%">
---

## 9. Supplementary scripts for de novo pangenome construction

In addition to the downstream analysis tutorials above, the complete pipeline for building a graph pangenome from scratch is available in the `pangenome_scripts/` directory of the MaizeGP GitHub repository. Each subdirectory addresses one stage of the construction process:

| Directory                    | Description                                                  |
| ---------------------------- | ------------------------------------------------------------ |
| `01.GenomeAssessment`        | Genome quality evaluation using LAI (LTR Assembly Index) and BUSCO (Benchmarking Universal Single-Copy Orthologs) |
| `02.OrthologIdentification`  | Identification of orthologous gene groups (OGGs) across genomes |
| `03.SVCalling`               | Whole-genome alignment with MUMmer, followed by SV calling with SyRI |
| `04.GraphConstruction`       | Graph pangenome construction from the SV calls and genome assemblies |
| `05.RNA-seqanalysis`         | RNA-seq data processing and differential expression analysis pipeline |
| `06.ATAC-seqanalysis`        | ATAC-seq data processing and peak-calling pipeline           |
| `07.panNLR`                  | NLR (nucleotide-binding leucine-rich repeat) gene annotation and pan-NLRome analysis |
| `08.heritability_estimation` | Heritability estimation using LDAK (Linkage Disequilibrium Adjusted Kinships) |

> **Repository:** The MaizeGP code and scripts are available at [https://github.com/Jia-nianhua/Maize-Graph-Pangenome](https://github.com/Jia-nianhua/Maize-Graph-Pangenome).

> **Note:** The graph files, SV catalogs, and TE libraries themselves are too large for GitHub hosting. Download these from the [MaizeGP Portal](https://www.maizepan.cn) Data Downloads module.

If you have any questions, feel free to contact me. For more data resources, please visit [maizepan.cn](http://maizepan.cn)

Please cite the following if you use MaizeGP, the portal, or downloaded datasets in your work:

*MaizeGP: An open graph-based maize pangenome platform for high-resolution, multi-scale exploration of structural variation.* (submitted,2026)

# Maize-Graph-Pangenome

# MaizeGP downstream analysis tutorial

## Introduction

This tutorial provides example scripts for downstream analyses using MaizeGP datasets. These examples are designed to help users reproduce core analyses and extend MaizeGP data to their own research questions.

## Software requirements

The following tools are required:

- Graph-based analysis
  - `vg` (>=1.40)
  - `minigraph`
  - `k8` (for mgutils scripts)
- Variant processing
  - `bcftools`
  - `samtools`
  - `bedtools`
- SV genotyping
  - `Paragraph`
- TE annotation
  - `RepeatMasker`
- GWAS
  - `PLINK`
  - `GEMMA` or `GAPIT`

Recommended environment:

```
conda create -n maizegp -c bioconda vg minigraph bcftools samtools bedtools repeatmasker plink gemma
conda activate maizegp
```

## 1.SV genotyping in resequencing populations using the graph pangenome

This section describes how to use the MaizeGP graph (`48pan.gfa.gz`) to genotype structural variants (SVs) in resequencing populations.

```
mkdir -p index gam pack vcf
vg autoindex --workflow giraffe --gfa graph/48pan.gfa.gz --prefix index/maizegp_48pan --threads 32
vg snarls index/maizegp_48pan.xg > index/maizegp_48pan.snarls

vg giraffe -Z index/maizegp_48pan.gbz -m index/maizegp_48pan.min -d index/maizegp_48pan.dist -f reads/sample1_R1.fq.gz -f reads/sample1_R2.fq.gz -t 32 > gam/sample1.gam
vg stats -a gam/sample1.gam > gam/sample1.gam.stat

vg pack -x index/maizegp_48pan.xg -g gam/sample1.gam -Q 5 -o pack/sample1.pack -t 32

vg call -s index/maizegp_48pan.xg -r index/maizegp_48pan.snarls -k pack/sample1.pack -s sample1 -a -A -t 32 > vcf/sample1.vcf

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

## 2. SV genotyping of new assemblies using the graph pangenome

This section shows how to map a new assembly to the pangenome graph and extract SV information.

```
minigraph -cxasm --call -t 30 pangenome.sv.gfa.gz sample.fa > sample.bed

paste *.bed | ./k8 mgutils.js merge -s <(./agc listset maizepan) - | gzip > maizepan.sv.bed.gz
./k8 mgutils-es6.js merge2vcf -r0 maizepan.sv.bed.gz > maizepan.sv.vcf
```

## 3. Adding a new assembly to the MaizeGP graph

This section describes how to incorporate a new genome assembly into the existing MaizeGP graph to expand the pangenome.

```
minigraph -cxggs -t 16 48pan.gfa NEW.fa > NEW.gfa
```

## 4. SV genotyping using a known SV panel (Paragraph)

This section demonstrates how to genotype known SVs using short-read data.

```
idxdepth -b sample.bam -r B73NAMV5.fa -o sample.txt

/opt/paragraph/bin/multigrmpy.py -i 48pansv.vcf.gz -m sample.txt -r B73NAMV5.fa -o sample --threads 128

ls *.vcf.gz > vcf.list
bcftools merge -l vcf.list -Oz -o maizegp.population.sv.vcf.gz
bcftools index -c maizegp.population.sv.vcf.gz
```

## 5.Overlap of SVs with genes, promoters, and ATAC-seq peaks

```
bedtools intersect -a sv.bed -b genes.gff3 -wa -wb > sv.gene.overlap.txt
bedtools intersect -a sv.bed -b promoter_2kb.bed -wa -wb > sv.promoter.overlap.txt
bedtools intersect -a sv.bed -b atac_peaks.bed -wa -wb > sv.atac.overlap.txt
```

## 6. Annotation of SVs or user-defined sequences using the 48-pan-TE library

```
RepeatMasker -e rmblast -pa 60 -qq -lib TElib.clean.fa SV.fa 
```

## 7. SV-based GWAS

GWAS can be performed using GEMMA or GAPIT.

```
plink --bfile SNP --pca 3 --out pca

awk 'NR==FNR { pc[$2]=$3" "$4" "$5; next } { print 1, pc[$2] }' pca.eigenvec PLINK.fam > covariates.txt

gemma -bfile SNP -gk 2 -p multi_trait.gemma.pheno
gemma -bfile SNP -k output/result.sXX.txt -lmm 2 -p multi_trait.gemma.pheno -c covariates.txt -n 1
```
![Main Figure](https://github.com/Jnhcau/Maize-Graph-Pangenome/blob/main/image/mainFig.jpg)

If you have any questions, feel free to contact me. For more data resources, please visit [maizepan.cn](http://maizepan.cn)

Please cite the following if you use MaizeGP, the portal, or downloaded datasets in your work:

*MaizeGP: A Graph-Based Pangenome for High-Resolution Structural Variation and Exploratory Trait Analysis in Maize.* (submitted)

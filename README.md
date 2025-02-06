# GATK4SCNA
Somatic copy number calling with GATK

This pipeline was built for the purpose of calling somatic copy number
variations (CNVs) from whole exome sequence (WES) data and whole genome
sequence (WGS) data. To do this it uses the GATK4 package and follows
the steps outlined in [part
1](https://gatk.broadinstitute.org/hc/en-us/articles/360035531092) and
[part 2](https://gatk.broadinstitute.org/hc/en-us/articles/360035890011)
of the GATK4 how to guide. It was designed for use on compute1 and uses
the GRCh38.d1.vd1.fa reference genome. It can be used to call CNVs on
tumor samples with or without paired normals as it relies upon a panel
of normals generated from normal samples generated via the same
sequencing methodology.

## Dependencies

This pipeline has two version, one for WGS and one for WES. The pipeline
relies upon a number of software dependencies and reference files. While
these versions use the same software dependencies they do not all use
the same reference files and scripts. They should not be run
interchangeably. Do not use the WES version of the pipeline with WGS
data. Do not run the WGS version of the pipeline with WES data. All of
these dependencies are listed in the config files for the respective
pipelines:

    WES:
    /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    or
    /code/config/config.gatk4scna.compute1.ini

    WGS:
    /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini

### Software

A docker image containing all of the necessary software pre-installed as
well as the scripts to run the different steps is available on [docker
hub
(austins2/gatk4scna)](https://hub.docker.com/repository/docker/austins2/gatk4scna).
This docker image contains and relies on the following software and
packages. There are no differences between the WES and WGS versions of
the pipeline regarding software used in the docker image. The scripts
however are different.

-   GATK4 [version
    4.1.9.0](https://hub.docker.com/r/broadinstitute/gatk/tags)
-   python3 version 3.6 or higher
    -   numpy
    -   pandas
    -   seaborn
-   Java (jre-8u152-linux-x64)

### Reference Files

##### Genome: [GRCh38.d1.vd1](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)

The WES and WGS versions of this pipeline rely upon the same genome.

    /storage1/fs1/dinglab/Active/Projects/austins2/references/GRCh38.d1.vd1/GRCh38.d1.vd1.fa

##### Genome dictionary

The WES and WGS versions of this pipeline rely upon the same genome
dictionary.

    /storage1/fs1/dinglab/Active/Projects/austins2/references/GRCh38.d1.vd1/GRCh38.d1.vd1.dict

##### Target Interval List

There are three versions of the target interval list file that I have
generated so far. The first two are compatible with the WES pipeline and
are present within the docker image. The third config that currently
exists is not present within the docker image and the location must be
separately specified when running on compute1. **If a different
reference is ever to be used then these files should be remade using the
new reference genome.**

The target file that is used as the default config file for this docker
image was made from the bed file for IDT_xGen_Exome_Research_Panel_v1
for GRCh38. The path to this file can be found in
`/code/config/config.gatk4scna.compute1.ini` inside the docker image. If
your WES data was not generated using this target panel then you will
need to generate and use your own config file. Copy the current config
file to a location outside of the docker image and then change the
`TargetIntervalList=` variable to the path to your target interval_list.
Another TargetIntervalList that comes with this docker image is
`nexterarapidcapture_exome_targetedregions_v1.2.hg38.bed.target.preprocessed.exome.interval_list`.
This is for WES data from Nextera Rapid Capture Targeted regions v1.2
using GRCh38. The location of the third target interval file, used for
the WGS version of the pipeline, can be found within the WGS config file
`/storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini`

If your sequencing data was not generated via one of these two target
panels then you will need to generate and make you own interval list. In
such a case there are 3 options:

1.  Repeat these steps with your target bed file
    1.  convert the bed file to an interval list using
        [Picard](https://broadinstitute.github.io/picard/)

<!-- -->

    java -jar picard.jar BedToIntervalList
        I=input.bed
        O=list.interval_list
        SD=reference_sequence.dict

::#

<li value="2">

Preprocess the intervals to prepare for coverage collection with GATK4
[1](https://gatk.broadinstitute.org/hc/en-us/articles/360050815612-PreprocessIntervals)

</li>

    gatk PreprocessIntervals \
        -L targets_C.interval_list \
        -R /gatk/ref/Homo_sapiens_assembly38.fasta \
        --bin-length 0 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O sandbox/targets_C.preprocessed.interval_list

:::#

<li value="1">

The commands used for the IDT WES bed file in the docker were:

</li>

    $ bsub -q research-hpc -n 1 -R "select[mem>30000] rusage[mem=30000]" -M 30000000 -a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
    -oo /gscmnt/gc2521/dinglab/austins2/tools/GATK4SCNA/db/logs/IDTxGEN_bed_to_interval_list.log \
    /gscmnt/gc2521/dinglab/austins2/software/jre1.8.0_152/bin/java -jar \
    /gscmnt/gc2521/dinglab/austins2/software/picard.jar BedToIntervalList \
        I=IDT_xGen_Exome_Research_Panel_v1.removed_alt_chr.bed \
        O=IDT_xGen_Exome_Research_Panel_v1.removed_alt_chr.bed.interval_list \
        SD=/gscmnt/gc7202/dinglab/common/Reference/A_Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.dict

    $ bsub -q research-hpc -n 1 -R "select[mem>30000] rusage[mem=30000]" -M 30000000 -a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
    -oo /gscmnt/gc2521/dinglab/austins2/tools/GATK4SCNA/db/logs/IDTxGEN_preprocessed_interval_list.log \
    /gscmnt/gc2521/dinglab/austins2/software/jre1.8.0_152/bin/java -jar \
    /gscmnt/gc2521/dinglab/austins2/software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar PreprocessIntervals \
        -L IDT_xGen_Exome_Research_Panel_v1.removed_alt_chr.bed.interval_list  \
        -R /gscmnt/gc7202/dinglab/common/Reference/A_Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.fa \
        --bin-length 0 \
        --interval-merging-rule OVERLAPPING_ONLY \
        --padding 250 \
        -O IDT_xGen_Exome_Research_Panel_v1.removed_alt_chr.bed.target.preprocessed.exome.interval_list

1.  <li value="2">

    Or for WES data launch the docker image using a mapped volume and
    run the command
    `bash /code/src/0.run.gatk4scna.preprocessIntervals.sh -C /code/config/config.gatk4scna.compute1.ini -T /path/to/your/WES/target.bed`

    </li>

2.  Or for WGS data (no target bed file) data run this command:
    `bash /code/src/0.run.gatk4scna.preprocessIntervals.sh -C /code/config/config.gatk4scna.compute1.ini`
    and then filter out blacklisted regions of the genome. The commands
    used for the WGS interval list file in the docker were as follows:

::#

<li value="1">

Generate an initial Picard style interval_list file that can be used to
filtered to generate the final interval list

</li>

    $ bsub -G compute-dinglab -q dinglab -N -R 'select[mem>4GB] rusage[mem=4GB]' -M 4GB -J intervals-GRCh38 -a 'docker(austins2/gatk4scna:test5)' \
        "bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/src/0.run.gatk4scna.preprocessIntervals.sh \
        -C /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini \
        -O /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/db"

::#

<li value="2">

The ENCODE blacklist was downloaded from
[github](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz)
and unzipped for use. These are regions in the human genome with "that
have anomalous, unstructured, or high signal in next-generation
sequencing experiments independent of cell line or experiment. The
removal of the ENCODE blacklist is an essential quality measure when
analyzing functional genomics data." for more information on the
blacklist please see the following publication Amemiya, H.M., Kundaje,
A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic
Regions of the Genome. Sci Rep 9, 9354 (2019).
[2](https://doi.org/10.1038/s41598-019-45839-z)

</li>

::# Convert the blacklist bedfile to an interval list so it can be used
for filtering

    $ bsub -G compute-dinglab -q dinglab -N -R 'select[mem>4GB] rusage[mem=4GB]' -M 4GB -J blacklist \
        -oo /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/db/hg38_blacklist_bed_to_interval_list.log \
        -a 'docker(austins2/gatk4scna:test5)' \
        "/code/jre1.8.0_152/bin/java -jar /storage1/fs1/dinglab/Active/Projects/austins2/software/picard.jar \
        BedToIntervalList I=hg38-blacklist.v2.bed O=hg38-blacklist.v2.bed.interval_list \
        SD=/storage1/fs1/dinglab/Active/Projects/austins2/references/GRCh38.d1.vd1/GRCh38.d1.vd1.dict"

::# Filter out the blacklisted regions from the interval list

    $ bsub -G compute-dinglab -q dinglab -N -R 'select[mem>4GB] rusage[mem=4GB]' -M 4GB -J blacklist_filter \
        -oo /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/db/hg38_blacklist_filter.log \
        -a 'docker(austins2/gatk4scna:test5)' \
        "/code/jre1.8.0_152/bin/java -jar /gatk/gatk-package-4.1.9.0-SNAPSHOT-local.jar \
        FilterIntervals -L GRCh38.d1.vd1.fa.targets.preprocessed.1k.interval_list -XL hg38-blacklist.v2.bed.interval_list \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O GRCh38.d1.vd1.fa.targets.preprocessed.blacklist_filtered.1k.interval_list"

##### Common Biallelic SNP database

This input file contains a whitelist of SNPs that are commonly present
in the population. This file is the same across all versions of the
pipeline. The SNPs documented in it are used to "tabulate counts of the
reference allele and counts of the dominant alternate allele for each
site in a given genomic intervals list". Once done these counts can be
passed to ModelSegments function where the counts will be used to aid in
modeling segments. In this case the original SNP database used is the
af-only-gnomad.hg38.filtered.vsf downloaded from the [broad
gatk-best-practices google cloud
bucket](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).
It was then filtered to
`remove entries with AF <0.01 and AF > 0.2 and entries not marked with 'PASS'`.
It is listed as `'COMMON_BIALLELIC='` in the config file.

##### Protein coding genes bed file

The original genes bed file was been updated to Gencode v34 when I
updated the pipeline to work on compute1. This file is not necessary to
run steps 0-4 which generate the initial segment calls as well as QC
plots for the GATK4SCNA pipeline. It is used in steps 5-6 to map the
copy ratios from the segments to individual gene level calls to easily
investigate CNVs at the individual gene levels. This file was generated
from the [gencode.v34.annotated.gtf
file](https://www.gencodegenes.org/human/release_34.html) as follows.

    awk '{ if ($3 == "gene") { print } }' gencode.v34.annotation.gtf > gencode.v34.annotation.gene_filtered.gtf
    awk '{ if ($12 == "\"protein_coding\"\;") { print } }' gencode.v34.annotation.gene_filtered.gtf > gencode.v34.annotation.gene_filtered.protein_coding.gtf
    awk '{printf("%s\t%s\t%s\t%s\t%s\n",$1,$4,$5,$14,$10);}' gencode.v34.annotation.gene_filtered.protein_coding.gtf > gencode.v34.annotation.gene_filtered.protein_coding.pre-bed
    awk '{gsub(/;/, "", $4); printf("%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5)}' gencode.v34.annotation.gene_filtered.protein_coding.pre-bed > tmp.bed
    awk '{gsub(/;/, "", $5); printf("%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5)}' tmp.bed > tmp2.bed
    awk '{gsub(/"/, "", $4); printf("%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5)}' tmp2.bed > tmp3.bed
    awk -F '[.]' '{printf("%s\n",$1)}' gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID.bed > gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.bed

The intermediate files were filtered to include only gene symbols also
recognized by HGNC and listed as protein-coding in the
hgnc_complete_set_2020-07-01.txt file that was downloaded from the
[https://www.genenames.org/download/archive/ HGNC complete set
Archive](https://www.genenames.org/download/archive/_HGNC_complete_set_Archive "wikilink").

    curl http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/quarterly/tsv/hgnc_complete_set_2020-07-01.txt --output hgnc_complete_set_2020-07-01.txt
    grep "protein-coding gene" hgnc_complete_set_2020-07-01.txt > hgnc_2020-07-01_protein_coding.txt

A python script was then written and used to filter the gencode
annotation file further to only those genes also listed in the HGNC
annotation file. The python script used can be found at
/code/db/remaking_protein-coding_v3/filtering_by_gene_symbol.py inside
of the docker image.

    conda activate py3.9
    python filtering_by_hgnc_gene_name.py hgnc_2020-07-01_protein_coding.txt gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.bed

^This step generates two files:

    gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.ensembl_100_missing_from_hgnc.txt
    and
    gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.protein-coding_hgnc_filtered.txt

Genes names that are listed more than once in the
gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.protein-coding_hgnc_filtered.txt
were then identified using unique_values.py and saved to the file
Gencode_unique.txt.

    python unique_values.py gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.protein-coding_hgnc_filtered.txt > Gencode_unique.txt 

For each instance of a duplicate gene_name that is listed in the
Gencode_unique.txt file only one of them was retained in the final
output file. This was done by manually checking which emsembl gene_ID
was listed in the hgnc_2020-07-01_protein_coding.txt and retaining only
the one present there. Duplicate gene names were not removed if the
gene(s) were present on both the X and Y chromosomes but not on the
autosomes. Then I removed the last column which contains the Ensembl
IDs.

    awk '{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4);}' gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.protein-coding_hgnc_filtered.duplicates_removed.txt > gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.protein-coding_hgnc_filtered.duplicates_removed.ensembl_ID_removed.txt

What is left is only the chromosome, the start and stop positions, and
the gene name in a bed file format.

## Input Files

There are two input files that will vary from run to run. These input
files are tsv formatted tables that contain the locations of Tumor and
normal samples bams as well as sample IDs. **All input files should be
tab-separated!** The first input file is the bam.list file and is
formatted as such:

| CaseID | NormalBam | TumorBam | Disease |
|----|----|----|----|
| Case1 | /path/to/Case1_normal.bam | /path/to/Case1_tumor.bam | PDAC |
| Case2 | NA | /path/to/Case2_tumor.bam | PDAC |
| Case3-sample1 | /path/to/Case3-sample1-by-a-different-name_normal.bam | /path/to/Case3-sample1-by-a-different-name_tumor.bam | PDAC |
| Case3-sample2 | /path/to/Case3-sample2_normal.bam | /path/to/Case3-sample2-by-a-different-name_tumor.bam | PDAC |

It is important to note that for the first input file the disease column
is not actually used. However, it is generally a good idea to include
the disease at the end as depending on the project it is not always
clear where the sampleID, caseIDs, or bam file names where the samples
originate from. This column can thus be used to trace the samples back
to which project they originated from. Additionally all of the IDs
listed in the TumorID column must be unique. So if you are planning to
look at multiple samples from the same case then you will need to
distinguish samples from the same case in this column. The NormalBam
column contains the full path to the normal sample bam file. The values
here do not necessarily need to be unique. The TumorBam column contains
the full path to the tumor bam file (while it does not necessarily need
to be unique I recommend that it should be).

The second input file is the normalBam.list file and is formatted as
such:

| TumorID | NormalBam | NormalID |
|----|----|----|
| Case1 | /path/to/Case1_normal.bam | Normal_ID_for_Case1 |
| Case3-sample1 | /path/to/Case3-sample1-by-a-different-name_normal.bam | Normal_ID_for_Case3 |
| Case3-sample2 | /path/to/Case3-sample2_normal.bam | Normal_ID_for_Case3 |

The TumorID file for the normalBam.list file must be identical to the
TumorID provided in the bam.list file that is provided. In the case
where multiple samples rely on the same normal bam file then only one
entry needs to be provided. The normalBam.list file is mostly used for
QC purposes by generating CNV calls on the normal bams. In the event
that a tumor sample does not have an associated normal then it can be
excluded from the normalBam.list file.

There is a final optional input file for if you would like to plot the
segment call and copy ratios on a per sample basis, while highlighting
the positions of certain genes. This option lets you plot either the
entire genome or specific contigs (chromosomes). This table is the same
as the bam.list file that is described above but it has two additional
column at end of each row. The 5th column contains a comma separated
list of genes that should be highlighted for each sample.

| CaseID | NormalBam | TumorBam | Disease | Gene_csv | Congtigs |
|----|----|----|----|----|----|
| Case1 | /path/to/Case1_normal.bam | /path/to/Case1_tumor.bam | PDAC | NA | NA |
| Case2 | NA | /path/to/Case2_tumor.bam | PDAC | KRAS | NA |
| Case3-sample1 | /path/to/Case3-sample1-by-a-different-name_normal.bam | /path/to/Case3-sample1-by-a-different-name_tumor.bam | PDAC | KRAS,TP53 |  |
| Case3-sample2 | /path/to/Case3-sample2_normal.bam | /path/to/Case3-sample2-by-a-different-name_tumor.bam | PDAC | Not_a_gene | chr1 |
| Case4 | NA | /path/to/Case4_tumor.bam | PDAC | NA | chr1,chr17 |
| Case5 | NA | /path/to/Case5_tumor.bam | PDAC |  | chr12 |
| Case6 | NA | /path/to/Case6_tumor.bam | PDAC |  |  |

The `Gene_csv` column takes a comma separated list of genes like shown
above. If no genes are to be marked on the plot then either leave this
column blank or include the value `NA` in this column. The Contigs
column includes a comma separated list of chromosome names as they
appear in the `Case*.T.modelFinal.seg` output file. These also match the
contig names included in the sequence dictionary file that is associated
with the reference genome. If you want to plot all contigs (aka
chromosomes) in the output file then either leave this column blank for
that sample/case or include the value `NA` for that sample. The
contig(s) can be provided in any order. If one or more contigs are
specified then only those contigs will be included in the output plot.

## Example Commands

<strong> WES un-paired (aka tumor-only) CNV calling example commands for
compute1: </strong>

    # Step 1 (precall)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p precall -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly/WES_allTumorBAM.list -M /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly/WES_normalBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 2 (PON)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p pon -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 3 (CNV calls all using tumor-only version) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p callcn_tumor_only -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly/WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 3.5 (CNV calls normal) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p callnormal -M /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly/WES_normalBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 4 (plot all Tumor) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p plot -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly/WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 4.5 (plot normals) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p plotNormal -M /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly/WES_normalBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 5 (Calls gene-level) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p geneLevel -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly/WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 6 (Merge gene-level files to one file) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p merge -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 5.5 (Calls chr_arm-level) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p chrarmLevel -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly/WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 6.5 (Merge arm-level files to one file) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p chrarmmerge -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini

<strong> WES paired (aka tumor-normal) CNV calling example commands for
compute1: </strong>

    # Step 1 (precall)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p precall -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair/WES_allTumorBAM.list -M /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair/WES_normalBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 2 (PON)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p pon -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 3 (CNV calls all using tumor-only version) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p callcn -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair/WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 3.5 (CNV calls normal) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p callnormal -M /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair/WES_normalBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 4 (plot all Tumor) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p plot -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair/WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 4.5 (plot normals) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p plotNormal -M /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair/WES_normalBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 5 (Calls gene-level) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p geneLevel -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair/WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 6 (Merge gene-level files to one file) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p merge -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 5.5 (Calls chr_arm-level) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p chrarmLevel -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair/WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
    # Step 6.5 (Merge arm-level files to one file) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p chrarmmerge -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wes_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini

<strong> WGS un-paired (aka tumor-only) CNV calling example commands for
compute1: </strong>

    # Step 1 (precall)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p precall -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly/WGS_CNV_gatk_allTumorBam.list -M /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly/WGS_CNV_gatk_normalBam.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 2 (PON)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p pon -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 3 (CNV calls all using tumor-only version) # done
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p callcn_tumor_only -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly/WGS_CNV_gatk_allTumorBam.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 3.5 (CNV calls normal) # run
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p callnormal -M /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly/WGS_CNV_gatk_normalBam.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 4 (plot all Tumor)  # run
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p plot -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly/WGS_CNV_gatk_allTumorBam.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 4.5 (plot normals)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p plotNormal -M /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly/WGS_CNV_gatk_normalBam.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 5 (Calls gene-level) # run
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p geneLevel -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly/WGS_CNV_gatk_allTumorBam.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 6 (Merge gene-level files to one file)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p merge -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 5.5 (Calls chr_arm-level) # run
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p chrarmLevel -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly/WGS_CNV_gatk_allTumorBam.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 6.5 (Merge arm-level files to one file) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p chrarmmerge -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_tonly -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini

<strong> WGS paired (aka tumor-normal) CNV calling example commands for
compute1: </strong>

    # Step 1 (precall)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p precall -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair/WGS_CNV_gatk_pairedBam.list -M /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair/WGS_CNV_gatk_normalBam.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 2 (PON)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p pon -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 3 (CNV calls all using tumor-only version)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p callcn -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair/WGS_CNV_gatk_pairedBam.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 3.5 (CNV calls normal)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p callnormal -M /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair/WGS_CNV_gatk_normalBam.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 4 (plot all Tumor) 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p plot -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair/WGS_CNV_gatk_pairedBam.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 4.5 (plot normals)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p plotNormal -M /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair/WGS_CNV_gatk_normalBam.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 5 (Calls gene-level)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p geneLevel -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair/WGS_CNV_gatk_pairedBam.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 6 (Merge gene-level files to one file)
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.wgs.cnv.compute1.sh -p merge -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 5.5 (Calls chr_arm-level) # run
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p chrarmLevel -t /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair/WGS_CNV_gatk_pairedBam.list -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini
    # Step 6.5 (Merge arm-level files to one file) # 
    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p chrarmmerge -o /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/test_run_v1.1/wgs_pair -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini

## Output Files

**Area under construction! Hard hats required!**</br> This pipeline will
also generate a directory `/path/to/your/run/directory/logs/` where all
log files will be saved. It is recommended to check all log files to
ensure that each step of the pipeline was successfully completed. If on
sample fails but the rest are successfuly, then the pipeline can
continue to be run provided that the problem sample is removed from the
bam.list and normalBam.list file. The output file
merged.gene_level.from_seg.hg38.log2ratio.tsv will contain a matrix
(tsv) of tumor cnv calls. The columns are samples and the rows are
genes. The output file merged.segment_level.hg38.log2ratio.tsv contains
the all of the segment calls for each sample. It is formatted as
follows:

| Sample  | Chromosome | Start   | End     | Num_Probes | Segment_Mean | Call |
|---------|------------|---------|---------|------------|--------------|------|
| Sample1 | chr1       | 794002  | 4117000 | 2859       | 0.103157     | 0    |
| Sample1 | chr1       | 4117002 | 4455000 | 333        | 1.479743     | \+   |
| Sample1 | chr1       | 794002  | 4117000 | 37         | -0.510732    | \-   |

The Num_Probes column lists the number of probes (single nucleotide
polymorphisms) from the input Common Biallelic SNP database that were
appropriately detected within the bounds of the segment. These probes
were then used to calculate the segment mean copy ratio that is reported
in the Segment_Mean column. The values in the Segment_Mean column are
log2ratios from copy numbers that are then centered at 0 using the
calculated mean of neutral segments. Neutral segments are regions
without any detected copy number variations. These are reported as 0 in
the Call column. Amplified segments are regions where more than 2 copies
of the genome have been detected and are reported as a +. Deleted
segments are regions where fewer than 2 copies of the genome have been
detected and are reported as a -.

In addition to this for each individual sample there will be a sample
folder that is generated. Inside of this folder will the the temporary
outputs from steps 1, 3 and 4. Step 2 will generate an additional folder
for the generation of the panel of normals file that is used in the
later steps. Each sample folder will also contain the plots that are
generated at steps 4 and 4.1.

Optional command for final copy number call plotting (this works with
the tumor-only and tumor-normal versions of the pipeline as well as the
WGS and WES version of the pipeline. The final plot is saved to the pdf
called `CaseID_PlotModeledSegments_complete.pdf`. Each part of the
complete plot is also saved to a separate pdf corresponding to that
part. This is done to ensure that the legend is available since it is
sometimes cutoff in the `CaseID_PlotModeledSegments_complete.pdf`. That
said if you open the `CaseID_PlotModeledSegments_complete.pdf` file in
adobe illustrator you will find that the legend is still there. It just
gets cut off by the boundaries of the pdf artboard.

    bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p plotFinal -t /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v2_2024-07-23/v2_WES_geneTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v2_2024-07-23/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini

## Troubleshooting

If you have had issues using this pipeline, then please indicate below
what the problem was as well as a solution if you or others were able to
find one. If linking to an outside source then please copy the relevant
information from that source to this page as sometimes links can change
and posts can go down on outside forums. While troubleshooting I
recommend making a copy on the config file, the wrapper script (either
the WES or WGS one depending on the pipeline you are trying to use), and
the src folder that is listed in the scriptDir of the wrapper script to
your working directory. You can then edit the scripts in that copied
folder.

<strong>Step 1 fail exit code 137 on compute1</strong>  
After running each step it is a good idea to use
`grep "Exited with" logs/*` to make sure that all prior steps
successfully completed. There was one instance where step 1 threw an
error code 137.

    [a.n.southard-smith@compute1-client-4 CNV]$ grep "Exited with" logs/*
    logs/gatk4cn.s1.C3N-01366.WGS.T.log:Exited with exit code 137.
    logs/gatk4cn.s1.C3N-02256.WGS.T.log:Exited with exit code 137.
    logs/gatk4cn.s1.C3N-02599.WGS.T.log:Exited with exit code 137.

According to the [IBM
documentation](https://www.ibm.com/docs/en/spectrum-lsf/10.1.0?topic=SSWRJV_10.1.0/best_practices/Enforcing%20job%20memory%20and%20swap%20with%20Linux%20cgroups.html)
this means that the job was killed due to it reaching a memory limit
`TERM_MEMLIMIT: job killed after reaching LSF memory usage limit.`. I
increased the memory threshold for step 1 from 4GB to 6GB for all jobs.
This was done by changing the beginning of line 105 from
`bash $submitJob 4 1 gatk4cn.s1.${CaseID}` to
`bash $submitJob 6 1 gatk4cn.s1.${CaseID}`. The same was done on line
111. The jobs for this step will now request 6GB of mem regardless of if
all of it is used or not. Additionally I changed the script that is
called `$scriptDir/src/1.run.gatk4scna.collectReadCounts.Normal.sh` so
that Java now requests 6G of memory instead of 4G. For this line 57 was
changed from `-Xmx4g -jar ${GATK} CollectReadCounts \` to
`-Xmx6g -jar ${GATK} CollectReadCounts \`.

<strong> Permission denied for BAM file at step 5 </strong>  
Step 5 does not need the BAM file(s) for each case and should not be
trying to access them. It only neads the CaseID from each line. It has
been found that step 5 in particular does not appropriately separate out
the fields in the file when spaces are used. I don't have time to fix
this right now. In the mean time you should make sure that the fields in
the input files are tab separated as this solves the problem. This can
be quickly done by
`awk '{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4);}' /path/to/your/current/input/bam.list > /path/to/new/tab/separated/bam.list`

<strong> Step 3/3.5 does not generate output files for all samples
</strong>  
There are a number of possible reasons for this to occur. The first
thing to check is that all of the input file paths are correct. If you
are running the tumor-only version of the pipeline you need to check
that there are no NA values present in the `normalBam.list` input file.
This file should only have one entry per normal sample. If there is no
matching normal sample for a case then leave that case out of the
normalBam.list file. These `NA` values for the `bam.list` file can only
be in the second column which is used to specify the matching normal Bam
file for a sample. The `NA` vales should only be present when no
matching normal bam file is available. If there are multiple tumor bam
files for a single case (example: a primary bam file and a matched
metastasis bam file) from the same case then each is to be listed as
it’s own entry in the `bam.list` file. If this case with two tumor bam
files has a single normal bam file then each of those entries in the
`bam.list` file will have the same Normal bam file in the second column
(literally just use the same file path). However, in the normalBam.list
file there should only be a single entry for that multiple tumor samples
case where the single normal bam file is listed once and only once. Do
not repeat normal bam files in the `normalBam.list` input file. This can
cause problems when the output from step3.5 is written to storage1 or
scratch1 and result in either crashes or corrupted output files for
steps 3.5 and 4.5. The only way it is possible to list the same normal
bam multiple times in the `normalBam.list` file is if it is listed under
different file names. If this is going to be done then that normal bam
file should be soft linked somewhere else under a different name and
then it may be listed again under that new file name. I strongly
recommend against this option as it will result in confusion for whoever
has to dig through your results later after they are published.  

### Maintainers

\* Austin Southard-Smith

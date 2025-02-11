There are three versions of the target interval list file that I have generated so far. The first two are compatible with the WES pipeline and are present within the docker image. The third config that currently exists is not present within the docker image and the location must be separately specified when running on compute1.

The target file that is used as the default config file for this docker image was made from the bed file for IDT_xGen_Exome_Research_Panel_v1 for GRCh38. The path to this file can be found in /code/config/config.gatk4scna.compute1.ini inside the docker image. If your WES data was not generated using this target panel then you will need to generate and use your own config file. Copy the current config file to a location outside of the docker image and then change the TargetIntervalList= variable to the path to your target interval_list. Another TargetIntervalList that comes with this docker image is nexterarapidcapture_exome_targetedregions_v1.2.hg38.bed.target.preprocessed.exome.interval_list. This is for WES data from Nextera Rapid Capture Targeted regions v1.2 using GRCh38. The location of the third target interval file, used for the WGS version of the pipeline, can be found within the WGS config file /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/config/config.gatk4scna.wgs.compute1.ini

If your sequencing data was not generated via one of these two target panels then you will need to generate and make you own interval list. In such a case there are 3 options:

    Repeat these steps with your target bed file
        convert the bed file to an interval list using Picard [https://broadinstitute.github.io/picard/]

java -jar picard.jar BedToIntervalList
    I=input.bed
    O=list.interval_list
    SD=reference_sequence.dict

            Preprocess the intervals to prepare for coverage collection with GATK4 [https://gatk.broadinstitute.org/hc/en-us/articles/360050815612-PreprocessIntervals] 

gatk PreprocessIntervals \
    -L targets_C.interval_list \
    -R /gatk/ref/Homo_sapiens_assembly38.fasta \
    --bin-length 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O sandbox/targets_C.preprocessed.interval_list

                The commands used for the IDT WES bed file in the docker were: 

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

    Or for WES data launch the docker image using a mapped volume and run the command bash /code/src/0.run.gatk4scna.preprocessIntervals.sh -C /code/config/config.gatk4scna.compute1.ini -T /path/to/your/WES/target.bed
    Or for WGS data (no target bed file) data run this command: bash /code/src/0.run.gatk4scna.preprocessIntervals.sh -C /code/config/config.gatk4scna.compute1.ini and then filter out blacklisted regions of the genome. The commands used for the WGS interval list file in the docker were as follows:

            Generate an initial Picard style interval_list file that can be used to filtered to generate the final interval list 

$ bsub -G compute-dinglab -q dinglab -N -R 'select[mem>4GB] rusage[mem=4GB]' -M 4GB -J intervals-GRCh38 -a 'docker(austins2/gatk4scna:test5)' \
    "bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/src/0.run.gatk4scna.preprocessIntervals.sh \
    -C /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini \
    -O /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/db"

            The ENCODE blacklist was downloaded from github (https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz )and unzipped for use. These are regions in the human genome with "that have anomalous, unstructured, or high signal in next-generation sequencing experiments independent of cell line or experiment. The removal of the ENCODE blacklist is an essential quality measure when analyzing functional genomics data." for more information on the blacklist please see the following publication Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019) [https://doi.org/10.1038/s41598-019-45839-z ]. 

            Convert the blacklist bedfile to an interval list so it can be used for filtering

$ bsub -G compute-dinglab -q dinglab -N -R 'select[mem>4GB] rusage[mem=4GB]' -M 4GB -J blacklist \
    -oo /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/db/hg38_blacklist_bed_to_interval_list.log \
    -a 'docker(austins2/gatk4scna:test5)' \
    "/code/jre1.8.0_152/bin/java -jar /storage1/fs1/dinglab/Active/Projects/austins2/software/picard.jar \
    BedToIntervalList I=hg38-blacklist.v2.bed O=hg38-blacklist.v2.bed.interval_list \
    SD=/storage1/fs1/dinglab/Active/Projects/austins2/references/GRCh38.d1.vd1/GRCh38.d1.vd1.dict"

            Filter out the blacklisted regions from the interval list

$ bsub -G compute-dinglab -q dinglab -N -R 'select[mem>4GB] rusage[mem=4GB]' -M 4GB -J blacklist_filter \
    -oo /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/wgs/db/hg38_blacklist_filter.log \
    -a 'docker(austins2/gatk4scna:test5)' \
    "/code/jre1.8.0_152/bin/java -jar /gatk/gatk-package-4.1.9.0-SNAPSHOT-local.jar \
    FilterIntervals -L GRCh38.d1.vd1.fa.targets.preprocessed.1k.interval_list -XL hg38-blacklist.v2.bed.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O GRCh38.d1.vd1.fa.targets.preprocessed.blacklist_filtered.1k.interval_list"

Common Biallelic SNP database:
This input file contains a whitelist of SNPs that are commonly present in the population. This file is the same across all versions of the pipeline. The SNPs documented in it are used to "tabulate counts of the reference allele and counts of the dominant alternate allele for each site in a given genomic intervals list". Once done these counts can be passed to ModelSegments function where the counts will be used to aid in modeling segments. In this case the original SNP database used is the af-only-gnomad.hg38.filtered.vsf downloaded from the broad gatk-best-practices google cloud bucket (https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?prefix=&forceOnObjectsSortingFiltering=false ). It was then filtered to remove entries with AF <0.01 and AF > 0.2 and entries not marked with 'PASS'. It is listed as 'COMMON_BIALLELIC=' in the config file. 


#############################
Protein coding genes bed file:

This file is not necessary to run steps 0-4 which generate the initial segment calls as well as QC plots for the GATK4SCNA pipeline. It is used in steps 5-6 to map the copy ratios from the segments to individual gene level calls to easily investigate CNVs at the individual gene levels. This file was generated from the gencode.v34.annotated.gtf file as follows. If you want to use one other than the file that is provided here and listed in the current config file then make a copy of the provided config file and change the path for the ProteinCodingGene variable to the other file you would like to use. The file needs to be in bed format with the gene IDs in column 4 (1-based).

awk '{ if ($3 == "gene") { print } }' gencode.v34.annotation.gtf > gencode.v34.annotation.gene_filtered.gtf
awk '{ if ($12 == "\"protein_coding\"\;") { print } }' gencode.v34.annotation.gene_filtered.gtf > gencode.v34.annotation.gene_filtered.protein_coding.gtf
awk '{printf("%s\t%s\t%s\t%s\t%s\n",$1,$4,$5,$14,$10);}' gencode.v34.annotation.gene_filtered.protein_coding.gtf > gencode.v34.annotation.gene_filtered.protein_coding.pre-bed
awk '{gsub(/;/, "", $4); printf("%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5)}' gencode.v34.annotation.gene_filtered.protein_coding.pre-bed > tmp.bed
awk '{gsub(/;/, "", $5); printf("%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5)}' tmp.bed > tmp2.bed
awk '{gsub(/"/, "", $4); printf("%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5)}' tmp2.bed > tmp3.bed
awk -F '[.]' '{printf("%s\n",$1)}' gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID.bed > gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.bed

The intermediate files were filtered to include only gene symbols also recognized by HGNC and listed as protein-coding in the hgnc_complete_set_2020-07-01.txt file that was downloaded from the HGNC complete set Archive website https://www.genenames.org/download/archive/ . 

curl http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/quarterly/tsv/hgnc_complete_set_2020-07-01.txt --output hgnc_complete_set_2020-07-01.txt
grep "protein-coding gene" hgnc_complete_set_2020-07-01.txt > hgnc_2020-07-01_protein_coding.txt

A python script was then written and used to filter the gencode annotation file further to only those genes also listed in the HGNC annotation file. This brings the number of entries down from 60669 to 19135 and ensures the file follows the bed format. The python script used can be found at /code/db/remaking_protein-coding_v3/filtering_by_gene_symbol.py inside of the docker image. 

python filtering_by_hgnc_gene_name.py hgnc_2020-07-01_protein_coding.txt gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.bed

^this script produces the following two output files.
gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.ensembl_100_missing_from_hgnc.txt
and
gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.protein-coding_hgnc_filtered.txt

Genes names that are listed more than once in the gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.protein-coding_hgnc_filtered.txt were then identified using unique_values.py and saved to the file Gencode_unique.txt.

python unique_values.py gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.protein-coding_hgnc_filtered.txt > Gencode_unique.txt 

For each instance of a duplicate gene_name  that is listed in the Gencode_unique.txt file only one of them was retained in the final output file. This was done by manually checking which emsembl gene_ID was listed in the hgnc_2020-07-01_protein_coding.txt and retaining only the one present there. Duplicate gene names were not removed if the gene(s) were present on both the X and Y chromosomes but not on the autosomes. Then I removed the last column which contains the Ensembl IDs.

awk '{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4);}' gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.protein-coding_hgnc_filtered.duplicates_removed.txt > gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.protein-coding_hgnc_filtered.duplicates_removed.ensembl_ID_removed.txt

What is left is only the chromosome, the start and stop positions, and the gene name in a bed file format.

#!/bin/bash

# Austin Southard-Smith
# a.n.southard-smith@wustl.edu

# 2023-02-27 v2
# 2021-02-02 v0.2
# 2021-01-30 v0.1


## USAGE:
# Note: Please set config.ini before running pipeline

# sh gatk_somatic.cnv.sh -c <config.ini> -p <programname> -t <table> -o <outdir>



submitJob=/storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/lsf_submit.sh
outdir=`pwd`

# set script dir
#scriptDir="/code/"
scriptDir="/storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/"

# set config.ini
config=${scriptDir}/config/config.gatk4scna.mgi.ini
config=${scriptDir}/config/config.gatk4scna.compute1.ini

CNVPON=''
hdf5Dir=''
while getopts "c:p:t:M:d:n:o:" opt; do
  case $opt in
    c)
      config=$OPTARG
      ;;
    p)
      program=$OPTARG
      ;;
    t)
      table=$OPTARG
      ;;
    M)
      normalTable=$OPTARG
      ;;
    d)
      hdf5Dir=$OPTARG
      ;;
    n)
      CNVPON=$OPTARG
      ;;
    o)
      outdir=$OPTARG
      ;;
    \?)
      echo "script usage: $(basename $0) [-t] [-n] " >&2
      exit 1
      ;;
  esac
done


source $config



###############################
##  STEP-1   Pre-call for normal
###############################

#if [[ $program == "precall" ]] || [[ $program == "s1" ]]; then
#
#    if [ ! -e $table ]; then
#        echo "[ERROR] The Table $table does not exist !" 1>&2
#        exit
#    fi
#    
#    if [ ! -d $outdir ]; then
#        echo "[ERROR] The OutDir $outdir does not exist !" 1>&2
#        exit
#    fi
#    
#    
#    sed '1d' $table | cut -f 1-3 | while read caseID normalBam tumorBam
#    do
#        bash $submitJob 4 1 gatk4cn.s1.${caseID} "bash $scriptDir/src/1.run.gatk4scna.collectReadCounts.Normal.sh -C ${config} -N ${caseID} -T Normal -B ${normalBam} -O ${outdir}"
#    done
#
#fi

if [[ $program == "precall" ]] || [[ $program == "s1" ]]; then

    if [ ! -e $table ]; then
        echo "[ERROR] The Table $table does not exist !" 1>&2
        exit
    fi
    
    if [ ! -d $outdir ]; then
        echo "[ERROR] The OutDir $outdir does not exist !" 1>&2
        exit
    fi
    # The normalTable version is meant to be used in the case that there are multiple tumor samples from a single case. Otherwise the calls from  normal sample might be represented more than once in the panel of normals 
    if [[ $normalTable != '' ]]; then
        echo $normalTable
        echo "using normal table for normal calling"
        sed '1d' $normalTable | cut -f 1-3 | while read CaseID normalBam normalID
        do
            bash $submitJob 4 1 gatk4cn.s1.${CaseID} "bash $scriptDir/src/1.run.gatk4scna.collectReadCounts.Normal.sh -C ${config} -N ${CaseID} -T Normal -B ${normalBam} -O ${outdir}"
        done
    else
        echo "using table for normal calling"
        sed '1d' $table | cut -f 1-3 | while read CaseID normalBam tumorBam
        do
            bash $submitJob 4 1 gatk4cn.s1.${CaseID} "bash $scriptDir/src/1.run.gatk4scna.collectReadCounts.Normal.sh -C ${config} -N ${CaseID} -T Normal -B ${normalBam} -O ${outdir}"
        done
    fi

fi



###############################
##  STEP-2   Pool Normal
###############################

if [[ $program == "pon" ]] || [[ $program == "s2" ]]; then
    
    if [[ $hdf5Dir == '' ]]; then
        hdf5Dir=$outdir/PON
    fi
    
    bash $submitJob 8 1 gatk4cn.s2.pon "bash $scriptDir/src/2.run.gatk4scna.createPON.sh -C ${config} -D ${hdf5Dir}"

fi




###############################
##  STEP-3   Call total cn
###############################


if [[ $program == "callcn" ]] || [[ $program == "s3" ]]; then
    
    if [[ $CNVPON == '' ]]; then
        CNVPON="$outdir/PON/gatk4scnaPON.Normal.hdf5"
    fi
    
    sed '1d' $table | cut -f 1-3 | while read caseID normalBam tumorBam
    do    
        bash $submitJob 16 1 gatk4cn.s3.${caseID} "bash $scriptDir/src/3.run.gatk4scna.callCNA.pair.sh -C ${config} -P ${CNVPON} -S ${caseID} -N ${normalBam} -T ${tumorBam} -O ${outdir}"
    done
    
fi

if [[ $program == "callcn_tumor_only" ]] || [[ $program == "s3.1" ]]; then

    if [[ $CNVPON == '' ]]; then
        CNVPON="$outdir/PON/gatk4scnaPON.Normal.hdf5"
    fi

    sed '1d' $table | cut -f 1-3 | while read caseID normalBam tumorBam
    do
        bash $submitJob 16 1 gatk4cn.s3.${caseID} "bash $scriptDir/src/3.run.gatk4scna.callCNA.tumor_only.sh -C ${config} -P ${CNVPON} -S ${caseID} -N ${normalBam} -T ${tumorBam} -O ${outdir}"
    done

fi


###############################
##  STEP-3.5   Call/plot normal cn
###############################


if [[ $program == "callnormal" ]] || [[ $program == "s3.5" ]]; then
    
    if [[ $CNVPON == '' ]]; then
        CNVPON="$outdir/PON/gatk4scnaPON.Normal.hdf5"
    fi
    
    sed '1d' $normalTable | cut -f 1-3 | while read TumorID normalBam normalID
    do    
        bash $submitJob 16 1 gatk4cn.s3.5.${normalID} "bash $scriptDir/src/3.5.run.gatk4scna.callCNA.normal.sh -C ${config} -P ${CNVPON} -S ${TumorID} -B ${normalBam} -N ${normalID} -O ${outdir}"
    done
    
fi

###############################
##  STEP-4   plot tumor initial QC
###############################


if [[ $program == "plot" ]] || [[ $program == "s4" ]]; then
    
    
    sed '1d' $table | cut -f 1-3 | while read CaseID normalBam tumorBam 
    do    
        bash $submitJob 8 1 gatk4cn.s4.${CaseID} "bash $scriptDir/src/4.run.gatk4scna.plot.sh -C ${config} -S ${CaseID} -O ${outdir}"
    done
    
fi

###############################
##  STEP-4.5   plot normal initial QC
###############################


if [[ $program == "plotNormal" ]] || [[ $program == "s4.5" ]]; then
    

    sed '1d' $normalTable | cut -f 1-3 | while read TumorID normalBam NormalID
    do
	bash $submitJob 8 1 gatk4cn.s4.${NormalID} "bash $scriptDir/src/4.5.run.gatk4scna.plot_normal.sh -C ${config} -S ${TumorID} -N ${NormalID} -O ${outdir}"
    done
    
fi

###############################
##  STEP-5   gene-level
###############################

if [[ $program == "geneLevel" ]] || [[ $program == "s5" ]]; then
    
    
    sed '1d' $table | cut -f 1-3 | while read CaseID normalBam tumorBam 
    do  
        DIR=${outdir}/${CaseID}
	bash $submitJob 1 1 gatk4cn.s5.${CaseID} "${PYTHON3} $scriptDir/src/segment_to_geneLevel_v4.py --prefix ${CaseID}.T --name ${CaseID} --seg ${DIR}/${CaseID}.T.called.igv.seg --gene ${ProteinCodingGene} -o ${DIR}"
    done
    
fi

###############################
##  STEP-5.5   chromosome-arm-level
###############################

if [[ $program == "chrarmLevel" ]] || [[ $program == "s5.5" ]]; then


    sed '1d' $table | cut -f 1-3 | while read CaseID normalBam tumorBam 
    do
        DIR=${outdir}/${CaseID}
        bash $submitJob 1 1 gatk4cn.chrarm.s5.5.${CaseID} "${PYTHON3} $scriptDir/src/segment_to_chr_arm_level_v4.py  --prefix ${CaseID}.T --name ${CaseID} --seg ${DIR}/${CaseID}.T.called.igv.seg --band $scriptDir/db/cytoBand.txt -o ${DIR}"
    done

fi


###############################
##  STEP-6  merge gene-level
###############################

if [[ $program == "merge" ]] || [[ $program == "s6" ]]; then
    
    bash $submitJob 2 1 gatk4cn.s6.${caseID} "${PYTHON3} $scriptDir/src/mergeMultipleFilesToOne_v2.py ${outdir}"
    
fi

###############################
##  STEP-6  merge chromosome-arm-level
###############################

if [[ $program == "chrarmmerge" ]] || [[ $program == "s6.5" ]]; then

    bash $submitJob 2 1 gatk4cn.chrarm.s6.5.${caseID} "${PYTHON3} $scriptDir/src/mergeMultiple_chr_arm_FilesToOne.py ${outdir}"

fi

###############################
##  STEP-7  plot final
###############################

if [[ $program == "plotFinal" ]] || [[ $program == "s7" ]]; then

    sed '1d' $table | cut -f 1-6 | while read CaseID normalBam tumorBam disease gene_csv plot_contigs
    do
        DIR=${outdir}/${CaseID}
	mkdir -p ${DIR}
	mkdir -p ${DIR}/plots
        if [[ ${gene_csv} == "NA" ]] || [[ ${gene_csv} == "" ]]; then
            gene_csv="no_gene_provided"
	fi
	if [[ ${plot_contigs} != "NA" ]] && [[ ${plot_contigs} != "" ]]; 
	then
	    bsub -G compute-dinglab -q dinglab -oo ${outdir}/logs/gatk4cn.s7.${CaseID}.log -N -a 'docker(austins2/ggplot_gatk4scna:v.2024.08.19)' \
            "Rscript $scriptDir/src/plot_modeled_CNV_segments_to_pdf.R \
            --sample_name ${CaseID}.T \
            --denoised_copy_ratios_file ${DIR}/${CaseID}.T.denoisedCR.tsv \
            --copy_ratio_call_file ${DIR}/${CaseID}.T.called.seg \
            --modeled_segments_file ${DIR}/${CaseID}.T.modelFinal.seg \
            --allelic_counts_file ${DIR}/${CaseID}.T.hets.tsv \
            --sequence_dictionary ${GENOME_DICT} \
            --cytoband_file /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/db/cytoBand.txt \
            --gene_location_file /storage1/fs1/dinglab/Active/Projects/austins2/tools/inferCNV/genes.Cellranger-2020-A.refdata-gex-GRCh38-2020-A.gene_filterd.remove_excess_columns.trimmed.no_duplicates.txt \
            --maximum_copy_ratio 4 \
            --point_size_copy_ratio 0.2 \
            --point_size_allele_fraction 0.4 \
            --output_folder ${DIR}/plots/ \
            --plot_genes ${gene_csv} \
            --retain_contigs ${plot_contigs}"
        else
            bsub -G compute-dinglab -q dinglab -oo ${outdir}/logs/gatk4cn.s7.${CaseID}.log -N -a 'docker(austins2/ggplot_gatk4scna:v.2024.08.19)' \
            "Rscript $scriptDir/src/plot_modeled_CNV_segments_to_pdf.R \
            --sample_name ${CaseID}.T \
            --denoised_copy_ratios_file ${DIR}/${CaseID}.T.denoisedCR.tsv \
            --copy_ratio_call_file ${DIR}/${CaseID}.T.called.seg \
            --modeled_segments_file ${DIR}/${CaseID}.T.modelFinal.seg \
            --allelic_counts_file ${DIR}/${CaseID}.T.hets.tsv \
            --sequence_dictionary ${GENOME_DICT} \
            --cytoband_file /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/db/cytoBand.txt \
            --gene_location_file /storage1/fs1/dinglab/Active/Projects/austins2/tools/inferCNV/genes.Cellranger-2020-A.refdata-gex-GRCh38-2020-A.gene_filterd.remove_excess_columns.trimmed.no_duplicates.txt \
            --maximum_copy_ratio 4 \
            --point_size_copy_ratio 0.2 \
            --point_size_allele_fraction 0.4 \
            --output_folder ${DIR}/plots/ \
            --plot_genes ${gene_csv}"
        fi
    done
fi


# Austin Southard-Smith

# 2021-01-21

# USAGE

# set 16Gb memory for MGI-server

CONFIG='/code/config/config.gatk4scna.mgi.ini'
OUTDIR='.'

min_contigLen=10000000  # based on this len to filter chr

while getopts "C:S:N:T:P:L:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    S)
      caseName=$OPTARG
      ;;
    N)
      BAM_Normal=$OPTARG
      ;;
    T)
      BAM_Tumor=$OPTARG
      ;;
    P)
      CNVPON=$OPTARG
      ;;
    L)
      min_contigLen=$OPTARG
      ;;    
    O)
      OUTDIR=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


source $CONFIG



# make outdir folder
OUT=$OUTDIR/$caseName

mkdir -p $OUT


##---------- TotalCNA Tumor BAM
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx16g -jar ${GATK} CollectReadCounts \
      -I ${BAM_Tumor} \
      -L ${TargetIntervalList} \
      --interval-merging-rule OVERLAPPING_ONLY \
      -O ${OUT}/${caseName}.T.counts.hdf5



${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx16g -jar ${GATK} DenoiseReadCounts \
      -I ${OUT}/${caseName}.T.counts.hdf5 \
      --count-panel-of-normals ${CNVPON} \
      --standardized-copy-ratios ${OUT}/${caseName}.T.standardizedCR.tsv \
      --denoised-copy-ratios ${OUT}/${caseName}.T.denoisedCR.tsv


# Need large memory
##---------- AlleleCNA Normal BAMs
#${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
#      -Dsamjdk.use_async_io_write_samtools=true \
#      -Dsamjdk.use_async_io_write_tribble=false \
#      -Dsamjdk.compression_level=2 \
#      -Xmx16g -jar ${GATK} CollectAllelicCounts \
#      -L ${COMMON_BIALLELIC} \
#      -I ${BAM_Normal} \
#      -R ${GENOME} \
#      -O ${OUT}/${caseName}.N.allelicCounts.tsv



##---------- AlleleCNA Tumor BAMs
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx16g -jar ${GATK} CollectAllelicCounts \
      -L ${COMMON_BIALLELIC} \
      -I ${BAM_Tumor} \
      -R ${GENOME} \
      -O ${OUT}/${caseName}.T.allelicCounts.tsv

      

##---------- AlleleCNA Nor & Tum 
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx16g -jar ${GATK} ModelSegments \
      --denoised-copy-ratios ${OUT}/${caseName}.T.denoisedCR.tsv \
      --allelic-counts ${OUT}/${caseName}.T.allelicCounts.tsv \
      --output ${OUT} \
      --output-prefix ${caseName}.T
    

    
##---------- Segment 
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx16g -jar ${GATK} CallCopyRatioSegments \
      --input ${OUT}/${caseName}.T.cr.seg \
      --output ${OUT}/${caseName}.T.called.seg




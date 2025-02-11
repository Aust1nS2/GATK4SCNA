
# Austin Southard-Smith

# 2021-12-08

# USAGE

# set 16Gb memory for compute1-server

CONFIG='/code/config/config.gatk4scna.mgi.ini'
OUTDIR='.'

#min_contigLen=10000000  # plotting has been moved to script 4 so this is no longer needed

while getopts "C:S:B:N:P:L:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    S)
      caseName=$OPTARG
      ;;
    B)
      BAM_Normal=$OPTARG
      ;;
    N)
      normalID=$OPTARG
      ;;
    P)
      CNVPON=$OPTARG
      ;;
#    L)
#      min_contigLen=$OPTARG # plotting has been moved to script 4 so this is no longer needed
#      ;;    
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
mkdir -p ${OUT}/normal

##---------- TotalCNA Tumor BAM
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx16g -jar ${GATK} CollectReadCounts \
      -I ${BAM_Normal} \
      -L ${TargetIntervalList} \
      --interval-merging-rule OVERLAPPING_ONLY \
      -O ${OUT}/normal/${caseName}.N.counts.hdf5

${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx16g -jar ${GATK} DenoiseReadCounts \
      -I ${OUT}/normal/${caseName}.N.counts.hdf5 \
      --count-panel-of-normals ${CNVPON} \
      --standardized-copy-ratios ${OUT}/normal/${caseName}.N.standardizedCR.tsv \
      --denoised-copy-ratios ${OUT}/normal/${caseName}.N.denoisedCR.tsv


# Need large memory
##---------- AlleleCNA Normal BAMs
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx16g -jar ${GATK} CollectAllelicCounts \
      -L ${COMMON_BIALLELIC} \
      -I ${BAM_Normal} \
      -R ${GENOME} \
      -O ${OUT}/normal/${caseName}.N.allelicCounts.tsv

##---------- AlleleCNA Normal (pretend it's tumor) BAMs
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx16g -jar ${GATK} CollectAllelicCounts \
      -L ${COMMON_BIALLELIC} \
      -I ${BAM_Normal} \
      -R ${GENOME} \
      -O ${OUT}/normal/${normalID}.N.allelicCounts.tsv

      

##---------- AlleleCNA Nor & Tum 
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx16g -jar ${GATK} ModelSegments \
      --denoised-copy-ratios ${OUT}/normal/${caseName}.N.denoisedCR.tsv \
      --allelic-counts ${OUT}/normal/${normalID}.N.allelicCounts.tsv \
      --normal-allelic-counts ${OUT}/normal/${caseName}.N.allelicCounts.tsv \
      --output ${OUT}/normal \
      --output-prefix ${caseName}.N
    

    
##---------- Segment 
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx16g -jar ${GATK} CallCopyRatioSegments \
      --input ${OUT}/normal/${caseName}.N.cr.seg \
      --output ${OUT}/normal/${caseName}.N.called.seg


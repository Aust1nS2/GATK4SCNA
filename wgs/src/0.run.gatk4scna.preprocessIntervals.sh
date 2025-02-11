
# originally written by Hua Sun. Modified by Austin Southard-Smith

# 2021-01-18

# USAGE

# Target (no header)
# <chr> <start> <stop>

# set 4Gb memory for compute1-server

CONFIG='/code/config/config.gatk4scna.mgi.ini'
TARGET=''
OUTDIR='.'

while getopts "C:T:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    T)
      TARGET=$OPTARG
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
mkdir -p $OUTDIR


# If there is no exome target file (i.e. if the data is whole genome)
if [[ $TARGET == "" ]]; then

    ${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx4g -jar ${GATK} PreprocessIntervals \
      -R ${GENOME} \
      --bin-length 1000 \
      --interval-merging-rule OVERLAPPING_ONLY \
      -O ${OUTDIR}/targets.preprocessed.1k.interval_list
fi
# the recommended default bin length is 1000 for WGS (https://gatk.broadinstitute.org/hc/en-us/articles/360035531092--How-to-part-I-Sensitively-detect-copy-ratio-alterations-and-allelic-segments). (Hua had this set to 5000 previously)



# If there is a exome target file
if [[ $TARGET != "" ]]; then

    ${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx4g -jar ${GATK} PreprocessIntervals \
      -R ${GENOME} \
      -L ${TARGET} \
      --bin-length 0 \
      --interval-merging-rule OVERLAPPING_ONLY \
      -O ${OUTDIR}/targets.preprocessed.exome.interval_list
fi



${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx4g -jar ${GATK} AnnotateIntervals \
      -R ${GENOME} \
      -L ${OUTDIR}/targets.preprocessed.exome.interval_list \
      --interval-merging-rule OVERLAPPING_ONLY \
      -O ${OUTDIR}/targets.preprocessed.exome.annotated.interval_list


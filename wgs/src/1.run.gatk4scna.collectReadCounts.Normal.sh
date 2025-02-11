
# Austin Southard-Smith

# 2021-01-18

# USAGE
# bash run.sh -C contig.ini -N samplename -T Normal -B /path/bam -O /path/outfolder

# set 4Gb memory for compute1-server

# The "TYPE" should be Tumor/Normal

CONFIG='/code/config/config.gatk4scna.mgi.ini'
OUTDIR='.'

while getopts "C:N:T:B:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    N)
      NAME=$OPTARG
      ;;
    T)
      TYPE=$OPTARG
      ;;
    B)
      BAM=$OPTARG
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

OUT=$OUTDIR/PON
mkdir -p $OUT


if [[ $TYPE == "Normal" ]]; then

    ${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx6g -jar ${GATK} CollectReadCounts \
      -I ${BAM} \
      -L ${TargetIntervalList} \
      --interval-merging-rule OVERLAPPING_ONLY \
      -O ${OUT}/${NAME}.N.counts.hdf5
 
fi 




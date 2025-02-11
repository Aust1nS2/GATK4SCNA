
# Austin Southard-Smith

# 2021-01-18

# USAGE
# bassh run.sh -C config.inf -D path/hdf5_dir
# set 8Gb memory for MGI-server

CONFIG='/code/config/config.gatk4scna.mgi.ini'

while getopts "C:D:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    D)
      DIR=$OPTARG
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


all_normal_hdf5=`ls ${DIR}/*.N.counts.hdf5 | while read file; do printf -- "-I $file "; done`

${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx128g -jar ${GATK} CreateReadCountPanelOfNormals \
      --minimum-interval-median-percentile 5.0 \
      -O ${DIR}/gatk4scnaPON.Normal.hdf5 \
      $all_normal_hdf5





# Austin Southard-Smith

# 2021-01-21

# USAGE

# set 8Gb memory for compute1-server

CONFIG='/code/config/config.gatk4scna.mgi.ini'
OUTDIR='.'

min_contigLen=10000000  # based on this len to filter chr

while getopts "C:S:N:L:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    S)
      caseName=$OPTARG
      ;;
    N)
      normalID=$OPTARG
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
OUT=${OUTDIR}/${caseName}

mkdir -p ${OUT}


##---------- Plot-TotalCN
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx8g -jar ${GATK} PlotDenoisedCopyRatios \
      --standardized-copy-ratios ${OUT}/normal/${caseName}.N.standardizedCR.tsv \
      --denoised-copy-ratios ${OUT}/normal/${caseName}.N.denoisedCR.tsv \
      --sequence-dictionary ${GENOME_DICT} \
      --minimum-contig-length ${min_contigLen} \
      --output ${OUT}/plots \
      --output-prefix ${normalID}.N.PlotDenoisedCopyRatios \
      --tmp-dir ${OUT}


##---------- Plot-alleleCN
${JAVA}  -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -jar ${GATK} PlotModeledSegments \
      --denoised-copy-ratios ${OUT}/normal/${caseName}.N.denoisedCR.tsv \
      --allelic-counts ${OUT}/normal/${caseName}.N.hets.tsv \
      --segments ${OUT}/normal/${caseName}.N.modelFinal.seg \
      --sequence-dictionary ${GENOME_DICT} \
      --minimum-contig-length ${min_contigLen} \
      --output ${OUT}/plots \
      --output-prefix ${normalID}.N.PlotModeledSegments \
      --tmp-dir ${OUT}






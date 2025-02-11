
# Austin Southard-Smith

# 2021-01-21

# USAGE

# set 8Gb memory for compute1-server

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

OUT_STANDARD=${OUTDIR}/${caseName}/${caseName}.T.standardizedCR.tsv
OUT_DENOISED=${OUTDIR}/${caseName}/${caseName}.T.denoisedCR.tsv
OUT_PLOTS=${OUTDIR}/${caseName}/plots
OUT_PREFIX=${caseName}.T.PlotDenoisedCopyRatios

##---------- Plot-TotalCN
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx8g -jar ${GATK} PlotDenoisedCopyRatios \
      --standardized-copy-ratios ${OUT_STANDARD} \
      --denoised-copy-ratios ${OUT_DENOISED} \
      --sequence-dictionary ${GENOME_DICT} \
      --minimum-contig-length ${min_contigLen} \
      --output ${OUT_PLOTS} \
      --output-prefix ${OUT_PREFIX} \
      --tmp-dir ${OUT}



##---------- Plot-alleleCN
${JAVA}  -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -jar ${GATK} PlotModeledSegments \
      --denoised-copy-ratios ${OUT}/${caseName}.T.denoisedCR.tsv \
      --allelic-counts ${OUT}/${caseName}.T.hets.tsv \
      --segments ${OUT}/${caseName}.T.modelFinal.seg \
      --sequence-dictionary ${GENOME_DICT} \
      --minimum-contig-length ${min_contigLen} \
      --output ${OUT}/plots \
      --output-prefix ${caseName}.T.PlotModeledSegments \
      --tmp-dir ${OUT}
    




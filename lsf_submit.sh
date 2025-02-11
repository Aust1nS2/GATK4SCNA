#!/bin/sh
# originally written by Hua Sun. Modified by Austin Southard-Smith
## USAGE
## sh lsf_submit.sh <memory int> <threads int> <any name> <any command>

MEM=$1;shift
THREADS=$1;shift
NAME=$1;shift

DIR=`pwd`
mkdir -p $DIR/logs

export LSF_DOCKER_PRESERVE_ENVIRONMENT=false
bsub -G compute-dinglab -q dinglab -M ${MEM}000000 -n ${THREADS} -R "select[mem>${MEM}000] rusage[mem=${MEM}000]" -oo $DIR/logs/${NAME}.log -eo $DIR/logs/${NAME}.err -a "docker(austins2/gatk4scna:v1.1)" "$@"

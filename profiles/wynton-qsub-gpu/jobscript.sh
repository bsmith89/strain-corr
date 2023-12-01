#!/bin/sh
#
#$ -clear
# properties = {properties}
#
# ## Pre-job summary, if running as a job
# [[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID" >&2
echo $HOSTNAME $JOB_ID $JOB_NAME
#
# # # This *should* prevent CUDA from
# # # using GPUs that were not assigned.
# # # FIXME: How to be sure?
# # # FIXME: Why is SGE_GPU not set??
#
export CUDA_VISIBLE_DEVICES=${{SGE_GPU:-0}}
echo CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES SGE_GPU=$SGE_GPU

# source ./env
export PROJECT_ROOT=$PWD
#
# mkdir -p /scratch/bsmith
#
{exec_job}
_status=$?

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID" >&2
exit $_status

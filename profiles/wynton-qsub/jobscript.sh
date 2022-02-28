#!/bin/sh
# properties = {properties}
#
# ## Pre-job summary, if running as a job
# [[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID" >&2
#
# # # This *should* prevent CUDA from
# # # using GPUs that were not assigned.
# # # FIXME: How to be sure?
# # # FIXME: Why is SGE_GPU not set??
#
# source ./env
#
# mkdir -p /scratch/bsmith
#
{exec_job}
# _status=$?
#
# ## End-of-job summary, if running as a job
# [[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID" >&2
#
# env
#
# exit $_status

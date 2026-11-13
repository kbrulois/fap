## fap on the Sherlock cluster: the one place that says where its R library is.
##
## Every fap script starts with library(fap) and nothing else; R finds fap (and
## ligandFinder, which fap imports) through R_LIBS. Shell lines and the sbatch
## commands fap writes use $R_LIBS/fap/... for installed scripts, and sbatch
## jobs inherit R_LIBS from the shell that submits them.
##
## Set it up once, from a Sherlock login shell:
##   echo 'source /home/groups/ebutcher/programs/pipeline/R_libs4.1/fap/cluster/sherlock_env.sh' >> ~/.bashrc
##
## A different library: change it here, reinstall fap, and nothing else changes.
## It must be a single directory: the scripts build paths as $R_LIBS/fap/...
export R_LIBS=/home/groups/ebutcher/programs/pipeline/R_libs4.1

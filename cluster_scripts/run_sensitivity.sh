# Recording date and time of inference
INF=$LSB_JOBNAME
echo "Performing inference $INF"
START=$( date '+%F_%H:%M:%S' )
echo $START

echo "Number: $1"
echo "Condition: $2"

Rscript ./simulations/scripts/sensitivity_mcmc.R $1 $2

# Recoding end date and time of inference
echo "Finished inference $INF"
END=$( date '+%F_%H:%M:%S' )
echo $END

T=$(printf '\t')
echo "$INF$T$START$T$END" >> /storage1/fs1/michael.landis/Active/empire/joblogs/run_log.txt

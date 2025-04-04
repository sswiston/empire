# Set up volumes
LSF_DOCKER_VOLUMES="/storage1/fs1/michael.landis/Active/ellipses:/storage1/fs1/michael.landis/Active/ellipses"
JOBDIR="/storage1/fs1/michael.landis/Active/ellipses/joblogs"

# Set up analysis
if [ $# == 2 ]; then
	CONDITION=$2
	NUMBER=$1
	NAME="$CONDITION.rase.$NUMBER"
fi

# Submit job
bsub -G compute-michael.landis \
-g /k.swiston/ellipses \
-cwd /storage1/fs1/michael.landis/Active/ellipses/ \
-o $JOBDIR/$NAME.stdout.txt \
-J $NAME \
-q general \
-n 1 -M 4GB -R "rusage [mem=4GB] span[hosts=1]" \
-a 'docker(sswiston/ellipses:2)' /bin/bash /storage1/fs1/michael.landis/Active/ellipses/cluster_scripts/run_rase.sh $NUMBER $CONDITION

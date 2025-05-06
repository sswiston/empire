RUN_LIST=$(seq -w 001 100)

for i in ${RUN_LIST[@]}
do
	source run_noisy_job.sh ${i} "full"
done

# OPTIONS:
# "under_prior" > experiment on a small tree under the prior
# "no_autmentation" > coverage experiment on small trees without data augmentation moves
# "small" > coverage experiment on small trees, including data augmentation
# "full" > full coverage experiment, including large trees

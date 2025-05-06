RUN_LIST_1=$(seq -w 001 010)
RUN_LIST_2=$(seq -w 001 005)

for i in ${RUN_LIST_1[@]}
do
	source run_skinks_job.sh ${i} standard
done

for i in ${RUN_LIST_2[@]}
do
	source run_skinks_job.sh ${i} shuffled
	source run_skinks_job.sh ${i} under_prior
done

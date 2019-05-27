for dfile in {100..103}
do
    for model in {1..4}
    do
        python automatic_jobs.py $dfile $model
        # python edit_job.py $dfile $model
        if [ $? -lt 0 ]
        then
            for iter in {1..$?}
            do
                sbatch ciecc
            done
        fi
    done
done
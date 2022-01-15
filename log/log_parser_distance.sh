#!/bin/bash
# rm all_logs_Dense seq_logs_Dense.csv MPI_logs_Dense.csv
# rm all_logs_Sparse seq_logs_Sparse.csv MPI_logs_Sparse.csv

logfile=strong/big_job/best_params_nx8096.csv
rm $logfile

echo "np,time,ite,resi" >> $logfile

for np in 1 2 4 6 8 12 16 20 24 28 32 36 40 48 56 60 64 70 72 74 78 80 82 84
do
  filename=strong/big_job/log-param${np}

  time=$(grep "took" ${filename} | cut -d " " -f3)
  ite=$(grep "Iterations =" ${filename} | cut -d "=" -f2 | cut -d " " -f2)
  resi=$(grep "Norm =" ${filename} | cut -d "=" -f2 | cut -d " " -f2)

  echo "${np},${time},${ite},${resi}" >> $logfile
done

echo "Done !"

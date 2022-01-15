#!/bin/bash
# rm all_logs_Dense seq_logs_Dense.csv MPI_logs_Dense.csv
# rm all_logs_Sparse seq_logs_Sparse.csv MPI_logs_Sparse.csv

logfile=strong/params/best_params_nx8096.csv
rm $logfile

echo "paramvalue,time,ite,resi" >> $logfile

for param in 1 2 4 6 8 12 16 20 24 28
do
  filename=strong/best_params/log-param${param}

  time=$(grep "took" ${filename} | cut -d " " -f3)
  ite=$(grep "Iterations =" ${filename} | cut -d "=" -f2 | cut -d " " -f2)
  resi=$(grep "Norm =" ${filename} | cut -d "=" -f2 | cut -d " " -f2)

  echo "${param},${time},${ite},${resi}" >> $logfile
done

echo "Done !"

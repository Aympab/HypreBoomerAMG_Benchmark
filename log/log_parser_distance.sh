#!/bin/bash
# rm all_logs_Dense seq_logs_Dense.csv MPI_logs_Dense.csv
# rm all_logs_Sparse seq_logs_Sparse.csv MPI_logs_Sparse.csv

logfile=strong/aggressive_level/agglevel_np28_nx8096.csv
rm $logfile

echo "paramvalue,time,ite,resi" >> $logfile

for param in 0 1 2 3 4 5 7 8 9 10 11 12 16 20
do
  filename=strong/aggressive_level/log-agglevel${param}

  time=$(grep "took" ${filename} | cut -d " " -f3)
  ite=$(grep "Iterations =" ${filename} | cut -d "=" -f2 | cut -d " " -f2)
  resi=$(grep "Norm =" ${filename} | cut -d "=" -f2 | cut -d " " -f2)

  echo "${param},${time},${ite},${resi}" >> $logfile
done

echo "Done !"

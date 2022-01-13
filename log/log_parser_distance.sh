#!/bin/bash
# rm all_logs_Dense seq_logs_Dense.csv MPI_logs_Dense.csv
# rm all_logs_Sparse seq_logs_Sparse.csv MPI_logs_Sparse.csv

logfile=node_distance/dist_np2_nx1024.csv
rm $logfile

echo "dist,time,ite,resi" >> $logfile

for dist in 0 1 2 3
do
  filename=node_distance/log-distance-${dist}

  time=$(grep "took" ${filename} | cut -d " " -f3)
  ite=$(grep "Iterations =" ${filename} | cut -d "=" -f2 | cut -d " " -f2)
  resi=$(grep "Norm =" ${filename} | cut -d "=" -f2 | cut -d " " -f2)

  echo "${dist},${time},${ite},${resi}" >> $logfile
done

echo "Done !"

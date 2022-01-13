#!/bin/bash
# rm all_logs_Dense seq_logs_Dense.csv MPI_logs_Dense.csv
# rm all_logs_Sparse seq_logs_Sparse.csv MPI_logs_Sparse.csv
rm weak_logs.csv

echo "np,n,time,ite,resi" >> weak_logs.csv
for np in 1 2 4 8 16 32 64 128
do
  for n in 128 256 512 1024 2048 4096 8192 16384
  do
    time=$(grep "took" log-emilia-NPI${np} | cut -d " " -f3)
    ite=$(grep "Iterations =" log-emilia-NPI${np} | cut -d "=" -f2 | cut -d " " -f2)
    resi=$(grep "Norm =" log-emilia-NPI${np} | cut -d "=" -f2 | cut -d " " -f2)

    echo "${np},${n},${time},${ite},${resi}" >> weak_logs.csv
  done
done

echo "Done !"

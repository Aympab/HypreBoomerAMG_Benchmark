#!/bin/bash
# rm all_logs_Dense seq_logs_Dense.csv MPI_logs_Dense.csv
# rm all_logs_Sparse seq_logs_Sparse.csv MPI_logs_Sparse.csv
rm emilia_strong_logs.csv

echo "np,time,ite,resi" >> emilia_strong_logs.csv
for np in 1 2 4 6 8 12 16 20 24 28 32 36 40 48 56
do
  time=$(grep "took" log-emilia-NPI${np} | cut -d " " -f3)
  ite=$(grep "Iterations =" log-emilia-NPI${np} | cut -d "=" -f2 | cut -d " " -f2)
  resi=$(grep "Norm =" log-emilia-NPI${np} | cut -d "=" -f2 | cut -d " " -f2)

  echo "${np},${time},${ite},${resi}" >> emilia_strong_logs.csv
done

echo "Done !"

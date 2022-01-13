#!/bin/bash
# rm all_logs_Dense seq_logs_Dense.csv MPI_logs_Dense.csv
# rm all_logs_Sparse seq_logs_Sparse.csv MPI_logs_Sparse.csv

filename=weak/weak_start256.csv
rm $filename

echo "np,nx,time,ite,resi" >> $filename

#array=( 1 2 4 8 16 28 32 56 64 84 )
#array2=( 100 200 400 800 1600 2800 3200 5600 6400 8400)

#array=( 1 2 4 8 16 28 32 56 64 84 )
#array2=( 33 66 132 264 528 924 1056 1848 2112 2772)

array=( 1 2 4 8 16 28 32 56 64 84 )
array2=( 256 512 1024 2048 4096 7168 8192 14336 16752 21504)

for i in "${!array[@]}"; do
  np=${array[i]}
  nx=${array2[i]}

  time=$(grep "took" weak/np${np}_nx${nx}.logs | cut -d " " -f3)
  ite=$(grep "Iterations =" weak/np${np}_nx${nx}.logs | cut -d "=" -f2 | cut -d " " -f2)
  resi=$(grep "Norm =" weak/np${np}_nx${nx}.logs | cut -d "=" -f2 | cut -d " " -f2)

  echo "${np},${nx},${time},${ite},${resi}" >> $filename
done

echo "Done !"

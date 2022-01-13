#!/bin/bash
# rm all_logs_Dense seq_logs_Dense.csv MPI_logs_Dense.csv
# rm all_logs_Sparse seq_logs_Sparse.csv MPI_logs_Sparse.csv
for np in 1 2 4 8 16 24
do
    for nx in 128 256 512 1024 2048 4096
    do
        grep DenseMV_ DenseMV/log-NX${nx}-NPI${np} | awk '{ printf ""'${nx}'"_"'${np}'""; print }' | sed 's/DenseMV//' | sed 's/:/_/g' >> all_logs_Dense
        grep SpMV_ SpMV/log-NX${nx}-NPI${np} | awk '{ printf ""'${nx}'"_"'${np}'""; print }' | sed 's/SpMV//' | sed 's/:/_/g' >> all_logs_Sparse
    done
done


nps=( 1 2 4 8 16 24 )
nxs=( 128 256 512 1024 2048 4096 )

$MPIRUN -np ${np}  $MONPROGCSR --nx ${nx} > logs/SpMV/log-NX${nx}-NPI${np}

for i in "${!nps[@]}"; do
    printf "%s is in %s\n" "${nps[i]}" "${nxs[i]}"
done

grep DenseMV_ DenseMV/log-NX${nx}-NPI${np} | awk '{ printf ""'${nx}'"_"'${np}'""; print }' | sed 's/DenseMV//' | sed 's/:/_/g' >> all_logs_Dense

for np in 1 2 4 8 16 24
do
    for nx in 128 256 512 1024 2048 4096
    do
        grep n${nps}_np${nxs}
    done
done


sed -i -e 's/_/;/g' all_logs_Dense
grep seq all_logs_Dense > seq_logs_Dense.csv
grep MPI all_logs_Dense > MPI_logs_Dense.csv
sed -i -e 's/;MPI//' MPI_logs_Dense.csv
sed -i -e 's/;seq//' seq_logs_Dense.csv
sed -i -e '1i nx;np;time' seq_logs_Dense.csv
sed -i -e '1i nx;np;num_p;time' MPI_logs_Dense.csv

sed -i -e 's/_/;/g' all_logs_Sparse
grep seq all_logs_Sparse > seq_logs_Sparse.csv
grep MPI all_logs_Sparse > MPI_logs_Sparse.csv
sed -i -e 's/;MPI//' MPI_logs_Sparse.csv
sed -i -e 's/;seq//' seq_logs_Sparse.csv
sed -i -e '1i nx;np;time' seq_logs_Sparse.csv
sed -i -e '1i nx;np;num_p;time' MPI_logs_Sparse.csv

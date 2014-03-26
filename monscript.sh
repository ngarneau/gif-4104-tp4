#PBS -A colosse-users
#PBS -l walltime=600
#PBS -l nodes=1:ppn=8

module load compilers/gcc/4.8.0
cd "${PBS_O_WORKDIR}"

./tp4

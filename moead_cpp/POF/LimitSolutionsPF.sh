# $1 -> number of objective functions  2
# $2 -> limit of solutions to return   200
# $3 -> input file                     fichero salida

# ALL FUNCTIONS ARE FOR MINIMIZATION
# IN CASE OF MAXIMIZATION, just multiply by -1

if [ $# -ne 3 ];
then
	echo " ********** ERROR: $0 <nobj> <LimitOfSolutions> <input-file>"
	exit
fi

nobj=$1
nSolutions=$2
file=$3

if [ ! -s $file ]; then
  echo "[ERROR] $file not found"
  exit
fi

rank=`expr $nobj + 1`
cd=`expr $nobj + 2`
./RANKandCD $nobj $file | awk -v rank=$rank '{if($rank==1) print $0}' | LC_ALL=C sort -g -rk$cd | uniq | head -$nSolutions | cut -f1-$nobj


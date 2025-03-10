#! /bin/sh

bm="cube-fractional-imbalance"
tries=1
ppn=24
benchmark="${bm}"
hemocell=~/HemoCell-mpi-2/

nodes=16
config="config.xml"


#for value in "${imbalances[@]}"
#do
#done

size="600,400,200"
iters=500

imbalance=0
benchmarkName=case-particled-balanced-0
uniqdir=$benchmarkName$((1 + $RANDOM % 10000))
setupdir=setups/$uniqdir
python3 $hemocell/Hemocell-Performance-Benchmarks/scripts/setup-experiment.py -i ./ -o setups/ -n $uniqdir -t $iters -s $size -p $(($ppn * $nodes)) --fli_fluid 0 --fli_part $imbalance --fli_part_base 65 --fli_part_stack 1
sbatch --array=1-$tries --exclusive -N $nodes --ntasks-per-node=${ppn} experiment.job ${config} ${benchmark} ${setupdir} $benchmarkName

imbalance=1.5
benchmarkName=case-particled-imbalanced-1
uniqdir=$benchmarkName$((1 + $RANDOM % 10000))
setupdir=setups/$uniqdir
python3 $hemocell/Hemocell-Performance-Benchmarks/scripts/setup-experiment.py -i ./ -o setups/ -n $uniqdir -t $iters -s $size -p $(($ppn * $nodes)) --fli_fluid 0 --fli_part $imbalance --fli_part_base 65 --fli_part_stack 1
sbatch --array=1-$tries --exclusive -N $nodes --ntasks-per-node=${ppn} experiment.job ${config} ${benchmark} ${setupdir} $benchmarkName

imbalance=1.5
benchmarkName=case-particled-2.8-imbalanced-2
uniqdir=$benchmarkName$((1 + $RANDOM % 10000))
setupdir=setups/$uniqdir
python3 $hemocell/Hemocell-Performance-Benchmarks/scripts/setup-experiment.py -i ./ -o setups/ -n $uniqdir -t $iters -s $size -p $(($ppn * $nodes)) --fli_fluid 0 --fli_part $imbalance --fli_part_base 65 --fli_part_stack 1
sbatch --array=1-$tries --exclusive -N $nodes --ntasks-per-node=${ppn} experiment.job ${config} ${benchmark} ${setupdir} $benchmarkName

imbalance=1.5
benchmarkName=case-particled-energy-3
uniqdir=$benchmarkName$((1 + $RANDOM % 10000))
setupdir=setups/$uniqdir
python3 $hemocell/Hemocell-Performance-Benchmarks/scripts/setup-experiment.py -i ./ -o setups/ -n $uniqdir -t $iters -s $size -p $(($ppn * $nodes)) --fli_fluid 0 --fli_part $imbalance --fli_part_base 65 --fli_part_stack 1
touch ${setupdir}/freq.txt
echo 1.5   >> ${setupdir}/freq.txt
echo 0 2.8 >> ${setupdir}/freq.txt
echo 1 2.8 >> ${setupdir}/freq.txt
sbatch --array=1-$tries --exclusive -N $nodes --ntasks-per-node=${ppn} experiment.job ${config} ${benchmark} ${setupdir} $benchmarkName

imbalance=1.5
benchmarkName=case-particled-energy-4
uniqdir=$benchmarkName$((1 + $RANDOM % 10000))
setupdir=setups/$uniqdir
python3 $hemocell/Hemocell-Performance-Benchmarks/scripts/setup-experiment.py -i ./ -o setups/ -n $uniqdir -t $iters -s $size  -p $(($ppn * $nodes)) --fli_fluid 0 --fli_part $imbalance --fli_part_base 65 --fli_part_stack 1
touch ${setupdir}/freq.txt
echo 2.4   >> ${setupdir}/freq.txt
echo 0 2.8 >> ${setupdir}/freq.txt
echo 1 2.8 >> ${setupdir}/freq.txt
sbatch --array=1-$tries --exclusive -N $nodes --ntasks-per-node=${ppn} experiment.job ${config} ${benchmark} ${setupdir} $benchmarkName

imbalance=1.5
benchmarkName=case-particled-energy-5
uniqdir=$benchmarkName$((1 + $RANDOM % 10000))
setupdir=setups/$uniqdir
python3 $hemocell/Hemocell-Performance-Benchmarks/scripts/setup-experiment.py -i ./ -o setups/ -n $uniqdir -t $iters -s $size  -p $(($ppn * $nodes)) --fli_fluid 0 --fli_part $imbalance --fli_part_base 65 --fli_part_stack 1
touch ${setupdir}/freq.txt
echo 1.5   >> ${setupdir}/freq.txt
echo 0 2.4 >> ${setupdir}/freq.txt
echo 1 2.4 >> ${setupdir}/freq.txt
sbatch --array=1-$tries --exclusive -N $nodes --ntasks-per-node=${ppn} experiment.job ${config} ${benchmark} ${setupdir} $benchmarkName

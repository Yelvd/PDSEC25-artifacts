#! /bin/sh

bm="cube-fractional-imbalance"
tries=1
ppn=24
benchmark="${bm}"
hemocell=~/HemoCell-mpi-2/

nodes=16
fli=1
config="config.xml"

iters=500

imbalance=0
benchmarkName=case-fluid-balanced-0
uniqdir=$benchmarkName$((1 + $RANDOM % 10000))
setupdir=setups/$uniqdir
python3 $hemocell/Hemocell-Performance-Benchmarks/scripts/setup-experiment.py -i ./ -o setups/ -n $uniqdir -t $iters -s 600,400,200 -p $(($ppn * $nodes)) --fli_fluid $imbalance --fli_part 0 --fli_part_base 38 --fli_part_stack 1 --ELI_case 3
sbatch --array=1-$tries --exclusive -N $nodes --ntasks-per-node=${ppn} experiment.job ${config} ${benchmark} ${setupdir} $benchmarkName

imbalance=$fli
benchmarkName=case-fluid-imbalanced-1
uniqdir=$benchmarkName$((1 + $RANDOM % 10000))
setupdir=setups/$uniqdir
python3 $hemocell/Hemocell-Performance-Benchmarks/scripts/setup-experiment.py -i ./ -o setups/ -n $uniqdir -t $iters -s 600,400,200 -p $(($ppn * $nodes)) --fli_fluid $imbalance --fli_part 0 --fli_part_base 38 --fli_part_stack 1 --ELI_case 3
sbatch --array=1-$tries --exclusive -N $nodes --ntasks-per-node=${ppn} experiment.job ${config} ${benchmark} ${setupdir} $benchmarkName

imbalance=$fli
benchmarkName=case-fluid-high-imbalanced-2
uniqdir=$benchmarkName$((1 + $RANDOM % 10000))
setupdir=setups/$uniqdir
python3 $hemocell/Hemocell-Performance-Benchmarks/scripts/setup-experiment.py -i ./ -o setups/ -n $uniqdir -t $iters -s 600,400,200 -p $(($ppn * $nodes)) --fli_fluid $imbalance --fli_part 0 --fli_part_base 38 --fli_part_stack 1 --ELI_case 3
touch ${setupdir}/freq.txt
echo 2.8   >> ${setupdir}/freq.txt
sbatch --array=1-$tries --exclusive -N $nodes --ntasks-per-node=${ppn} experiment.job ${config} ${benchmark} ${setupdir} $benchmarkName


imbalance=$fli
benchmarkName=case-fluid-energy-3
uniqdir=$benchmarkName$((1 + $RANDOM % 10000))
setupdir=setups/$uniqdir
python3 $hemocell/Hemocell-Performance-Benchmarks/scripts/setup-experiment.py -i ./ -o setups/ -n $uniqdir -t $iters -s 600,400,200 -p $(($ppn * $nodes)) --fli_fluid $imbalance --fli_part 0 --fli_part_base 38 --fli_part_stack 1 --ELI_case 3
touch ${setupdir}/freq.txt
echo 1.5   >> ${setupdir}/freq.txt
echo 0 2.8 >> ${setupdir}/freq.txt
echo 1 2.8 >> ${setupdir}/freq.txt
sbatch --array=1-$tries --exclusive -N $nodes --ntasks-per-node=${ppn} experiment.job ${config} ${benchmark} ${setupdir} $benchmarkName

imbalance=$fli
benchmarkName=case-fluid-energy-4
uniqdir=$benchmarkName$((1 + $RANDOM % 10000))
setupdir=setups/$uniqdir
python3 $hemocell/Hemocell-Performance-Benchmarks/scripts/setup-experiment.py -i ./ -o setups/ -n $uniqdir -t $iters -s 600,400,200 -p $(($ppn * $nodes)) --fli_fluid $imbalance --fli_part 0 --fli_part_base 38 --fli_part_stack 1 --ELI_case 3
touch ${setupdir}/freq.txt
echo 2.4   >> ${setupdir}/freq.txt
echo 0 2.8 >> ${setupdir}/freq.txt
echo 1 2.8 >> ${setupdir}/freq.txt
sbatch --array=1-$tries --exclusive -N $nodes --ntasks-per-node=${ppn} experiment.job ${config} ${benchmark} ${setupdir} $benchmarkName

imbalance=$fli
benchmarkName=case-fluid-energy-5
uniqdir=$benchmarkName$((1 + $RANDOM % 10000))
setupdir=setups/$uniqdir
python3 $hemocell/Hemocell-Performance-Benchmarks/scripts/setup-experiment.py -i ./ -o setups/ -n $uniqdir -t $iters -s 600,400,200 -p $(($ppn * $nodes)) --fli_fluid $imbalance --fli_part 0 --fli_part_base 38 --fli_part_stack 1 --ELI_case 3
touch ${setupdir}/freq.txt
echo 1.5   >> ${setupdir}/freq.txt
echo 0 2.4 >> ${setupdir}/freq.txt
echo 1 2.4 >> ${setupdir}/freq.txt
sbatch --array=1-$tries --exclusive -N $nodes --ntasks-per-node=${ppn} experiment.job ${config} ${benchmark} ${setupdir} $benchmarkName

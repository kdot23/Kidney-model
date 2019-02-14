declare -a pops
for i in $(seq 1 500)
do
	pops[$(($i-1))]=pop50_100/data$i.dat
done
echo "${pops[@]}"
python greedyMatching2.py --inputFiles ${pops[@]} -o graph_stuff/OAES_count.csv --agents graph_stuff/OAES_count_A.csv


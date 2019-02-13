
declare -a data
for i in $(seq 1 500)
do
	data[$(($i-1))]=pop50_100/data$i.dat
done
python demoAndAgent.py --data ${data[@]} --agent graph_stuff/OAES_count_A.csv -T 50 -K 100 -o graph_stuff/OAES_count_DA.csv



'''
declare -a train
for i in $(seq 1 100)
do
	train[$(($i-1))]=pop50_100/data$i.dat
done

python generateBeta.py --inputFiles ${train[@]} -o train50_100_2_7.csv --quality --graph_state 
'''
declare -a pops
for i in $(seq 102 150)
do
	pops[$(($i-1))]=pop50_100/data$i.dat
done

#python onlineMatching.py --trainFiles train.csv --testFiles ${pops[@]} --quality -o odasse.csv --agents odasse_agents.csv --graph_state --lpRepeat

python onlineMatching.py --trainFiles train50_100.csv --testFiles ${pops[@]} --quality -o results2_7/odasse3_50_100_10_2.csv --agents results2_7/odasse3_50_100_10_2agents.csv --fwd_proj 10 --fwd_proj_repeat

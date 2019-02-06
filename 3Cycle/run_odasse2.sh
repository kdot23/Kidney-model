'''
declare -a train
for i in $(seq 1 100)
do
	train[$(($i-1))]=pop50_100/data$i.dat
done

python generateBeta.py --inputFiles ${train[@]} -o train50_100.csv --quality --graph_state 
'''

declare -a pops
for i in $(seq 101 150)
do
	pops[$(($i-1))]=pop50_100/data$i.dat
done
#python onlineMatching.py --trainFiles train.csv --testFiles ${pops[@]} --quality -o odasse.csv --agents odasse_agents.csv --graph_state --lpRepeat

python onlineMatching.py --trainFiles train50_100.csv --testFiles ${pops[@]} --quality -o odasse3_50_100.csv --agents odasse3_50_100_agents.csv --fwd_proj 50 --fwd_proj_repeat


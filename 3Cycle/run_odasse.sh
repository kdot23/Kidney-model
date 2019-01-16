: '
declare -a train
for i in $(seq 1 100)
do
	train[$(($i-1))]=pop50_100/data$i.dat
done

python generateBeta.py --inputFiles $train -o train.csv --quality --graph_state 
'

declare -a pops
for i in $(seq 101 500)
do
	pops[$(($i-1))]=pop50_100/data$i.dat
done
python onlineMatching.py --trainFiles train.csv --testFiles ${pops[@]} --quality -o odasse.csv --agents odasse_agents.csv --graph_state --lpRepeat

python onlineMatching.py --trainFiles train.csv --testFiles ${pops[@]} --quality -o odasse.csv --agents odasse_agents.csv --graph_state --lpRepeat


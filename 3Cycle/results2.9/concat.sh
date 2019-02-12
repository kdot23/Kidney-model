: '
declare -a train
for i in $(seq 1 100)
do
	train[$(($i-1))]=pop50_100/data$i.dat
done

python generateBeta.py --inputFiles ${train[@]} -o train.csv --quality --graph_state 
'
for i in $(seq 101 50 451)
do
	j=$(($i+49))
	echo $j
	cat odasseQ_fwdProj_$i-$j*agents.csv >> odasse_fwd_A.csv 
done

#python onlineMatching.py --trainFiles train.csv --testFiles ${pops[@]} --quality -o odasse.csv --agents odasse_agents.csv --graph_state --lpRepeat



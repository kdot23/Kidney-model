for i in {1..500}
do
	echo $i
	python KidneyDataGen.py -K 50 -T 50 -o pop50_50/data$i.dat
	python KidneyDataGen.py -K 100 -T 100 -o pop100_100/data$i.dat
done


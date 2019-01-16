for i in {1..500}
do
	echo $i
	python KidneyDataGen.py -K 100 -T 50 -o pop50_100/data$i.dat
done


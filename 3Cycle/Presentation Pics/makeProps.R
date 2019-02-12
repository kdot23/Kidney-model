greedyA = read.csv('fancy_greedy_DA.csv', header=FALSE, sep='\t')
odasseA = read.csv('fwd_rep_odasse_DA.csv', header=FALSE, sep='\t')

greedyO = subset(greedyA, grepl('I', V1) & V7 == 1)
totalGreedyO = nrow(greedyO)
propG_O_1 = nrow(subset(greedyO, V2 <=50))/totalGreedyO
propG_O_2 = nrow(subset(greedyO, V2 ==51))/totalGreedyO
propG_O_3 = nrow(subset(greedyO, V2 >51))/totalGreedyO

greedyNO = subset(greedyA, grepl('I', V1) & V7 == 1)
totalGreedyNO = nrow(greedyNO)
propG_NO_1 = nrow(subset(greedyNO, V2 <=50))/totalGreedyNO
propG_NO_2 = nrow(subset(greedyNO, V2 ==51))/totalGreedyNO
propG_NO_3 = nrow(subset(greedyNO, V2 >51))/totalGreedyNO

odasseO = subset(odasseA, grepl('I', V1) & V8 == 1)
totalOdasseO = nrow(odasseO)
propO_O_1 = nrow(subset(odasseO, V2 <=50))/totalOdasseO
propO_O_2 = nrow(subset(odasseO, V2 ==51))/totalOdasseO
propO_O_3 = nrow(subset(odasseO, V2 >51))/totalOdasseO

odasseNO = subset(odasseA, grepl('I', V1) & V8 == 1)
totalOdasseNO = nrow(odasseNO)
propO_NO_1 = nrow(subset(odasseNO, V2 <=50))/totalOdasseNO
propO_NO_2 = nrow(subset(odasseNO, V2 ==51))/totalOdasseNO
propO_NO_3 = nrow(subset(odasseNO, V2 >51))/totalOdasseNO

props = cbind(c(propG_O_1, propG_O_2, propG_O_3),  
              c(propO_O_1, propO_O_2, propO_O_3), 
              c(propG_NO_1, propG_NO_2, propG_NO_3),
              c(propO_NO_1, propO_NO_2, propO_NO_3))


par(mar=c(10,5,2,2), xpd=TRUE)
barplot(props, col=c('ivory4', 'ivory3', 'ivory1'), names=c('Fancy Greedy O', 'fwd rep ODASSE O', 'Fancy Greedy NO', 'fwd rep ODASSE NO'), ylab='Proportion of Matches', main='Proportions of Matches by Algorithm and Recipient Blood Type')
legend('bottom',legend=c('Online', 'End', 'Unmatched'), fill=c('ivory4', 'ivory3', 'ivory2'), inset=c(-.4,-.8), ncol=3)

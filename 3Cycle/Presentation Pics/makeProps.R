greedyA = read.csv('greedyDA.csv', header=FALSE, sep='\t')
odasseA = read.csv('onlineLinearGraphAgents.csv', header=FALSE, sep='\t')

greedyA = subset(greedyA, grepl('I', V1))
totalGreedy = nrow(greedyA)
propG1 = nrow(subset(greedyA, V2 <=50))/totalGreedy
propG2 = nrow(subset(greedyA, V2 ==51))/totalGreedy
propG3 = nrow(subset(greedyA, V2 >51))/totalGreedy


odasseA = subset(odasseA, grepl('I', V1))
totalOdasse = nrow(odasseA)
propO1 = nrow(subset(odasseA, V2 <=50))/totalOdasse
propO2 = nrow(subset(odasseA, V2 ==51))/totalOdasse
propO3 = nrow(subset(odasseA, V2 >51))/totalOdasse

props = cbind(c(propG1, propG2, propG3), c(propO1, propO2, propO3))


par(mar=c(10,5,2,2), xpd=TRUE)
barplot(props, col=c('lightblue2', 'palegreen', 'palegoldenrod'), names=c('Greedy', 'ODASSE'), ylab='Proportion of Matches', main='Proportions of Matches by Algorithm')
legend('bottom',legend=c('Received Kidney at End','Received Kidney During Online Part','Did not Receive a Kidney'), fill=c('lightblue2', 'palegreen', 'palegoldenrod'), inset=c(-.4,-.4))

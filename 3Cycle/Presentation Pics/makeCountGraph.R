incompat_only = read.csv('incompat_only.csv', header=FALSE, sep='\t')
fancy_greedy = read.csv('fancy_greedy.csv', header=FALSE, sep='\t')
odasse = read.csv('odasse.csv', header=FALSE, sep='\t')
fwd_rep_odasse = read.csv('fwd_rep_odasse.csv', header=FALSE, sep='\t')
oracle = read.csv('oracle.csv', header=FALSE, sep='\t')

incompat_only$Cat = c(rep('incompat only', nrow(incompat_only)))
fancy_greedy$Cat = c(rep('fancy greedy', nrow(fancy_greedy)))
odasse$Cat = c(rep('odasse', nrow(odasse)))
fwd_rep_odasse$Cat = c(rep('fwd rep odasse', nrow(fwd_rep_odasse)))
oracle$Cat = c(rep('oracle', nrow(oracle)))

total = rbind(incompat_only, fancy_greedy, odasse, fwd_rep_odasse, oracle)
total$Cat = factor(total$Cat, c('incompat only', 'fancy greedy', 'odasse', 'fwd rep odasse', 'oracle'))
par(cex.main=1.5, cex.lab=1, cex.axis=.7)
a = boxplot(V2~Cat, total, col='powderblue', names=c('Baseline (Compatible and \n Incompatible Separate)', 
                                                 'OAES', expression(paste('ODASSE \n (ML estimates of ', beta, ')')),
                                                 expression(paste('ODASSE \n(Simulated estimates of ', beta, ')')), 
                                                 expression(paste('Oracle \n(ODASSE with perfect ', beta, ')'))), 
        main='Total Expected Graft Survival by Algorithm', ylab='Total Expected Graft Survival', boxlwd=2, las=2, srt=45)

ggplot(total, aes(x=Cat, y=V1)) + geom_boxplot(fill="gray")+
  labs(title="Total Expected Graft Survival by Algorithm",x="", y = "Total Expected Graft Survival")+
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) 

par(cex.main=1.25, cex.lab=1.5, cex.axis=.6)
boxplot(V2~Cat, total, col='powderblue', names=c('incompat only', 'fancy greedy', 'odasse', 'fwd odasse', 'fwd rep odasse', 'oracle'), main='Total Expected Graft Survival by Algorithm', ylab='Total Expected Graft Survival', boxlwd=2)


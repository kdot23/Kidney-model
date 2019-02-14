incompat_only = read.csv('incompat_only_DA.csv', header=FALSE, sep='\t')
greedy_count = read.csv('OAES_count_DA.csv', header=FALSE, sep='\t') #CHANGE TO GREEDY COUNT AFTER SIMULATION

#table w/ prop of incompat matched, prop of type O matched and prop high PRA matched by greedy count and baseline
incompat_baseline = subset(incompat_only, grepl('I', V1))
type_O_baseline = subset(incompat_baseline, V7 == 1)
high_pra_baseline = subset(incompat_baseline, V25 == .9)
baseline_stats = c(nrow(subset(incompat_baseline, V2 <= 51)) / nrow(incompat_baseline),
                   nrow(subset(type_O_baseline, V2 <= 51)) / nrow(type_O_baseline),
                   nrow(subset(high_pra_baseline, V2 <= 51)) / nrow(high_pra_baseline)) 
  
incompat_greedy = subset(greedy_count, grepl('I', V1))
type_O_greedy = subset(greedy_count, V7 == 1)
high_pra_greedy = subset(greedy_count, V25 == .9)
greedy_stats = c(nrow(subset(incompat_greedy, V2 <= 51)) / nrow(incompat_greedy),
                   nrow(subset(type_O_greedy, V2 <= 51)) / nrow(type_O_greedy),
                   nrow(subset(high_pra_greedy, V2 <= 51)) / nrow(high_pra_greedy)) 

a = rbind(baseline_stats, greedy_stats)
par(mar=c(8,5,2,2), xpd=TRUE, cex.lab=1, cex.main=1.25)
b= barplot(a, beside=TRUE, names.arg = c("Incompatible Pairs", 
                                      "Pairs with recipients\n with type O blood", "Pairs with high PRA"),
        main = "Proportion of Pairs Matched by Algorithm", col=c('seagreen4', 'seagreen2'), ylim=c(0,1),
        ylab="Proportion Matched")

legend("bottom",
       c("Baseline (Compatible and Incompatible Separate)", "OAES optimized for number of matches"),
       fill = c('seagreen4', 'seagreen2'), inset=c(-.4,-.5), cex=.8)
text(b, 0, round(a, 3),cex=.8,pos=3) 

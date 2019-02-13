baseline = read.csv('incompat_only_A.csv', header=FALSE, sep='\t')
oaes = read.csv('fancy_greedy_A.csv', header=FALSE, sep='\t')
odasse = read.csv('odasse_A.csv', header=FALSE, sep='\t')
oracle = read.csv('oracle_A.csv', header=FALSE, sep='\t')

baseline_incompat = subset(baseline, grepl('I', V1))
baseline_stats = c(nrow(subset(baseline_incompat, V2 <= 51)) / nrow(baseline_incompat),
                   mean(subset(baseline, grepl('C', V1))$V3),
                   mean(subset(baseline_incompat, V2 <= 51)$V3))
oaes_incompat = subset(oaes, grepl('I', V1))
oaes_stats = c(nrow(subset(oaes_incompat, V2 <= 51)) / nrow(oaes_incompat),
                   mean(subset(oaes, grepl('C', V1))$V3),
               mean(subset(oaes_incompat, V2 <= 51)$V3))
odasse_incompat = subset(odasse, grepl('I', V1))
odasse_stats = c(nrow(subset(odasse_incompat, V2 <= 51)) / nrow(odasse_incompat),
                   mean(subset(odasse, grepl('C', V1))$V3), 
                 mean(subset(odasse_incompat, V2 <= 51)$V3))
oracle_incompat = subset(oracle, grepl('I', V1))
oracle_stats = c(nrow(subset(oracle_incompat, V2 <= 51)) / nrow(oracle_incompat),
                   mean(subset(oracle, grepl('C', V1))$V3),
                 mean(subset(oracle_incompat, V2 <= 51)$V3))
a = cbind(baseline_stats, oaes_stats, odasse_stats, oracle_stats)

#for compatible pairs only, what is the expected improvement in graft survival 
#_conditional_ on matching with an incompatible pair?
#average EGS received by compatible pairs when matching with themselves in baseline
self_match = mean(subset(baseline, grepl('C', V1))$V3)

#average EGS received by compatible pairs when matching with incompatible in oaes
exchange = mean(subset(oaes, grepl('C', V1) & V4 == "I")$V3)

exchange - self_match

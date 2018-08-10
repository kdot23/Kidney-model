library(readr)
#Data from pool of 50 compatible, 100 incompatible. Optimized with up to 3 cycles over 500 populations
greedy <- read_delim("~/Documents/Kidney-model/3Cycle/results50_100/greedyDA2.csv", 
                             "\t", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE)
odasse <- read_delim("~/Documents/Kidney-model/3Cycle/results50_100/odasseDA.csv", 
                     "\t", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)

odasse$algorithm = c(rep('odasse'))
greedy$algorithm = c(rep('greedy'))

agents = rbind(odasse, greedy)
colnames(agents) = c("id","time","get_quality","get_type","give_quality", "give_type", "beta", 
                     "blood_patient_o", "blood_patient_a", "blood_patient_b", "blood_patient_ab", 
                     "blood_donor_o", "blood_donor_a", "blood_donor_b", "blood_donor_ab", "donor_afam", 
                     "donor_age", "donor_sex", "donor_cig", "rec_sex", "donor_weight", "rec_weight", 
                     "donor_bmi", "donor_eGFR", "donor_sbp","rec_pra", "algorithm")

agents$betaCat = cut(agents$beta, c(-100,5,10,15,20,25), labels = c('<5','5-10','10-15','15-20','>20'))
agents$time_type = agents$time_type = ifelse(agents$time <= 50, 'during', ifelse(agents$time == 51, 'after', 'never'))
agents$timeCat = cut(agents$time, c(0,10,20,30,40,50,51), labels = c('0-10','10-20','20-30','30-40','40-50', 'End'))

#proportion of time_type by beta & algorithm (poster & prez)
t = prop.table(cbind(table(subset(agents, betaCat == "<5" & algorithm == "greedy" & grepl("I",id))$time_type), 
                     table(subset(agents, betaCat == "<5" & algorithm == "odasse" & grepl("I",id))$time_type), 
                     table(subset(agents, betaCat == "5-10" & algorithm == "greedy" & grepl("I",id))$time_type), 
                     table(subset(agents, betaCat == "5-10" & algorithm == "odasse" & grepl("I",id))$time_type), 
                     table(subset(agents, betaCat == "10-15" & algorithm == "greedy" & grepl("I",id))$time_type), 
                     table(subset(agents, betaCat == "10-15" & algorithm == "odasse" & grepl("I",id))$time_type), 
                     table(subset(agents, betaCat == "15-20" & algorithm == "greedy" & grepl("I",id))$time_type), 
                     table(subset(agents, betaCat == "15-20" & algorithm == "odasse" & grepl("I",id))$time_type), 
                     c(0,table(subset(agents, betaCat == ">20" & algorithm == "greedy" & grepl("I",id))$time_type)), 
                     table(subset(agents, betaCat == ">20" & algorithm == "odasse" & grepl("I",id))$time_type)),2)
par(mar = c(11.5,4,2,4),xpd=T)
barplot(t, col = c("lightblue2", "palegreen", "palegoldenrod"), #main = "Time Incompatible Pairs\n Match by Beta and Algorithm", 
        names = c("Greedy β < 5", "ODASSE β<5", 
                 "Greedy 5<β<10", "ODASSE 5<β<10", 
                  "Greedy 10<β<15","ODASSE 10<β<15", 
                  "Greedy 15<β<20", "ODASSE 15<β<20", 
                  "Greedy β>20", "ODASSE β>20"), las=2)
legend("bottom", inset = c(-1.6,-1.6), legend = 
         c("Received kidney at the end", "Received kidney during the online part", "Did not receive a kidney"), 
       fill = c("lightblue2", "palegreen", "palegoldenrod"), cex = .75)

#time_type by afam and algorithm (poster & prez)
t = prop.table(cbind(table(subset(agents, grepl("I", id) & donor_afam == 1 & algorithm == 'greedy')$time_type),
                     table(subset(agents, grepl("I", id) & donor_afam == 1 & algorithm == 'odasse')$time_type),
                     table(subset(agents, grepl("I", id) & donor_afam == 0 & algorithm == 'greedy')$time_type),
                     table(subset(agents, grepl("I", id) & donor_afam == 0 & algorithm == 'odasse')$time_type)),2)
par(mar = c(11.5,4,2,4), xpd=T)
barplot(t, names = c("Greedy- African\nAmerican Donor", "ODASSE- African\nAmerican Donor", "Greedy- non-African\nAmerican Donor", "ODASSE-non-African\nAmerican Donor"),
        col = c("lightblue2", "palegreen", "palegoldenrod"), 
        #main = "Time Incompatible Pairs with\nAfrican American Donors Match", 
        las=2)
legend("bottom", inset = c(-1.55,-1.55), legend = 
         c("Received kidney at the end", "Received kidney during the online part", "Did not receive a kidney"), 
       fill = c("lightblue2", "palegreen", "palegoldenrod"), cex = .65)

#time type by blood type & algorithm (poster)
par(mar = c(9,4,2,4), xpd=T)
t = prop.table(cbind(table(subset(agents, grepl("I", id) & blood_patient_o == 1 & algorithm == 'greedy')$time_type),
                     table(subset(agents, grepl("I", id) & blood_patient_o == 1 & algorithm == 'odasse')$time_type),
                     table(subset(agents, grepl("I", id) & blood_patient_a == 1 & algorithm == 'greedy')$time_type),
                     table(subset(agents, grepl("I", id) & blood_patient_a == 1 & algorithm == 'odasse')$time_type),
                     table(subset(agents, grepl("I", id) & blood_patient_b == 1 & algorithm == 'greedy')$time_type),
                     table(subset(agents, grepl("I", id) & blood_patient_b == 1 & algorithm == 'odasse')$time_type),
                     table(subset(agents, grepl("I", id) & blood_patient_ab == 1 & algorithm == 'greedy')$time_type),
                     table(subset(agents, grepl("I", id) & blood_patient_ab == 1 & algorithm == 'odasse')$time_type)),2)
barplot(t, col = c("lightblue2", "palegreen", "palegoldenrod"), #main = "Time Incompatible Pairs Match\nby Patient Blood Type", 
        names = c("Greedy-O", "ODASSE-O", "Greedy-A", "ODASSE-A", "Greedy-B", "ODASSE-B", 
                  "Greedy-AB", "ODASSE-AB"), las=2)
legend("bottom", inset = c(-.9,-.9), legend = 
         c("Received kidney at the end", "Received kidney during the online part", "Did not receive a kidney"), 
       fill = c("lightblue2", "palegreen", "palegoldenrod"), cex = .75)

#time type by blood type & pra & algorithm (prez)
par(mar = c(12,4,2,4), xpd=T)
t = prop.table(cbind(table(subset(agents, grepl("I", id) & blood_patient_o == 1 & algorithm == 'greedy')$time_type),
                     table(subset(agents, grepl("I", id) & blood_patient_o == 1 & algorithm == 'odasse')$time_type),
                     table(subset(agents, grepl("I", id) & blood_patient_o == 0 & algorithm == 'greedy')$time_type),
                     table(subset(agents, grepl("I", id) & blood_patient_o == 0 & algorithm == 'odasse')$time_type),
                     table(subset(agents, grepl("I", id) & rec_pra == .90 & algorithm == 'greedy')$time_type),
                     table(subset(agents, grepl("I", id) & rec_pra == .90 & algorithm == 'odasse')$time_type),
                     table(subset(agents, grepl("I", id) & rec_pra != .90 & algorithm == 'greedy')$time_type),
                     table(subset(agents, grepl("I", id) & rec_pra != .90 & algorithm == 'odasse')$time_type)),2)
barplot(t, col = c("lightblue2", "palegreen", "palegoldenrod"), #main = "Time Incompatible Pairs Match\nby Patient Blood Type", 
        names = c("Greedy-O", "ODASSE-O", "Greedy-Not O", "ODASSE-Not O", "Greedy-High PRA", "ODASSE-High PRA", 
                  "Greedy-Med or Low PRA", "ODASSE-Med or Low PRA"), las=2)
legend("bottomleft", inset = c(0,-1.7), legend = 
         c("Received kidney at the end", "Received kidney during the online part", "Did not receive a kidney"), 
       fill = c("lightblue2", "palegreen", "palegoldenrod"), cex = .75)

#get quality by time matched (prez)
par(mar=c(7,4,2,3), xpd = T)
boxplot(get_quality~timeCat*algorithm, data = subset(agents, get_type != "N" & (algorithm == "greedy" | algorithm == "odasse")),
        las=2, col = c(rep("chocolate", 6), rep("cadetblue", 6)), ylab = "Quality Received (in EGS)", #main = "Quality Received by Time Matched of Incompatible Pairs",
        names = c("0<t<10", "10<t<20","20<t<30", "30<t<40", "40<t<50", "End", "0<t<10", "10<t<20","20<t<30", "30<t<40", "40<t<50", "End"))
legend("bottom", inset = c(-.55,-.55), legend = c("Greedy", "ODASSE"), fill = c("chocolate", "cadetblue"), 
       cex = .75, horiz = TRUE)


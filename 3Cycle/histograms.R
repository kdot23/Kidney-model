online$betaCat = cut(online$beta, c(-100,5,10,15,20,25), labels = c('<5','5-10','10-15','15-20','>20'))
greedy$betaCat = cut(greedy$beta, c(-100,5,10,15,20,25), labels = c('<5','5-10','10-15','15-20','>20'))

#histogram of all agents
greedyh = hist(agentsOnline50$get_quality, col= rgb(0,0,1,1/4), , main = "Histogram of Quality for All Pairs", 
               xlab="Quality (Expected Graft Survival)", breaks = c(0:25), ylim = c(0,8000))
onlineh = hist(agentsGreedy50$get_quality, col=rgb(1,0,0,1/4), breaks = c(0:25), add=T)
legend("topright", c("Online", "Greedy"), col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), fill=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))

#barplot of Received type by algorithm
#table gets counts of each category
r = cbind(table(onlinec$get_type),table(greedyc$get_type),table(onlinei$get_type),table(greedyi$get_type))
r[3,1:2] = 0
#mar gives margins of c(bottom,left,top,right), xpd means things can be written in margin
par(mar=c(6.1, 3.1, 2.1, 3.1), xpd=TRUE)
barplot(r, ylim = c(0, 25000), col = c("lightblue2", "palegreen", "palegoldenrod"), 
        names = c("C-Online", "C-Greedy", "I-Online", "I-Greedy"), main = "Matches By Algorithm")
#inset is how far from the plot
legend("bottom", inset = c(-.45,-.45), legend = 
         c("Received kidney from a compatible pair", "Received kidney from an incompatible pair", "No match"), 
       fill = c("lightblue2", "palegreen", "palegoldenrod"), cex = .75)

#barplot of what happens by beta (may use proportions)
onlinei = subset(agentsOnline50, grepl("I", id))
t = cbind(table(onlinei[which(onlinei$betaCat == '<5'), ]$get_type),
      table(onlinei[which(onlinei$betaCat == '5-10'), ]$get_type),
      table(onlinei[which(onlinei$betaCat == '10-15'), ]$get_type),
      table(onlinei[which(onlinei$betaCat == '15-20'), ]$get_type),
      table(onlinei[which(onlinei$betaCat == '>20'), ]$get_type))
#proportions
tp = cbind(t[,1]/nrow(onlinei[which(onlinei$betaCat == '<5'),]),
                t[,2]/nrow(onlinei[which(onlinei$betaCat == '5-10'),]),
                t[,3]/nrow(onlinei[which(onlinei$betaCat == '10-15'),]),
                t[,4]/nrow(onlinei[which(onlinei$betaCat == '15-20'),]),
                t[,5]/nrow(onlinei[which(onlinei$betaCat == '>20'),]))
par(mar=c(6.1, 3.1, 3.1, 3.1), xpd=TRUE)
barplot(tp, ylim = c(0, 1), col = c("lightblue2", "palegreen", "palegoldenrod"), 
        names = c('<5','5-10','10-15','15-20','>20'), main = "Proportion of Match Types By Beta")
legend("bottom", inset = c(-.45,-.45), legend = 
         c("Received kidney from a compatible pair", "Received kidney from an incompatible pair", "No match"), 
       fill = c("lightblue2", "palegreen", "palegoldenrod"), cex = .75)

#boxplot of quality by beta if matched with compatible
par(mar=c(5,4,4,4))
fmatched = subset(online, grepl("I",id) & online$get_type == 'C')
gmatched = subset(greedy, grepl("I",id) & greedy$get_type == 'C')
match = rbind(fmathced,gmatched)
boxplot(fmatched[which(fmatched$betaCat == '<5'), ]$get_quality, gmatched[which(gmatched$betaCat == '<5'), ]$get_quality,
        fmatched[which(fmatched$betaCat == '5-10'), ]$get_quality, gmatched[which(gmatched$betaCat == '5-10'), ]$get_quality,
        fmatched[which(fmatched$betaCat == '10-15'), ]$get_quality, gmatched[which(gmatched$betaCat == '10-15'), ]$get_quality,
        fmatched[which(fmatched$betaCat == '15-20'), ]$get_quality, gmatched[which(gmatched$betaCat == '15-20'), ]$get_quality,
        fmatched[which(fmatched$betaCat == '>20'), ]$get_quality, gmatched[which(gmatched$betaCat == '>20'), ]$get_quality,
        main = "Quality Received When Matched With Compatible", ylim = c(0,25), lab = "Beta Value", ylab = "Quality",
        names= c("<5: Online","5-10: Online","10-15: Online","15-20: Online", ">20: Online",
                 "<5: Greedy","5-10: Greedy","10-15: Greedy","15-20: Greedy", ">20: Greedy"))

#boxplot of avg time by beta if matched with compatible
boxplot(fmatched[which(fmatched$betaCat == '<5'), ]$time, 
        fmatched[which(fmatched$betaCat == '5-10'), ]$time,
        fmatched[which(fmatched$betaCat == '10-15'), ]$time, 
        fmatched[which(fmatched$betaCat == '15-20'), ]$time,
        fmatched[which(fmatched$betaCat == '>20'), ]$time, 
        main = "Time Matched With Compatible", xlab = "Beta Value", ylab = "Time",
        names= c("<5","5-10","10-15","15-20", ">20"))
library(readr)
greedy <- read_delim("~/Documents/Kidney-model/3Cycle/results/greedyAgents.csv", 
                     "\t", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
online <- read_delim("~/Documents/Kidney-model/3Cycle/results/demoAndAgentsOnline.csv", 
                     "\t", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
colnames(greedy)= c("id","time","get_quality","get_type","give_quality", "give_type", "beta", 
                    "blood_patient_o", "blood_patient_a", "blood_patient_b", "blood_patient_ab", 
                    "blood_donor_o", "blood_donor_a", "blood_donor_b", "blood_donor_ab", "donor_afam", 
                    "donor_age", "donor_sex", "donor_cig", "rec_sex", "donor_weight", "rec_weight", 
                    "donor_bmi", "donor_eGFR", "donor_sbp","rec_pra")





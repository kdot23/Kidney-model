incompat = subset(agents, grepl("I",id))
par(mar=c(7,4,4,4),xpd=T)
b =boxplot(get_quality~algorithm*betaCat, data=incompat,
        ylab="Quality Received", las=2) #,  col=(c("gold","darkgreen","blue")), main = "Quality Received by Beta Value of Incompatible Pairs")
compat = subset(agents, grepl("C",id))
boxplot(compat$get_quality, add=TRUE)
#get_type of afam
afam = subset(agents, agents$donor_afam == 1)
a = subset(agents, agents$algorithm=='Greedy' | agents$algorithm=='Online')
boxplot(get_quality~algorithm*donor_afam, data = agents,  ylab="Quality Received", #names=c("Greedy-Other donor", "Online-other donor", "Greedy-Afam donor", "Online-Afam donor"),
        las=2,  col=c("gold","darkgreen", "blue"), main = "Received Quality for African American Donors by Algorithm")
bartable = table(afam$get_type, afam$algorithm)
bartable2 = table(agents[which(agents$donor_afam == 0),]$get_type, agents[which(agents$donor_afam == 0),]$algorithm)
barplot(prop.table(cbind(table(agents[which(agents$donor_afam == 1),]$get_type, agents[which(agents$donor_afam == 1),]$algorithm), 
                         table(agents[which(agents$donor_afam == 0),]$get_type, agents[which(agents$donor_afam == 0),]$algorithm)),2),
        names = c("Greedy-Afam", "OnlineLP-Afam", "Oracle-Afam", "Greedy-Other","OnlineLP-Other", "Oracle-Other"),
        col = c("lightblue2", "palegreen", "palegoldenrod"), main="Pairs with African American donors", las=2)
legend("bottom", inset = c(-.65,-.65), legend = 
         c("Received kidney from a compatible pair", "Received kidney from an incompatible pair", "No match"), 
       fill = c("lightblue2", "palegreen", "palegoldenrod"), cex = .75)
barplot(get_type~algorithm*donor_afam)

barplot(prop.table(cbind(table(incompat[which(incompat$algorithm == "Online3_50"),]$get_type), 
                         table(incompat[which(incompat$algorithm == "Online3_50_100"),]$get_type)),2),
        names = c("50-C, 50-I", "50-C, 100-I"), col = c("lightblue2", "palegreen", "palegoldenrod"), 
        main = "Proportion of Type Received by Incompatible Pairs")

par(mar=c(7,4,3,4))
boxplot(get_quality~algorithm*donor_blood, data = agents,  ylab="Quality Received", 
        #names=c("Greedy-O", "Online-O", "Greedy-A", "Online-A", "Greedy-B", "Online-B", "Greedy-AB", "Online-AB"),
        las=2,  col=c("gold","darkgreen","blue","orange"), main = "Received Quality by Donor Blood Type and Algorithm")
boxplot(get_quality~algorithm*rec_blood, data = agents,  ylab="Quality Received", 
        #names=c("Greedy-O", "Online-O", "Greedy-A", "Online-A", "Greedy-B", "Online-B", "Greedy-AB", "Online-AB"),
        las=2,  col=c("gold","darkgreen","blue","orange"), main = "Received Quality by Recipient Blood Type and Algorithm")
boxplot(get_quality~algorithm*ageCat, data = agents,  ylab="Quality Received", 
        las=2,  col=c("gold","darkgreen","blue","orange"), main = "Received Quality by Donor Age and Algorithm")
boxplot(get_quality~algorithm*donor_cig, data = agents,  ylab="Quality Received", 
       # names=c("Greedy-No", "Online-No", "Greedy-Yes", "Online-Yes"),
        las=2,  col=c("gold","darkgreen","blue","orange"), main = "Received Quality by Donor Cigarette Use and Algorithm")
boxplot(get_quality~algorithm*rec_pra, data = agents,  ylab="Quality Received", 
        las=2,  col=c("gold","darkgreen","blue", "orange"), main = "Received Quality by Recipient PRA and Algorithm")

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
onlinei = subset(online, grepl("I", id))
greedyi = subset(greedy, grepl("I", id))
t = cbind(table(onlinei[which(onlinei$betaCat == '<5'), ]$get_type),
      table(onlinei[which(onlinei$betaCat == '5-10'), ]$get_type),
      table(onlinei[which(onlinei$betaCat == '10-15'), ]$get_type),
      table(onlinei[which(onlinei$betaCat == '15-20'), ]$get_type),
      table(onlinei[which(onlinei$betaCat == '>20'), ]$get_type),
      table(greedyi[which(greedyi$betaCat == '<5'), ]$get_type),
      table(greedyi[which(greedyi$betaCat == '5-10'), ]$get_type),
      table(greedyi[which(greedyi$betaCat == '10-15'), ]$get_type),
      table(greedyi[which(greedyi$betaCat == '15-20'), ]$get_type),
      table(greedyi[which(greedyi$betaCat == '>20'), ]$get_type))
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
par(mar=c(7,4,2,4))
fmatched = subset(online, grepl("I",id) & online$get_type == 'C')
fmatched$type = c(rep("Online"))
gmatched = subset(greedy, grepl("I",id) & greedy$get_type == 'C')
gmatched$type = c(rep("Greedy"))
match = rbind(fmatched,gmatched)
boxplot(get_quality~type*betaCat, data=match, notch=TRUE, 
        col=(c("gold","darkgreen")), las =2, ylab = "Quality Received")
boxplot(time~type*betaCat, data=match, notch=TRUE, 
        col=(c("gold","darkgreen")), las =2, ylab = "Time of Match")




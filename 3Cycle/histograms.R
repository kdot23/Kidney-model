agentsOnline50$betaCat = cut(agentsOnline50$beta, c(-100,5,10,15,20,25), labels = c('<5','5-10','10-15','15-20','>20'))

#histogram of all agents
greedyh = hist(agentsOnline50$get_quality, col= rgb(0,0,1,1/4), , main = "Histogram of Quality for All Pairs", 
               xlab="Quality (Expected Graft Survival)", breaks = c(0:25), ylim = c(0,8000))
onlineh = hist(agentsGreedy50$get_quality, col=rgb(1,0,0,1/4), breaks = c(0:25), add=T)
legend("topright", c("Online", "Greedy"), col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), fill=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))

#barplot of recieved type by algorithm
#table gets counts of each category
r = cbind(table(onlinec$get_type),table(greedyc$get_type),table(onlinei$get_type),table(greedyi$get_type))
r[3,1:2] = 0
#mar gives margins of c(bottom,left,top,right), xpd means things can be written in margin
par(mar=c(6.1, 3.1, 2.1, 3.1), xpd=TRUE)
barplot(r, ylim = c(0, 25000), col = c("lightblue2", "palegreen", "palegoldenrod"), 
        names = c("C-Online", "C-Greedy", "I-Online", "I-Greedy"), main = "Matches By Algorithm")
#inset is how far from the plot
legend("bottom", inset = c(-.45,-.45), legend = 
         c("Recieved kidney from a compatible pair", "Recieved kidney from an incompatible pair", "No match"), 
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
         c("Recieved kidney from a compatible pair", "Recieved kidney from an incompatible pair", "No match"), 
       fill = c("lightblue2", "palegreen", "palegoldenrod"), cex = .75)

#boxplot of quality by beta if matched with compatible
par(mar=c(5,4,4,4))
fmatched = onlinei[which(onlinei$get_type == 'C'),]
boxplot(fmatched[which(fmatched$betaCat == '<5'), ]$get_quality, 
        fmatched[which(fmatched$betaCat == '5-10'), ]$get_quality,
        fmatched[which(fmatched$betaCat == '10-15'), ]$get_quality, 
        fmatched[which(fmatched$betaCat == '15-20'), ]$get_quality,
        fmatched[which(fmatched$betaCat == '>20'), ]$get_quality, 
        main = "Quality Recieved When Matched With Compatible", xlab = "Beta Value", ylab = "Quality",
        names= c("<5","5-10","10-15","15-20", ">20"))

#boxplot of avg time by beta if matched with compatible
boxplot(fmatched[which(fmatched$betaCat == '<5'), ]$time, 
        fmatched[which(fmatched$betaCat == '5-10'), ]$time,
        fmatched[which(fmatched$betaCat == '10-15'), ]$time, 
        fmatched[which(fmatched$betaCat == '15-20'), ]$time,
        fmatched[which(fmatched$betaCat == '>20'), ]$time, 
        main = "Time Matched With Compatible", xlab = "Beta Value", ylab = "Time",
        names= c("<5","5-10","10-15","15-20", ">20"))

"""
greedyc = subset(agentsGreedy50, grepl("C", id))
greedyi = subset(agentsGreedy50, grepl("I", id))

onlinec = subset(agentsOnline50, grepl("C", id))
onlinei = subset(agentsOnline50, grepl("I", id))

counts = table(onlinec$get_type, greedyc$get_type, onlinei$get_type, greedyi$get_type)
barplot(counts)
counts = table(greedyc$get_type, onlinec$get_type)
counts = table(greedyi$get_type, onlinei$get_type )
#histogram of compatibles
greedych = hist(greedyc$get_quality, col=rgb(0,0,1,1/4), , main = "Histogram of Quality for Compatible Pairs", 
                xlab="Quality (Expected Graft Survival)", breaks = c(1:25))
onlinech = hist(onlinec$get_quality, col=rgb(1,0,0,1/4), breaks = c(1:25), add=T)
legend("topright", c("Greedy", "Online"), col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), fill=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))

#histogram of incompatibles
greedyih = hist(greedyi$get_quality, col=rgb(0,0,1,1/4), , main = "Histogram of Quality for Incompatible Pairs", 
                xlab="Quality (Expected Graft Survival)", breaks = c(1:25))
onlineih = hist(onlinec$get_quality, col=rgb(1,0,0,1/4), breaks = c(1:25), add=T)
legend("topright", c("Greedy", "Online"), col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), fill=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))

#quality gotten by beta value
fmatched = onlinei[which(onlinei$get_quality > 0),]
percentMatched = c(length(intersect(matched$id, beta1$id)) / nrow(beta1), 
                   length(intersect(matched$id, beta2$id)) / nrow(beta2),
                   length(intersect(matched$id, beta3$id)) / nrow(beta3),
                   length(intersect(matched$id, beta4$id)) / nrow(beta4),
                   length(intersect(matched$id, beta5$id)) / nrow(beta5))
boxplot(onlinei[which(onlinei$beta <= 5), ]$get_quality, 
        onlinei[which(onlinei$beta > 5 & onlinei$beta <= 10), ]$get_quality,
        onlinei[which(onlinei$beta > 10 & onlinei$beta <= 15), ]$get_quality, 
        onlinei[which(onlinei$beta > 15 & onlinei$beta <= 20), ]$get_quality,
        onlinei[which(onlinei$beta > 20), ]$get_quality, 
        main = "Quality Recieved by Beta Value", xlab = "Beta Value", ylab = "Quality",
        names= c("0-5","5-10","10-15","15-20", ">20"))
#quality given by beta
boxplot(onlinei[which(onlinei$beta <= 5), ]$give_quality, 
        onlinei[which(onlinei$beta > 5 & onlinei$beta <= 10), ]$give_quality,
        onlinei[which(onlinei$beta > 10 & onlinei$beta <= 15), ]$give_quality, 
        onlinei[which(onlinei$beta > 15 & onlinei$beta <= 20), ]$give_quality,
        onlinei[which(onlinei$beta > 20), ]$give_quality, 
        main = "Quality Given by Beta Value", xlab = "Beta Value", ylab = "Quality",
        names= c("0-5","5-10","10-15","15-20", ">20"))
compat_match = subset(onlinei, onlinei$time <= 50)
incompat_match = subset(onlinei, onlinei$time == 51)
no_match = subset(onlinei, onlinei$time == 52)
beta1 = onlinei[which(onlinei$beta <= 5), ]
beta2 = onlinei[which(onlinei$beta > 5 & onlinei$beta <= 10), ]
beta3 = onlinei[which(onlinei$beta > 10 & onlinei$beta <= 15), ]
beta4 = onlinei[which(onlinei$beta > 15 & onlinei$beta <= 20), ]
beta5 = onlinei[which(onlinei$beta > 20), ]

#percentage of matched incompatible pairs that match with another incompatible pair
p = c(length(intersect(incompat_match$id, beta1$id)) / nrow(beta1), 
      length(intersect(incompat_match$id, beta2$id)) / nrow(beta2),
      length(intersect(incompat_match$id, beta3$id)) / nrow(beta3),
      length(intersect(incompat_match$id, beta4$id)) / nrow(beta4),
      length(intersect(incompat_match$id, beta5$id)) / nrow(beta5))
barplot(p*100, xlab = "Beta Value", ylab = "Percentage", 
        names = c("0-5","5-10","10-15","15-20", ">20"), ylim = c(0,100), 
        main = "Percentage of Incompatible Matches \n After Online Portion")
#time that incompatible pairs matched with compatible pairs by beta value 
boxplot(compat_match[which(compat_match$beta <= 5), ]$time, 
        compat_match[which(compat_match$beta > 5 & compat_match$beta <= 10), ]$time,
        compat_match[which(compat_match$beta > 10 & compat_match$beta <= 15), ]$time, 
        compat_match[which(compat_match$beta > 15 & compat_match$beta <= 20), ]$time,
        compat_match[which(compat_match$beta > 20), ]$time, 
        main = "Time Matched by Beta Value", xlab = "Beta Value", ylab = "Time Matched",
        names= c("0-5","5-10","10-15","15-20", ">20"))
boxplot(greedyc$get_quality, onlinec$get_quality)
boxplot(greedyi$get_quality, onlinei$get_quality)



setwd("~/desktop")
library(ggplot2)
library(gridExtra)
df <- data.frame(cases=c(1150, 1294, 879, 756), lower = c(935, 943, 879,756), upper=c(3146,2484,879,756),year=c(2015,2016, 2015, 2016),Cases=c("Estimated","Estimated", "Observed", "Observed"))
df.S <- data.frame(cases=c(823,327,1104,190), lower = c(636, 299,773,170), upper=c(1464,1682,2183,301),year=c(2015,2015,2016,2016),Sex=c("Female","Male","Female","Male"))
df.A <- data.frame(cases=c(495,655,692,602), lower = c(365,570,460,483), upper=c(938,2208,1393,1091),year=c(2015,2015,2016,2016),Age=c("Minor","Adult","Minor","Adult"))
df.E <- data.frame(cases=c(702,448,1047,247), lower = c(548,387,717,226), upper=c(1230,1916,2110,374),year=c(2015,2015,2016,2016),Exploitation=c("Sexual","Other","Sexual","Other"))
df.D <- data.frame(cases=c(567,583,338,956), lower = c(498,437,321,622), upper=c(2037,1109,431,2053),year=c(2015,2015,2016,2016),Destination=c("External","Internal","External","Internal"))

# png("TOTAL.png", width=540, height=380, res=500)
ggplot(data=df, aes(x=year, y=cases, color=Cases)) + geom_point() + geom_line() + ylim(0, 3500) + ggtitle("Number of Human Trafficking Cases, Romania") + 
  scale_x_continuous(breaks = seq(2015, 2016, by = 1)) + scale_color_manual(values=c("#F8766D", "#000000")) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = Cases),  linetype=2, alpha = 0.2) + scale_fill_manual(values=c("#F8766D", "#000000")) + 
  theme(plot.title = element_text(size=12, hjust = 0.5), axis.text=element_text(size=10), axis.title=element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10)) + 
  xlab("Year") + ylab("Cases")
ggsave(filename = "TOTAL.png", units = "mm", path = "~/desktop", width = 110, height = 75, device='png', dpi=500)

SEX <- ggplot(data=df.S, aes(x=year, y=cases, colour=Sex)) + geom_point() + geom_line() + ggtitle("Sex") + ylim(0, 2500) +
    scale_x_continuous(breaks = seq(2015, 2016, by = 1)) + 
    geom_ribbon(aes(ymin=lower, ymax=upper, fill = Sex),linetype=2, alpha = 0.2) + 
  theme(plot.title = element_text(size=12, hjust = 0.5), axis.text=element_text(size=10), axis.title=element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10)) + 
  xlab("Year") + ylab("Estimated Cases")

AGE <- ggplot(data=df.A, aes(x=year, y=cases, colour=Age)) + geom_point() + geom_line() + ggtitle("Age") + ylim(0, 2500) +
  scale_x_continuous(breaks = seq(2015, 2016, by = 1)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = Age), linetype=2, alpha = 0.2) + 
  theme(plot.title = element_text(size=12, hjust = 0.5), axis.text=element_text(size=10), axis.title=element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10)) + 
  xlab("Year") + ylab("Estimated Cases")

EXP <- ggplot(data=df.E, aes(x=year, y=cases, colour=Exploitation)) + geom_point() + geom_line() + ggtitle("Exploitation") + ylim(0, 2500) +
  scale_x_continuous(breaks = seq(2015, 2016, by = 1)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = Exploitation), linetype=2, alpha = 0.2) + 
  theme(plot.title = element_text(size=12, hjust = 0.5), axis.text=element_text(size=10), axis.title=element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10)) + 
  xlab("Year") + ylab("Estimated Cases")

DES <- ggplot(data=df.D, aes(x=year, y=cases, colour=Destination)) + geom_point() + geom_line() + ggtitle("Destination") + ylim(0, 2500) +
  scale_x_continuous(breaks = seq(2015, 2016, by = 1)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = Destination), linetype=2, alpha = 0.2) + 
  theme(plot.title = element_text(size=12, hjust = 0.5), axis.text=element_text(size=10), axis.title=element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10)) + 
  xlab("Year") + ylab("Estimated Cases")

p <- grid.arrange(SEX, AGE, ncol=2)
ggsave(filename = "SEXAGE.png", plot = p, units = "mm", path = "~/desktop", width = 220, height = 75, device='png', dpi=500)
q <- grid.arrange(EXP, DES, ncol=2)
ggsave(filename = "EXPDES.png", plot = q, units = "mm", path = "~/desktop", width = 220, height = 75, device='png', dpi=500)


X = rnorm(50, mean=0, sd=1) 
w <- which(X <= 0) 
X[w] <- 10*X[w]
print(X) 
wmin <- which.min(X) 
wmax <- which.max(X)  
which.min( abs(X - mean(X)) )
m <- Inf; wm <- 0 
for(i in 1:length(X))
{
   if( abs(X[i] - mean(X)) < m ) 
   {
      m <- abs(X[i] - mean(X))
      wm <- i
   }
}
Z <- ( (X <= 1) & (X >= 0) )  
Z
mean(Z) 
count <- 0 
for(i in 1:length(X))
{
   if( (X[i] <= 1) & (X[i] >= 0) ) 
   { 
      count <- count + 1
   }
}
count/length(X)
Y <- read.table('chips.dt', header=T, row.names=1) 
Y$Data 
Y$Data <- Y$Data/median(Y$Data) 
attach(Y)
Data <- rep(0, length(Data))
Data 
Y$Data
p<-.7
rbinom(1, 1, p)
(runif(1) < p)
p_l <- .41 
person <- "immune"
if( runif(1) < p )  person <- "susceptible" 
people <- rep("immune", 10000)
for(i in 1:10000)
{
 if( runif(1) < p_l ) people[i] <- "susceptible"
}
mean( people == "susceptible" ) 
A = matrix(0, 50, 50) 
for(i in 1:50){
for(j in 1:50){
if( (i*j) > 1000 ){
A[i,j] <- i/j } else A[i,j] <- i*j}}
for(i in 1:50){
   for(j in 1:50){
      if( (i*j) > 1000 ){
         A[i,j] <- i/j
      } else{
         A[i,j] <- i*j}}}
for(i in 1:50)
{
   for(j in 1:50)
   {
      if( (i*j) > 1000 ) 
      {
         A[i,j] <- i/j 
      }  else
      {
         A[i,j] <- i*j
      }
   }
}
install.packages("RSQLite", dependencies=T) 
library(RSQLite)
drive <- dbDriver("SQLite") 
connect <- dbConnect(drive, "baseball.db")
dbListTables(connect)
dbListFields(connect, "Salaries")
Salaries = dbReadTable(connect, "Salaries")
Years <- unique( Salaries$yearID ) 
MedSal <- matrix(0, length(Years), 2)
for(year in Years)
{
   Sal.year <- Salaries[which(Salaries$yearID == year), ] 
   if( length( unique(Sal.year$lgID) ) == 2 )
   {
     whichAL <- which(Sal.year$lgID == "AL")
     sal.AL <- Sal.year$salary[whichAL]
     sal.NL <- Sal.year$salary[-whichAL]
     MedSal[which(Years==year), ] <- c( median(sal.AL), median(sal.NL) ) 
  }
}
ix <- which(MedSal[,1] != 0)
Years <- Years[ix] 
MedSal <- MedSal[ix,]
IX <- sort( Years, index=T ) 
Years <- IX$x
MedSal <- MedSal[IX$ix, ]
plot(Years, MedSal[,1], xlab="Year", ylab="Median salary in hundred thousands", 
      main="Median Salary vs. Year", axes=F, col=3, ylim=c(1e5, 1.4e6)) 
lines(Years, MedSal[,1], col=3)
points(Years, MedSal[,2], col=4) 
lines(Years, MedSal[,2], col=4)
axis(1, seq(1980, 2005, by=5), seq(1980, 2005, by=5))
axis(2, seq(2e5, 1.4e6, by=2e5), seq(2, 14, by=2))
box()
dbListTables(connect)
dbListFields(connect, "Allstar")
D <- dbGetQuery(connect, "SELECT playerID from Allstar")
D <- unique(D)
D.AL <- dbGetQuery(connect, "SELECT playerID from Allstar WHERE lgID in ('AL') ")
D.AL <- unique(D.AL)
D.AL1972 <- dbGetQuery(connect,
    "SELECT playerID from Allstar WHERE lgID in ('AL') AND yearID>1972")
D.AL1972 <- unique(D.AL1972)
nrow(D.AL1972)/nrow(D)
nrow(D.AL1972)/nrow(D.AL)
AS <- dbGetQuery(connect, "SELECT playerID from Allstar")
HOF <- dbGetQuery(connect, "SELECT hofID, inducted from HallOfFame")
is.same <- function(pID, hID) (hID == paste(pID,"h",sep=""))
AllstarIDs <- unique(AS$playerID) 
HOFer <- unique( HOF$hofID )
got.in <- rep(0, length(HOFer))
for(i in 1:length(got.in))
{
    induct <- HOF$inducted[which(HOF$hofID==HOFer[i])]
    got.in[i] <- ( sum(induct=="Y") > 0 )    
}
HOFer <- HOFer[which(got.in == 1)]
MATCH <- rep(0, length(AllstarIDs))
for(i in 1:length(AllstarIDs))
{ 
   for(j in 1:length(HOFer))
   {  
      if( is.same( AllstarIDs[i], HOFer[j] ) )
      {     
         MATCH[i] <- 1
         break
      }      
   }
}
mean(MATCH) 
dbListFields(connect, "Master")
A = A[-which(is.na(A)),]








### Exercise 1.4

### The vector I'll use to illustrate order and rank contains times for the 2012 Olympic Men's 100 Meters listed in order of lane
athlete=c("Thompson", "Powell", "Gay", "Blake", "Gatlin", "Bolt", "Bailey", "Martina")
lane=c(2:9)
x=c(9.98, 11.99, 9.80, 9.75, 9.79, 9.63, 9.88, 9.94 ) #Final time
h=c(10.02, 9.94, 9.90, 9.85, 9.82, 9.87, 9.96, 9.91) #Heat time

### To compare order and rank, first check out the help menu
help(rank)
help(order)

rank(x)
## rank does what we expect; it tells us that the first value in x is the 7th when listed in ascending order, the second value in x is the 8th when listed in ascending order, and so forth
order(x)
## order is less intuitive; it is telling us that the smallest x value is in position 6, the next smallest in position 4, and so forth.
## order is an intermediary step in sorting.  Note that the following two commands yeild the same result
x[order(x)]
sort(x)

## Suppose we want to report the name of the winning athlete or the lane of the winning athlete
lane[rank(x)==1]
## We can get the same pointer by asking for the first position of the order vector
athlete[order(x)[1]]

##Other commands related to rank and order are which.min and which.max.  They give the position of the min and max
which.min(x)
order(x)[1]

### rep repeats its first argument.  The argument can be a single value, a vector, or even a list
help(rep)
rep(1,20)  ## gives a vector of 20 1's
rep(c(1:9), 20)  ## repeats the sequence 1,2,..,9  20 times for a vector of length 180
rep(c(1:9),length.out=20)  ## repeats the whole sequence twice, then a partial version to get to length 20
rep(c(1:9), each=20)  ## each element gets repeated 20 times
sort(rep(c(1:9),20)) ### note that this gives the same result as the previous command
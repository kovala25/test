setwd("~/Spring 2022 Research")

source("test.r")

# A LUAD single-cell dataset
# clone_sizes = c(4,128,537)
# mutations_LC18 = read.csv("LC18_mutation_counts.csv")
# non_misisng_mutations_LC18 = read.csv("LC18_nonmissing_mutation_counts.csv")

library(xlsx)
# truth = read.xlsx(file = "table_parameters_2_clones.xlsx", sheetIndex = 1)


# New 1000 cell simulated data
dat = read.csv("mutational_data_500_700_simulated_SFS_1_TRUTH.csv")
# dat = read.csv("mutational_data_500_700_simulated_SFS_14_TRUTH.csv")
# dat = read.csv("mutational_data_500_700_simulated_SFS_15_TRUTH.csv")
# dat = read.csv("mutational_data_500_700_simulated_SFS_19_TRUTH.csv")
# dat = read.csv("mutational_data_500_700_simulated_SFS_20_TRUTH.csv")
cell.status = as.vector(t(dat[1,-1]))
table(cell.status)
table(cell.status)/sum(table(cell.status))

rownames(dat) = dat[,1]
dat = dat[-1,-1]


temp = rowSums(dat[,-1])


#### Clean the data
singletons = c()
count = 1
cutoff = 1

for(i in 1:length(temp)){
 if(temp[i] <= cutoff){
   singletons[count] = i
   count = count + 1
 } 
}


non.singleton.dat = dat[-singletons,]


mutation.freqs = rowSums(non.singleton.dat)
hist(mutation.freqs, breaks = 200)
table(cell.status)/sum(table(cell.status))


# check = hist(mutation.freqs, breaks = 200)
# check.dat = as.data.frame(cbind(check$breaks, check$counts))
# check.dat2 = check.dat[which(check.dat$V1 > 300),]
# colSums(check.dat2) # 181 mutations of frequency > 300.
# table(cell.status)/sum(table(cell.status))


neutral = mutation.freqs[which(mutation.freqs <= 50)]
m = seq(2, 50, 1)
# For each "-tons" 
A = 16000
formula = (1/m)/(m-1)*A # times big A (some constant)
# formula2 = (1/(m*(m-1)) )*A # times big A (some constant)
# singleton.rate = (1000*18*log(1000*0.01271/2))/(2*1000*0.01271/2)
singleton.rate = (1000*1/0.01271)*log(1000*0.01271)
formula = c(singleton.rate, formula)

check = hist(neutral, breaks = length(unique(neutral)))
check$counts
hist(neutral, breaks = length(unique(neutral)))
lines(formula)


#####################################################
#####################################################
#####################################################
### 2.) Where do those humps at .6 and .46 come from?
#####################################################
#####################################################
#####################################################

check2 = hist(mutation.freqs, breaks = length(unique(neutral)))

## Get the breaks/counts together.
comp.dat = as.data.frame(cbind(check2$breaks, check2$counts))

## Divide dataset by subclone status.
clone0 = dat[, which(cell.status=="0")]
clone1 = dat[, which(cell.status=="1")]
clone2 = dat[, which(cell.status=="2")]

clones  = list(clone0, clone1, clone2)

## Clean the data for each subclone.

#### Clean the data

non.singleton.clones = list()

for(j in 1:3){
  singletons = c()
  count = 1
  cutoff = 1
  temp = rowSums(clones[[j]])
  temp.dat = clones[[j]]
  for(i in 1:length(temp)){
    if(temp[i] <= cutoff){
      singletons[count] = i
      count = count + 1
    } 
  }
  # Line is redundant, but for my own clarity
  temp.non.singleton.dat = temp.dat[-singletons,]
  non.singleton.clones[[j]] = temp.non.singleton.dat
  
}

mutation.freqs = rowSums(non.singleton.clones[[1]])
hist(mutation.freqs, breaks = 100)
mutation.freqs = rowSums(non.singleton.clones[[2]])
hist(mutation.freqs, breaks = 100)
mutation.freqs = rowSums(non.singleton.clones[[3]])
hist(mutation.freqs, breaks = 100)



# which mutations have a frequency between 600-620?
temp = rowSums(dat)
odd.hump1.mutations = temp[which(temp > 590 & temp < 630)]

# So some of the cells in subclone 1 have all of these mutations, and some do not. 
# ... But all the cells in subclone 2 have these mutations.
odd.hump1.placements = rbind(cell.status, dat[which(temp > 590 & temp < 630),])
write.csv(odd.hump1.placements, file= "590-630_hump.csv")
# Same as above for the grouping below.
odd.hump2.placements = rbind(cell.status, dat[which(temp > 400 & temp < 500),])
write.csv(odd.hump2.placements, file= "400-500_hump.csv")

# 
odd.hump3.placements = rbind(cell.status, dat[which(temp > 150 & temp < 200),])
#
odd.hump4.placements = rbind(cell.status, dat[which(temp > 800 & temp < 820),])


odd.hump5.placements = rbind(cell.status, dat[which(temp > 350 & temp < 400),])
write.csv(odd.hump5.placements, file= "Expected_350-400_hump.csv")


#### 
#### 
#### 
#### 
table(cell.status)
clone2.shared.mutations = rbind(cell.status, dat[which(temp > 350 & temp < 400),])
temp.sums = colSums(clone2.shared.mutations[-1,])

# Precisely 351 cells, all with 18 shared mutations.
temp1c = temp.sums[which(temp.sums == dim(clone2.shared.mutations[-1,])[1])]

# 14 cells, with two sets of shared mutations (5 and 13).
temp1d = temp.sums[which(temp.sums < dim(clone2.shared.mutations[-1,])[1] & temp.sums > 0)]

temp1a = temp.sums[which(temp.sums > 0)]
temp1b = which(temp.sums > 0)
dat1 = as.data.frame(cbind(temp1b, cell.status[unique(temp1b)], temp.sums[unique(temp1b)]))
table(dat1$V2)

#####
#####
#####
#####
table(cell.status)
clone1.and.2.shared.mutations = rbind(cell.status, dat[which(temp > 750 & temp < 800),])
temp.sums = colSums(clone1.and.2.shared.mutations[-1,])

# Precisely 768cells with 74 shared mutations.
temp1c = temp.sums[which(temp.sums == dim(clone1.and.2.shared.mutations[-1,])[1])]

# 19 cells with different shared mutation amounts (62, 52, 3, 47, and 56).
temp1d = temp.sums[which(temp.sums < dim(clone1.and.2.shared.mutations[-1,])[1] & temp.sums > 0)]
unique(temp1d)

temp1b = which(temp.sums > 0)
dat2 = as.data.frame(cbind(temp1b, cell.status[unique(temp1b)]))
table(dat2$V2)


############# Ask Khanh ###############
############# Ask Khanh ###############
############# Ask Khanh ###############
############# Ask Khanh ###############

# Ask about Newick form (or whatever else he is using)
# Ask about Newick form (or whatever else he is using)
# Ask about Newick form (or whatever else he is using)
# Ask about Newick form (or whatever else he is using)
# Check if possible to put branch length and # of mutations
# Check if possible to put branch length and # of mutations
# Check if possible to put branch length and # of mutations

# More spread in growth rates.
# More spread in growth rates.
# More spread in growth rates.
# More spread in growth rates.
# Draw tree 200-300 cells, low mutation rates, mark mutations.
# Draw tree 200-300 cells, low mutation rates, mark mutations.
# Draw tree 200-300 cells, low mutation rates, mark mutations.
# Draw tree 200-300 cells, low mutation rates, mark mutations.
# Draw tree 200-300 cells, low mutation rates, mark mutations.

############# Ask Khanh ###############
############# Ask Khanh ###############
############# Ask Khanh ###############
############# Ask Khanh ###############

####
####
####
####
table(cell.status)
odd.hump1.placements = rbind(cell.status, dat[which(temp > 590 & temp < 630),])
temp.sums = colSums(odd.hump1.placements[-1,])

# 612 cells with 21 shared mutations.
temp1c = temp.sums[which(temp.sums == dim(odd.hump1.placements[-1,])[1])]
unique(temp1c)

# none.
temp1d = temp.sums[which(temp.sums < dim(odd.hump1.placements[-1,])[1] & temp.sums > 0)]
unique(temp1d)

# the rest are 0.
test = temp.sums[which( temp.sums ==0)]

temp1b = which(temp.sums > 0)
dat2 = as.data.frame(cbind(temp1b, cell.status[unique(temp1b)]))
table(dat2$V2)
# 261 cells are from clone 1 and 351 cells are from clone 2.

##### 
##### 
##### 
#####
table(cell.status)
odd.hump2.placements = rbind(cell.status, dat[which(temp > 400 & temp < 500),])

temp.sums = colSums(odd.hump2.placements[-1,])


# 457 cells have 28 shared mutations.
temp1c = temp.sums[which(temp.sums == dim(odd.hump2.placements[-1,])[1])]
unique(temp1c)

# 5 cells have 15 shared mutations.
temp1d = temp.sums[which(temp.sums < dim(odd.hump2.placements[-1,])[1] & temp.sums > 0)]
unique(temp1d)

temp1b = which(temp.sums > 0)
dat2 = as.data.frame(cbind(temp1b, cell.status[unique(temp1b)]))
table(dat2$V2)
# 111 clone1 cells and 351 clone 2 cells.


#####
#####
#####


# How many mutations are shared by ALL clone 2 cells?
# clone2.dat = non.singleton.clones[[3]]
clone2.labels = which(cell.status=="2")
clone2.dat = dat[,clone2.labels]

shared.by.all.clone2 = c()
count = 1

for(i in 1:dim(clone2.dat)[1]){
    temp = sum(clone2.dat[i,])
    
    if(temp==dim(clone2.dat)[2]){
      shared.by.all.clone2[count] = i
      count = count + 1
    }
}

# 157 shared mutations by all cells in clone 2.
shared.by.all.clone2
length(shared.by.all.clone2)

###
###
###
###

# clone1.dat = non.singleton.clones[[2]]
clone1.labels = which(cell.status=="1")
clone1.dat = dat[,clone1.labels]

shared.by.all.clone1 = c()
count = 1

for(i in 1:dim(clone1.dat)[1]){
  temp = sum(clone1.dat[i,])
  
  if(temp==dim(clone1.dat)[2]){
    shared.by.all.clone1[count] = i
    count = count + 1
  }
}


# 90 mutations shared by all cells in clone 1.
shared.by.all.clone1
length(shared.by.all.clone1)

####
####
####
####
# clone0.dat = non.singleton.clones[[1]]
clone0.labels = which(cell.status=="0")
clone0.dat = dat[,clone0.labels]

shared.by.all.clone0 = c()
count = 1

for(i in 1:dim(clone0.dat)[1]){
  temp = sum(clone0.dat[i,])
  
  if(temp==dim(clone0.dat)[2]){
    shared.by.all.clone0[count] = i
    count = count + 1
  }
}


# 10 mutations shared by all cells in clone 0.
shared.by.all.clone0

shared.by.all.clone1
shared.by.all.clone2

matches = c()
match.index = c()
count = 1

for(i in 1:length(shared.by.all.clone1)){
  for(j in 1:length(shared.by.all.clone2)){
    if(shared.by.all.clone1[i]==shared.by.all.clone2[j]){
      matches[count] = shared.by.all.clone1[i]
      match.index[count] = i
      count = count + 1
    }
  }
}

matches
match.index

# 67 mutations unique to clone 2
unique.clone2.mutations = shared.by.all.clone2[-match.index]
# 80 other mutations from clone 1
unique.clone1.mutations = shared.by.all.clone1[-(1:10)]

total = length(unique.clone2.mutations) + length(unique.clone1.mutations) + 10
total

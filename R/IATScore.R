#' Package to score Implicit Association Test (IAT) output
#' @description This minimalist package developed by Daniel Storage is designed to quickly score raw data outputted from an Implicit Association Test (IAT; Greenwald, McGhee, & Schwartz, 1998). IAT scores are calculated as specified by Greenwald, Nosek, and Banaji (2003). Outputted values can be interpreted as effect sizes. Refer to the DESCRIPTION file for more detailed help information (or contact me at sabbystorage@gmail.com).
#' @param IAT The name of the dataset to be analyzed.
#' @param Trials The number of trials across your entire IAT. The default is set to 220, which is typical of most IATs.
#' @param First Whether participants first sorted Congruent or Incongruent trials. The default is set to Congruent.
#' @keywords IAT
#' @examples
#' IATScore(IAT)
#' IATScore(IAT, Trials=220)
#' IATScore(TooFastIAT, Trials=220)
#' IATScore(BriefIAT, Trials=140)
#' IATScore(IAT, Trials=220, First="Congruent")
#' IATScore(IAT, Trials=220, First="Incongruent")
#' @export

IATScore <- function(IAT, Trials, First){

start.time <- Sys.time()

colnames(IAT) <- c("Block", "Trial", "Category", "Cat_Item", "Correct", "RT")

RT=Block=Correct=sd=NULL

# Set default number of trials to 221, which is standard for the IAT
if(missing(Trials)){
  Trials = 221
} else {
  Trials = Trials + 1
}

# Step 1: Delete any reaction times > 10,000 ms ####
i <- nrow(IAT) # define i counting variable for while loop
while (i < Trials) { # define while loop for Step 1
  if (IAT$RT[i] > 10000) {IAT$RT[i] <- 0}
  i = i + 1
}
IAT2 <- subset(IAT, RT!=0) # new data frame, excluding trials over 10,000 ms

# Step 2: Check for exclusion based on response speed (10% trials < 300 ms) ####
SpeedCount <- length(which(IAT2$RT<300)) # count number of RTs under 300
SpeedCount # display the number of RTs under 300
SpeedProp <- SpeedCount/nrow(IAT2) # calculate proportion of RTs under 300
SpeedProp # display proportion of RTs under 300

# Step 3: Compute means of correct trials in blocks 2, 3, 5, and 6 ####
Block2trials <- subset(IAT2, Block==2) # subset data frame for only Block 2 trials
Block2correct <- subset(Block2trials, Correct==0) # subset data frame for only correct trials
Block2correctMean <- mean(Block2correct$RT) # mean of Block 2 correct trials
Block3trials <- subset(IAT2, Block==3) # subset data frame for only Block 3 trials
Block3correct <- subset(Block3trials, Correct==0) # subset data frame for only correct trials
Block3correctMean <- mean(Block3correct$RT) # mean of Block 3 correct trials
Block5trials <- subset(IAT2, Block==5) # subset data frame for only Block 5 trials
Block5correct <- subset(Block5trials, Correct==0) # subset data frame for only correct trials
Block5correctMean <- mean(Block5correct$RT) # mean of Block 5 correct trials
Block6trials <- subset(IAT2, Block==6) # subset data frame for only Block 6 trials
Block6correct <- subset(Block6trials, Correct==0) # subset data frame for only correct trials
Block6correctMean <- mean(Block6correct$RT) # mean of Block 6 correct trials

# Step 4: Replace incorrect trials with avg RT by block + 600 ####
newBlock2 <- Block2correctMean + 600
newBlock3 <- Block3correctMean + 600
newBlock5 <- Block5correctMean + 600
newBlock6 <- Block6correctMean + 600

i <- 1 # define i counting variable for while loop
while (i < nrow(IAT2) + 1) { # create while loop for Block 2 incorrect trial replacement
  if (IAT2$Block[i] == 2 && IAT2$Correct[i] == 1) {
  IAT2$RT[i] <- newBlock2 }
  i <- i + 1
}

i <- 1 # define i counting variable for while loop
while (i < nrow(IAT2) + 1) { # create while loop for Block 3 incorrect trial replacement
  if (IAT2$Block[i] == 3 && IAT2$Correct[i] == 1) {
    IAT2$RT[i] <- newBlock3 }
  i <- i + 1
}

i <- 1 # define i counting variable for while loop
while (i < nrow(IAT2) + 1) { # create while loop for Block 5 incorrect trial replacement
  if (IAT2$Block[i] == 5 && IAT2$Correct[i] == 1) {
    IAT2$RT[i] <- newBlock5 }
  i <- i + 1
}

i <- 1 # define i counting variable for while loop
while (i < nrow(IAT2) + 1) { # create while loop for Block 6 incorrect trial replacement
  if (IAT2$Block[i] == 6 && IAT2$Correct[i] == 1) {
    IAT2$RT[i] <- newBlock6 }
  i <- i + 1
}

# Step 5: Calculate stdevs for all trials in blocks 2&5 and 3&6 ####
IAT3 <- subset(IAT2, Block==2 | Block==5) # pool all trials in blocks 2&5
sd25 <- sd(IAT3$RT) # stdev of all trials in blocks 2&5
IAT4 <- subset(IAT2, Block==3 | Block==6) # pool all trials in blocks 3&6
sd36 <- sd(IAT4$RT) # stdev of all trials in blocks 3&6

# Step 6: Calculate means for all trials in blocks 2, 3, 5, and 6 ####
IAT5 <- subset(IAT2, Block==2) # only Block 2 trials
block2mean <- mean(IAT5$RT) # block 2 mean
IAT6 <- subset(IAT2, Block==3) # only Block 3 trials
block3mean <- mean(IAT6$RT) # block 3 mean
IAT7 <- subset(IAT2, Block==5) # only Block 5 trials
block5mean <- mean(IAT7$RT) # block 5 mean
IAT8 <- subset(IAT2, Block==6) # only Block 6 trials
block6mean <- mean(IAT8$RT) # block 6 mean

# Step 7: Compute mean differences in test trials ####
meandiff1 <- block5mean - block2mean # first mean difference (blocks 5 - 2)
meandiff2 <- block6mean - block3mean # second mean difference (blocks 6 - 3)

# Step 8: Divide each difference by its associated pooled stdev ####
value1 <- meandiff1/sd25
value2 <- meandiff2/sd36

# Step 9: Average these two values to get your IAT effect size ####
IATeffect <- (value1+value2)/2 # calculate IAT effect size
IATeffect # print IAT effect size

# Sets the default value for "First" to "Congruent", so the function can run without specifying this argument
if(missing(First)){
  First = "Congruent"
}

# Flips the sign of the IAT effect if the IAT dataset comes from an "incongruent first" condition
if(First == "Incongruent"){
  IATeffect = IATeffect*(-1)
  message("NOTE: IAT effect size sign flipped for this participant (Incongruent First)")
}

end.time <- Sys.time()
time.taken <- end.time - start.time
message("Seconds taken to run this code:")
message(round(time.taken, 3))
message(" ")

# Tell the researcher to exclude this participant if they went too fast on over 10% of the trials (< 300 ms)
if (SpeedProp > 0.10) {
  message("NOTE: THE PARTICIPANT WENT TOO FAST ON THE IAT -- EXCLUDE")
  message("Proportion of trials with responses faster than 300 ms:")
  message(round(SpeedProp, 3))
  message(" ")
  message("The IAT effect size for this participant is:")
  return(IATeffect)
}

# Otherwise, simply give them the appropriate IAT effect size for this participant
if (SpeedProp < 0.10) {
  message("The IAT effect size for this participant is:")
  return(IATeffect)
}

}





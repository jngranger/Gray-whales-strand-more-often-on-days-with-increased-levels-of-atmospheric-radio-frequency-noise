####In this file you will find functions to do the following:
#Remove rows with NAs in specific columns
#Merge the null and stranding dataframes in the desired way
#Merge a bootstrapped version of the null and stranding dataframes in the desired way
#Bin RF, Ap and SS


#A function to remove rows with NA values in the desired column

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}


#Create functions to create a combined dataset with the strandings and the null data

create=function(LivedataUse,nulldata)
{
  
  #Add a column to both data sets that is a boolean for whether or not there was a "stranding"
  Yes=rep(1,nrow(LivedataUse)) #Yes, stranding for the data
  LivedataUse$Stranded=Yes
  No=rep(0, nrow(nulldata))    #No, no stranding for the null
  nulldata$Stranded=No
  
  #Trim both datasets to only the desired columns
  LivedataUse=LivedataUse[c("Date","AP", "SS", "PDO",
                            "RadioFlux", "Stranded", 
                            "DateAsCont", "RF.scale", "PDO.scale",
                            "AP.scale", "SS.scale","Year")]
  nulldata=nulldata[c("Date","AP", "SS", "PDO", 
                      "RadioFlux", "Stranded",
                      "DateAsCont", "RF.scale", "PDO.scale",
                      "AP.scale", "SS.scale", "Year")]
  
  #Choose a random sample from the null data, the same size as the stranding data frame
  randomsample=nulldata[sample(nrow(nulldata), nrow(LivedataUse)), ]
  rbind(LivedataUse, randomsample)  #Bind them together
  
}

#Do the same thing as before, but choose 600 samples from both datasets, with replacement
create.bootstrap=function(LivedataUse,nulldata)
{
    #Randomly sample 600 with replacement
  randomsample=nulldata[sample(nrow(nulldata), 600, replace=TRUE), ]
  
  #Add a column to both data sets that is a boolean for whether or not there was a "stranding"
  Yes=rep(1,nrow(LivedataUse))  #Yes, stranding for the data
  LivedataUse$Stranded=Yes
  No=rep(0, nrow(randomsample))  #No, no stranding for the null
  randomsample$Stranded=No
  
  #Trim both datasets to only the desired columns
  LivedataUse=LivedataUse[c("Date","AP", "SS", "PDO", "RadioFlux", "Season", "Stranded",
                            "DateAsCont", "RF.scale", "PDO.scale","AP.scale", "SS.scale","Year")]
  randomsample=randomsample[c("Date","AP", "SS", "PDO", "RadioFlux", "Season", "Stranded",
                              "DateAsCont", "RF.scale", "PDO.scale","AP.scale", "SS.scale","Year")]
  
  #Randomly sample 600 with replacement
  sampleLivedataUse=LivedataUse[sample(nrow(LivedataUse), 600, replace=TRUE),]
  #Bind together
  rbind(sampleLivedataUse, randomsample)
  
}


#Create Binning Functions for Sunspots (SS), RF, and Ap (K)



#Sunspot quintiles to know where the bins should start and end
#quantile(nulldata$SS, prob = seq(0, 1, length = 7), type = 5)

binSS=function(datause){
  #Group one =0-12
  zero=datause[datause$SS>=0,]  #Use only the data greater than 0, and less than 12
  zero=datause[datause$SS<12,]
  zero.dead=sum(zero$Stranded==0)  #How many were "No Stranding"?
  zero.alive=sum(zero$Stranded==1) #How many were "Stranding"?
  G1=zero.alive/zero.dead          #Take the proportion
  #group two =12-27                #Repeat for all bins
  one=datause[datause$SS>=12,]
  one=one[one$SS<27,]
  one.dead=sum(one$Stranded==0)
  one.alive=sum(one$Stranded==1)
  G2=one.alive/one.dead
  #Group three=27-56
  two=datause[datause$SS<56,]
  two=two[two$SS>=27,]
  two.dead=sum(two$Stranded==0)
  two.alive=sum(two$Stranded==1)
  G3=two.alive/two.dead
  #Group four= 56-95
  three=datause[datause$SS>=56,]
  three=three[three$SS<95,]
  three.dead=sum(three$Stranded==0)
  three.alive=sum(three$Stranded==1)
  G4=three.alive/three.dead
  #Group five= 95-151
  four=datause[datause$SS<151,]
  four=four[four$SS>=95,]
  four.dead=sum(four$Stranded==0)
  four.alive=sum(four$Stranded==1)
  G5=four.alive/four.dead
  #Group Six=151+
  five=datause[datause$SS>=151,]
  five.dead=sum(five$Stranded==0)
  five.alive=sum(five$Stranded==1)
  G6=five.alive/five.dead
  #Hist
  groups=c(G1, G2, G3, G4, G5, G6) #Combine the proportions
  groups
}


#Repeat the above code but for RF
#RF quintiles
  #nulldata.RF=completeFun(nulldata,"RadioFlux")
  #quantile(nulldata.RF$RadioFlux, prob = seq(0, 1, length = 7), type = 5)

binRF=function(datause){
  #Group one =0-72
  zero=datause[datause$RadioFlux>=0,]
  zero=datause[datause$RadioFlux<72,]
  zero.dead=sum(zero$Stranded==0)
  zero.alive=sum(zero$Stranded==1)
  G1=zero.alive/zero.dead
  #group two =72-82
  one=datause[datause$RadioFlux>=72,]
  one=one[one$RadioFlux<82,]
  one.dead=sum(one$Stranded==0)
  one.alive=sum(one$Stranded==1)
  G2=one.alive/one.dead
  #Group three=82-100
  two=datause[datause$RadioFlux<100,]
  two=two[two$RadioFlux>=82,]
  two.dead=sum(two$Stranded==0)
  two.alive=sum(two$Stranded==1)
  G3=two.alive/two.dead
  #Group four= 56-95
  three=datause[datause$RadioFlux>=100,]
  three=three[three$RadioFlux<126,]
  three.dead=sum(three$Stranded==0)
  three.alive=sum(three$Stranded==1)
  G4=three.alive/three.dead
  #Group five= 95-151
  four=datause[datause$RadioFlux<165,]
  four=four[four$RadioFlux>=126,]
  four.dead=sum(four$Stranded==0)
  four.alive=sum(four$Stranded==1)
  G5=four.alive/four.dead
  #Group Six=151+
  five=datause[datause$RadioFlux>=165,]
  five.dead=sum(five$Stranded==0)
  five.alive=sum(five$Stranded==1)
  G6=five.alive/five.dead
  #Hist
  groups=c(G1, G2, G3, G4, G5, G6)
  groups
}


#Repeat the above code but for Ap

#quantile((nulldata$AP), prob = seq(0, 1, length = 7), type = 5)

binK=function(datause)
{
  #Group one =0-9
  zero=datause[datause$AP>=(0),]
  zero=zero[zero$AP<(9),]
  zero.dead=sum(zero$Stranded==0)
  zero.alive=sum(zero$Stranded==1)
  G1=zero.alive/zero.dead
  #group two =9-15.25
  one=datause[datause$AP>=(9),]
  one=one[one$AP<(15.25),]
  one.dead=sum(one$Stranded==0)
  one.alive=sum(one$Stranded==1)
  G2=one.alive/one.dead
  #Group three=15.25-23.875
  two=datause[datause$AP>=(15.25),]
  two=two[two$AP<(23.875),]
  two.dead=sum(two$Stranded==0)
  two.alive=sum(two$Stranded==1)
  G3=two.alive/two.dead
  #Group four= 23.875-35.75
  three=datause[datause$AP>=(23.875),]
  three=three[three$AP<(35.75),]
  three.dead=sum(three$Stranded==0)
  three.alive=sum(three$Stranded==1)
  G4=three.alive/three.dead
  #Group five= 35.75-54.625
  four=datause[datause$AP>=(35.75),]
  four=four[four$AP<(54.625),]
  four.dead=sum(four$Stranded==0)
  four.alive=sum(four$Stranded==1)
  G5=four.alive/four.dead
  #Group Six=54.625+
  five=datause[datause$AP>=54.625,]
  five.dead=sum(five$Stranded==0)
  five.alive=sum(five$Stranded==1)
  G6=five.alive/five.dead
  #Hist
  groups=c(G1, G2, G3, G4, G5, G6)
  groups
}

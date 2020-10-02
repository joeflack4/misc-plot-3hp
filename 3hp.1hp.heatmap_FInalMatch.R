library("RColorBrewer")
library(gplots)

source("3hp.1hp.model_0116_FinalMatch.R")

#Change the setup to produce a heatmap image 

heatmap.params = data.frame(model.param)


#Modified Parameters
heatmap.params$C_rif = 0.231 / 2
heatmap.params$Pre_LTBI = 0.5

#n.values<-70

y.values = 27 #How many values along the y axis do we want
y.axis.label = "Completion for 1HP"
y.axis.min = 0.74 #Completion 1HP
y.axis.max = 1.0 
x.values = 21 #How many values along the x axis do we want
x.axis.label = "Efficacy of 1HP"
x.axis.min = 0.9 #0.63 #Efficacy 1HP
x.axis.max = 1.0

#Progressive, Iterative
y.axis.seq = rev(seq( from = y.axis.min, to = y.axis.max, length.out = y.values ))
x.axis.seq = seq( from = x.axis.min, to = x.axis.max, length.out = x.values )

sample.matrix  = matrix ( data = rep(heatmap.params, times=x.values * y.values),
                          ncol = ncol(heatmap.params),
                          nrow = x.values * y.values,
                          byrow = TRUE,
                          dimnames = list(c(),names(heatmap.params)))

sample.matrix [,'E_1HP'] = rep(x.axis.seq, each = y.values)
sample.matrix [,'P_1HPtreatment'] = rep(y.axis.seq, times = x.values)

#Perform your simulations using a combination of apply() and do.call()
#apply() will iteratively perform the same function to each row and/or column of a matrix argument
#The matrix of sampled input values will serve as our matrix argument X
#We will set MARGIN=1 to indicate that we want to iterate over the rows (not columns) of the matrix
#We will create a custom function to apply over each row (r)
#The key of the custom function is do.call() 
#do.call() calls a function of your choice and passes user-defined arguments to it
#we will use it to call our previously-defined function treatment.sim()
#We will supply it our row of sampled input values as user-defined arguments
#arguments for do.call() must be in the form of a list: as.list(r)
#apply() returns the outputs of each iteration as a column of a matrix
#for ease of use, we will transpose - t() - the matrix so that each iteration will be returned as a row

result.sim = t(apply(
  X=sample.matrix, 
  MARGIN=1, 
  FUN=function(r) do.call(what=treatment.sim, args=list(r))))

#Becase treatment.sim() returns a vector of 4 results 
#sim.out contains 4 columns, each corresponding to one of the outputs of treatment.sim()
#The first row of sim.out corresponds to the results of the simulation created by the first row of sampled input values in sampled.value
ICER = matrix(ncol=3, nrow=x.values*y.values, dimnames = list(c(), c("Incremental.Cost","Incremental.Effectiveness","Ratio")))


#Perform the ICER calculations
ICER[,"Incremental.Cost"] <- result.sim[,"Cost.1HP.total"]- result.sim[,"Cost.3HP.total"]
ICER[,"Incremental.Effectiveness"] <- result.sim[,"DALY.3HP.total"] - result.sim[,"DALY.1HP.total"]
ICER[,"Ratio"] <- ICER[,"Incremental.Cost"] / ICER[,"Incremental.Effectiveness"]
#sample.array <- cbind (sample.array,ICER[,"Ratio"])

limited.ICER = ICER[,"Ratio"]

ICER.upper.limit = 5000 #Preferred value 5000
ICER.lower.limit = 0 

limited.ICER[ limited.ICER > ICER.upper.limit ] = ICER.upper.limit
limited.ICER[ limited.ICER < ICER.lower.limit ] = ICER.lower.limit

#x.label.mask = c(FALSE,TRUE,TRUE,FALSE,TRUE,TRUE)
#y.label.mask = c(FALSE,TRUE,TRUE,FALSE,TRUE,TRUE)

y.axis.tick.label = round(y.axis.seq,3)
#y.axis.tick.label[y.label.mask] = ""

x.axis.tick.label = round(x.axis.seq,3)
#x.axis.tick.label[x.label.mask] = ""

heat.matrix =  t(matrix(data=limited.ICER, ncol=y.values, nrow=x.values, byrow = TRUE,
                        dimnames = list(
                          x.axis.tick.label,
                          y.axis.tick.label
                        )))

#Plot the heatmap

par( mar = c(5,2,2,3) )

heatmap.vals = heatmap.2(heat.matrix, 
                         dendrogram = 'none', 
                         #breaks = 12,
                         breaks = c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000),
                         Rowv=FALSE, 
                         Colv=FALSE, 
                         trace='none', 
                         scale='none',
                         xlab=x.axis.label,
                         ylab=y.axis.label,
                         cexRow=1.0,
                         cexCol=1.0,
                         symkey = FALSE,
                         #symbreaks = FALSE,
                         key = FALSE,
                         margins = c(5,5),
                         col=rev(brewer.pal(10,"RdBu"))
)

# icer_trunc = format(ICER.limit,big.mark=",", trim=TRUE)
# trunc_string = sprintf ("ICER values above %s have been set to %s", icer_trunc, icer_trunc)
# text (x=1, y=1, labels=trunc_string, cex=.8)

heatmap.legend.labels = array(data=NA,dim=c(nrow(heatmap.vals$colorTable), 1))

for (i in 1:nrow(heatmap.legend.labels)) {
  heatmap.legend.labels[i] = 
    sprintf("%.0f - %.0f", 
            round(heatmap.vals$colorTable[i,1],0), 
            round(heatmap.vals$colorTable[i,2],0)
    )
}

legend("left",      
       legend = heatmap.legend.labels,
       col=heatmap.vals$col,
       title="ICER",
       bty="n",
       lty= 1,             
       lwd = 8,
       cex=.9
)

mtext("Incremental Cost-Effectiveness of 1HP vs 3HP (2019 USD per DALY averted)\nLTBI prev = 0.5, Cost of rif = 0.1155", line = -3)
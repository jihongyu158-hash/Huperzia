library(phytools)
library(geiger)
library(nlme)
library(hisse)
library(devtools)
mydata <- read.table("trait.xls",row.names=1,header = TRUE)
mytree <- read.tree("23species_cds_beast_time_tree_onlyHuperzia.tree")
comparison <- name.check(phy=mytree,data=mydata)
comparison
name.check(phy=mytree,data=mydata)
comparison <- name.check(phy=mytree,data=mydata)

cat("Running Full HiSSE Model 1 \n")
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
pp1 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual,output.type="raw")
pp1

cat("Bissee model all free, run 2 \n")
turnover.anc = c(1,2,0,0)
eps.anc = c(1,2,0,0)
#if one wanted to run a Bisse model in Hisse
trans.rates.bisse = TransMatMaker.old(hidden.states=FALSE)
trans.rates.bisse
pp2 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse)
pp2

cat("Bissee model, extinction 0=1, run 3 \n")
turnover.anc = c(1,2,0,0)
eps.anc = c(1,1,0,0)
#if one wanted to run a Bisse model in Hisse
trans.rates.bisse = TransMatMaker.old(hidden.states=FALSE)
trans.rates.bisse
pp3 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse)
pp3

cat("Bissee model with equal q's, run 4 \n")
turnover.anc = c(1,2,0,0)
eps.anc = c(1,2,0,0)
#if one wanted to run a Bisse model in Hisse
trans.rates.bisse = TransMatMaker.old(hidden.states=FALSE)
trans.rates.bisse
trans.rates.bisse.equal = ParEqual(trans.rates.bisse, c(1,2))
trans.rates.bisse.equal
pp4 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse.equal)
pp4

cat("Bissee model with equal q's and e0=e1, run 5 \n")
turnover.anc = c(1,2,0,0)
eps.anc = c(1,1,0,0)
#if one wanted to run a Bisse model in Hisse
trans.rates.bisse = TransMatMaker.old(hidden.states=FALSE)
trans.rates.bisse
trans.rates.bisse.equal = ParEqual(trans.rates.bisse, c(1,2))
trans.rates.bisse.equal
pp5 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse.equal)
pp5

cat("2-state character independent CID-2 model, equal q's but different transition rates, run 6 \n")
turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,2,2)
#full 8 transition model
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
pp6 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")
pp6

cat("2-state character independent CID-2 model, equal q's and e's but different transition rates, run 7 \n")
turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,1,1)
#full 8 transition model
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
pp7 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")
pp7

cat("full 8 transition model \n")
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
pp8.hisse.null4 <- hisse.old(mytree, mydata, f=c(0.5,0.5), turnover.anc=rep(c(1,2,3,4),2),eps.anc=rep(c(1,2,3,4),2), trans.rate=trans.rates.nodual.allequal)
pp8.hisse.null4

cat("CD-4 model, q's and e's equal, run9 \n")
eps.anc = c(1,1,1,1)
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
pp9.hisse.null4 = hisse.old(mytree, mydata, f=c(0.5,0.5), turnover.anc=rep(c(1,2,3,4),2),eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal)
pp9.hisse.null4

cat("run the full model with transition rates equal, run 10 \n")
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)
#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates
#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
pp10 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")
pp10

cat("run the full model with transition rates and e's equal, run 11 \n")
turnover.anc = c(1,2,3,4)
eps.anc = c(1,1,1,1)
#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates
#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
pp11 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")
pp11


###################################################################################
#####run the full model with transition rates equal, and turnover 0a=1a=0b, and extinction 0a=1a=0b run 12
cat("run the full model with transition rates equal, and turnover 0a=1a=0b, and extinction 0a=1a=0b run 12 \n")

turnover.anc = c(1,1,1,2)
eps.anc = c(1,1,1,2)

#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp12 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp12
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=1a=0b, and extinction equal run 13
cat("run the full model with transition rates equal, and turnover 0a=1a=0b, and extinction equal run 13 \n")

turnover.anc = c(1,1,1,2)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp13 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp13
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=0b, and e0a=e0b, run 14
cat("run the full model with transition rates equal, and turnover 0a=0b, and e0a=e0b, run 14 \n")

turnover.anc = c(1,2,1,3)
eps.anc = c(1,2,1,3)

#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp14 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp14
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=0b, and extinction equal, run 15
cat("run the full model with transition rates equal, and turnover 0a=0b, and extinction equal, run 15 \n")

turnover.anc = c(1,2,1,3)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp15 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp15
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=1a, and e0a=e1a, run 16
cat("run the full model with transition rates equal, and turnover 0a=1a, and e0a=e1a, run 16 \n")

turnover.anc = c(1,1,2,3)
eps.anc = c(1,1,2,3)

#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp16 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp16
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=1a, and extinction equal, run 17
cat("run the full model with transition rates equal, and turnover 0a=1a, and extinction equal, run 17 \n")

turnover.anc = c(1,1,2,3)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp17 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp17
###################################################################################
#####run the full model with different turnover.anc, eps.anc, and but rates q0b1b=0,q1b0b=0, all other equals, run 18
cat("run the full model with different turnover.anc, eps.anc, and but rates q0b1b=0,q1b0b=0, all other equals, run 18 \n")
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)

#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp18 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp18
###################################################################################
#####run the full model with different e's equal and but rates q0b1b=0,q1b0b=0, all other equals, run 19
cat("run the full model with different e's equal and but rates q0b1b=0,q1b0b=0, all other equals, run 19 \n")

turnover.anc = c(1,2,3,4)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp19 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp19
###################################################################################
#####run the full model with r0a=r1a=r0b, e0a=e0b=e0b and rates q0b1b=0,q1b0b=0, all other equals, run 20
cat("run the full model with r0a=r1a=r0b, e0a=e0b=e0b and rates q0b1b=0,q1b0b=0, all other equals, run 20 \n")

turnover.anc = c(1,1,1,2)
eps.anc = c(1,1,1,2)

#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp20 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp20
###################################################################################
#####run the full model with r0a=r1a=r0b, extinction equals and rates q0b1b=0,q1b0b=0, all other equals, run 21
cat("run the full model with r0a=r1a=r0b, extinction equals and rates q0b1b=0,q1b0b=0, all other equals, run 21 \n")

turnover.anc = c(1,1,1,2)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp21 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp21
###################################################################################
#####run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 22
cat("run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 22 \n")

turnover.anc = c(1,2,1,3)
eps.anc = c(1,2,1,3)

#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp22 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp22
###################################################################################
#####run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 23
cat("run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 23 \n")

turnover.anc = c(1,2,1,3)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp23 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp23
###################################################################################
#####run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 24
cat("run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 24 \n")

turnover.anc = c(1,1,2,3)
eps.anc = c(1,1,2,3)

#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp24 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp24
###################################################################################
#####run the full model withr0a=r0b, extinction equal, and rates q0b1b=0,q1b0b=0, all other equals, run 25
cat("run the full model withr0a=r0b, extinction equal, and rates q0b1b=0,q1b0b=0, all other equals, run 25 \n")

turnover.anc = c(1,1,2,3)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp25 = hisse.old(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp25


###################################################################################
##Write a dataframe for comparison of modles

Model_test_leave_col<- data.frame(loglik= c(pp1[["loglik"]],pp2[["loglik"]],pp3[["loglik"]],pp4[["loglik"]],pp5[["loglik"]],pp6[["loglik"]],
                                            pp7[["loglik"]],pp8.hisse.null4[["loglik"]],pp9.hisse.null4[["loglik"]],pp10[["loglik"]],pp11[["loglik"]],pp12[["loglik"]],
                                            pp13[["loglik"]],pp14[["loglik"]],pp15[["loglik"]],pp16[["loglik"]],pp17[["loglik"]],pp18[["loglik"]],
                                            pp19[["loglik"]],pp20[["loglik"]],pp21[["loglik"]],pp22[["loglik"]],pp23[["loglik"]],pp24[["loglik"]],pp25[["loglik"]]),
                                  AIC= c(pp1[["AIC"]],pp2[["AIC"]],pp3[["AIC"]],pp4[["AIC"]],pp5[["AIC"]],pp6[["AIC"]],
                                         pp7[["AIC"]],pp8.hisse.null4[["AIC"]],pp9.hisse.null4[["AIC"]],pp10[["AIC"]],pp11[["AIC"]],pp12[["AIC"]],
                                         pp13[["AIC"]],pp14[["AIC"]],pp15[["AIC"]],pp16[["AIC"]],pp17[["AIC"]],pp18[["AIC"]],
                                         pp19[["AIC"]],pp20[["AIC"]],pp21[["AIC"]],pp22[["AIC"]],pp23[["AIC"]],pp24[["AIC"]],pp25[["AIC"]]))

write.csv(Model_test_leave_col,"Model_test_leave_col.csv")
setwd("C:/Users/84289/Desktop/paper/11.性状与多样化/叶脉")
getwd()
pp.recon_leave_col <- MarginRecon.old(mytree, mydata, f=c(0.5,0.5), pars=pp21$solution, hidden.states=TRUE)
save(pp.recon_leave_col, file="full_HiSSE_recon_leave_col.RData")
rates_leave_col<- GetModelAveRates(pp.recon_leave_col, type = "tips")
save(rates_leave_col, file="Tip_rates_fullHiSSE_model_leave_col.RData")
write.table(rates_leave_col, file="Tip_rates_fullHiSSE_model_leave_col.txt",quote=FALSE,sep="\t",row.names=FALSE)

#after combinging state calls into one new column called states
rates <- read.table("Tip_rates_fullHiSSE_model_leave_col.txt", header=TRUE)

pdf(file="leave_col_recon.pdf",20,20)
plot.hisse.states(pp.recon_leave_col, type = "phylogram",rate.param="net.div", show.tip.label=TRUE,
                  fsize=3,do.observed.only=TRUE, rate.colors=c("steelblue","tomato"))

dev.off()
legend("topleft", legend=c("0","1"), fill=c("steelblue","tomato"),cex=2,title="States")

#
pdf(file="speciation_leave_col.pdf",20,20)
plot.hisse.states(pp.recon_leave_col, type = "phylogram",rate.param="speciation", show.tip.label=TRUE,
                  fsize=3,do.observed.only=TRUE, rate.colors=c("steelblue","tomato"))

dev.off()
pdf(file="extinction_leave_col.pdf",20,20)
plot.hisse.states(pp.recon_leave_col, type = "phylogram",rate.param="extinction", show.tip.label=TRUE,
                  fsize=3,do.observed.only=TRUE, rate.colors=c("steelblue","tomato"))

dev.off()
#
pdf(file="Boxpot_netdiv_leave_col.pdf")
boxplot(rates$net.div~rates$state, boxwex=0.5, 
        notch = TRUE,main="Net Diversification",
        col = c("steelblue","tomato"), 
        xlab="Same vs. Different", ylab="Diversification rate")
dev.off()

pdf(file="Boxplot_speciation_leave_col.pdf")
boxplot(rates$speciation~rates$state, boxwex=0.5,
        notch = TRUE,main="Speciation",
        col = c("steelblue","tomato"), 
        xlab="Same vs. Different", ylab="Speciation rate")
dev.off()

pdf(file="Boxplot_extinction_leave_col.pdf")
boxplot(rates$extinction~rates$state,boxwex=0.5, 
        notch = TRUE,main="Extinction", 
        col = c("steelblue","tomato"),
        xlab="Same vs. Different", ylab="Extinction rate")
dev.off()


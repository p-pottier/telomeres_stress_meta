## 2019/08/06
## telomeres-stress meta-analysis - code for Chatelain et al. 2020
## Szymek Drobniak (szymek.drobniak_at_gmail.com)


# load necessary packages ----
####


library(ape)
# library(asreml)
# library(MCMCglmm)
library(metafor)
library(MuMIn)
eval(metafor:::.MuMIn)
# commented-out packages left here for compatibility with earlier version and in case they're needed


# load data ----
####

data <- read.table("./data/data_stressor exposure_telomeres.csv", sep=";", head=T)
data <- read.table("./data/190801_FINAL_DATA_cleaned_SZD.csv", sep=";", head=T)

# fix specific levels in data
# for species where no phylogenetic information is available
# here we assign a sister species wherever possible
levels(data$species)[1] <- "Acipenser fulvescens"
levels(data$species)[25] <- "Myoxus glis"
levels(data$species)[37] <- "Oncorhynchus nerka mykiss"
levels(data$species)[c(65,66)] <- "Taeniopygia guttata"

# add underscore to all names
levels(data$species) <- gsub(" ", "_", levels(data$species), fixed=T)

# relevel factors with non-optimal (e.g. containing spaces) factor names or redundant levels
levels(data$measure) <- c("tel_len", "tel_short")
levels(data$corr) <- c("N", "N", "Y")

# calculate sampling variances and weights
data$mvar <- 1/(data$n_tot-3)
data$weights <- 1/data$mvar

# this is the tau-corrected V_m, tau = 0.1098 in our case
# tau comes from initial models - you have to estimate it yourself
# as overall between-study variance 
data$weights_tau <- data$n_tot-3+(1/0.1098)

# remove invalid data if any
data.analysis <- data[!is.na(data$weights),]

# simplify the age variable to single-out studies with juveniles
data.analysis$age2 <- data.analysis$age
levels(data.analysis$age2) <- c("AD","AD","JUV")

# produce duplicate species column
data.analysis$species2 <- data.analysis$species




# load phylogeny ----
####

# the trees where combined root to root from individual trees manually directly in Newick files
# all branch length set to unity due to nonuniform units in original phylogenies

all_taxa <- read.tree("_mamm_bir_rep_amph_fish.txt")
plot(all_taxa, root.edge = T, cex = 0.6)

all_taxa$edge.length <- NULL
all_taxa$node.label <- NULL
all_taxa <- compute.brlen(all_taxa, method="Grafen")



# pre-process data ----
####

# CURRENTLY REDUNDANT
# compute the inverse of phylogenetic correlation matrix
IA <- inverseA(all_taxa, nodes = "TIPS")
IA.asreml <- sm2asreml(IA$Ainv, IA$node.names)

# CURRENTLY REDUNDANT
# calculate sigma_m
sigma_m <- sum(data.analysis$weights*nrow(data.analysis))/(sum(data.analysis$weights)^2-sum(data.analysis$weights^2))




# generate full and intercept models ----
####

# generate phylogenetic VCV matrix
vcv.phylo(all_taxa, cor = T) -> taxa_cor


# full and intercept models without phylogenetic effect (1a and 2a)
model1a <- rma.mv(yi = Final_Z,
                mods = ~age2+sex+phylum+wild_captive+tissue_lifespan+stressor_cat+exper_correl + measure+method_T+corr+accounting_for_age +
                  measure:age2 + measure:stressor_cat + tissue_lifespan:stressor_cat + age2:stressor_cat,
                V = mvar, data = data.analysis, test = "t", random = list(~1|species, ~1|DOI, ~1|ES_id))
model2a <- rma.mv(yi = Final_Z, V = mvar, data = data.analysis, test = "t", random = list(~1|species, ~1|DOI, ~1|ES_id))
summary(model1a)
summary(model2a)

# full and intercept models with both species and phylogenetic effects (1t and 2t)
model1t <- rma.mv(yi = Final_Z,
                  mods = ~age2+sex+phylum+wild_captive+tissue_lifespan+stressor_cat+exper_correl + measure+method_T+corr+accounting_for_age +
                    measure:age2 + measure:stressor_cat + tissue_lifespan:stressor_cat + age2:stressor_cat,
                  V = mvar, data = data.analysis, test = "t", random = list(~1|species, ~1|species2, ~1|DOI, ~1|ES_id), R = list(species = taxa_cor))
model2t <- rma.mv(yi = Final_Z, V = mvar, data = data.analysis, test = "t", random = list(~1|species, ~1|species2, ~1|DOI, ~1|ES_id), R = list(species = taxa_cor))
summary(model1t)
summary(model2t)

# phylogenetic full and intercept models (species only linked to phylogeny, no non-phylogenetic species effect) (1p and 2p)
model1p <- rma.mv(yi = Final_Z,
                  mods = ~age2+sex+phylum+wild_captive+tissue_lifespan+stressor_cat+exper_correl + measure+method_T+corr+accounting_for_age +
                    measure:age2 + measure:stressor_cat + tissue_lifespan:stressor_cat + age2:stressor_cat,
                  V = mvar, data = data.analysis, test = "t", random = list(~1|species, ~1|DOI, ~1|ES_id), R = list(species = taxa_cor))
model2p <- rma.mv(yi = Final_Z, V = mvar, data = data.analysis, test = "t", random = list(~1|species, ~1|DOI, ~1|ES_id), R = list(species = taxa_cor))
summary(model1p)
summary(model2p)




# random terms tests - intercept model ----
####

# these are the radnom tests mentioned in the paper

# phylogeny vs taxonomic
as.numeric(1 - pchisq(2*(logLik(model2t) - logLik(model2p)), 1))

# phylogeny vs regular
as.numeric(1 - pchisq(2*(logLik(model2p) - logLik(model2a)), 1))

# test the species term
model2a.sp <- rma.mv(yi = Final_Z,
                     V = mvar, data = data.analysis, test = "t", random = list(~1|DOI, ~1|ES_id))
as.numeric(1 - pchisq(2*(logLik(model2a) - logLik(model2a.sp)), 1))

# test the DOI term
model2a.doi <- rma.mv(yi = Final_Z,
                      V = mvar, data = data.analysis, test = "t", random = list(~1|species, ~1|ES_id))
as.numeric(1 - pchisq(2*(logLik(model2a) - logLik(model2a.doi)), 1))



# random terms tests - full models ----
####

# performed to complement the previous batch and check the results hold

# # phylogeny vs taxonomic
as.numeric(1 - pchisq(2*(logLik(model1t) - logLik(model1p)), 1))

# phylogeny vs regular
as.numeric(1 - pchisq(2*(logLik(model1p) - logLik(model1a)), 1))

# test the species term alone
model1a.sp <- rma.mv(yi = Final_Z,
                  mods = ~age2+sex+phylum+wild_captive+tissue_lifespan+stressor_cat+exper_correl + measure+method_T+corr+accounting_for_age +
                    measure:age2 + measure:stressor_cat + tissue_lifespan:stressor_cat + age2:stressor_cat,
                  V = mvar, data = data.analysis, test = "t", random = list(~1|DOI, ~1|ES_id))
as.numeric(1 - pchisq(2*(logLik(model1a) - logLik(model1a.sp)), 1))

# test the DOI term alone
model1a.doi <- rma.mv(yi = Final_Z,
                     mods = ~age2+sex+phylum+wild_captive+tissue_lifespan+stressor_cat+exper_correl + measure+method_T+corr+accounting_for_age +
                       measure:age2 + measure:stressor_cat + tissue_lifespan:stressor_cat + age2:stressor_cat,
                     V = mvar, data = data.analysis, test = "t", random = list(~1|species, ~1|ES_id))
as.numeric(1 - pchisq(2*(logLik(model1a) - logLik(model1a.doi)), 1))



####
# model selection of full model, method = ML ----

model1a.ML <- rma.mv(yi = Final_Z,
                  mods = ~age2+sex+phylum+wild_captive+tissue_lifespan+stressor_cat+exper_correl + measure+method_T+corr+accounting_for_age +
                    measure:age2 + measure:stressor_cat + tissue_lifespan:stressor_cat + age2:stressor_cat,
                  V = mvar, data = data.analysis, test = "t", random = list(~1|species, ~1|ES_id), method = "ML")




# # !!! NOTE - takes v. long time to run !!!
# # the below lines are commented-out to avoid accidentally running the model selection
# # which takes a long time to run
# # to explore the results please load the final_selection.Rdata saved object by running:
# load("190807_final_selection.Rdata")

# # proper model selection code starts below
# model1a.ML.selection <- dredge(model1a.ML, trace = T)
# 
# # save the resulting selection table to save on time if you want to return to results later
# save(model1a.ML.selection, file = "DATEYOURFILEHERE_model1a.ML.selection.Rdata")

# lookup the selection table - data in Table 2
model1a.ML.selection


# model average estimates ----
####

# (prior to moving here please load the .Rdata object referenced above)

model1a.avg <- model.avg(model1a.ML.selection, subset = delta <= 2)
summary(model1a.avg)
importance(model1a.avg)

# the below produces data in Table 3
summary(model1a.avg)
importance(model1a.avg)


# publication bias and related analyses ----
####

# run simple version of the best model - metafor cannot handle here hierarchical models
model1a.best1.uni <- rma(yi = Final_Z,
                        vi = mvar, data = data.analysis, test = "t")
funnel(model1a.best1.uni)
# test for publication bias using Egger's regression
regtest(model1a.best1.uni)

# estimate trim-and-fill profiles using 3 available estimators
trimfill(model1a.best1.uni, estimator = "L0")
trimfill(model1a.best1.uni, estimator = "R0")
trimfill(model1a.best1.uni, estimator = "Q0")

# check for year trend
model1a.yr <- rma.mv(yi = Final_Z,
                        mods = ~ year_publication,
                        V = mvar, data = data, test = "t", random = list(~1|species, ~1|ES_id))
summary(model1a.yr)

# as a supplement - see how effect sizes partition between animal phyla
model1a.ph <- rma.mv(yi = Final_Z,
                        mods = ~ phylum-1,
                        V = mvar, data = data.analysis, test = "t", random = list(~1|species, ~1|ES_id))
summary(model1a.ph)



# oxidative stress vs telomeres analysis ----
####

# load data and create sampling variances for both effect size sets
ox_table <- read.csv("./data/data_stressor exposure_telomeres_oxidative stress.csv", sep = ";", head = T)
ox_table$mev_o <- 1/(ox_table$n_OS-3)
ox_table$mev_t <- 1/(ox_table$n_T-3)

summary(ox_table)

# factorise variable indexing individual rows
ox_table$ES_id <- as.factor(ox_table$ES_id)

# fit overall model to look at general patterns of dependence (note - only approximately valid)
modelox <- rma.mv(yi = ox_table$Z_OS, V = ox_table$mev_o, data = ox_table, random = list(~1|species, ~1|ES_id),
                  mods = ~OS_category_pooled)
summary(modelox)


# look at simple (non-metanalytical) dependence
library(lme4)
library(lmerTest)
modeloxtel <- lmer(Z_T~Z_OS+(1|species), data = ox_table)
summary(modeloxtel)

# check assumptions; slight variance heterogeneity visible
plot(modeloxtel)

# prepare area for plotting and add raw telomere ~ ox_stress data for the first subset of data
par(mfrow=c(1,2))
plot(Z_T~Z_OS, data = ox_table, pch = 19, cex = 1.5, col = "gray",
     xlab = "Fisher's Z: oxidative stress", ylab = "Fisher's Z: telomere length/attrition")

# add iteratively sampling errors for each pair of effect sizes
for(i in 1:nrow(ox_table)) {
  cen.point.x <- ox_table[i, "Z_OS"]
  cen.point.y <- ox_table[i, "Z_T"]
  
  #plot horizontals
  segments(cen.point.x, cen.point.y, cen.point.x-ox_table[i, "mev_o"], cen.point.y)
  segments(cen.point.x, cen.point.y, cen.point.x+ox_table[i, "mev_o"], cen.point.y)
  #plot verticals
  segments(cen.point.x, cen.point.y, cen.point.x, cen.point.y-ox_table[i, "mev_t"])
  segments(cen.point.x, cen.point.y, cen.point.x, cen.point.y+ox_table[i, "mev_t"])
}


# randomisation procedure for OX - TEL analysis (ALL data) ----
####

# form empty columns to store sampled data
ox_table$sampled_O <- numeric(nrow(ox_table))
ox_table$sampled_T <- numeric(nrow(ox_table))

# set simulation size
N <- 5000
# create object to store simulated outcomes
out <- numeric(N)

## !!!NOTE - takes some time to run
for (j in 1:N) {
  # in each of N simulations ...
  for(i in 1:nrow(ox_table)) {
    # ...go through the rows of original data table...
    
    # ... for each effect size sample T and OS effect size from distribution
    # suggested by sampling variances of the give effect-sizes pair ...
    sample.O <- rnorm(n = 1, mean = ox_table[i, "Z_OS"], sd = sqrt(ox_table[i, "mev_o"]))
    sample.T <- rnorm(n = 1, mean = ox_table[i, "Z_T"], sd = sqrt(ox_table[i, "mev_t"]))
    
    # ...store sampled values in the pre-prepared table...
    sample.O -> ox_table[i, "sampled_O"]
    sample.T -> ox_table[i, "sampled_T"]
  }
  
  # ...calculate correlation of sampled effect sizes in the Nth simulation ...
  cor(ox_table$sampled_O, ox_table$sampled_T) -> out[j]
  
  # ...and draw the otcome of the Nth simulation as additional line on the plot
  abline(lm(ox_table$sampled_T ~ ox_table$sampled_O), col = "#DDDDDD")
}

# add the line signifying oroginal relationship between effect sizes
abline(lm(ox_table$Z_T ~ ox_table$Z_OS), col = "orangered", lwd = 2)
axis(2) # fiz the Y axis

# plot the histogram of sampled correlations - with average sampled r added as red line
hist(out, 100, col = "gray", border = NA, xlab = "Correlation of effect sizes", main = NA)
abline(v = mean(out), col = "orangered", lwd = 2)
sort(out)[c(250,4750)] # approx. CI interval for r_sampled
mean(out) # average r_sampled




# randomisation procedure for OX - TEL analysis only for antioxidants levels  ----
####

# (the logic is the same as described above)

ox_table <- read.csv("stress_tel_SZD190807.csv", sep = ";", head = T)
ox_table$mev_o <- 1/(ox_table$n_OS-3)
ox_table$mev_t <- 1/(ox_table$n_T-3)
#ox_table <- gdata::drop.levels(ox_table[1:183,])
summary(ox_table)
names(ox_table)[1] <- "ES_id"
ox_table$ES_id <- as.factor(ox_table$ES_id)

ox_table <- subset(ox_table, OS_category_pooled == "antioxidants")
ox_table$sampled_O <- numeric(nrow(ox_table))
ox_table$sampled_T <- numeric(nrow(ox_table))
N <- 5000
out <- numeric(N)

par(mfrow=c(1,2))
plot(Z_T~Z_OS, data = ox_table, pch = 19, cex = 1.5, col = "gray", bty = "n",
     xlab = "Fisher's Z: oxidative stress", ylab = "Fisher's Z: telomere length/attrition", type = "n")

for(i in 1:nrow(ox_table)) {
  cen.point.x <- ox_table[i, "Z_OS"]
  cen.point.y <- ox_table[i, "Z_T"]
  points(cen.point.x, cen.point.y, pch = 19, cex = 1.5, col = "gray")
  #plot horizontals
  segments(cen.point.x, cen.point.y, cen.point.x - ox_table[i, "mev_o"], cen.point.y)
  segments(cen.point.x, cen.point.y, cen.point.x + ox_table[i, "mev_o"], cen.point.y)
  #plot verticals
  segments(cen.point.x, cen.point.y, cen.point.x, cen.point.y - ox_table[i, "mev_t"])
  segments(cen.point.x, cen.point.y, cen.point.x, cen.point.y + ox_table[i, "mev_t"])
}

for (j in 1:N) {
  for(i in 1:nrow(ox_table)) {
    sample.O <- rnorm(n = 1, mean = ox_table[i, "Z_OS"], sd = sqrt(ox_table[i, "mev_o"]))
    sample.T <- rnorm(n = 1, mean = ox_table[i, "Z_T"], sd = sqrt(ox_table[i, "mev_t"]))
    
    sample.O -> ox_table[i, "sampled_O"]
    sample.T -> ox_table[i, "sampled_T"]
  }
  cor(ox_table$sampled_O, ox_table$sampled_T) -> out[j]
  abline(lm(ox_table$sampled_T ~ ox_table$sampled_O), col = "#DDDDDD")
}
abline(lm(ox_table$Z_T ~ ox_table$Z_OS), col = "orangered", lwd = 2)
axis(2)

hist(out, 100, col = "gray", border = NA, xlab = "Correlation of effect sizes", main = NA)
abline(v = mean(out), col = "orangered", lwd = 2)

sort(out)[c(250,4750)]
mean(out)



# randomisation procedure for OX - TEL analysis only for level of oxidative damages  ----
# the logic is the same as in the first example
ox_table <- read.csv("stress_tel_SZD190807.csv", sep = ";", head = T)
ox_table$mev_o <- 1/(ox_table$n_OS-3)
ox_table$mev_t <- 1/(ox_table$n_T-3)
#ox_table <- gdata::drop.levels(ox_table[1:183,])
summary(ox_table)
names(ox_table)[1] <- "ES_id"
ox_table$ES_id <- as.factor(ox_table$ES_id)

ox_table <- subset(ox_table, OS_category_pooled == "oxidative damage")
ox_table$sampled_O <- numeric(nrow(ox_table))
ox_table$sampled_T <- numeric(nrow(ox_table))
N <- 5000
out <- numeric(N)

par(mfrow=c(1,2))
plot(Z_T~Z_OS, data = ox_table, pch = 19, cex = 1.5, col = "gray", bty = "n",
     xlab = "Fisher's Z: oxidative stress", ylab = "Fisher's Z: telomere length/attrition", type = "n")

for(i in 1:nrow(ox_table)) {
  cen.point.x <- ox_table[i, "Z_OS"]
  cen.point.y <- ox_table[i, "Z_T"]
  points(cen.point.x, cen.point.y, pch = 19, cex = 1.5, col = "gray")
  #plot horizontals
  segments(cen.point.x, cen.point.y, cen.point.x - ox_table[i, "mev_o"], cen.point.y)
  segments(cen.point.x, cen.point.y, cen.point.x + ox_table[i, "mev_o"], cen.point.y)
  #plot verticals
  segments(cen.point.x, cen.point.y, cen.point.x, cen.point.y - ox_table[i, "mev_t"])
  segments(cen.point.x, cen.point.y, cen.point.x, cen.point.y + ox_table[i, "mev_t"])
}

## !!!NOTE - takes some time to run
for (j in 1:N) {
  for(i in 1:nrow(ox_table)) {
    sample.O <- rnorm(n = 1, mean = ox_table[i, "Z_OS"], sd = sqrt(ox_table[i, "mev_o"]))
    sample.T <- rnorm(n = 1, mean = ox_table[i, "Z_T"], sd = sqrt(ox_table[i, "mev_t"]))
    
    sample.O -> ox_table[i, "sampled_O"]
    sample.T -> ox_table[i, "sampled_T"]
  }
  cor(ox_table$sampled_O, ox_table$sampled_T) -> out[j]
  abline(lm(ox_table$sampled_T ~ ox_table$sampled_O), col = "#DDDDDD")
}
abline(lm(ox_table$Z_T ~ ox_table$Z_OS), col = "orangered", lwd = 2)
axis(2)

hist(out, 100, col = "gray", border = NA, xlab = "Correlation of effect sizes", main = NA)
abline(v = mean(out), col = "orangered", lwd = 2)

sort(out)[c(250,4750)]
mean(out)
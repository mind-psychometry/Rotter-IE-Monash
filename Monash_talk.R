# install.packages("devtools")
# devtools::install_github("masiraji/tabledown",force = TRUE)
# install.packages("pacman")
# pacman::p_load(tidyverse,psych,mirt,lavaan,gtsummary, VIM)
# devtools::install_github("moshagen/semPower")

# Load the Libraries
library(tidyverse)
library(psych)
library(mirt)
library(lavaan)
library(gtsummary)
library(VIM)
library(tabledown)
library(semPower)
library(labelled)
library(MOTE)
library(correlation)
library(qgraph)
library(semTools)
library(semTools)
# Sample Size Estimation


## df calculation
df <-  (p(p+1)/2)-q

lavmodel <-  "LOC =~ item18 + item25 +item11 + item15 + item16 + item10 +item9 + 
item23 + item13 + item5 + item28 + item4" 
df <- semPower::semPower.getDf(lavmodel)

ap <- semPower::semPower.aPriori(effect = .05, effect.measure = 'RMSEA', alpha = .05, 
                       power = .80, df = df) 
summary(ap)

# ph <- semPower.postHoc(effect = .05, effect.measure = 'RMSEA', alpha = .05, N = 178, 
#                  df = 54) 
# summary(ph)


#Load the data
data_0 <- readRDS("WRotter.rds")
data <- labelled::remove_labels(data_0)
##EFA data
Rotter.EFA <- data[1:300,11:33]
##EFA Descriptivs
efa.descriptives <- data[1:300,c(3,4,5,6,7,8)]
##CFA Data
Rotter.CFA <- data[301:478, 11:33]
##CFA Descriptives
cfa.descriptives <- data[301:478,c(3,4,5,6,7,8)]
##IRT data
WRotter.IRT <- data[,11:33]
##IRT Descriptives
irt.descriptives <- data[,c(3,4,5,6,7,8)]


## Data for descriptives (Publication)
Descriptives.gender <- data %>% 
  dplyr::select("Gender", "Age", "Education Years", "Social Stance","Marital Status", 
                "sample")



#Check Missing Values
missing.EFA <- VIM::aggr(Rotter.EFA, plot =T)
missing.CFA<- VIM::aggr(Rotter.CFA, plot =T)
missing.IRT <- VIM::aggr(WRotter.IRT, plot =T)

#Descriptives
#Base form
data_summary <- psych::describeBy(Descriptives.gender, group ="sample")

# Beautiful Table
desc_table <-  gtsummary::tbl_summary(Descriptives.gender,
                                      by = sample,
                                      statistic = list(
                                        all_continuous() ~ "{mean} ({sd})",
                                        all_categorical() ~ "{n} ({p}%)"),
                                      type = list(`Education Years` ~ 'continuous'),
                                      digits = all_continuous() ~ 2,
                                      missing = "no") %>%add_overall() %>% add_p()
#Item Summary
des.tab <- tabledown::des.tab(Rotter.EFA)

#multivariate normality
mardia <- mardia(Rotter.EFA, na.rm = T, plot =T)

#Correlation

##Correlation matrix
correlations <- psych::tetrachoric(Rotter.EFA,smooth=TRUE)
correlations_con <- as.data.frame(cor(Rotter.EFA))
upper <- correlations$rho[upper.tri(correlations$rho, diag = F)]

min.cor <- MOTE::apa(min((upper)),2,F) #minimum cor.coefficient of the matrix
max.cor <- MOTE::apa(max((upper)),2,F) # max. correlation coefficient of the matrix
##Quality of the matrix
bartlet <- psych::cortest.bartlett(Rotter.EFA, n =300)

##Correlation Plot
#qgraph
qgraph::qgraph(correlations$rho, cut =.30, details =T, posCol = "green", negCol = "orangered3", 
               labels = names(Rotter.EFA))

#Psych
psych::corPlot(correlations$rho, diag =F, scale=T,
               upper =F, numbers = T, cex.axis =.5, stars=T, cex=.8 )



#ggcorplot
p.mat <- correlation::cor_to_p(correlations$rho, 300, method = "tetrachoric")

corplot <- ggcorrplot::ggcorrplot(correlations$rho, p.mat=p.mat$p, insig = "pch", hc.order = FALSE, outline.col = "white", type = "lower", legend.title = "Correlation",sig.level = 0.05, pch.cex=1.8,
                                  # insig = "blank",
                                  ggtheme = ggplot2::theme_minimal(),tl.srt = 90,colors = c("#E64B35FF", "white", "#3C5488FF") )+
  theme( panel.grid.major = element_blank(),
         axis.text.x = element_text(size = 8),
         axis.text.y = element_text(size = 8),
         panel.border = element_rect(color = "black",
                                     fill = NA,
                                     size = 1))

#CTT Item Analysis
Item_analysis <- psych::alpha(Rotter.EFA)
R.cor.sum <- summary(Item_analysis$item.stats$r.cor)


#Factor Identification
#Scree
Scree <- scree(correlations$rho,factors=TRUE,pc=TRUE,main="Scree Plot",
               hline=NULL,add=F)
#HUll method
EFA.MRFA::hullEFA(Rotter.EFA,extr = "ULS", index_hull = "CAF", display = TRUE, 
                  graph = T,
                  details = TRUE)

#Factor analysis
fa.2F.1 <- fa(r=correlations$rho, nfactors = 2, fm= "wls", rotate = "varimax",
              residuals = TRUE, SMC = TRUE, n.obs =300)
print(fa.2F.1, cut = .3, digits = 3, sort = T)

#Refitting the model
reduced.model.2F.1 <- dplyr::select(Rotter.EFA, -item22, -item6,-item29,-item25, -item9)
correlations.redu.2F.1 <-  tetrachoric(reduced.model.2F.1)
fa.2F.2 <- fa(r=correlations.redu.2F.1$rho, nfactors = 2, fm= "wls", rotate = "varimax",
              residuals = TRUE, SMC = TRUE, n.obs =300)
print(fa.2F.2, cut = .3, digits = 3, sort = TRUE)

Factor_analysis_tab <- FAc <-tabledown::fac.tab( fa.2F.2, .3, complexity = F)

# ---------------------------5mins break-----------------------------------------


Rotter.one <- "LOC =~ item18 + item25 +item11 + item15 + item16 + item10 +item9 + 
item23 + item13 + item5 + item28 + item4"

Rotter.one.fit <- cfa(Rotter.one,  Rotter.CFA,  mimic = "Mplus",
                      ordered = names(Rotter.CFA),estimator = "WLSMV")

One.factor.sum <- summary(Rotter.one.fit, fit.measures =TRUE,standardized = TRUE,
                          rsq =TRUE)

One.factor.Cfa <- fitmeasures (Rotter.one.fit,c("gfi", "agfi", "nfi","rfi", 
                                                "cfi.scaled","tli.scaled",
                                                "rmsea", 
                                                "srmr","aic"))

One.factor.reliability <- semTools::reliability(Rotter.one.fit)
One.factor.reliability[5]
modfit.Cor.one <- modindices(Rotter.one.fit, sort. = TRUE) 
modfit.Cor.one[modfit.Cor.one$mi>3.84,]


Rotter.one.2 <- "LOC =~ item18 + item25 +item11 + item15 + item16 + item10 +item9 + item13 + item5 + item28 + item4"

Rotter.one.fit.2 <- cfa(Rotter.one.2,  Rotter.CFA,  mimic = "Mplus",
                        ordered = names(Rotter.CFA),estimator = "WLSMV")

One.factor.Cfa.mod.sum <- summary(Rotter.one.fit.2, fit.measures =TRUE,standardized = TRUE, rsq =TRUE)

One.factor.Cfa.mod.fit <- fitmeasures (Rotter.one.fit.2,c("gfi", "agfi", "nfi","rfi", 
                                                          "cfi.scaled","tli.scaled",
                                                          "rmsea", 
                                                          "srmr","aic"))

One.factor.reliability.mod <- semTools::reliability(Rotter.one.fit.2)

CFA_plot1 <- semPlot::semPaths (Rotter.one.fit.2 , 
                                  what= "std",
                                  intercepts = F,
                                  style ="OpenMx",
                                  theme = "colorblind",
                                  nCharNodes = 0,
                                  reorder =T,
                                  rotation =2,
                                  layout ="tree",
                                  cardinal = T,
                                  curvePivot =T,
                                  sizeMan =6,#items length
                                  sizeMan2 = 2,#items
                                  sizeLat = 12,#factors
                                  thresholds = FALSE,
                                  equalizeManifests =F,
                                  fade = FALSE)


#IRT
Rotter.one.IRT <- WRotter.IRT %>% 
  dplyr::select(item18, item25, item11, item15, item16, item10, item9, item13,  
                item5,  item28, item4)


unidimensional_fit <- mirt::mirt(Rotter.one.IRT, model = 1, 
                           itemtype = '2PL', 
                           SE = TRUE, 
                           Se.type = 'MHRM',
                           technical = list(NCYCLES = 10000))
##Model Fit
unidim_fit_indices <- mirt::M2(unidimensional_fit)

##Item fit
unidimensional_itemfit <- mirt::itemfit(unidimensional_fit,
                                  fit_stats = c("S_X2", "G2","Zh", "infit"))

unidimensional_params <- mirt::coef(unidimensional_fit, IRTpars = TRUE, simplify = TRUE, rawug = FALSE) 

unidimensional_items <- data.frame(unidimensional_params$items)

##Person-fit
unidimensional.personfit <- mirt::personfit(unidimensional_fit) 

unidimensional.personfit_model_misfits <- subset(unidimensional.personfit, Zh < -2)

rownames(unidimensional.personfit_model_misfits)

nrow(unidimensional.personfit_model_misfits)

personfit <- ggplot(unidimensional.personfit, aes(x=Zh)) + geom_histogram(binwidth=.5,col=I("black"),fill="#00A08799")+geom_vline(xintercept = -2)+
  labs( y="Number of Participants", x = "Zh Statistics")


#IRT Curves
## Conventionals
item_info <- mirt::plot(unidimensional_fit, type= "infotrace")
test_info <- mirt::plot(unidimensional_fit, type= "info")
occ_info <- mirt::plot(unidimensional_fit, type= "trace")



tic <- tabledown::ggtestinfo(Rotter.one.IRT, unidimensional_fit)+
  geom_area()+aes(fill =  "#E64B3599")+theme(legend.position = "none")+
  scale_fill_identity()+xlab(expression(theta))+ylab(expression(I(theta)))

item_info <- tabledown::ggiteminfo(unidimensional_fit, 1,6)

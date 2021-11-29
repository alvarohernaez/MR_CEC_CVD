rm(list=ls())

# Load necessary packages #

library(devtools)
library(googleAuthR)
library(MendelianRandomization)
library(mr.raps)
library(meta)
library(MRPRESSO)
library(MRInstruments)
library(MRMix)
library(RadialMR)
library(ieugwasr)
library(TwoSampleMR)


# Data aspect functions #

guapa<-function(x)
{
  redondeo<-ifelse(abs(x)<0.0001,signif(x,1),
                   ifelse(abs(x)<0.001,signif(x,1),
                          ifelse(abs(x)<0.1,round(x,3),
                                 ifelse(abs(x)<1,round(x,2),signif(x,3)))))
  return(redondeo)
}

ic_guapa<-function(x,y,z)
{
  ic<-paste(x," [",y,"; ",z,"]",sep="")
  return(ic)
}

pval_guapa<-function(x)
{
  pval<-ifelse(x<0.00001,"<0.00001",
               ifelse(x<0.001,"<0.001",round(x,3)))
  return(pval)
}

mean_ic_guapa <- function(x, na.rm=FALSE) 
{
  if (na.rm) x <- na.omit(x)
  se<-sqrt(var(x)/length(x))
  z<-qnorm(1-0.05/2)
  media<-mean(x)
  ic95a<-guapa(media-(z*se))
  ic95b<-guapa(media+(z*se))
  media<-guapa(media)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

mean_sd_guapa <- function(x) 
{
  media<-guapa(mean(x, na.rm=TRUE))
  sd<-guapa(sd(x, na.rm=TRUE))
  end<-paste(media," (",sd,")",sep="")
  return(end)
}

beta_se_ic_guapa <- function(x, y) 
{
  z<-qnorm(1-0.05/2)
  ic95a<-guapa(x-(z*y))
  ic95b<-guapa(x+(z*y))
  media<-guapa(x)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa2 <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-round(exp(x),5)
  ic95a<-round(exp(x-(z*y)),5)
  ic95b<-round(exp(x+(z*y)),5)
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

pval_ic_guapa <- function(ratio, lo, up) 
{
  z<-qnorm(1-0.05/2)
  x<-log(ratio)/((log(up)-log(lo))/(2*z))
  pval<-pval_guapa(exp(-0.717*x-(0.416*x^2)))
  return(pval)
}

header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))] }


dir.create("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Data/cec_cvd")
dir.create("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Data/cec_cvd/results")
setwd("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Data/cec_cvd")


######################################################################
### EXPOSURE: CEC - All: https://pubmed.ncbi.nlm.nih.gov/30369316/ ###
######################################################################

cec<-read.csv2("./cec_lowkam_all_2018.csv",
                header=TRUE,sep=";",dec=".")
cec$Phenotype<-c("Cholesterol efflux capacity")
cec$units<-c("units")
cec$id<-c("Low-Kam C, J Am Heart Assoc, 2018")
cec$samplesize<-5293
cec<-rename.vars(cec,
                  from=c("ï..rsid","ref_allele","m1_b","m1_se","m1_pval"),
                  to=c("SNP","other_allele","beta","se","pval"))

cec_mod1<-cec[,c("Phenotype","SNP","beta","se","eaf","effect_allele","other_allele","pval",
                "units","samplesize","id")]
cec_mod1<-subset2(cec_mod1,"cec_mod1$pval<5e-08")
write.table(cec_mod1,"./cec_mod1.txt",sep="\t",row.names=FALSE, quote=FALSE)
cec_mod1_harm<-read_exposure_data(filename = "./cec_mod1.txt",
                                 clump = FALSE,
                                 sep="\t",
                                 phenotype_col = "Phenotype",
                                 snp_col = "SNP",
                                 beta_col = "beta",
                                 se_col = "se",
                                 eaf_col = "eaf",
                                 effect_allele_col = "effect_allele",
                                 other_allele_col = "other_allele",
                                 pval_col = "pval",
                                 units_col = "units",
                                 samplesize_col = "samplesize",
                                 id_col = "id",
                                 min_pval = 1e-200,
                                 log_pval = FALSE)
write.table(cec_mod1_harm,"./cec_mod1_harm.txt",sep="\t",row.names=FALSE, quote=FALSE)
save(cec_mod1_harm,file="./cec_mod1_harm.RData")


### OUTCOME: CHD UKBB + CARDIoGRAMplusC4D ('ebi-a-GCST005195'): https://pubmed.ncbi.nlm.nih.gov/29212778/ ###
### OUTCOME: Stroke meta-analysis ('ebi-a-GCST005838'): https://pubmed.ncbi.nlm.nih.gov/29531354/ ###
### OUTCOME: Only ischemic stroke meta-analysis ('ebi-a-GCST005843'): https://pubmed.ncbi.nlm.nih.gov/29531354/ ###
### OUTCOME: Heart failure meta-analysis ('ebi-a-GCST009541'): https://pubmed.ncbi.nlm.nih.gov/31919418/ ###


chd_harm<-extract_outcome_data(snps=cec_mod1_harm$SNP,outcomes='ebi-a-GCST005195')
dat<-harmonise_data(cec_mod1_harm, chd_harm, action = 2)
save(dat,file="./cec_mod1_chd.RData")

str_harm<-extract_outcome_data(snps=cec_mod1_harm$SNP,outcomes='ebi-a-GCST005838')
dat<-harmonise_data(cec_mod1_harm, str_harm, action = 2)
save(dat,file="./cec_mod1_str.RData")

is_harm<-extract_outcome_data(snps=cec_mod1_harm$SNP,outcomes='ebi-a-GCST005843')
dat<-harmonise_data(cec_mod1_harm, is_harm, action = 2)
save(dat,file="./cec_mod1_is.RData")

hf_harm<-extract_outcome_data(snps=cec_mod1_harm$SNP,outcomes='ebi-a-GCST009541')
dat<-harmonise_data(cec_mod1_harm, hf_harm, action = 2)
save(dat,file="./cec_mod1_hf.RData")


#########################
### STEIGER FILTERING ###
#########################

### HDL-C (Klarin et al., Nat Genet, 2018) ###

hdlc<-as.data.frame(read.delim("./hdlc_summary.txt",
                               header=TRUE,sep="\t"))
names(hdlc)<-tolower(names(hdlc))
hdlc$a1<-toupper(hdlc$a1)
hdlc$a2<-toupper(hdlc$a2)
hdlc<-rename.vars(hdlc,
                  from=c("rsid","a1","a2","p.value","n","freq.a1.1000g.eur"),
                  to=c("SNP","effect_allele","other_allele","pval","samplesize","eaf"))
hdlc<-hdlc[,c("SNP","effect_allele","other_allele","beta","se","eaf","pval","samplesize")]
hdlc$Units<-c("Units")
hdlc$Phenotype<-c("hdlc")
write.table(hdlc,"./hdlc_summary_clean.txt",sep="\t",row.names=FALSE, quote=FALSE)

hdlc_harm <- read_outcome_data(
  snps = cec_mod1_harm$SNP,
  filename = "./hdlc_summary_clean.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

cec_mod1_hdlc<-harmonise_data(cec_mod1_harm, hdlc_harm, action = 2)
cec_mod1_hdlc_steiger<-steiger_filtering(cec_mod1_hdlc)
cec_mod1_hdlc_steiger<-subset2(cec_mod1_hdlc_steiger,"cec_mod1_hdlc_steiger$steiger_dir==TRUE")
dim(cec_mod1_hdlc)[1]-dim(cec_mod1_hdlc_steiger)[1]
cec_mod1_hdlc_steiger$SNP


# TRIGLYCERIDES (Klarin et al., Nat Genet, 2018) #

tg<-as.data.frame(read.delim("./tg_summary.txt",
                             header=TRUE,sep="\t"))
names(tg)<-tolower(names(tg))
tg$a1<-toupper(tg$a1)
tg$a2<-toupper(tg$a2)
tg<-rename.vars(tg,
                from=c("rsid","a1","a2","p.value","n","freq.a1.1000g.eur"),
                to=c("SNP","effect_allele","other_allele","pval","samplesize","eaf"))
tg<-tg[,c("SNP","effect_allele","other_allele","beta","se","eaf","pval","samplesize")]
tg$Units<-c("Units")
tg$Phenotype<-c("tg")
write.table(tg,"./tg_summary_clean.txt",sep="\t",row.names=FALSE, quote=FALSE)

tg_harm <- read_outcome_data(
  snps = cec_mod1_harm$SNP,
  filename = "./tg_summary_clean.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

cec_mod1_tg<-harmonise_data(cec_mod1_harm, tg_harm, action = 2)
cec_mod1_tg_steiger<-steiger_filtering(cec_mod1_tg)
cec_mod1_tg_steiger<-subset2(cec_mod1_tg_steiger,"cec_mod1_tg_steiger$steiger_dir==TRUE")
dim(cec_mod1_tg)[1]-dim(cec_mod1_tg_steiger)[1]
cec_mod1_tg_steiger$SNP


# LDL-C (Klarin et al., Nat Genet, 2018) #

ldlc<-as.data.frame(read.delim("./ldlc_summary.txt",
                               header=TRUE,sep="\t"))
names(ldlc)<-tolower(names(ldlc))
ldlc$a1<-toupper(ldlc$a1)
ldlc$a2<-toupper(ldlc$a2)
ldlc<-rename.vars(ldlc,
                  from=c("rsid","a1","a2","p.value","n","freq.a1.1000g.eur"),
                  to=c("SNP","effect_allele","other_allele","pval","samplesize","eaf"))
ldlc<-ldlc[,c("SNP","effect_allele","other_allele","beta","se","eaf","pval","samplesize")]
ldlc$Units<-c("Units")
ldlc$Phenotype<-c("ldlc")
write.table(ldlc,"./ldlc_summary_clean.txt",sep="\t",row.names=FALSE, quote=FALSE)

ldlc_harm <- read_outcome_data(
  snps = cec_mod1_harm$SNP,
  filename = "./ldlc_summary_clean.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

cec_mod1_ldlc<-harmonise_data(cec_mod1_harm, ldlc_harm, action = 2)
cec_mod1_ldlc_steiger<-steiger_filtering(cec_mod1_ldlc)
cec_mod1_ldlc_steiger<-subset2(cec_mod1_ldlc_steiger,"cec_mod1_ldlc_steiger$steiger_dir==TRUE")
dim(cec_mod1_ldlc)[1]-dim(cec_mod1_ldlc_steiger)[1]
cec_mod1_ldlc_steiger$SNP


# SMALL HDL PARTICLES #

shdlp_harm <- extract_outcome_data(
  snps = cec_mod1_harm$SNP,
  outcomes = 'met-c-922'
)

cec_mod1_shdlp<-harmonise_data(cec_mod1_harm, shdlp_harm, action = 2)
cec_mod1_shdlp_steiger<-steiger_filtering(cec_mod1_shdlp)
cec_mod1_shdlp_steiger<-subset2(cec_mod1_shdlp_steiger,"cec_mod1_shdlp_steiger$steiger_dir==TRUE")
dim(cec_mod1_shdlp)[1]-dim(cec_mod1_shdlp_steiger)[1]
cec_mod1_shdlp_steiger$SNP


# MEDIUM-SIZED HDL PARTICLES #

mhdlp_harm <- extract_outcome_data(
  snps = cec_mod1_harm$SNP,
  outcomes = 'met-c-902'
)

cec_mod1_mhdlp<-harmonise_data(cec_mod1_harm, mhdlp_harm, action = 2)
cec_mod1_mhdlp_steiger<-steiger_filtering(cec_mod1_mhdlp)
cec_mod1_mhdlp_steiger<-subset2(cec_mod1_mhdlp_steiger,"cec_mod1_mhdlp_steiger$steiger_dir==TRUE")
dim(cec_mod1_mhdlp)[1]-dim(cec_mod1_mhdlp_steiger)[1]
cec_mod1_mhdlp_steiger$SNP


# LARGE HDL PARTICLES #

lhdlp_harm <- extract_outcome_data(
  snps = cec_mod1_harm$SNP,
  outcomes = 'met-c-878'
)

cec_mod1_lhdlp<-harmonise_data(cec_mod1_harm, lhdlp_harm, action = 2)
cec_mod1_lhdlp_steiger<-steiger_filtering(cec_mod1_lhdlp)
cec_mod1_lhdlp_steiger<-subset2(cec_mod1_lhdlp_steiger,"cec_mod1_lhdlp_steiger$steiger_dir==TRUE")
dim(cec_mod1_lhdlp)[1]-dim(cec_mod1_lhdlp_steiger)[1]
cec_mod1_lhdlp_steiger$SNP


# VERY-LARGE HDL PARTICLES #

xlhdlp_harm <- extract_outcome_data(
  snps = cec_mod1_harm$SNP,
  outcomes = 'met-c-946'
)

cec_mod1_xlhdlp<-harmonise_data(cec_mod1_harm, xlhdlp_harm, action = 2)
cec_mod1_xlhdlp_steiger<-steiger_filtering(cec_mod1_xlhdlp)
cec_mod1_xlhdlp_steiger<-subset2(cec_mod1_xlhdlp_steiger,"cec_mod1_xlhdlp_steiger$steiger_dir==TRUE")
dim(cec_mod1_xlhdlp)[1]-dim(cec_mod1_xlhdlp_steiger)[1]
cec_mod1_xlhdlp_steiger$SNP


snp_steiger<-intersect(cec_mod1_hdlc_steiger$SNP,cec_mod1_ldlc_steiger$SNP)
snp_steiger01<-intersect(snp_steiger,cec_mod1_tg_steiger$SNP)
snp_steiger02<-intersect(snp_steiger01,cec_mod1_shdlp_steiger$SNP)
snp_steiger02<-intersect(snp_steiger02,cec_mod1_mhdlp_steiger$SNP)
snp_steiger02<-intersect(snp_steiger02,cec_mod1_lhdlp_steiger$SNP)
snp_steiger02<-intersect(snp_steiger02,cec_mod1_xlhdlp_steiger$SNP)


### HARMONIZED DATASETS ###

chd_harm<-extract_outcome_data(snps=snp_steiger01,outcomes='ebi-a-GCST005195')
dat<-harmonise_data(cec_mod1_harm, chd_harm, action = 2)
save(dat,file="./cec_mod1_chd_steiger01.RData")

str_harm<-extract_outcome_data(snps=snp_steiger01,outcomes='ebi-a-GCST005838')
dat<-harmonise_data(cec_mod1_harm, str_harm, action = 2)
save(dat,file="./cec_mod1_str_steiger01.RData")

is_harm<-extract_outcome_data(snps=snp_steiger01,outcomes='ebi-a-GCST005843')
dat<-harmonise_data(cec_mod1_harm, is_harm, action = 2)
save(dat,file="./cec_mod1_is_steiger01.RData")

hf_harm<-extract_outcome_data(snps=snp_steiger01,outcomes='ebi-a-GCST009541')
dat<-harmonise_data(cec_mod1_harm, hf_harm, action = 2)
save(dat,file="./cec_mod1_hf_steiger01.RData")


chd_harm<-extract_outcome_data(snps=snp_steiger02,outcomes='ebi-a-GCST005195')
dat<-harmonise_data(cec_mod1_harm, chd_harm, action = 2)
save(dat,file="./cec_mod1_chd_steiger02.RData")

str_harm<-extract_outcome_data(snps=snp_steiger02,outcomes='ebi-a-GCST005838')
dat<-harmonise_data(cec_mod1_harm, str_harm, action = 2)
save(dat,file="./cec_mod1_str_steiger02.RData")

is_harm<-extract_outcome_data(snps=snp_steiger02,outcomes='ebi-a-GCST005843')
dat<-harmonise_data(cec_mod1_harm, is_harm, action = 2)
save(dat,file="./cec_mod1_is_steiger02.RData")

hf_harm<-extract_outcome_data(snps=snp_steiger02,outcomes='ebi-a-GCST009541')
dat<-harmonise_data(cec_mod1_harm, hf_harm, action = 2)
save(dat,file="./cec_mod1_hf_steiger02.RData")


#####################
### TWO-SAMPLE MR ###
#####################

vars01<-c("cec_mod1","cec_mod1","cec_mod1","cec_mod1",
          "cec_mod1","cec_mod1","cec_mod1","cec_mod1",
          "cec_mod1","cec_mod1","cec_mod1","cec_mod1")
vars02<-c("chd","str","is","hf",
          "chd_steiger01","str_steiger01","is_steiger01","hf_steiger01",
          "chd_steiger02","str_steiger02","is_steiger02","hf_steiger02")

z<-qnorm(1-0.05/2)
tab<-NULL
for(i in 1:length(vars01))
  
{
  namedat<-paste("./",vars01[i],"_",vars02[i],".RData",sep="")
  load(namedat)
  
  mr_results<-mr(dat,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
  mr_results$pval<-pval_guapa(mr_results$pval)
  mr_results$beta<-paste(guapa(mr_results$b)," (",guapa(mr_results$se),")",sep="")
  mr_results$or<-risk_se_ic_guapa(mr_results$b,mr_results$se)
  mr_raps<-mr.raps(b_exp=dat$beta.exposure,b_out=dat$beta.outcome,se_exp=dat$se.exposure,se_out=dat$se.outcome,diagnosis=FALSE)
  mr_raps_beta<-paste(guapa(mr_raps$beta.hat)," (",guapa(mr_raps$beta.se),")",sep="")
  mr_raps_or<-risk_se_ic_guapa(mr_raps$beta.hat,mr_raps$beta.se)
  mr_raps_pval<-pval_guapa(mr_raps$beta.p.value)
  mr_ivw_forest01<-round(exp(mr_results$b[1]),5)
  mr_ivw_forest02<-round(exp(mr_results$b[1]-(z*mr_results$se[1])),5)
  mr_ivw_forest03<-round(exp(mr_results$b[1]+(z*mr_results$se[1])),5)
  mr_egger_forest01<-round(exp(mr_results$b[2]),5)
  mr_egger_forest02<-round(exp(mr_results$b[2]-(z*mr_results$se[2])),5)
  mr_egger_forest03<-round(exp(mr_results$b[2]+(z*mr_results$se[2])),5)
  mr_wme_forest01<-round(exp(mr_results$b[3]),5)
  mr_wme_forest02<-round(exp(mr_results$b[3]-(z*mr_results$se[3])),5)
  mr_wme_forest03<-round(exp(mr_results$b[3]+(z*mr_results$se[3])),5)
  mr_wmo_forest01<-round(exp(mr_results$b[4]),5)
  mr_wmo_forest02<-round(exp(mr_results$b[4]-(z*mr_results$se[4])),5)
  mr_wmo_forest03<-round(exp(mr_results$b[4]+(z*mr_results$se[4])),5)
  mr_raps_forest01<-round(exp(mr_raps$beta.hat),5)
  mr_raps_forest02<-round(exp(mr_raps$beta.hat-(z*mr_raps$beta.se)),5)
  mr_raps_forest03<-round(exp(mr_raps$beta.hat+(z*mr_raps$beta.se)),5)
  q1<-paste(round(mr_heterogeneity(dat)$Q[2],2)," (P=",pval_guapa(mr_heterogeneity(dat)$Q_pval[2]),")",sep="")
  q2<-paste(round(mr_heterogeneity(dat)$Q[1],2)," (P=",pval_guapa(mr_heterogeneity(dat)$Q_pval[1]),")",sep="")
  plei<-pval_guapa(mr_pleiotropy_test(dat)$pval)
  
  tab<-rbind(tab,cbind(mr_results$nsnp[1],
                       mr_results$beta[1],mr_results$or[1],mr_results$pval[1],
                       mr_results$beta[2],mr_results$or[2],mr_results$pval[2],
                       mr_results$beta[3],mr_results$or[3],mr_results$pval[3],
                       mr_results$beta[4],mr_results$or[4],mr_results$pval[4],
                       mr_raps_beta,mr_raps_or,mr_raps_pval,plei,q1,q2,
                       mr_ivw_forest01,mr_ivw_forest02,mr_ivw_forest03,
                       mr_egger_forest01,mr_egger_forest02,mr_egger_forest03,
                       mr_wme_forest01,mr_wme_forest02,mr_wme_forest03,
                       mr_wmo_forest01,mr_wmo_forest02,mr_wmo_forest03,
                       mr_raps_forest01,mr_raps_forest02,mr_raps_forest03))
}

colnames(tab)<-c("SNPs","IVW_b","IVW_or","IVW_p","Egger_b","Egger_or","Egger_p",
                 "WMe_b","WMe_or","Wme_p","WMo_b","WMo_or","WMo_p",
                 "RAPS_b","RAPS_or","RAPS_p","Egger_pleio","Cochran_Q","Rucker_Q",
                 "IVW_beta","IVW_lo","IVW_hi",
                 "Egger_beta","Egger_lo","Egger_hi",
                 "W-Median_beta","W-Median_lo","W-Median_hi",
                 "W-Mode_beta","W-Mode_lo","W-Mode_hi",
                 "RAPS_beta","RAPS_lo","RAPS_hi")
rownames(tab)<-paste(vars01,"_",vars02,sep="")
write.table(tab,file="./results/cec_cvd_mr.csv",sep=";",col.names=NA)


####################
### FOREST PLOTS ###
####################

library(forestplot)
tab<-read.csv2("./results/cec_cvd_mr.csv",header=TRUE,sep=";",dec=".")

level00<-c(" ")
level01<-c(" Coronary heart disease")
level02<-c("      Main result (MR-RAPS method)")
level03<-paste("      MR-Egger (intercept test: p-value = ",tab[1,18],")",sep="")
level04<-c("      Weighted median method")
level05<-c("      Weighted mode method")
level06<-c("      MR-RAPS: SNPs unrelated to HDL-C, LDL-C and TGs")
level07<-c("      MR-RAPS: SNPs unrelated to lipids and HDL particles")
level08<-c(" Stroke")
level09<-c("      Main result (MR-RAPS method)")
level10<-paste("      MR-Egger (intercept test: p-value = ",tab[2,18],")",sep="")
level11<-c("      Weighted median method")
level12<-c("      Weighted mode method")
level15<-c(" Ischemic stroke only")
level16<-c("      Main result (MR-RAPS method)")
level17<-paste("      MR-Egger (intercept test: p-value = ",tab[3,18],")",sep="")
level18<-c("      Weighted median method")
level19<-c("      Weighted mode method")
level22<-c(" Heart failure")
level23<-c("      Main result (MR-RAPS method)")
level24<-paste("      MR-Egger (intercept test: p-value = ",tab[4,18],")",sep="")
level25<-c("      Weighted median method")
level26<-c("      Weighted mode method")

level00<-c(" ")
level01<-c(" ")
level02<-paste(tab[1,16],sep="")
level03<-paste(tab[1,7],sep="")
level04<-paste(tab[1,10],sep="")
level05<-paste(tab[1,13],sep="")
level06<-paste(tab[5,16],sep="")
level07<-paste(tab[9,16],sep="")
level08<-c(" ")
level09<-paste(tab[2,16],sep="")
level10<-paste(tab[2,7],sep="")
level11<-paste(tab[2,10],sep="")
level12<-paste(tab[2,13],sep="")
level15<-c(" ")
level16<-paste(tab[3,16],sep="")
level17<-paste(tab[3,7],sep="")
level18<-paste(tab[3,10],sep="")
level19<-paste(tab[3,13],sep="")
level22<-c(" ")
level23<-paste(tab[4,16],sep="")
level24<-paste(tab[4,7],sep="")
level25<-paste(tab[4,10],sep="")
level26<-paste(tab[4,13],sep="")

jpeg("./results/cec_cvd_mr_ors.jpg", width = 10000, height = 7000, res=1200)
row_names<-c(level00,level01,level02,level03,level04,level05,level06,level07,level08,level09,level10,
             level11,level12,level15,level16,level17,level18,level19,level22,level23,level24,level25,level26)
forestplot(labeltext=row_names,
           c(NA,NA,tab[1,33],tab[1,24],tab[1,27],tab[1,30],tab[5,33],tab[9,33],
             NA,tab[2,33],tab[2,24],tab[2,27],tab[2,30],
             NA,tab[3,33],tab[3,24],tab[3,27],tab[3,30],
             NA,tab[4,33],tab[4,24],tab[4,27],tab[4,30]),
           c(NA,NA,tab[1,34],tab[1,25],tab[1,28],tab[1,31],tab[5,34],tab[9,34],
             NA,tab[2,34],tab[2,25],tab[2,28],tab[2,31],
             NA,tab[3,34],tab[3,25],tab[3,28],tab[3,31],
             NA,tab[4,34],tab[4,25],tab[4,28],tab[4,31]),
           c(NA,NA,tab[1,35],tab[1,26],tab[1,29],tab[1,32],tab[5,35],tab[9,35],
             NA,tab[2,35],tab[2,26],tab[2,29],tab[2,32],
             NA,tab[3,35],tab[3,26],tab[3,29],tab[3,32],
             NA,tab[4,35],tab[4,26],tab[4,29],tab[4,32]),
           lwd.zero=1,
           boxsize = 0.35,
           cex=1,
           txt_gp = fpTxtGp(ticks=gpar(cex=0.8),xlab=gpar(cex=0.8),label=gpar(cex=0.9),title=gpar(cex=1.1)),
           fn.ci_norm = c(fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI,
                          fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI,
                          fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI,
                          fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI,
                          fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI),
           col = fpColors(all.elements="black"),
           xlog=TRUE,
           #graphwidt=unit(0.3,"npc"),
           lwd.ci=2,
           lwd.xaxis=2,
           ci.vertices=T,
           ci.vertices.height = 0.18,
           xticks = c(0.25,0.5,1,1.5),
           mar=unit(c(-0.1,0.1,0,0.5),"cm"),
           hrzl_lines=list("9"=gpar(lty=2),"14"=gpar(lty=2),"19"=gpar(lty=2)),
           clip=c(0.18,1.51),
           xlab = "Disease risk (odds ratio [95% CI])")
dev.off()


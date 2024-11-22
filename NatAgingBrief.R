
library(tidyr)
library(tidyverse)
library(ggplot2)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(Rmisc)
library(boot)
library(table1)
library(lme4)
library(lmerTest)
library(car)
library(psycho)
library(gridExtra)
library(gtable)
library(dplyr)
library(sensemakr)
library(flextable)
library(rstatix)
library(reshape2)
library(ggrepel)
library(ggpubr)

# Load MVPA Data
mvpa_enc_long = read.csv('mvpa_encoding_VTC_ANG.csv')

mvpa_enc <- mvpa_enc_long %>%
  mutate(category=if_else(category=='WF','face','place')) %>%
  mutate(category=factor(category), roi=factor(roi)) %>% 
  group_by(subject,roi) %>%
  dplyr::summarise(logit = mean(logit),
                   acc = mean(acc))

mvpa_reinst_long = read.csv('mvpa_reinstatement_VTC_ANG.csv')

mvpa_reinst <- mvpa_reinst_long %>%
  group_by(subject,roi,trial_type) %>%
  dplyr::summarise(logit = mean(logit),
                   acc = mean(acc),
                   RT = mean(RT),
                   num_trials = n()) %>%
  filter(num_trials>=5)

sub_list <- unique(mvpa_enc$subject) ## N=159 with MVPA encoding data

# Load biomarkers and demographics
sams_biomarkers = read.csv('SAMS_Biomarker_Summary_AcrossWaves_102424.csv')

sams_biomarkers <- sams_biomarkers %>% filter(Wave_Ref=='Baseline') %>% dplyr::select(-APOE, -Education)

# load memory data
sams_mem = read.csv('SAMS_Baseline_Level3_0224.csv')

sams_biomarkers <- sams_biomarkers %>% left_join(sams_mem[, c('pidn','subject','age.cni', 'educ', 'APOE', 'APOE_e4', 'assoc_dp', 'assoc_hr', 'assoc_far','mmse', 'delayed_recall_comp')], by='pidn') 

# combine and filter for fmri subs
sams_fmri <- sams_biomarkers %>% 
  filter(subject%in%sub_list) %>% 
  mutate(amyloid_status = factor(amy_status_gold_trans),
         amy_status_csf = factor(amy_status_csf),
         amy_status_pet = factor(amy_status_pet),
         amy_status_gold=factor(amy_status_gold)) %>%
  mutate(sex = factor(Gender, levels=c("Male","Female"))) %>%
  dplyr::rename(pTau181=CSF_P_Tau181,
                AB42=CSF.AB42,
                AB42_AB40=CSF_AB42_AB40) %>%
  mutate(logpTau181=log(pTau181),
         ptau181_AB42 = pTau181/AB42,
         logptau181_AB42 = log(ptau181_AB42)) %>%
  mutate(APOE_e4 = factor(APOE_e4),
         APOE=factor(APOE))


## prepare encoding data
mvpa_enc_wide <- mvpa_enc %>%
  group_by(subject,roi) %>%
  dplyr::summarise(logit=mean(logit), acc=mean(acc)) %>%
  pivot_wider(id_cols = c('subject'), names_from = roi, values_from = c('logit','acc')) 

# prepare mvpa reinstatement data
reinst_ah <- mvpa_reinst %>% 
  filter(trial_type == 'AH') %>% 
  pivot_wider(id_cols = c('subject'), names_from = roi, values_from = c('logit','acc'))

# join all together
sams_fmri <- sams_fmri %>% left_join(mvpa_enc_wide, by='subject') %>% left_join(reinst_ah, by='subject')


## create residualized reinstatement variables
m<-lm(logit_VTC.y~logit_VTC.x, data = sams_fmri)
sams_fmri$VTCres <- rstandard(m) ## this yields a standardized residual 

m<-lm(logit_ANG.y~logit_ANG.x, data = sams_fmri)
sams_fmri$ANGres <- rstandard(m) ## this yields a standardized residual 
sams_fmri$reinst_avg <-rowMeans(sams_fmri[c('VTCres','ANGres')])

sams_fmri$VTC.enc.z <- scale(sams_fmri$logit_VTC.x) # z-score to create encoding composite (for control analyses)
sams_fmri$ANG.enc.z <- scale(sams_fmri$logit_ANG.x)
sams_fmri$encoding_avg <-rowMeans(sams_fmri[c('VTC.enc.z','ANG.enc.z')])



# Create Table 1

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.0f %%)", FREQ, PCT))))
}

mytable<- table1::table1(~age.cni + Gender + educ + Race + Ethnicity + APOE + APOE_e4 + amyloid_status + AB42_AB40 + fs7_centiloid + pTau181 + mmse, data = sams_fmri, overall=c(left="Total"),render.categorical=my.render.cat)

t1flex(mytable) %>%
  save_as_docx(path="NatAgingTable1.docx")


# Models predicting biomarkers

## age & APOE predict amyloid;  age, APOE, amyloid explain unique variance in pTau
m1<-lm(AB42_AB40~age.cni+sex+APOE_e4, data = sams_fmri)
summary(m1)

m2<-lm(logpTau181~age.cni+sex+APOE_e4+AB42_AB40, data = sams_fmri)
summary(m2)

m3<-lm(logpTau181~age.cni+sex+APOE_e4+amyloid_status, data = sams_fmri)
summary(m3)

tab_model(m1,m2,m3,collapse.ci = TRUE, show.std = TRUE)


# Models predicting reinstatement by ROI

reinst_stack <- sams_fmri %>% 
  dplyr::select(subject,age.cni, sex, educ, assoc_dp, APOE_e4, amyloid_status, AB42_AB40, logpTau181, PlasmaABRatio, PlasmaPTau181, VTCres, ANGres) %>%
  pivot_longer(cols=VTCres:ANGres, names_to='roi', values_to ='reinst') %>%
  mutate(roi=factor(roi))

m1 = lmer(reinst ~ roi*(age.cni + sex) + (1|subject), data=reinst_stack)
summary(m1)

m1 = lmer(reinst ~ roi*(age.cni + sex + APOE_e4) + (1|subject), data=reinst_stack)
summary(m1)

m2 = lmer(reinst ~ roi*(age.cni + sex + amyloid_status) + (1|subject), data=reinst_stack)
summary(m2)

m3 = lmer(reinst ~ roi*(age.cni + sex + logpTau181) + (1|subject), data=reinst_stack)
summary(m3)

tab_model(m1,m2,m3,collapse.ci = TRUE, show.std = TRUE)


# Models predicting mean reinstatement score

## Age effects

m1 = lm(reinst_avg ~ age.cni + sex , data=sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m1)

# Main effects of individual AD factors
m1 = lm(reinst_avg ~ age.cni + sex + APOE_e4, data=sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m1)

m2 = lm(reinst_avg ~ age.cni + sex + amyloid_status, data=sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m2)

m3 = lm(reinst_avg ~ age.cni + sex + logpTau181, data=sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m3)


## APOE effects controlling for amyloid and tau
m1 = lm(reinst_avg ~ age.cni + sex + APOE_e4 + amyloid_status, data=sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m1)

m2= lm(reinst_avg ~ age.cni + sex + APOE_e4 + AB42_AB40, data=sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m2)

m3 = lm(reinst_avg ~ age.cni + sex + APOE_e4 + logpTau181, data=sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m3)


## ptau effects controlling for amyloid

m1 = lm(reinst_avg ~ age.cni + sex + logpTau181, data=sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m1)

m2 = lm(reinst_avg ~ age.cni + sex + logpTau181 + amyloid_status , data=sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m2)

m3 = lm(reinst_avg ~ age.cni + sex + logpTau181 + AB42_AB40, data=sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m3)

## APOE x biomarker interactions

m1 = lm(reinst_avg ~ age.cni + sex + APOE_e4*amyloid_status, data=sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m1)

m2 = lm(reinst_avg ~ age.cni + sex + APOE_e4*logpTau181, data=sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m2)

m3= lm(reinst_avg ~ age.cni + sex + APOE_e4*logpTau181+amyloid_status, data=sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m3)

## APOE x Amyloid

# e4 effect is evident in amyloid positive but not amy negative 
m1 = lm(reinst_avg ~ age.cni + sex + APOE_e4*amyloid_status, data=sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m1)

m2 = lm(reinst_avg ~ age.cni + sex + APOE_e4, data=sams_fmri  %>% filter(!subject%in%subs_exclude) %>% filter(amyloid_status=="A+"))
summary(m1)

m3 = lm(reinst_avg ~ age.cni + sex + APOE_e4, data=sams_fmri %>% filter(!subject%in%subs_exclude) %>% filter(amyloid_status=="A-"))
summary(m3)

tab_model(m1,m2,m3,collapse.ci = TRUE, show.std = TRUE)


## APOE x pTau
m1 = lm(reinst_avg ~ age.cni + sex + APOE_e4*logpTau181, data=sams_fmri %>% filter(!is.na(logpTau181)))
summary(m1)

## no effect in non-carriers, sig effect in carriers 
m2 = lm(reinst_avg ~ age.cni + sex + logpTau181, data=sams_fmri %>% filter(!is.na(logpTau181)) %>% filter(APOE_e4==0))
summary(m2)

m3 = lm(reinst_avg ~ age.cni + sex + logpTau181, data=sams_fmri %>% filter(!is.na(logpTau181)) %>% filter(APOE_e4==1))
summary(m3)

# remains significant controlling for amyloid status
m4 = lm(reinst_avg ~ age.cni + sex + logpTau181 + amyloid_status, data=sams_fmri %>% filter(!is.na(logpTau181)) %>% filter(APOE_e4==1))
summary(m4)

# Models predicting memory

## Reinstatement 
m1<-lm(reinst_avg~age.cni+sex, data = sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m1)

m2<-lm(assoc_dp~age.cni+sex+educ+reinst_avg, data = sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m2)

m3<-lm(delayed_recall_comp~age.cni+sex+educ+reinst_avg, data = sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m3)

## AD markers
m1<-lm(assoc_dp~age.cni+sex+educ+APOE_e4, data = sams_fmri)
summary(m1)

m2<-lm(assoc_dp~age.cni+sex+educ+AB42_AB40, data = sams_fmri)
summary(m2)

m3<-lm(assoc_dp~age.cni+sex+educ+amyloid_status, data = sams_fmri)
summary(m3)

m4<-lm(assoc_dp~age.cni+sex+educ+logpTau181, data = sams_fmri)
summary(m4)

## Multiple AD markers
m1<-lm(assoc_dp~age.cni+sex+educ+APOE_e4+logpTau181, data = sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m1)

m2<-lm(assoc_dp~age.cni+sex+educ+APOE_e4+AB42_AB40, data = sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m2)

m3<-lm(assoc_dp~age.cni+sex+educ+APOE_e4+amyloid_status, data = sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m3)

## APOE x pTau
m1<-lm(assoc_dp~age.cni+sex+educ+APOE_e4*logpTau181, data = sams_fmri)
summary(m1)

m2<-lm(assoc_dp~age.cni+sex+educ+logpTau181, data = sams_fmri %>% filter(APOE_e4==0))
summary(m2)

m3<-lm(assoc_dp~age.cni+sex+educ+logpTau181, data = sams_fmri %>% filter(APOE_e4==1))
summary(m3)

m4<-lm(assoc_dp~age.cni+sex+educ+logpTau181+amyloid_status, data = sams_fmri %>% filter(APOE_e4==1))
summary(m4)

## APOE x Amyloid
m1<-lm(assoc_dp~age.cni+sex+educ+APOE_e4*amyloid_status, data = sams_fmri %>% filter(!subject%in%subs_exclude) )
summary(m1)


## Multimodal with interaction term
m1<-lm(assoc_dp~age.cni+sex+educ+reinst_avg, data = sams_fmri %>% filter(!subject%in%subs_exclude) %>% filter(!is.na(APOE_e4)) %>% filter(!is.na(logpTau181)))
summary(m1)

m2<-lm(assoc_dp~age.cni+sex+educ+APOE_e4*logpTau181, data = sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m2)

m3<-lm(assoc_dp~age.cni+sex+educ+APOE_e4*logpTau181+reinst_avg, data = sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m3)

## pTau + reinst 

# in e4 non-carriers
m1<-lm(assoc_dp~age.cni+sex+educ+reinst_avg+logpTau181, data = sams_fmri %>% filter(!subject%in%subs_exclude) %>% filter(APOE_e4==0))
summary(m1)

# in e4 carriers
m2<-lm(assoc_dp~age.cni+sex+educ+reinst_avg+logpTau181, data = sams_fmri %>% filter(!subject%in%subs_exclude) %>% filter(APOE_e4==1))
summary(m2)

# controlling for amyloid
m3<-lm(assoc_dp~age.cni+sex+educ+reinst_avg+logpTau181+amyloid_status, data = sams_fmri %>% filter(!subject%in%subs_exclude) %>% filter(APOE_e4==1))
summary(m3)


# Delayed recall

## Memory by AD markers 
m1<-lm(delayed_recall_comp~age.cni+sex+educ+APOE_e4, data = sams_fmri)
summary(m1)

m2<-lm(delayed_recall_comp~age.cni+sex+educ+AB42_AB40, data = sams_fmri)
summary(m2)

m3<-lm(delayed_recall_comp~age.cni+sex+educ+amyloid_status, data = sams_fmri)
summary(m3)

m4<-lm(delayed_recall_comp~age.cni+sex+educ+logpTau181, data = sams_fmri)
summary(m4)

# APOE + biomarkers
m1<-lm(delayed_recall_comp~age.cni+sex+educ+APOE_e4+amyloid_status, data = sams_fmri)
summary(m1)

m2<-lm(delayed_recall_comp~age.cni+sex+educ+APOE_e4+logpTau181, data = sams_fmri)
summary(m2)

# Multimodal models with interaction term

m1<-lm(delayed_recall_comp~age.cni+sex+educ+reinst_avg, data = sams_fmri %>% filter(!subject%in%subs_exclude) %>% filter(!is.na(APOE_e4)) %>% filter(!is.na(logpTau181)))
summary(m1)

m2<-lm(delayed_recall_comp~age.cni+sex+educ+APOE_e4*logpTau181, data = sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m2)

m3<-lm(delayed_recall_comp~age.cni+sex+educ+APOE_e4*logpTau181 + reinst_avg, data = sams_fmri %>% filter(!subject%in%subs_exclude))
summary(m3)

# Stratified by APOE
m1<-lm(delayed_recall_comp~age.cni+sex+educ+logpTau181, data = sams_fmri %>% filter(APOE_e4==0))
summary(m1)

m2<-lm(delayed_recall_comp~age.cni+sex+educ+logpTau181, data = sams_fmri %>% filter(APOE_e4==1))
summary(m2)

# Figures

sams_plot <- sams_fmri %>% filter(!subject%in%subs_exclude)
sams_plot$APOE_e4 <- factor(sams_plot$APOE_e4, levels=c(0,1), labels=c("non-carrier","carrier"))
sams_plot$amyloid_status <- factor(sams_plot$amyloid_status, levels=c('A-','A+'), labels=c("Aβ-","Aβ+"))

## Figure 1
p1<-ggplot(sams_plot %>% filter(!is.na(APOE)), aes(x=APOE_e4, y=reinst_avg,color=APOE_e4)) + 
  geom_boxplot(outlier.shape = NA) + scale_color_manual(name=expression(paste(italic("APOE4"))),values=c("#56B4E9","#E69F00")) +
  geom_jitter(alpha=0.7,size=1,width = 0.2) + scale_x_discrete(labels = c("non-carrier","carrier")) +
  labs(title="", x=expression(paste(italic("APOE4"))), y="Reinstatement Strength (z)",tag="a")+guides(color="none")

p2<-ggplot(sams_plot %>% filter(!is.na(amyloid_status)), aes(x=amyloid_status, y=reinst_avg,color=amyloid_status)) + 
  geom_boxplot(outlier.shape = NA) + scale_color_manual(values=c("#00BFC4","#F8766D"))+
  geom_jitter(alpha=0.7,size=1,width = 0.2) + 
  labs(title="", x="Aβ Status", y="Reinstatement Strength (z)",tag="b")+guides(color="none")

p3<-ggplot(sams_plot %>% filter(!is.na(APOE)) %>% filter(!is.na(amyloid_status)), aes(x=amyloid_status, y=reinst_avg, color=APOE_e4)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(alpha=.7, size=1,position=position_jitterdodge()) +
  scale_color_manual(name=expression(paste(italic("APOE4"))),values=c("#56B4E9","#E69F00")) +
  labs(title="", x="Aβ Status", y="Reinstatement Strength (z)",tag="d")

upper <- grid.arrange(p1, p2, nrow=1,ncol=2)

p <- grid.arrange(upper, p3, ncol=1, nrow=2)


p4<-ggplot(sams_plot %>% filter(!is.na(logpTau181)), aes(x=logpTau181, y=reinst_avg, color=amyloid_status)) +
  geom_point(alpha=0.7,size=1.5,shape=16) +
  geom_smooth(formula = y ~ x, method=lm, color='black') + scale_color_manual(name="",values=c("#00BFC4","#F8766D", "grey69")) +
  labs(x="CSF pTau181 (log)", y="Reinstatement Strength (z)", tag="c") +
  theme(legend.justification=c(1,1), legend.position=c(0.9,1.1))

p5<-ggplot(sams_plot %>% filter(!is.na(APOE)) %>% filter(!is.na(logpTau181)), aes(x=logpTau181, y=reinst_avg, color=APOE_e4)) +
  geom_point(alpha=0.7) +
  geom_smooth(formula = y ~ x, method=lm) + 
  scale_color_manual(name="",values=c("#56B4E9","#E69F00"), labels=c(expression(paste(italic("APOE4"), " non-carrier")),expression(paste(italic("APOE4"), " carrier       ")))) +
  theme(legend.justification=c(1,1), legend.position=c(1.0,1.1))+
  labs(x="CSF pTau181 (log)", y="Reinstatement Strength (z)", tag="e")  

g<-grid.arrange(p4,p5,ncol=1)

fig<-grid.arrange(p,g,ncol=2)

ggsave("~/OneDrive - Stanford/AgingR01/manuscripts/Trelle_APOE_CSF_Reinstatement/figures/Figure1.png", width = 10, height = 8,fig)


## Figure 2
p1<-ggplot(sams_plot %>% filter(!is.na(APOE)), aes(x=APOE_e4, y=assoc_dp,color=APOE_e4)) + 
  geom_boxplot(outlier.shape = NA)+ scale_color_manual(name=expression(paste(italic("APOE4"))),values=c("#56B4E9","#E69F00")) +
  geom_jitter(alpha=.7, size=1,width = 0.2) + scale_x_discrete(labels = c("non-carrier","carrier")) +
  labs(title="", x=expression(paste(italic("APOE4"))), y=expression(paste("Associative Memory (", italic("d') "))),tag="a") + guides(color="none")

p2<-ggplot(sams_plot %>% filter(!is.na(amyloid_status)), aes(x=amyloid_status, y=assoc_dp,color=amyloid_status)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=.7, size=1,width = 0.2) + scale_color_manual(name="",values=c("#00BFC4","#F8766D")) +
  labs(title="", x="Aβ Status", y=expression(paste("Associative Memory (", italic("d') "))),tag="b")+guides(color="none")

p3<-ggplot(sams_plot %>% filter(!is.na(APOE)) %>% filter(!is.na(amyloid_status)), aes(x=amyloid_status, y=assoc_dp, color=APOE_e4)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(alpha=.7, size=1,position=position_jitterdodge()) + 
  scale_color_manual(name=expression(paste(italic("APOE4"))),values=c("#56B4E9","#E69F00")) +
  labs(title="", x="Aβ Status", y=expression(paste("Associative Memory (", italic("d') "))),tag="d")  

upper <- grid.arrange(p1, p2, nrow=1,ncol=2)

p <- grid.arrange(upper, p3, ncol=1, nrow=2)


p4<-ggplot(sams_plot %>% filter(!is.na(logpTau181)), aes(x=logpTau181, y=assoc_dp, color=amyloid_status)) +
  geom_point(alpha=0.7,size=1.5,shape=16) +
  geom_smooth(formula = y ~ x, method=lm, color='black') + scale_color_manual(name="",values=c("#00BFC4","#F8766D", "grey69")) +
  labs(x="CSF pTau181 (log)", y=expression(paste("Associative Memory (", italic("d') "))), tag="c") + 
  theme(legend.justification=c(1,1), legend.position=c(0.9,1.1))

# p5<-ggplot(sams_plot %>% filter(!is.na(APOE)) %>% filter(!is.na(logpTau181)), aes(x=logpTau181, y=assoc_dp, color=APOE_e4)) +
#   geom_point(alpha=.7,size=1.5,aes(shape=amyloid_status)) + scale_shape_manual(values=c(16, 17))+
#   geom_smooth(formula = y ~ x, method=lm) + scale_color_manual(name="APOE ε4",values=c("#56B4E9","#E69F00")) +
#   labs(x="CSF pTau181 (log)", y=expression(paste("Associative Memory (", italic("d') "))), tag="E") +guides(color="none",shape="none") 

p5<-ggplot(sams_plot %>% filter(!is.na(APOE)) %>% filter(!is.na(logpTau181)), aes(x=logpTau181, y=assoc_dp, color=APOE_e4)) +
  geom_point(alpha=0.7) +
  geom_smooth(formula = y ~ x, method=lm) + 
  scale_color_manual(name="",values=c("#56B4E9","#E69F00"), labels=c(expression(paste(italic("APOE4"), " non-carrier")),expression(paste(italic("APOE4"), " carrier       ")))) +
  theme(legend.justification=c(1,1), legend.position=c(1.0,1.1))+
  labs(x="CSF pTau181 (log)", y=expression(paste("Associative Memory (", italic("d') "))), tag="e")  


g<-grid.arrange(p4,p5,ncol=1)


fig<-grid.arrange(p,g,ncol=2)

ggsave("~/OneDrive - Stanford/AgingR01/manuscripts/Trelle_APOE_CSF_Reinstatement/figures/Figure2.png", width = 10, height = 8,fig)


## Extended Data Figure 1 
p1<-ggplot(sams_plot %>% filter(!is.na(APOE_e4)) %>% filter(!is.na(AB42_AB40)), aes(x=APOE_e4, y=AB42_AB40,color=APOE_e4)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.7,width = 0.2,size=1.5) +
  scale_x_discrete(labels = c("non-carrier","carrier")) +
  scale_color_manual(name=expression(paste(italic("APOE4"))),values=c("#56B4E9","#E69F00")) +
  labs(x=expression(paste(italic("APOE4"))), y=expression("CSF A"*beta*"42/40"), tag="a") +guides(color="none")

p2<-ggplot(sams_plot %>% filter(!is.na(APOE_e4)) %>% filter(!is.na(pTau181)), aes(x=amyloid_status, y=pTau181,color=APOE_e4)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(alpha=0.7,size=1.5,position = position_jitterdodge()) +
  # scale_shape_discrete(name="", solid=TRUE)+
  # scale_color_manual(name="", values=c("#00BFC4","#F8766D"))+
  scale_color_manual(name="",values=c("#56B4E9","#E69F00"), labels = c("non-carrier","carrier"))+
  # scale_x_discrete(labels = c("non-carrier","carrier")) +
  labs(x="Aβ Status", y="CSF pTau181 (pg/ml)", tag="b") + 
  theme(legend.justification=c(1,1), legend.position=c(0.50,1.1)) 


p<-grid.arrange(p1,p2,ncol=2)

ggsave("~/OneDrive - Stanford/AgingR01/manuscripts/Trelle_APOE_CSF_Reinstatement/figures/SuppFig1.png", width = 7, height = 4,p)

## Extended Data Figure 2 and 4
p1<- ggplot() + labs(x="",y="",tag="a")

p2<-ggplot(sams_plot, aes(x=age.cni, y=assoc_dp)) +
  geom_point(alpha=0.7) +   
  geom_smooth(formula = y ~ x, method=lm, color='black') +
  labs(x="Age",y=expression(paste("Associative Memory (", italic("d') "))), tag="b")

p3<-ggplot(sams_plot, aes(x=age.cni, y=reinst_avg)) +
  geom_point(alpha=0.7) +   
  geom_smooth(formula = y ~ x, method=lm, color='black') + 
  labs(x="Age",y="Reinstatement Strength (z)" , tag="a")

p4<-ggplot(sams_plot, aes(x=reinst_avg, y=assoc_dp)) +
  geom_point(alpha=0.7) +   
  geom_smooth(formula = y ~ x, method=lm, color='black') + 
  labs(x="Reinstatement Strength (z)",y=expression(paste("Associative Memory (", italic("d') "))) , tag="b")

p5<-ggplot(sams_plot, aes(x=reinst_avg, y=delayed_recall_comp)) +
  geom_point(alpha=0.7) +   
  geom_smooth(formula = y ~ x, method=lm, color='black') + 
  labs(x="Reinstatement Strength (z)",y="Delayed Recall Score (z)", tag="c")

p<-grid.arrange(p1,p2,ncol=2)

ggsave("~/OneDrive - Stanford/AgingR01/manuscripts/Trelle_APOE_CSF_Reinstatement/figures/SuppFig2.png", width = 10, height = 4,p)

p<- grid.arrange(p3,p4,p5, ncol=2)
ggsave("~/OneDrive - Stanford/AgingR01/manuscripts/Trelle_APOE_CSF_Reinstatement/figures/SuppFig4.png", width = 10, height = 8,p)


## Extended Data Figure 3

p1<-ggplot(reinst_stack, aes(x=age.cni, y=reinst, color=roi)) +
  geom_point(alpha=0.7) +
  geom_smooth(formula = y ~ x, method=lm) + scale_color_manual(name="",values=c("#1B9E77","#7570B3"), labels=c("ANG","VTC")) + 
  labs(x="Age", y="Reinstatement Strength (std res)", tag="") + theme(legend.position = 'top')

p2<-ggplot(reinst_stack, aes(x=reinst, y=assoc_dp, color=roi)) +
  geom_point(alpha=0.7) +
  geom_smooth(formula = y ~ x, method=lm) + scale_color_manual(name="",values=c("#1B9E77","#7570B3"), labels=c("ANG","VTC")) + 
  labs(x="Reinstatement Strength (std res)", y=expression(paste("Associative Memory (", italic("d') "))), tag="") + theme(legend.position = 'top')

p<-grid.arrange(p1,p2,ncol=2)

p3<-ggplot(reinst_stack %>% filter(!is.na(APOE_e4)), aes(x=roi, y=reinst, color=APOE_e4)) + 
  geom_boxplot(outlier.shape = NA) + scale_color_manual(name="",values=c("#56B4E9","#E69F00"), labels=c('ε4-','ε4+')) +
  scale_x_discrete(labels=c("ANG","VTC"))+theme(legend.position = 'top')+
  geom_point(alpha=.7, size=1,position=position_jitterdodge()) +  
  labs(x="", y="Reinstatement Strength (std res)", tag="a")

p4<-ggplot(reinst_stack %>% filter(!is.na(amyloid_status)), aes(x=roi, y=reinst, color=amyloid_status)) + 
  geom_boxplot(outlier.shape = NA) + scale_color_manual(name="",values=c("#00BFC4","#F8766D"), labels=c("Aβ-","Aβ+")) +
  scale_x_discrete(labels=c("ANG","VTC"))+theme(legend.position = 'top')+
  geom_point(alpha=.7, size=1,position=position_jitterdodge()) +  
  labs(x="", y="Reinstatement Strength (std res)", tag="b")

pp<-grid.arrange(p3,p4, ncol=2)

p5<-ggplot(reinst_stack %>% filter(!is.na(logpTau181)), aes(x=logpTau181, y=reinst, color=roi)) +
  geom_point(alpha=0.7) +
  geom_smooth(formula = y ~ x, method=lm) + scale_color_manual(name="",values=c("#1B9E77","#7570B3"), labels=c("ANG","VTC")) + 
  labs(x="CSF pTau181 (log)", y="Reinstatement Strength (std res)", tag="c") + theme(legend.position = 'top')

ppp <- grid.arrange(pp, p5, ncol=2)


ggsave("~/OneDrive - Stanford/AgingR01/manuscripts/Trelle_APOE_CSF_Reinstatement/figures/SuppFig3.png", width = 10, height = 4.5, ppp)


## Extended Data Figure 5
p1<-ggplot(sams_plot %>% filter(!is.na(APOE)), aes(x=APOE_e4, y=delayed_recall_comp,color=APOE_e4)) + 
  geom_boxplot(outlier.shape = NA) + scale_color_manual(name=expression(paste(italic("APOE4"))),values=c("#56B4E9","#E69F00")) +
  geom_jitter(alpha=0.7,size=1,width = 0.2) + scale_x_discrete(labels = c("non-carrier","carrier")) +
  labs(title="", x=expression(paste(italic("APOE4"))), y="Delayed Recall Score (z)",tag="a")+guides(color="none")

p2<-ggplot(sams_plot %>% filter(!is.na(amyloid_status)), aes(x=amyloid_status, y=delayed_recall_comp,color=amyloid_status)) + 
  geom_boxplot(outlier.shape = NA) + scale_color_manual(values=c("#00BFC4","#F8766D"))+
  geom_jitter(alpha=0.7,size=1,width = 0.2) + 
  labs(title="", x="Aβ Status", y="Delayed Recall Score (z)",tag="b")+guides(color="none")

p3<-ggplot(sams_plot %>% filter(!is.na(APOE)) %>% filter(!is.na(amyloid_status)), aes(x=amyloid_status, y=delayed_recall_comp, color=APOE_e4)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(alpha=.7, size=1,position=position_jitterdodge()) + 
  scale_color_manual(name=expression(paste(italic("APOE4"))),values=c("#56B4E9","#E69F00")) +
  labs(title="", x="Aβ Status", y="Delayed Recall Score (z)",tag="d")  

upper <- grid.arrange(p1, p2, nrow=1,ncol=2)

p <- grid.arrange(upper, p3, ncol=1, nrow=2)


p4<-ggplot(sams_plot %>% filter(!is.na(logpTau181)), aes(x=logpTau181, y=delayed_recall_comp, color=amyloid_status)) +
  geom_point(alpha=0.7,size=1.5,shape=16) +
  geom_smooth(formula = y ~ x, method=lm, color='black') + scale_color_manual(name="",values=c("#00BFC4","#F8766D", "grey69")) +
  labs(x="CSF pTau181 (log)", y="Delayed Recall Score (z)", tag="c") + 
  theme(legend.justification=c(1,1), legend.position=c(0.9,1.1))

p5<-ggplot(sams_plot %>% filter(!is.na(APOE)) %>% filter(!is.na(logpTau181)), aes(x=logpTau181, y=reinst_avg, color=APOE_e4)) +
  geom_point(alpha=0.7) +
  geom_smooth(formula = y ~ x, method=lm) + 
  scale_color_manual(name="",values=c("#56B4E9","#E69F00"), labels=c(expression(paste(italic("APOE4"), " non-carrier")),expression(paste(italic("APOE4"), " carrier       ")))) +
  theme(legend.justification=c(1,1), legend.position=c(1.0,1.1))+
  labs(x="CSF pTau181 (log)", y="Delayed Recall Score (z)", tag="e") 

g<-grid.arrange(p4,p5,ncol=1)


fig<-grid.arrange(p,g,ncol=2)

ggsave("~/OneDrive - Stanford/AgingR01/manuscripts/Trelle_APOE_CSF_Reinstatement/figures/SuppFig5_legend_edit.png", width = 10, height = 8,fig)


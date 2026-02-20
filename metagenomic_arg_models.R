# Author: Noel Waters Feb 2026
# Purpose: generalized linear mixed effects models for modelling arg carriage rates
# in bacteria in municipal wastewater treatment plant influents vs hospitals and biofilms vs wastewaters


{
  library(tidyverse)
  library(DHARMa)
  library(emmeans)
  library(glmmTMB)
}




# Preprocess input data. --------------------------------------------------


resfinder_counts<-read.csv("../../metagenomics/ResFinder_Counts_Final.csv",header = TRUE) 
metadata<-readxl::read_xlsx("../../metagenomics/metadata.xlsx") |> rename(Sample_Type=Type,
                                                                          Site_Type=Type_Site)

metadata$Sample_Type[metadata$Sample_Type=="Wastewater"]<-"Water"

# merge these two:
resfinder_counts<- resfinder_counts |> merge(metadata,by="Sample")

resfinder_counts<- resfinder_counts |> rename(Total_Reads=`Total Reads_After`) |> dplyr::select(-c(X))

resfinder_counts$CPM<-resfinder_counts$Count/resfinder_counts$Total_Reads*1e6



# Check that names are correctly spelled.

resfinder_counts$Site<-str_replace_all(str_to_title(str_replace_all(resfinder_counts$Site,"_"," "))," ","_")



# updated total reads, just 1/3 of the previous values..
resfinder_counts$Site_Type<-ifelse(resfinder_counts$Site_Type=="Hospital","Hospital",
                                   ifelse(grepl("influent|incoming",ignore.case = TRUE,x = resfinder_counts$Site),"Influent",
                                          ifelse(grepl("effluent",ignore.case = TRUE,x = resfinder_counts$Site),"Effluent","Fix")
                                   )
)
# sum across ALL arg types since we care about the overall
resfinder_totals<-resfinder_counts |>  group_by(Sample_ID,Sample,Sample_Type,Site_Type,Site,Total_Reads) |> 
  summarise(Total_CPM=sum(CPM),
            Total_Count=sum(Count)) |>
  ungroup() |>
  mutate(Total_Millions=Total_Reads/1e6) |> 
  mutate(Site_Type=factor(Site_Type)) |> 
  mutate(Site_Type=relevel(Site_Type,ref="Influent")) |> 
  mutate(Sample_Type=relevel(factor(Sample_Type),ref="Water")) |> ungroup() |> 
  mutate(log2CPM=log(Total_CPM,base=2)) |> 
  mutate(logTotalMillions=log(Total_Millions))



#saveRDS(resfinder_totals,file = "resfinder_totals_16_02_2026.rds")




# Getting an overview of the counts data by site type and water type
# We see that 1: mölndal biofilm stands out considerably and there are no effluent wastewater samples sequenced.
# We will, like for the e.coli rates, not 
ggplot(resfinder_totals,aes(x=Site_Type,y=log2CPM,
                              tal_CPM,fill=Sample_Type)) + geom_boxplot()+ geom_point(position=position_jitterdodge(jitter.width=0.1), pch=24,size=2,alpha=0.6)

# Here we see the distribution across site types
ggplot(resfinder_totals ,
       aes(y=log2CPM,x=Site,color=Sample_Type, group=Sample))+
  geom_point(position=position_dodge(width=0.4)) + theme(axis.text.x = element_text(angle=90)) + facet_wrap(~Site_Type, scales="free_x")

# We se consistently lower arg carraige in influent biofilms vs waters, and the opposite trend for hospitals.
# With effluents there are no water data so can no say.
# Now to check if differences are statistically significant, let's run a linear mixed effects model that can account 
# for the fact that some samples come from the same sites.


# For the modelling, we drop the effluent samples

table(resfinder_totals$Site_Type)

resfinder_totals_no_effluent<-resfinder_totals |> filter(Site_Type!="Effluent")

# the outlying mölndal biofilm sample
mölndal<-resfinder_totals |> filter(Site=="Mölndal" & Sample_Type=="Biofilm") |> pull(Sample)

no_mölndal<- resfinder_totals_no_effluent |> filter(Sample!=mölndal)




# First, A Comparison of effluent and influent biofilms ----------------------------
municipal_biofilms<-resfinder_totals |> filter(Sample_Type=="Biofilm" & Site_Type !="Hospital")
influents<-municipal_biofilms |> filter(Site_Type=="Influent") |> pull(log2CPM)
effluents<-municipal_biofilms |> filter(Site_Type=="Effluent") |> pull(log2CPM)

t.test(influents,effluents, var.equal = FALSE)
#unpaired, comparing all to all is highly significant. A paired test would also be.




# Modelling of arg counts -------------------------------------------------
# Idea is to Model the Total Counts and use the Counts per Million as offset, such that what
#is effectively modelled are the ARG counts per million reads
# This is anologous to the model used for ecoli where it was resistant ecoli per ecoli.
# Since we have counts data again,and some relatedness of samples (some come from the same site) 
# negative binomial regression with random effect for the sites comes in handy.

# 1: Fit one overall model for all site and sample types.
# If the interaction term is significant in an anova, proceed to split into two models, one for each site type.


nb_fitter<-function(df){
  
  
  fit1<-glmmTMB(
    Total_Count ~ Sample_Type * Site_Type + (1|Site) + offset(logTotalMillions) ,
    data = df,
    family = nbinom2)  
  
  fit1_no_interaction<-glmmTMB(
    Total_Count ~ Sample_Type + Site_Type + (1|Site) + offset(logTotalMillions), 
    data = df,
    family = nbinom2)  
  
  # For more complex random structure. Only use if it improved fit ( assess by AIC and DHARMA)
  fit2<-glmmTMB(
    Total_Count ~ Sample_Type * Site_Type + (1|Site/Sample) + offset(logTotalMillions) , 
    data = df,
    family = nbinom2)  
  
  
  fit2_no_interaction<-glmmTMB(
    Total_Count ~ Sample_Type + Site_Type + (1|Site/Sample) + offset(logTotalMillions), 
    data = df,
    family = nbinom2)  
  return(list(fit1=fit1,
              fit1_no_interaction=fit1_no_interaction,
              fit2=fit2,
              fit2_no_interaction=fit2_no_interaction))
}

nb_fit<-nb_fitter(df=resfinder_totals_no_effluent)

get_aic_and_anova<-function(four_fits){
print(AIC(four_fits$fit1,four_fits$fit1_no_interaction,
    four_fits$fit2,four_fits$fit2_no_interaction))
  
print(anova(four_fits$fit1_no_interaction,four_fits$fit1))

print(anova(four_fits$fit2_no_interaction,four_fits$fit2))  
}

get_aic_and_anova(nb_fit)

# The anova tells of a highly significant interaction term, and also that the simpler random effect structure was better. Proceed with that.




nb_model_by_site_type<-function(df,site_type, which_fit=1){
  indat<-df |> filter(Site_Type==site_type)
  
  fit1<-glmmTMB(
    Total_Count ~ Sample_Type + (1|Site) + offset(logTotalMillions) ,
    data = indat,
    family = nbinom2)  
  
  # Addition of deeper random structure, if needed
  fit2<-glmmTMB(
    Total_Count ~ Sample_Type + (1|Site/Sample) + offset(logTotalMillions) ,
    data = indat,
    family = nbinom2)  
  
  #AIC(fit1,fit2)
  
  if(which_fit==2){
    fit<-fit2
  }  else {
    fit<-fit1
  }
  
  emm <- emmeans(fit, ~ Sample_Type,
                 type = "link",
                 offset = 0)
  
  difference_over_sample_type<-contrast(emm, 
                                        method = "trt.vs.ctrl",
                                        ref="Water",
                                        adjust="none",
                                        type="link",
                                        infer=TRUE)
  
  return(list(fit=fit,
              fit1=fit1,
              fit2=fit2,
              Site_Type=site_type,
              key_stats=difference_over_sample_type |> data.frame() |> mutate(Site_Type=site_type),
              diffs_df=emm |> data.frame() |> mutate(Site_Type=site_type),
              emm=emm,
              preds=predict(fit,type="response") # predictions on the Total Count scale
              )
  )
}


site_types<-list("Hospital","Influent")

by_st<-purrr::map(site_types,\(x){nb_model_by_site_type(df=resfinder_totals_no_effluent,site_type = x)})

#Hospital
sims1<-by_st[[1]]$fit |> simulateResiduals() # SLightly sigmoidal residuals-No severe deviations
plot(sims1)

# Influent
sims2<-by_st[[2]]$fit |> simulateResiduals() # No severe deviations
plot(sims2)


key_stats_by_site_type<-purrr::map_df(by_st,pluck("key_stats"))
key_stats_by_site_type

emms_by_site_type<-purrr::map_df(by_st,pluck("diffs_df"))
emms_by_site_type

ggplot(emms_by_site_type , aes(y=emmean,x=Site_Type,color=Sample_Type))+
  geom_errorbar(aes(ymax=asymp.UCL,ymin=asymp.LCL), position=position_dodge(width=0.95))+
  geom_point(position=position_dodge(width=0.95))

# The overlap is large for the hospital, reflecting the large-ish pvalue of the biofilm effect

ggplot(resfinder_totals_no_effluent, aes(y = Total_Count/Total_Millions,x=Site,color=Sample_Type))+ geom_point() +
  theme(axis.text.x = element_text(angle=90))+ facet_wrap(~Site_Type)




# Same analysis mölndal biofilm dropped -------------------------------------------

nb_fit_no_mölndal<-nb_fitter(df=no_mölndal)
get_aic_and_anova(nb_fit_no_mölndal)
simulateResiduals(nb_fit_no_mölndal$fit1) |> plot() # no severe deviations

by_st_no_mölndal<-purrr::map(site_types,\(x){nb_model_by_site_type(df=no_mölndal,site_type = x)})

key_stats_by_site_type_no_mölndal<-purrr::map_df(by_st_no_mölndal,pluck("key_stats"))
key_stats_by_site_type_no_mölndal

emms_by_site_type_no_mölndal<-purrr::map_df(by_st_no_mölndal,pluck("diffs_df"))
emms_by_site_type_no_mölndal

# Illustration of model
ggplot(emms_by_site_type_no_mölndal , aes(y=emmean,x=Site_Type,color=Sample_Type))+
  geom_errorbar(aes(ymax=asymp.UCL,ymin=asymp.LCL), position=position_dodge(width=0.95))+
  geom_point(position=position_dodge(width=0.95))

# When mölndal is dropped, the p-value drops considerably for the hospital side and the interaction term in the overall model.



# Writing the results to file.

# Anova results
anova_res<-anova(nb_fit$fit1_no_interaction,nb_fit$fit1) 

anova_res_no_mölndal<-anova(nb_fit_no_mölndal$fit1_no_interaction,nb_fit_no_mölndal$fit1) 
anova_res_no_mölndal

writexl::write_xlsx(list(anova_for_interaction=anova_res,
                         metagenomics_stats=key_stats_by_site,
                         metagenomics_means=emms_by_site_type,
                         anova_for_interaction_no_mölndal=anova_res_no_mölndal,
                         metagenomics_stats_no_mölndal=key_stats_by_site_no_mölndal,
                         metagenomics_means_no_mölndal=emms_by_site_type_no_mölndal,
                         metagenomics_data=resfinder_totals), path = "./output/metagenomics_stats_18_02_2026.xlsx")




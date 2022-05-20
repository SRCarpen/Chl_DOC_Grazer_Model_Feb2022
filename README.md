# Chl_DOC_Grazer_Model_Feb2022
Chlorophyll response to DOC &amp; grazers: fits to data &amp; stability analysis

The R script Fit2015_GitHub_2022-02-02.R and two .Rdata files were used to fit and analyze the nonlinear dynamic model presented in 
S.I. Model of Phytoplankton Response to DOC and Enrichment and main text figures 2 and 3 of the manuscript

Carpenter, S.R., Pace M.L., and Wilkinson G.M. 2022. Organic color and resilience of phytoplankton to enrichment. 

The R script fits the model by maximum likelihood and conducts stability analyses for scenarios presented in the manuscript.

On 20 May 2022 the R script Fit2015_GitHub_2022-05-20.R was added to the repository.  This script generates additional panels for Figure 2 of the manuscript by Carpenter, Pace, and Wilkinson. These panels compare predicted and observed one-step rates of change in chlorophyll concentration. In addition the R script generates residual plots to be added to the online Supporting Information. 
The ReadMe.pdf was updated to include this new R script. 

The dataset of daily values Daily_PeterTuesday2015_all+ZB.Rdata was constructed from spreadsheets compiled during the original research summarized by 
Pace et al. 2017 and Wilkinson et al. 2018 and by computing daily means for high-frequency data from the EDI repository:

Pace, M., J. Cole, and S. Carpenter. 2020. Cascade project at North Temperate Lakes LTER - High Frequency Data for Whole Lake Nutrient Additions 2013-2015 ver 2. 
Environmental Data Initiative. https://doi.org/10.6073/pasta/cbe19041db41e720d84970f43156c042

The dataset ZmixModel_LakeYears_Hbird.Rdata was computed in https://github.com/SRCarpen/LakeYears_for_Zmix_Model.git
and explained in S.I. Estimating Thermocline Depth from Lake Area, DOC and P Load.  

References

Pace, M. L., Batt, R. D., Buelo, C. D., Carpenter, S. R., Cole, J. J., Kurtzweil, J. T., & Wilkinson, G. M. (2017). Reversal of a cyanobacterial bloom in response to early warnings. Proceedings of the National Academy of Sciences, 114(2), 352-357. doi:10.1073/pnas.1612424114

Wilkinson, G. M., Carpenter, S. R., Cole, J. J., Pace, M. L., Batt, R. D., Buelo, C. D., & Kurtzweil, J. T. (2018). Early warning signals precede cyanobacterial blooms in multiple whole-lake experiments. Ecological Monographs, 88(2), 188-203. doi:doi:10.1002/ecm.1286

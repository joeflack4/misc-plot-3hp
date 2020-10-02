source("3hp.1hp.model_0116_FinalMatch.R")

#You'll find the default parameters in the 3hp.1hp.model file, which contains the model.param dataframe with the default
#values.  We create a data.frame from it, which makes a copy, which we then use to run the model.
params = data.frame(model.param)

#If you want to modify parameters, you can do that here.  For a list, see the 3hp.1hp.model file model.param dataframe
#params['E_3HP'] <- .94

#You can either run the treatment.sim() function, which will give you a four element structure, eg:
#
#  Cost.1HP.total DALY.1HP.total Cost.3HP.total DALY.3HP.total 
#     1528168.164       7391.679    1522008.727       7396.724 

#Or you can run the treatment.details() function, which will collect some intermediate values and display them on completion, eg:
#
#           total_ly total_daly total_reg_complete total_deaths total_tb_death total_active_cases      ICER prev_therapy_cost
# 3HP     9378.52497  7396.7239           186.5732  638.7337729             NA           21.34555        NA         22721.861
# 1HP     9389.69791  7391.6795           236.9984  638.0300078             NA           11.13402        NA         30463.618
# 1HP-3HP   11.17294    -5.0444            50.4252   -0.7037651             NA          -10.21152 -1221.045          7741.758
#         actv_tb_cost     art_cost  total_cost
# 3HP         3894.302 1495392.5635 1522008.727
# 1HP         2027.951 1495676.5944 1528168.164
# 1HP-3HP    -1866.351     284.0309    6159.437

results = treatment.details ( params )
#results = treatment.sim( params )

#calculate the icer, primarily for treatment.sim() execution
ICER = (results["Cost.1HP.total"] - results["Cost.3HP.total"]) / (results["DALY.3HP.total"] - results["DALY.1HP.total"])
names(ICER) = c()
print(ICER)

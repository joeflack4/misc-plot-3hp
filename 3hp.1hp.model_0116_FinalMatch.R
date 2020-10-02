library(expm)

paramAddCosts = function ( mp ) {
  #Those that depend on the Stage1 values
  mp$C_3HP_complete = ((0.02*3)+(mp$C_rif*6))*12+(mp$C_outpatient*4)
  mp$C_3HP_incomplete = ((0.02*3)+(mp$C_rif*6))*4+(mp$C_outpatient*1)
  #One outpatient visit or two for 1HP complete?
  mp$C_1HP_complete = ((0.02*1)+(mp$C_rif*4))*28+(mp$C_outpatient*2)
  mp$C_1HP_incomplete = mp$C_1HP_complete/2
  mp
}

model.param = data.frame ( 
  # Epidemiologic and health system values
  Pre_LTBI = 0.261,   
  Disengagement_HIVcare = 0.108,
  P_3HPtreatment = 0.74,
  E_3HP = 0.9,   
  
  #TODO set these to sane starting values
  E_1HP = 0.90,
  P_1HPtreatment = 0.94,
  
  #Prevalence of nonlethal adverse event during prev. therapy
  Pre_nonlth_prev_ther = 0.034,
  # Mortality
  m_LTBI_onART = 0.0354,
  m_LTBI_offART = 0.1326,
  m_ATB_onART = 0.1,
  m_ATB_OffART = 0.81,
  
  # Morbidity
  R_reactivation_HIV_offART = 0.043,
  R_reactivation_onART = 0.35,
  P_adverse = 0.034,
  
  #Disability weights
  disability1 = 0.078, #on ART LTBI
  disability2 = 0.408, #ATB 
  disability3 = 0.582, #off ART LTBI
  disability4 = 0.408, #ATB
  
  # Costs
  C_rif = 0.231,
  C_outpatient = 2.76,
  C_ART= 192,
  C_ATB= 231
)  

base_transmat <- function ( mp ) {
  #Create the base Transmat matrix, using the values that are common
  #for each instance of the Transmat matrix
  Transmat<-array(0, dim=c(5,5))
  Transmat[1,3]<-mp$Disengagement_HIVcare
  Transmat[1,5]<-mp$m_LTBI_onART
  Transmat[2,1]<-1-mp$m_ATB_onART
  Transmat[2,5]<-mp$m_ATB_onART
  Transmat[2,2]<-1-(Transmat[2,1])-(Transmat[2,5])
  Transmat[3,5]<-mp$m_LTBI_offART
  Transmat[4,3]<-1-mp$m_ATB_OffART
  Transmat[4,5]<-mp$m_ATB_OffART
  Transmat[4,4]<-1-Transmat[4,3]-Transmat[4,5]
  Transmat[5,5]<-1
  Transmat
}

##1.  LY 3HP

transmat_1_LY_3HP <- function ( mp ) {
  #Create the specific Transmat matrix relevant to the 1. LY 3HP scenario
  base <- base_transmat( mp )
  base[1,2]<-mp$R_reactivation_onART*mp$R_reactivation_HIV_offART*(1-mp$E_3HP)
  base[1,1]<-1-(base[1,2])-(base[1,3])-(base[1,5])
  base[3,4]<-mp$R_reactivation_HIV_offART*(1-mp$E_3HP)
  base[3,3]<-1-(base[3,4])-(base[3,5])
  base
}


##1.  LY 1HP
transmat_1_LY_1HP <- function ( mp ) {
  #Create the specific Transmat matrix relevant to the 1. LY 1HP scenario
  base <- base_transmat( mp )
  base[1,2]<-mp$R_reactivation_onART*mp$R_reactivation_HIV_offART*(1-mp$E_1HP)
  base[1,1]<-1-(base[1,2])-(base[1,3])-(base[1,5])
  base[3,4]<-mp$R_reactivation_HIV_offART*(1-mp$E_1HP)
  base[3,3]<-1-(base[3,4])-(base[3,5])
  base
}

##2. LYs_3HP & 2. LYs_1HP

transmat_2_LY_3HP_1HP <- function ( mp ) {
  #Create the specific Transmat matrix relevant to the 2. LY 3HP & 2. LY 1HP scenario
  base <- base_transmat( mp )
  base[1,2]<-mp$R_reactivation_onART * mp$R_reactivation_HIV_offART
  base[1,1]<-1-(base[1,2])-(base[1,3])-(base[1,5])
  base[3,4]<-mp$R_reactivation_HIV_offART
  base[3,3]<-1-(base[3,4])-(base[3,5])
  base
}

##3&4 LYs_3HP & 3&4 LYs_1HP

transmat_3_4_LY_3HP_1HP <- function ( mp ) {
  #Create the specific Transmat matrix relevant to the 3&4. LY 3HP & 3&4. LY 1HP scenario  
  base <- base_transmat( mp )
  base[1,2]<-0
  base[1,1]<-1-(base[1,2])-(base[1,3])-(base[1,5])
  base[3,4]<-0
  base[3,3]<-1-(base[3,4])-(base[3,5])
  base 
}

scenario <- function ( initial_val, Transmat ) {
  ##Function now takes the initial value and the value of the Transmat matrix,
  ##which is generated using the above transmat_* functions.
  
  ## create a vector containing transition data
  statemat <- matrix(c(initial_val,0,0,0,0), byrow=T, ncol=5)
  
  #Add the first results that rely on disengagment from HIV care
  #The initial row of statemat ([1,]) is the initial multiplication base
  statemat <- rbind (statemat, 
                     statemat[1,]%*%Transmat, 
                     statemat[1,]%*%(Transmat%^%2), 
                     statemat[1,]%*%(Transmat%^%3))
  
  # Disengagement from HIV care is 0 from year 4
  Transmat[1,3] <- 0
  Transmat[1,1] <- 1-(Transmat[1,2])-(Transmat[1,3])-(Transmat[1,5]) 
  
  #The final state entered at this point is the new base for multiplication
  new_base <- nrow(statemat)
  
  #Calculate the remaining years
  for (i in 1:17) {
    statemat <- rbind (statemat, statemat[new_base,]%*%(Transmat%^%(i)))
  }

  statemat
}

DALY_calc <- function ( statemat_param, mp ) {
  
  #The first line of statemat_param is year 0, remove it

  statemat = statemat_param[-1,]
  
  #Total Life Years (Undiscounted)
  TotalLY<-array(NA, c(20,1))
  for (j in 1:20){
    TotalLY[j,1]<-(statemat[j,1]+statemat[j,2]+statemat[j,3]+statemat[j,4])
  }

  #Discount factor
  d_factor = 1.03^-(1:20)

  #Total life years discounted
  TotalLY_DC = sum(d_factor * TotalLY)
  
  #Total Deaths
  Total_Death = statemat[20,5]

  #Years of life with disability
  
  #DALY calculation
  YLD.state<-array(NA, c(20,4))
  for (j in 1:20){
    YLD.state[j,1]<-statemat[j,1]*mp$disability1
    YLD.state[j,2]<-statemat[j,2]*mp$disability2
    YLD.state[j,3]<-statemat[j,3]*mp$disability3
    YLD.state[j,4]<-statemat[j,4]*mp$disability4
  }
  
  #colnames(YLD.state)<-c("YLD.onARTLTB","YLD.onARTATB", "YLD.offARTLTB", "YLD.offARTATB")
  YLD.total<-sum(YLD.state[,1]+YLD.state[,2]+YLD.state[,3]+YLD.state[,4])

  #Death.state requires year 0, so we're using statemat_param instead of statemat  
  Death.state<-array(NA, c(20,4))
  for (j in 1:20){
    Death.state[j,1]<-statemat_param[j,1]*mp$m_LTBI_onART
    Death.state[j,2]<-statemat_param[j,2]*mp$m_ATB_onART
    Death.state[j,3]<-statemat_param[j,3]*mp$m_LTBI_offART
    Death.state[j,4]<-statemat_param[j,4]*mp$m_ATB_OffART
  }

  Annualized<-array(NA, c(20,1))
  for(j in 1:20){
    Annualized[j,1]<-(1+0.03)^-j
  }
  
  AnnualizedYLL<-array(NA, c(20,1))
  for(j in 1:20){
    AnnualizedYLL[j,1]<-sum(Annualized[j:20,1])
  }

  
  YLL.state<-array(NA, c(20,4))
  for (j in 1:20){
    YLL.state[,1]<-Death.state[,1]*AnnualizedYLL
    YLL.state[,2]<-Death.state[,2]*AnnualizedYLL
    YLL.state[,3]<-Death.state[,3]*AnnualizedYLL
    YLL.state[,4]<-Death.state[,4]*AnnualizedYLL
  }
  
    YLL.total<-sum(YLL.state[,1]+YLL.state[,2]+YLL.state[,3]+YLL.state[,4])

  DALY.state<-array(NA, c(20,1))
  for (j in 1:20){
    DALY.state[j,1]<-sum(YLL.state[j,1:4]+YLD.state[j,1:4])
  }

  c(sum(DALY.state), TotalLY_DC, Total_Death)
}

COST_calc <- function ( statemat_param, mp ) {
  
  #The first line of statemat_param is year 0, remove it
  statemat = statemat_param[-1,]
  
  #discounting costs, same as in DALY above
  d_factor = 1.03^-(1:20)

  #Cost calculation
  Cost.state<-array(NA, c(20,5))
  for (j in 1:20){
    Cost.state[j,1]<-statemat[j,1]*mp$C_ART
    Cost.state[j,2]<-statemat[j,2]*mp$C_ART
    Cost.state[j,3]<-statemat[j,2]*mp$C_ATB
    Cost.state[j,4]<-statemat[j,3]*0
    Cost.state[j,5]<-statemat[j,4]*mp$C_ATB
  }

  #Apply the discount factor to the results
  c(sum(Cost.state[,1] * d_factor),
    sum(Cost.state[,2] * d_factor),
    sum(Cost.state[,3] * d_factor),
    sum(Cost.state[,4] * d_factor),   
    sum(Cost.state[,5] * d_factor))
}

decision_tree <- function ( mp ) {
  total_population = 1000
  
  #Cells will be ordered left to right, top to bottom
  
  #Column C
  #Positive LTBI
  cell1 = total_population * mp$Pre_LTBI
  #Negative LTBI
  cell2 = total_population * ( 1 - mp$Pre_LTBI )
  
  #Column E, Top 2
  #Adverse Event <- Positive LTBI
  cell3 = cell1 * mp$Pre_nonlth_prev_ther
  #No Adverse Event <- Positive LTBI
  cell4 = cell1 * (1 - mp$Pre_nonlth_prev_ther)
  
  #Column E, Bottom 2
  #Adverse Event <- Negative LTBI
  cell5 = cell2 * mp$Pre_nonlth_prev_ther
  #No Adverse Event <- Negative LTBI
  cell6 = cell2 * (1 - mp$Pre_nonlth_prev_ther)
  
  #Column G, 1-2
  #Survive <- Adverse Event <- Positive LTBI
  cell7 = cell3 * (1 - mp$m_LTBI_onART)
  #Death <- Adverse Event <- Positive LTBI
  #Cell may never be used, worth calculating?
  cell8 = cell3 * mp$m_ATB_onART
  
  #Column G, 3-4
  #3HP Complete Therapy <- No Adverse Event <- Positive LTBI
  cell9a = cell4 * mp$P_3HPtreatment
  #1HP Complete Therapy <- No Adverse Event <- Positive LTBI
  cell9b = cell4 * mp$P_1HPtreatment
  #3HP Lost To Follow up <- No Adverse Event <- Positive LTBI
  cell10a = cell4 * (1 - mp$P_3HPtreatment)
  #1HP Lost To Follow up <- No Adverse Event <- Positive LTBI
  cell10b = cell4 * (1 - mp$P_1HPtreatment)
  
  #Column G, 5-6
  #Survive <- Adverse Event <- Negative LTBI
  cell11 = cell5 * (1 - mp$m_LTBI_onART)
  #Death <- Adverse Event <- Negative LTBI
  #Cell may never be used, worth calculating?
  cell12 = cell5 * mp$m_LTBI_onART
  
  #Column G, 7-8
  #3HP Complete Therapy <- No Adverse Event <- Negative LTBI
  cell13a = cell6 * mp$P_3HPtreatment
  #1HP Complete Therapy <- No Adverse Event <- Negative LTBI
  cell13b = cell6 * mp$P_1HPtreatment
  #3HP Lost to Follow up <- No Adverse Event <- Negative LTBI
  cell14a = cell6 * (1 - mp$P_3HPtreatment)
  #1HP Lost to Follow up <- No Adverse Event <- Negative LTBI
  cell14b = cell6 * (1 - mp$P_1HPtreatment)
  
  
  
  #1. LY 3HP = cell9a (3HP complete)
  #1. LY 1HP = cell9b (1HP complete)
  #2. LY 3HP = cell7 + cell10a (3HP incomplete)
  #2. LY 1HP = cell7 + cell10b (1HP incomplete)
  #3. LY 3HP = cell11 + cell14a (3HP negative incompelte)
  #3. LY 1HP = cell11 + cell14b (1HP negative incomplete)
  #4. LY 3HP = cell13a (3HP negative complete)
  #4. LY 1HP = cell13b (1HP negative complete)
  
  Cost.3HP.1 = cell9a * mp$C_3HP_complete  #(3HP complete) 
  Cost.1HP.1 = cell9b * mp$C_1HP_complete   #(1HP complete) 
  Cost.3HP.2 = (cell7 + cell10a)*mp$C_3HP_incomplete  #(3HP incompelte) 
  Cost.1HP.2 = (cell7 + cell10b)*mp$C_1HP_incomplete   #(1HP incomplete) 
  Cost.3HP.3 = (cell11 + cell14a)*mp$C_3HP_incomplete  #(3HP negative incompelte) 
  Cost.1HP.3 = (cell11 + cell14b)*mp$C_1HP_incomplete  #(1HP negative incompelte) 
  Cost.3HP.4 = cell13a*mp$C_3HP_complete #(3HP negative compelte)
  Cost.1HP.4 = cell13b*mp$C_1HP_complete  #(1HP negative complete)

  array( c( cell9a, cell9b, Cost.3HP.1, Cost.1HP.1,
            cell7 + cell10a, cell7 +cell10b, Cost.3HP.2, Cost.1HP.2,
            cell11 + cell14a, cell11 + cell14b, Cost.3HP.3, Cost.1HP.3,
            cell13a, cell13b, Cost.3HP.4, Cost.1HP.4), c(4,4))
  
}

treatment.details<-function( mp ) {

  #make sure the secondary costs are calculated before proceeding
  mp = paramAddCosts(mp)
  
  return_frame = data.frame( total_ly = c(NA,NA,NA), 
                             total_daly = c(NA,NA,NA), 
                             total_reg_complete = c(NA,NA,NA), 
                             total_deaths = c(NA,NA,NA), 
                             total_tb_death = c(NA,NA,NA), 
                             total_active_cases = c(NA,NA,NA), 
                             ICER = c(NA,NA,NA), 
                             #Costs
                             prev_therapy_cost = c(NA,NA,NA), 
                             actv_tb_cost = c(NA,NA,NA), 
                             art_cost = c(NA,NA,NA), 
                             total_cost = c(NA,NA,NA),
                             row.names=c("3HP","1HP","1HP-3HP"))
  
  DT_Results = decision_tree ( mp )
  
  #DT_Results is a array of 4 rows and 4 columns; the first row is
  #3HP and the second row is 1HP
  
  return_frame['3HP','total_reg_complete'] = DT_Results[1,1]
  return_frame['1HP','total_reg_complete'] = DT_Results[2,1]
  return_frame['1HP-3HP','total_reg_complete'] = DT_Results[2,1] - DT_Results[1,1]

  Model.3HP = DT_Results[1,]
  Model.1HP = DT_Results[2,]
  DTcost.3HP = DT_Results[3,]
  DTcost.1HP = DT_Results[4,]
  
  return_frame['3HP','prev_therapy_cost'] = sum(DTcost.3HP)
  return_frame['1HP','prev_therapy_cost'] = sum(DTcost.1HP)
  return_frame['1HP-3HP','prev_therapy_cost'] = sum(DTcost.1HP) - sum(DTcost.3HP)
  
  Transmat_3HP = list( transmat_1_LY_3HP( mp ), 
                       transmat_2_LY_3HP_1HP( mp ), 
                       transmat_3_4_LY_3HP_1HP( mp ), 
                       transmat_3_4_LY_3HP_1HP( mp ) )
  
  Transmat_1HP = list( transmat_1_LY_1HP( mp ),
                       transmat_2_LY_3HP_1HP( mp ), 
                       transmat_3_4_LY_3HP_1HP( mp ), 
                       transmat_3_4_LY_3HP_1HP( mp ) )
  DALY.1HP = c()
  DALY.3HP = c()
  COST.1HP = c()
  COST.3HP = c()
  
  total.active.cases.3HP = 0.0
  total.active.cases.1HP = 0.0
  cost.treating.active.3HP = 0.0
  cost.treating.active.1HP = 0.0
  cost.art.3HP = 0.0
  cost.art.1HP = 0.0
  tot.disc.ly.1HP = 0.0
  tot.disc.ly.3HP = 0.0
  tot.death.1HP = 0.0
  tot.death.3HP = 0.0
  
  for (i in 1:4) {

    statemat_3HP = scenario ( Model.3HP[i], Transmat_3HP[[i]] )
    
    statemat_1HP = scenario ( Model.1HP[i], Transmat_1HP[[i]] )
    
    total.active.cases.3HP = total.active.cases.3HP + 
      sum(statemat_3HP[,2]) + 
      sum(statemat_3HP[,4])
    total.active.cases.1HP = total.active.cases.1HP + 
      sum(statemat_1HP[,2]) + 
      sum(statemat_1HP[,4])

    daly.calc.1HP = DALY_calc(statemat_1HP, mp)
    daly.calc.3HP = DALY_calc(statemat_3HP, mp)
    
    tot.disc.ly.3HP = tot.disc.ly.3HP + daly.calc.3HP[2]
    tot.disc.ly.1HP = tot.disc.ly.1HP + daly.calc.1HP[2]
    
    tot.death.3HP = tot.death.3HP + daly.calc.3HP[3]
    tot.death.1HP = tot.death.1HP + daly.calc.1HP[3]
    
    cost.calc.1HP = COST_calc( statemat_1HP, mp )
    cost.calc.3HP = COST_calc( statemat_3HP, mp )
    
    cost.treating.active.3HP = cost.treating.active.3HP + cost.calc.3HP[3] + cost.calc.3HP[5]
    cost.treating.active.1HP = cost.treating.active.1HP + cost.calc.1HP[3] + cost.calc.1HP[5]
    
    cost.art.3HP = cost.art.3HP + cost.calc.3HP[1] + cost.calc.3HP[2]
    cost.art.1HP = cost.art.1HP + cost.calc.1HP[1] + cost.calc.1HP[2]
    
    DALY.1HP = c(DALY.1HP, daly.calc.1HP[1])
    DALY.3HP = c(DALY.3HP, daly.calc.3HP[1])
    COST.1HP = c(COST.1HP, sum(cost.calc.1HP))
    COST.3HP = c(COST.3HP, sum(cost.calc.3HP))
  }
  
  #Total active Cases
  return_frame['3HP','total_active_cases'] = total.active.cases.3HP
  return_frame['1HP','total_active_cases'] = total.active.cases.1HP
  return_frame['1HP-3HP','total_active_cases'] = total.active.cases.1HP - total.active.cases.3HP
  
  #Cost of active Cases
  return_frame['3HP','actv_tb_cost'] = cost.treating.active.3HP
  return_frame['1HP','actv_tb_cost'] = cost.treating.active.1HP
  return_frame['1HP-3HP','actv_tb_cost'] = cost.treating.active.1HP - cost.treating.active.3HP

  #Cost of ART
  return_frame['3HP','art_cost'] = cost.art.3HP
  return_frame['1HP','art_cost'] = cost.art.1HP
  return_frame['1HP-3HP','art_cost'] = cost.art.1HP - cost.art.3HP
  
  #Total Discounted Life Years
  return_frame['3HP','total_ly'] = tot.disc.ly.3HP
  return_frame['1HP','total_ly'] = tot.disc.ly.1HP
  return_frame['1HP-3HP','total_ly'] = tot.disc.ly.1HP - tot.disc.ly.3HP
  
  #Total Deaths
  return_frame['3HP','total_deaths'] = tot.death.3HP
  return_frame['1HP','total_deaths'] = tot.death.1HP
  return_frame['1HP-3HP','total_deaths'] = tot.death.1HP - tot.death.3HP
  
  #the probability of achieving cure after a single treatment course (default value = 50%) 
  DALY.1HP.total<-sum(DALY.1HP)
  #the cumulative probability of being cured after the second course of treatment, given that the first course of treatment was unsuccessful (default value = 75%)
  DALY.3HP.total<-sum(DALY.3HP)
  
  #the average cost-per-week of the first course of treatment
  Cost.1HP.total<-sum(COST.1HP)+sum(DTcost.1HP)
  Cost.3HP.total<-sum(COST.3HP)+sum(DTcost.3HP)
  
  return_frame['3HP','total_daly'] = DALY.3HP.total
  return_frame['1HP','total_daly'] = DALY.1HP.total
  return_frame['1HP-3HP','total_daly'] = DALY.1HP.total - DALY.3HP.total
  
  return_frame['3HP','total_cost'] = Cost.3HP.total
  return_frame['1HP','total_cost'] = Cost.1HP.total
  return_frame['1HP-3HP','total_cost'] = Cost.1HP.total - Cost.3HP.total
  
  return_frame['1HP-3HP','ICER'] = return_frame['1HP-3HP','total_cost'] /
                                   return_frame['1HP-3HP','total_daly']
  
  results<-c(Cost.1HP.total, DALY.1HP.total, Cost.3HP.total, DALY.3HP.total)
  names(results)<-c("Cost.1HP.total", "DALY.1HP.total","Cost.3HP.total", "DALY.3HP.total")
  
  print (return_frame)
  
  results
}

treatment.sim<-function( mp ){
  
  mp = paramAddCosts(mp)
  
  DT_Results = decision_tree ( mp )
  
  #DT_Results is a array of 4 rows and 4 columns; the first row is
  #3HP and the second row is 1HP
  
  Model.3HP = DT_Results[1,]
  Model.1HP = DT_Results[2,]
  DTcost.3HP = DT_Results[3,]
  DTcost.1HP = DT_Results[4,]
  
  
  Transmat_3HP = list( transmat_1_LY_3HP( mp ), 
                       transmat_2_LY_3HP_1HP( mp ), 
                       transmat_3_4_LY_3HP_1HP( mp ), 
                       transmat_3_4_LY_3HP_1HP( mp ) )
  
  Transmat_1HP = list( transmat_1_LY_1HP( mp ),
                       transmat_2_LY_3HP_1HP( mp ), 
                       transmat_3_4_LY_3HP_1HP( mp ), 
                       transmat_3_4_LY_3HP_1HP( mp ) )
  DALY.1HP = c()
  DALY.3HP = c()
  COST.1HP = c()
  COST.3HP = c()
  
  for (i in 1:4) {
  
    statemat_3HP = scenario ( Model.3HP[i], Transmat_3HP[[i]] )
    
    statemat_1HP = scenario ( Model.1HP[i], Transmat_1HP[[i]] )
    
    
    DALY.1HP = c(DALY.1HP, DALY_calc ( statemat_1HP, mp )[1])
    DALY.3HP = c(DALY.3HP, DALY_calc ( statemat_3HP, mp )[1])
    
    COST.1HP = c(COST.1HP, sum(COST_calc ( statemat_1HP, mp )))
    COST.3HP = c(COST.3HP, sum(COST_calc ( statemat_3HP, mp )))
  }
  
  #the probability of achieving cure after a single treatment course (default value = 50%) 
  DALY.1HP.total<-sum(DALY.1HP)
  #the cumulative probability of being cured after the second course of treatment, given that the first course of treatment was unsuccessful (default value = 75%)
  DALY.3HP.total<-sum(DALY.3HP)

  #the average cost-per-week of the first course of treatment
  Cost.1HP.total<-sum(COST.1HP)+sum(DTcost.1HP)
  Cost.3HP.total<-sum(COST.3HP)+sum(DTcost.3HP)

  results<-c(Cost.1HP.total, DALY.1HP.total, Cost.3HP.total, DALY.3HP.total)
  names(results)<-c("Cost.1HP.total", "DALY.1HP.total","Cost.3HP.total", "DALY.3HP.total")
  
  results
}

####END OF SIMULATION#####
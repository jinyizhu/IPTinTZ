#include <Rcpp.h>
#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <cerrno>
#include <cfenv>
#include <cstring>
using namespace Rcpp;

// [[Rcpp::export]]
// Function to convert rates to probabilities
double RateToProb(double rate, double t){
  return(1-exp(-rate*t));
}

// [[Rcpp::export]]
// Function to convert probabilities to rates
double ProbToRate(double prob, double t){
  return(-log(1-prob)/t);
}

// [[Rcpp::export]]
// Function to convert uniform distribution to truncated exponential
double UnifToExpo(double u, double lambda, double a){
  return( log(1+ u * exp(-lambda * a) -u ) / (-lambda) );
}

// [[Rcpp::export]]
List sim_pop_analysis(int nPop, 
             NumericVector age_vec, NumericVector male_vec, //Age and sex
             NumericVector cd4base_vec,  //Baseline CD4
             NumericVector prevTB_vec, //Baseline Active TB status 
             NumericVector enr_mo_vec, //Enrollment month
             NumericVector pdie_m, NumericVector pdie_f, // Life tables
             NumericVector LE_reference, //reference LE
             NumericVector RR_agemix_vec, //Age mixing matrix
             NumericVector agegroup_cali, NumericVector cd4group_cali, //Age category for calibration
             NumericVector IPTmo_vec, //IPT start month
             double HIVmort_vec_onART_par1, double HIVmort_vec_onART_par2, //Parameters for HIV mortality rates on ART
             double HIVmort_vec_offART_par1, double HIVmort_vec_offART_par2, //Parameters for HIV mortality rates off ART
             double prop_TB_HIVmort, //Proporiton of HIV mortality due to TB
             double mult_mort_3mo, //mortality multiplier for the first three months
             double LTFtheta1, double LTFtheta2, //LTFU parameters
             double Asym, double R0, double lrc, //CD4 rebound parameters
             double cd4_rate_postLTF, //CD4 declining rate post LTFU
             double TBmort, //TB mortality
             double phi_vec_par1, double phi_vec_par2, double phi_vec_par3, //TB reactivation
             double gamma_vec_par1, double gamma_vec_par2, //TB reactivation
             double TB_FoI, //TB force of infection
             double IRR_LTBI, //incidence rate ratio for infected
             double p_Detect, //probability of detecting TB each time
             double p_UndetectedTB, //Probability of undetected TB at CD4 = 0
             double p_UndetectedTB_scale, //Scaling factor for the probability of undetected TB
             double p_LatentTB, //Probability of latent TB at baseline
             double IPT_eff, //IPT efficacy
             double prob_IPTstart_new_post2019, double prob_IPTstart_old_pre2019, double prob_IPTstart_old_post2019, 
             double prob_IPTcomplete_2014, double prob_IPTcomplete_inc, double prob_IPTcomplete_2018,
             double p_Fail_pan, double p_Fail_INHR, double rate_INHR, double rate_Cure_IPT, double prev_INHR,
             double prob_trueLTF = 1,
             double tInf_scale = 0.8,
             int IPT_len = 6,
             double DW_HIV_onART = 0.078, double DW_HIV_offART_AIDS = 0.581, 
             double DW_HIV_offART_noAIDS = 0.274, double DW_TB = 0.4089,
             double c_ART_mo = 22.06, double c_TBtreat_mo = 83.33, double c_IPT_mo = 0.564,
             double disc_r = 0.03
){
  
  //Common Random Numbers system
  int nRN = nPop*12;
  NumericVector randomnums_NV = runif(nRN);
  int rini = 0; //initial position in the random number vector
  int tmax = 1000;
  //types of events
  int r_bgdeath = 0;
  int r_HIVdeath = 1;
  int r_TBdeath = 2;
  int r_LTF = 3;
  int r_TBinf = 4;
  int r_progression = 5;
  int r_detection = 6;
  int r_IPT = 7;
  int r_cure_IPT = 8;
  int r_cure_TBTx = 9;
  int r_INHR_IPT = 10;
  
  std::vector<double> randomnums(nRN);
  for(int i = 0; i < nRN; i++){randomnums[i] = randomnums_NV[i];}
  
  //Parameters
  double prop_HIVmortExTB = 1 - prop_TB_HIVmort; //proportion of HIV mortality NOT due to TB
  double clength = 1.0/12.0; //Cycle length
  double pdie_HIV = 0; //probability of dying from HIV
  int length_TbTx = 6; //length of TB treatment
  
  double pdie[2][100]; //lifetables -- input was monthly probability
  for(int i = 0; i < 100; i++){ 
    pdie[0][i] = pdie_f[i];
    pdie[1][i] = pdie_m[i];
  }

  //Output
  NumericVector TBcases_by_year(100);
  NumericVector TBdeaths_by_year(100);
  NumericVector deaths_by_year(100);
  NumericVector TBcases_by_cal(100);
  NumericVector TBdeaths_by_cal(100);
  NumericVector deaths_by_cal(100);
  NumericVector TBcases_by_year_INHR(100);
  std::vector<double> LE_vec(nPop, 0); //Life expectancy
  std::vector<double> tobs_vec(nPop, 0); //time observed in model
  std::vector<bool> lost_vec(nPop, 0); // lost to follow-up
  std::vector<bool> dead_vec(nPop, 0); //dead at the end of model -- should all be 1?
  std::vector<bool> dead_by_2018(nPop, 0); //person dead by 2018
  std::vector<bool> TBcase(nPop, 0); //1st TB case detection
  std::vector<int> FU_TBcase (nPop,-999); //TB case detection time
  std::vector<bool> TBdeath(nPop, 0); //TB death
  std::vector<int> FU_TBdeath (nPop,-999); //TB death FU time
  std::vector<bool> HIVdeath(nPop, 0); //HIV death
  std::vector<double> tInf(nPop, -999); //time since infection at baseline
  std::vector<bool> UndetectedTB_base(nPop, 0); //baseline undetected TB case
  std::vector<bool> LatentTB_base(nPop, 0); //baseline undetected latent TB
  std::vector<int> LTBI_time(nPop, 999); //LTBI time
  int n_newcases = 0; //number of new cases
  int n_recentInf = 0; //number of recent infection
  std::vector<double> DALY_vec(nPop, 0); //DALY
  std::vector<double> cART_vec(nPop, 0); //cost of ART per patient
  std::vector<double> cTB_vec(nPop, 0); //cost of TB per patient
  std::vector<double> cIPT_vec(nPop, 0); //cost of IPT per patient
  std::vector<int> tIPT_ART_vec(nPop, -999); //IPT time since ART
  std::vector<double> LE_disc_vec(nPop, 0); //LE disc
  std::vector<double> DALY_disc_vec(nPop, 0); //DALY disc
  std::vector<double> cART_disc_vec(nPop, 0); //cost of ART per patient disc
  std::vector<double> cTB_disc_vec(nPop, 0); //cost of TB per patient disc
  std::vector<double> cIPT_disc_vec(nPop, 0); //cost of IPT per patient disc
  std::vector<int> TBcase_true(nPop, 0); //True TB case
  std::vector<bool> TBdeath_true(nPop, 0); //True TB death

  
  //looping over population
  for(int id = 0; id < nPop; id++){
    
    //Read in individual characteristics
    double age = age_vec[id];
    int male = male_vec[id];
    double cd4base = cd4base_vec[id];
    bool prevTB = prevTB_vec[id];
    double p_LatentTB_i = p_LatentTB; //Latent TB
    int IPTmo = IPTmo_vec[id];
    int enr_mo = enr_mo_vec[id];
    
    //Health states
    bool ARTstate = 1; //on/off ART
    std::string TBstate("S");
    
    //IPT related variables
    bool IPT = 0;
    bool IPT_complete = 0;
    bool IPT_cured = 0;
    int  tIPT = 0;
    bool IPThist = 0;
    bool INHR = 0;
    double IPT_eff_id = IPT_eff;
    
    //Define trackers and output variables
    double modelage = age; //age in the simulation
    double modelcd4 = cd4base; //CD4 in the simulation
    double LTFprob = 0; //LTF prob
    bool lost = 0; //LTF
    int artmon = 1; //time since ART initiation
    bool dead_obs = 0; //observed death
    double tobs = 0; //observed time (until death or LTF)
    bool TBhist = 0; //TB history
    double t_TBinfection = 0; // time since last TB infection
    int t_TBtreat = 0;
    double mmult = mult_mort_3mo; //mortality multiplier
    double prob_IPTstart = 0; //prob IPT start per month
    
    //initiate cohort's TB status
    //INH resistance
    ++rini;
    if(randomnums[rini] <= prev_INHR){
      INHR = 1;
      IPT_eff_id = 1;
    }
    double prob_IPTcomplete = prob_IPTcomplete_2014 + enr_mo*prob_IPTcomplete_inc;
    if(prob_IPTcomplete > prob_IPTcomplete_2018){prob_IPTcomplete = prob_IPTcomplete_2018;}
    ++rini;
    if(randomnums[rini] <= prob_IPTcomplete){ //For those who can completed IPT
      IPT_complete = 1;
      ++rini;
      if(randomnums[rini] <= IPT_eff_id){
        IPT_cured = 1; //Cured by IPT or not
      }
    }else{
      ++rini;
      if(randomnums[rini] <= IPT_eff_id*IPTmo/6){
        IPT_cured = 1; //Cured by IPT or not
      }
    }
    
    //Prev TB
    if(prevTB == 1){
      TBstate = "I_treating";
      TBhist = 1;
      t_TBtreat = 1;
      LTBI_time[id] = 0;
      TBcase_true[id] = 1;
    }
    
    //Undetected TB
    if(TBstate == "S"){
      ++rini;
      if(randomnums[rini] <= p_UndetectedTB*exp(- p_UndetectedTB_scale * cd4base)){
        TBcase_true[id] = 1;
        TBstate = "I_undetected";
        UndetectedTB_base[id] = 1;
        LTBI_time[id] = 0;
      }
    }
    
    //Latent TB
    if(TBstate == "S"){
      ++rini;
      if(randomnums[rini] <= p_LatentTB_i){
        TBstate = "L";
        ++rini;
        t_TBinfection = tInf_scale * UnifToExpo(randomnums[rini], TB_FoI, modelage); //scaling factor
        LatentTB_base[id] = 1;
        LTBI_time[id] = 0;
      }
    }
    
    
    
    //run simulations
    while(modelage < 100){
      double DW_HIV_offART = DW_HIV_offART_noAIDS;
      if(modelcd4 < 200){
        DW_HIV_offART = DW_HIV_offART_AIDS;
      }
      
      //update age and mmult
      int ageint = (int) modelage;
      if(artmon > 3){
        mmult = 1;
      }
      
      //update IPT start prob
      if(enr_mo >= 61){ // post-2019
        prob_IPTstart = prob_IPTstart_new_post2019;
      }else{
        if(enr_mo + artmon - 1 >= 61){
          prob_IPTstart = prob_IPTstart_old_post2019;
        }else{
          prob_IPTstart = prob_IPTstart_old_pre2019;
        }
      }
      
      //Update IPT efficacy
      if(INHR == 1){
        IPT_eff_id = 1;
        IPT_cured = 0;}
      
      //Record cycle rewards
      LE_vec[id] += clength;
      LE_disc_vec[id] += clength * pow(1/(1 + disc_r), (artmon + enr_mo - 2)/12);

      if(ARTstate == 1){ //if on ART
        cART_vec[id] += c_ART_mo;
        cART_disc_vec[id] += c_ART_mo * pow(1/(1 + disc_r), (artmon + enr_mo - 2)/12);
        
        if(TBstate == "I_undetected"){
          DALY_vec[id] += DW_TB * clength;
          DALY_disc_vec[id] += DW_TB * clength * pow(1/(1 + disc_r), (artmon + enr_mo - 2)/12);
          if(IPT== 1){
            cIPT_vec[id] += c_IPT_mo;
            cIPT_disc_vec[id] += c_IPT_mo * pow(1/(1 + disc_r), (artmon + enr_mo - 2)/12);
          }
        }else if(TBstate == "I_treating" ){
          DALY_vec[id] += DW_TB * clength;
          DALY_disc_vec[id] += DW_TB * clength * pow(1/(1 + disc_r), (artmon + enr_mo - 2)/12);
          cTB_vec[id] += c_TBtreat_mo;
          cTB_disc_vec[id] += c_TBtreat_mo * pow(1/(1 + disc_r), (artmon + enr_mo - 2)/12);
        }else{
          DALY_vec[id] += DW_HIV_onART * clength;
          DALY_disc_vec[id] += DW_HIV_onART * clength * pow(1/(1 + disc_r), (artmon + enr_mo - 2)/12);
          if(IPT== 1){
            cIPT_vec[id] += c_IPT_mo;
            cIPT_disc_vec[id] += c_IPT_mo * pow(1/(1 + disc_r), (artmon + enr_mo - 2)/12);
          }
        }
      }else{ //if not on ART
        DALY_vec[id] += DW_HIV_offART * clength;
        DALY_disc_vec[id] += DW_HIV_offART * clength * pow(1/(1 + disc_r), (artmon + enr_mo - 2)/12);
        if(TBstate == "I_LTF_treating"){
          cTB_vec[id] += c_TBtreat_mo;
          cTB_disc_vec[id] += c_TBtreat_mo * pow(1/(1 + disc_r), (artmon + enr_mo - 2)/12);
        }
      }
      
      //1. Die from background mortality?
      if(randomnums[ nPop*3-1 + id + r_bgdeath*tmax + artmon ]  < 
        RateToProb(ProbToRate(pdie[male][ageint], 1)*mmult, 1)){
        if(ARTstate == 1){
          dead_obs = 1; if(artmon + enr_mo - 1 <= 59){
            dead_by_2018[id] = 1;
          }
          tobs = modelage;
        }
        ++deaths_by_year(artmon/12);
        ++deaths_by_cal((enr_mo + artmon - 1)/12);
        break;}
      
      //2. Die from HIV?
      if(ARTstate == 1){ //on ART
        pdie_HIV = RateToProb((exp(HIVmort_vec_onART_par1*modelcd4 + HIVmort_vec_onART_par2)*
          prop_HIVmortExTB*mmult),
          1.0/12.0);
        if(randomnums[ nPop*3-1 + id + r_HIVdeath*tmax + artmon ]  < pdie_HIV){
          dead_obs = 1; if(artmon + enr_mo - 1 <= 59){dead_by_2018[id] = 1;}
          tobs = modelage;
          HIVdeath[id] = 1;
          ++deaths_by_year(artmon/12);
          ++deaths_by_cal((enr_mo + artmon - 1)/12);
          break;}
      }else{ //Off ART
        pdie_HIV = RateToProb(exp(HIVmort_vec_offART_par1*modelcd4 + HIVmort_vec_offART_par2)*
          prop_HIVmortExTB*mmult, 
          1.0/12.0);
        if(randomnums[ nPop*3-1 + id + r_HIVdeath*tmax + artmon ]  < pdie_HIV){
          HIVdeath[id] = 1;
          ++deaths_by_year(artmon/12);
          ++deaths_by_cal((enr_mo + artmon - 1)/12);
          break;}
      }
      
      //3. Die from TB?
      if(TBstate == "I_undetected"){ //Active TB undetected
        if(IPT == 1 && INHR == 0){ // on IPT
          if(randomnums[ nPop*3-1 + id + r_TBdeath*tmax + artmon ]  < RateToProb((6-tIPT)/6*TBmort*mmult, 1.0/12.0) ){
            dead_obs = 1; if(artmon + enr_mo - 1 <= 59){dead_by_2018[id] = 1;}
            tobs = modelage; 
            TBdeath[id] = 1;
            TBdeath_true[id] = 1;
            FU_TBdeath[id] = artmon;
            ++deaths_by_year(artmon/12);
            ++deaths_by_cal((enr_mo + artmon - 1)/12);
            ++TBdeaths_by_year(artmon/12);
            ++TBdeaths_by_cal((enr_mo + artmon - 1)/12);
            break;}
        }else{ // not on IPT
          if(randomnums[ nPop*3-1 + id + r_TBdeath*tmax + artmon ]  < RateToProb(TBmort*mmult, 1.0/12.0) ){
            dead_obs = 1; if(artmon + enr_mo - 1 <= 59){dead_by_2018[id] = 1;}
            tobs = modelage; 
            TBdeath[id] = 1;
            TBdeath_true[id] = 1;
            FU_TBdeath[id] = artmon;
            ++deaths_by_year(artmon/12);
            ++deaths_by_cal((enr_mo + artmon - 1)/12);
            ++TBdeaths_by_year(artmon/12);
            ++TBdeaths_by_cal((enr_mo + artmon - 1)/12);
            break;}
        }
        
        
      }else if(TBstate == "I_treating" || TBstate == "I_LTF_treating"){ //active TB on treatment
        if(randomnums[ nPop*3-1 + id + r_TBdeath*tmax + artmon ]  < RateToProb((6-t_TBtreat)/6*TBmort*mmult, 1.0/12.0) ){
          dead_obs = 1; if(artmon + enr_mo - 1 <= 59){dead_by_2018[id] = 1;}
          tobs = modelage;
          TBdeath[id] = 1;
          TBdeath_true[id] = 1;
          FU_TBdeath[id] = artmon;
          ++deaths_by_year(artmon/12);
          ++deaths_by_cal((enr_mo + artmon - 1)/12);
          ++TBdeaths_by_year(artmon/12);
          ++TBdeaths_by_cal((enr_mo + artmon - 1)/12);
          break;}
      }else if(TBstate == "I_LTF"){ //active TB, LTFU
        if(randomnums[ nPop*3-1 + id + r_TBdeath*tmax + artmon ]  < RateToProb(TBmort*mmult, 1.0/12.0) ){
          TBdeath_true[id] = 1;
          ++deaths_by_year(artmon/12);
          ++deaths_by_cal((enr_mo + artmon - 1)/12);
          ++TBdeaths_by_year(artmon/12);
          ++TBdeaths_by_cal((enr_mo + artmon - 1)/12);
          break;}
      }
      
      if(ARTstate == 1){//On ART
        
        //4. Lost to follow-up?
        LTFprob = RateToProb(prob_trueLTF*pow(LTFtheta1, -LTFtheta2)*LTFtheta2*pow(artmon, (LTFtheta2 - 1)), 1);
        if(randomnums[ nPop*3-1 + id + r_LTF*tmax + artmon ] < LTFprob){
          ARTstate = 0;
          lost = 1;
          tobs = modelage;
          if(TBstate == "I_undetected" || TBstate == "I_treating"){
            TBstate = "I_LTF";
          }
        }
        
        //5. Check TB status
        if(TBstate == "S"){ //Susceptible 
          //Susceptible state -> can get infected
          if(randomnums[ nPop*3-1 + id + r_TBinf*tmax + artmon ]  < RateToProb(TB_FoI * RR_agemix_vec[modelage/5] * (1-(1-IPT_eff_id)*IPT), 1.0/12.0) ){
            TBstate = "L";
            t_TBinfection = 1; // first infection
            LTBI_time[id] = artmon;
          }
          //On IPT --> 
          if(IPT == 1){
            //time to stop IPT?
            if((IPT_complete == 1 && tIPT == 6)||(IPT_complete == 0 && tIPT == IPTmo)){ 
              IPT = 0;
              if(IPT_cured == 1){TBstate = "TP";}
            }
          }
          //Not on IPT yet --> time to get on IPT?
          if(IPThist == 0 && artmon > 1){
            if(randomnums[ nPop*3-1 + id + r_IPT*tmax + artmon ] < prob_IPTstart){
              IPT = 1;
              tIPT_ART_vec[id] = artmon - 1;
              IPThist = 1;
            }
          }
          
        }else if(TBstate == "L"){
          //Latent state -> can get reinfected or progress to active TB
          //Progression to active TB
          if(randomnums[ nPop*3-1 + id + r_progression*tmax + artmon ]  < RateToProb( (1-(1-IPT_eff_id)*IPT)*
             (exp(phi_vec_par1)+exp(phi_vec_par2)*(exp(-exp(phi_vec_par3)*t_TBinfection)))*
             (1 + exp(gamma_vec_par1)*(exp(-exp(gamma_vec_par2)*modelcd4))), 1)){
            
            TBcase_true[id] ++;
            TBcases_by_year[artmon/12] ++;
            if(INHR == 1){TBcases_by_year_INHR[artmon/12] ++;}
            TBcases_by_cal[(enr_mo + artmon - 1)/12]++;
            
            //DETECTED!!!!
            if(randomnums[ nPop*3-1 + id + r_detection*tmax + artmon ]  < p_Detect ){
              TBstate = "I_treating";
              IPT = 0;
              t_TBtreat = 0;
              
              //First TB detection
              if(TBhist == 0){
                TBcase[id] = 1;
                FU_TBcase[id] = artmon;
                tInf[id] = t_TBinfection;
                TBhist = 1;
                if(UndetectedTB_base[id] == 0){
                  ++ n_newcases;
                  if(t_TBinfection <= 24){++ n_recentInf; }
                }
              }
            }else{
              TBstate = "I_undetected";
              //Not on IPT yet --> time to get on IPT?
              if(IPThist == 0 && artmon > 1){
                if(randomnums[ nPop*3-1 + id + r_IPT*tmax + artmon ] < prob_IPTstart){
                  IPT = 1;
                  tIPT_ART_vec[id] = artmon - 1;
                  IPThist = 1;
                }
              }
            }
            //Age + 1/12 to make the loop running
            
            modelcd4 = cd4base-R0+R0*exp(-exp(lrc)*artmon*365.25/12);
            modelage += clength;
            ++artmon;
            tIPT += IPT;
            continue;
          }
          
          //Not on IPT yet --> time to get on IPT?
          if(IPThist == 0 && artmon > 1){
            if(randomnums[ nPop*3-1 + id + r_IPT*tmax + artmon ] < prob_IPTstart){
              IPT = 1;
              tIPT_ART_vec[id] = artmon - 1;
              IPThist = 1;
            }
          }
          
          if(randomnums[ nPop*3-1 + id + r_TBinf*tmax + artmon ]  < RateToProb(TB_FoI *RR_agemix_vec[modelage/5] * 
             IRR_LTBI * (1-(1-IPT_eff_id)*IPT), 1.0/12.0) ){
            //TB reinfection
            t_TBinfection = 1; // reset time since infection
            //Age + 1/12 to make the loop running
            modelcd4 = cd4base-R0+R0*exp(-exp(lrc)*artmon*365.25/12);
            modelage += clength;
            ++artmon;
            tIPT += IPT;
            continue;
          }
          
          if((IPT== 1) && ((IPT_complete == 1 && tIPT == 6)||(IPT_complete == 0 && tIPT == IPTmo))){
            //complete/stop IPT
            if(IPT_cured == 1){TBstate = "TP"; t_TBinfection = 0;}
            IPT = 0;
            //Age + 1/12 to make the loop running
            modelcd4 = cd4base-R0+R0*exp(-exp(lrc)*artmon*365.25/12);
            modelage += clength;
            ++artmon;
            tIPT += IPT;
            continue;
          }
          
          ++t_TBinfection;
          
        }else if(TBstate == "TP"){ //IPT treated, partially immune
          if(randomnums[ nPop*3-1 + id + r_TBinf*tmax + artmon ]  < RateToProb(TB_FoI *RR_agemix_vec[modelage/5] * IRR_LTBI, 1.0/12.0) ){
            TBstate = "L";
            t_TBinfection = 1; // first infection
            LTBI_time[id] = artmon;
          }
        }else if(TBstate == "I_undetected"){ //active TB, undetected
          //Cured by IPT
          if(IPT == 1 && INHR == 0){
            if(randomnums[ nPop*3-1 + id + r_cure_IPT*tmax + artmon ]  < RateToProb(rate_Cure_IPT, 1.0/12.0) ){
              TBstate = "R";
              t_TBinfection = 3; //reset TB infection time
            }else{
              //resistance?
              if(randomnums[ nPop*3-1 + id + r_INHR_IPT*tmax + artmon ]  < RateToProb(rate_INHR, 1.0/12.0) ){
                INHR = 1;
              }
            }
          }


          //Active TB but not detected yet -> can be detected to be put on treatment
          if(randomnums[ nPop*3-1 + id + r_detection*tmax + artmon ]  < p_Detect ){
            //DETECTED!!!!
            TBstate = "I_treating";
            IPT = 0;
            t_TBtreat = 0;
            
            //First TB detection
            if(TBhist == 0){
              TBcase[id] = 1;
              FU_TBcase[id] = artmon;
              tInf[id] = t_TBinfection;
              TBhist = 1;
              if(UndetectedTB_base[id] == 0){
                ++ n_newcases;
                if(t_TBinfection <= 24){++ n_recentInf; }
              }
            }
            
          }else{
            //On IPT --> 
            if(IPT == 1){
              //time to stop IPT?
              if((IPT_complete == 1 && tIPT == 6)||(IPT_complete == 0 && tIPT == IPTmo)){ 
                IPT = 0;
                if(IPT_cured == 1){TBstate = "TP";}
              }
            }
            //Not on IPT yet --> time to get on IPT?
            if(IPThist == 0 && artmon > 1){
              if(randomnums[ nPop*3-1 + id + r_IPT*tmax + artmon ] < prob_IPTstart){
                IPT = 1;
                tIPT_ART_vec[id] = artmon - 1;
                IPThist = 1;
              }
            }
          }

        }else if(TBstate == "I_treating"){ //active TB, on treatment
          if(t_TBtreat >= length_TbTx){ //finish treatment
            if(INHR == 0){ //no INH resistance
              if(randomnums[ nPop*3-1 + id + r_cure_TBTx*tmax + artmon ]  < p_Fail_pan){ //fail treatment
                TBstate = "I_undetected";
              }else{
                TBstate = "R";
              }
            }else{ //INH resistance
              if(randomnums[ nPop*3-1 + id + r_cure_TBTx*tmax + artmon ]  < p_Fail_INHR){ //fail treatment
                TBstate = "I_undetected";
              }else{
                TBstate = "R";
              }
            }
            t_TBinfection = 3; //reset TB infection time
          }else{
            t_TBtreat +=1;
          }
        }else if(TBstate == "R"){ //Recovered
          
          if(randomnums[ nPop*3-1 + id + r_progression*tmax + artmon ]  < RateToProb( 
              (exp(phi_vec_par1)+exp(phi_vec_par2)*(exp(-exp(phi_vec_par3)*t_TBinfection)))*
                (1 + exp(gamma_vec_par1) * (exp(-exp(gamma_vec_par2)*modelcd4))),1) ){
            //get TB again
            TBcase_true[id] ++;
            TBcases_by_year[artmon/12] ++;
            if(INHR == 1){TBcases_by_year_INHR[artmon/12] ++;}
            TBcases_by_cal[(enr_mo + artmon - 1)/12] ++;
            TBstate = "I_undetected";
            //Age + 1/12 to make the loop running
            modelcd4 = cd4base-R0+R0*exp(-exp(lrc)*artmon*365.25/12);
            modelage += clength;
            ++artmon;
            tIPT += IPT;
            continue;
          }
          
          if(randomnums[ nPop*3-1 + id + r_TBinf*tmax + artmon ]  < RateToProb(TB_FoI*RR_agemix_vec[modelage/5], 1.0/12.0) ){
            //TB reinfection
            t_TBinfection = 1; // reset time since infection
            //Age + 1/12 to make the loop running
            modelcd4 = cd4base-R0+R0*exp(-exp(lrc)*artmon*365.25/12);
            modelage += clength;
            ++artmon;
            tIPT += IPT;
            continue;
          }
          ++t_TBinfection;
        }
        
        //Update CD4
        modelcd4 = cd4base-R0+R0*exp(-exp(lrc)*artmon*365.25/12);
        
      }else{
        //4. Check TB status
        if(TBstate == "S"){
          //Susceptible state -> can get infected
          if(randomnums[ nPop*3-1 + id + r_TBinf*tmax + artmon ]  < RateToProb(TB_FoI*RR_agemix_vec[modelage/5], 1.0/12.0) ){
            TBstate = "L";
            t_TBinfection = 1; // first infection
            LTBI_time[id] = artmon;
          }
        }else if(TBstate == "L"){
          //Latent state -> can get reinfected or progress to active TB
          //Progression to active TB
          if(randomnums[ nPop*3-1 + id + r_progression*tmax + artmon ]  < RateToProb(
              (exp(phi_vec_par1)+exp(phi_vec_par2)*(exp(-exp(phi_vec_par3)*t_TBinfection)))*
                (1 + exp(gamma_vec_par1) * (exp(-exp(gamma_vec_par2)*modelcd4))),1) ){
            TBcase_true[id] ++;
            TBcases_by_year[artmon/12] ++;
            if(INHR == 1){TBcases_by_year_INHR[artmon/12] ++;}
            TBcases_by_cal[(enr_mo + artmon - 1)/12]++;
            
            //DETECTED!!!!
            if(randomnums[ nPop*3-1 + id + r_detection*tmax + artmon ]  < p_Detect/2 ){
              TBstate = "I_LTF_treating";
              t_TBtreat = 0;
            }else{
              TBstate = "I_LTF";;
            }
            
            //Update CD4
            modelcd4 += cd4_rate_postLTF;
            if(modelcd4 < 0){modelcd4 = 0;}
            //Age + 1/12 to make the loop running
            modelage += clength;
            ++artmon;
            tIPT += IPT;
            continue;
          }
          
          if(randomnums[ nPop*3-1 + id + r_TBinf*tmax + artmon ]  < RateToProb(TB_FoI*IRR_LTBI * RR_agemix_vec[modelage/5], 1.0/12.0) ){
            //TB reinfection
            t_TBinfection = 1; // reset time since infection
            //Update CD4
            modelcd4 += cd4_rate_postLTF;
            if(modelcd4 < 0){modelcd4 = 0;}
            //Age + 1/12 to make the loop running
            modelage += clength;
            ++artmon;
            tIPT += IPT;
            continue;
          }
          ++t_TBinfection;
        }else if(TBstate == "I_LTF"){
          //Active TB but not detected yet -> can be detected to be put on treatment
          if(randomnums[ nPop*3-1 + id + r_detection*tmax + artmon ]  < p_Detect/2 ){
            //DETECTED!!!!
            TBstate = "I_LTF_treating";
            t_TBtreat = 0;
          }
        }else if(TBstate == "I_LTF_treating"){
          if(t_TBtreat >= length_TbTx){ //finish treatment
            if(INHR == 0){ //no INH resistance
              if(randomnums[ nPop*3-1 + id + r_cure_TBTx*tmax + artmon ]  < p_Fail_pan){ //fail treatment
                TBstate = "I_LTF";
              }else{
                TBstate = "R";
              }
            }else{ //INH resistance
              if(randomnums[ nPop*3-1 + id + r_cure_TBTx*tmax + artmon ]  < p_Fail_INHR){ //fail treatment
                TBstate = "I_LTF";
              }else{
                TBstate = "R";
              }
            }
            t_TBinfection = 3; //reset TB infection time
          }else{
            t_TBtreat +=1;
          }
        }else if(TBstate == "R"){
          //Progression to active TB
          if(randomnums[ nPop*3-1 + id + r_progression*tmax + artmon ]  < RateToProb(
              (exp(phi_vec_par1)+exp(phi_vec_par2)*(exp(-exp(phi_vec_par3)*t_TBinfection)))*
                (1 + exp(gamma_vec_par1) * (exp(-exp(gamma_vec_par2)*modelcd4))),1) ){
            TBcase_true[id] ++;
            TBcases_by_year[artmon/12] ++;
            if(INHR == 1){TBcases_by_year_INHR[artmon/12] ++;}
            TBcases_by_cal[(enr_mo + artmon - 1)/12]++;
            TBstate = "I_LTF";
            //Update CD4
            modelcd4 += cd4_rate_postLTF;
            if(modelcd4 < 0){modelcd4 = 0;}
            //Age + 1/12 to make the loop running
            modelage += clength;
            ++artmon;
            tIPT += IPT;
            continue;
          }
          
          if(randomnums[ nPop*3-1 + id + r_TBinf*tmax + artmon ]  < RateToProb(TB_FoI * RR_agemix_vec[modelage/5], 1.0/12.0) ){
            //TB reinfection
            t_TBinfection = 1; // reset time since infection
            //Update CD4
            modelcd4 += cd4_rate_postLTF;
            if(modelcd4 < 0){modelcd4 = 0;}
            //Age + 1/12 to make the loop running
            modelage += clength;
            ++artmon;
            tIPT += IPT;
            continue;
          }
          ++t_TBinfection;
        }
        
        
        //Update CD4
        modelcd4 += cd4_rate_postLTF;
        if(modelcd4 < 0){modelcd4 = 0;}
      } 
      
      //Age + 1/12 to make the loop running
      modelage += clength;
      ++artmon;
      tIPT += IPT;
    }
    if(modelage >= 100){
      ++deaths_by_year(artmon/12);
      ++deaths_by_cal((enr_mo + artmon - 1)/12);
    }
    
    if (tobs == 0){tobs = modelage;}
    
    tobs_vec[id] = tobs - age + clength;
    lost_vec[id] = lost;
    dead_vec[id] = dead_obs;
    
    if (TBcase[id] == 0 && TBhist == 0){FU_TBcase[id] = (tobs - age)*12 + 1;}
    if (TBdeath[id] == 0){FU_TBdeath[id] = (tobs - age)*12 + 1;}
    
    DALY_vec[id] = DALY_vec[id] + LE_reference[modelage/5];
    double T1 = (artmon + enr_mo - 2)/12 ;
    double T2 = (artmon + enr_mo - 2)/12 + LE_reference[modelage/5];
    DALY_disc_vec[id] = DALY_disc_vec[id] + 1/pow((1 + disc_r), T1)*(1 - pow(1/(1+disc_r), T2- T1 + 1))/(1 - 1/(1 + disc_r));
  }

  
  //values to return

  NumericVector LE_vec_out(nPop); 
  NumericVector tobs_vec_out(nPop); 
  NumericVector lost_vec_out(nPop); 
  NumericVector dead_vec_out(nPop);
  NumericVector dead_by_2018_out(nPop);
  NumericVector TBcase_out(nPop);
  NumericVector FU_TBcase_out(nPop);
  NumericVector TBdeath_out(nPop);
  NumericVector FU_TBdeath_out(nPop);
  NumericVector HIVdeath_out(nPop);
  NumericVector tInf_out(nPop);
  NumericVector UndetectedTB_base_out(nPop);
  NumericVector LatentTB_base_out(nPop);
  NumericVector LTBI_time_out(nPop);
  double prop_recentInf = double(n_recentInf)/double(n_newcases);
  NumericVector DALY_vec_out(nPop);
  NumericVector cART_vec_out(nPop);
  NumericVector cTB_vec_out(nPop);
  NumericVector cIPT_vec_out(nPop);
  NumericVector tIPT_ART_vec_out(nPop);
  NumericVector LE_disc_vec_out(nPop); //LE disc
  NumericVector DALY_disc_vec_out(nPop); //DALY disc
  NumericVector cART_disc_vec_out(nPop); //cost of ART per patient disc
  NumericVector cTB_disc_vec_out(nPop); //cost of TB per patient disc
  NumericVector cIPT_disc_vec_out(nPop); //cost of IPT per patient disc
  NumericVector TBcase_true_out(nPop); //True TB case
  NumericVector TBdeath_true_out(nPop); //True TB death

  for(int i=0; i<nPop; i++) { 
    LE_vec_out(i)            = LE_vec[i];
    tobs_vec_out(i)          = tobs_vec[i];
    lost_vec_out(i)          = lost_vec[i];
    dead_vec_out(i)          = dead_vec[i];
    dead_by_2018_out(i)      = dead_by_2018[i];
    TBcase_out(i)            = TBcase[i];
    FU_TBcase_out(i)         = FU_TBcase[i];
    TBdeath_out(i)           = TBdeath[i];
    FU_TBdeath_out(i)        = FU_TBdeath[i];
    HIVdeath_out(i)          = HIVdeath[i];
    tInf_out(i)              = tInf[i];
    UndetectedTB_base_out(i) = UndetectedTB_base[i];
    LatentTB_base_out(i)     = LatentTB_base[i];
    LTBI_time_out(i)         = LTBI_time[i];
    DALY_vec_out(i)          = DALY_vec[i];
    cART_vec_out(i)          = cART_vec[i];
    cTB_vec_out(i)           = cTB_vec[i];
    cIPT_vec_out(i)          = cIPT_vec[i];
    tIPT_ART_vec_out(i)      = tIPT_ART_vec[i];
    LE_disc_vec_out(i)       = LE_disc_vec[i]; 
    DALY_disc_vec_out(i)     = DALY_disc_vec[i];
    cART_disc_vec_out(i)     = cART_disc_vec[i];
    cTB_disc_vec_out(i)      = cTB_disc_vec[i]; 
    cIPT_disc_vec_out(i)     = cIPT_disc_vec[i];
    TBcase_true_out(i)       = TBcase_true[i]; 
    TBdeath_true_out(i)      = TBdeath_true[i];
  }
  
  return(List::create(
      Named("LE")                   = LE_vec_out,
      Named("DALY")                 = DALY_vec_out,
      Named("cART")                 = cART_vec_out,
      Named("cTB")                  = cTB_vec_out,
      Named("cIPT")                 = cIPT_vec_out,
      Named("TBcase_true")          = TBcase_true_out,
      Named("TBdeath_true")         = TBdeath_true_out,
      Named("LE_disc")              = LE_disc_vec_out,
      Named("DALY_disc")            = DALY_disc_vec_out,
      Named("cART_disc")            = cART_disc_vec_out,
      Named("cTB_disc")             = cTB_disc_vec_out,
      Named("cIPT_disc")            = cIPT_disc_vec_out,
      Named("dead_onART")           = dead_vec_out,
      Named("TBcases_by_year")      = TBcases_by_year,
      Named("TBdeaths_by_year")     = TBdeaths_by_year,
      Named("deaths_by_year")       = deaths_by_year,
      Named("TBcases_by_cal")       = TBcases_by_cal,
      Named("TBdeaths_by_cal")      = TBdeaths_by_cal,
      Named("deaths_by_cal")        = deaths_by_cal,
      Named("TBcases_by_year_INHR") = TBcases_by_year_INHR
      )
  );
}

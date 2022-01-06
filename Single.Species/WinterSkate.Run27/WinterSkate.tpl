//For Winter skate, setting agesel phase = -4, and setting agesel == ficsel at beginning of parameter section

//SpZ, Run 32:
// Fitting to length proportions instead of age proportions; Did not test individual functions
// Adding small constant to zt calcs in FIC fx;  Mod. dat, FIC fxs (set FICtmp = TotFIC) and ofv fxs; TotC in thous mt, Incorporating FICfage, two trawl surveys, avg initial params not bounded vectors
//Removed //SSP Modified// statements (with exception of Pop dy function) as compared to Final SpZ 10-03-09 Version
//Parametering Yr1 as regular parameters, but leaving Age1 and Ft as means + deviations
//Function Details:
//  data_section = 01-22-10, SSP Version with revised sumFIC
//  parameter_section = 01-08-10a, SSP Version with 2 trawl surveys but avgs not bounded
//  procedure to top_of_main sections = 12-20-09, SSP Version
//  calc_initial_states = 01-08-10, SSP Version
//  calc_pop_dy = 12-16-09, SSP Version TotC in Thous mt
//  calc_survey_abundance = 10-19-10 Version
//  calc_biomass_penalty = 09-14-09 Version
//  calc_year1_constraint = 06-24-09 Version
//  calc_recruitment_penalty = 12-20-09 Version
//  objective_function through final section = 01-22-10a, SSP Version, with sumFIC and modified debug statements
//  report_section = 01-13-10, SSP Version, includes ofv_ideal, Age1CV, FICfage, 2 surveys, nFIC and FICmon

//Debug statements
  //  0's = Main body of program
      //  =  1:  Exits after data section
      //  =  2:  Exits after parameter section
      //  =  3:  Exits at end of procedure section
      //  =  4:  Prints checkpoints after each function in the procedure section, except those within the year loop
      //  =  5:  Prints checkpoints after each function within the year loop
      //  =  6:  Prints parameter estimates after each iteration
      //  =  7:  Prints pop dy variables after each iteration
      //  =  8:  Prints trophic matrices at end of year loop, then exits
      //  =  9:  Prints out predicted indices that go into the objective function after each iteration
  //10's = Initial states function
      //  = 12:  Outputs fishery and survey selectivity matrices at the end of the function and then exits
      //  = 13:  Outputs abundance and biomass arrays and then exits
      //  = 14:  Outputs Yr1, Age1 and iFt matrices to ensure 'means + devt'ns' parameterized correctly and then exits
  //40's = Population dynamics function
      //  = 40:  Outputs N, C_hat and Cprop_hat at end of function
      //  = 41:  Outputs N, C_hat and Cprop_hat at end of function and then exits
      //  = 42:  Outputs mortality components for each species in each year and exits at end of year loop after trophic =1
      //  = 43:  Outputs mortality components at end of function and then exits
  //50's = Survey abundance function
      //  = 50:  Prints intermediate arrays for species where survey data is one contiguous time series and then exits
      //  = 51:  Prints intermediate arrays for species where the survey data is split into multiple segments and then exits
      //  = 52:  Prints predicted q, FICs, FIC_hat and N for each species and then exits
      //  = 53:  Prints estimated q matrix at the end of each iteration
  //60's = Log likelihood function
      //  = 60: Prints checkpoints after each likelihood component
      //  = 61: Prints checkpoints for multinomial components within the year loop
      //  = 62: Prints predicted and log predicted indices for TotC and TotFIC
      //  = 63: Prints predicted and log predicted indices for Cprop
      //  = 64: Prints predicted and log predicted indices for Sprop
      //  = 65: FHprop, when added
      //  = 66: Prints summary of objective function components at end of each iteration
  //80's = Penalty functions
      // = 80: Biomass penalty function: Prints pre- and post- B for every species and year, regardless of whether a penalty is imposed
      // = 81: Biomass penalty function: Prints pre- and post- biomass and assorted arrays when biomass falls below threshold
      // = 82: Yr1 penalty function: Prints avgZ, thYr1, Yr1 and Ypen arrays and then exits at end of function
      // = 83: Recruitment penalty function: Prints Age1 parameter estimates, calculated CV and recruitment penalty


DATA_SECTION

  int debug;
  !!debug = 0;

  int i; int j; int t; int pd; int py; int b; int pdA; int a; int seg;
  int nf;
  !!nf = 1;
  number o;   //Tiny number for calculation of lognormal distributions
  !!o = 1.e-3;
  number p;   //Tiny number for calculation of multinomial residuals
  !!p = 1.e-30;

  init_int nsp;                        //Number of species
  int nsp2;                            //4darray index; init_int declarations (i.e. nsp) will not work
    !!nsp2 = nsp;
  init_int nFIC;                      //Number of trawl survey datasets
  int nFIC2;
    !!nFIC2 = nFIC;
  init_ivector fyr(1,nsp);        //First year in model; species-specific
  int minfyr;                          //Across all species, the earliest fyr
    !!minfyr = min(fyr);
  init_ivector lyr(1,nsp);         //Last year in model; species-specific
  int maxlyr;                         //Across all species, the latest lyr
    !!maxlyr = max(lyr);
  ivector nyr(1,nsp);              //Total number of years; species-specific
    !!nyr = 1 + lyr - fyr;

  int sage;                              //Total number of age classes summed across all species
  init_ivector nage(1,nsp);       //Number of age classes for each species
    !!sage = sum(nage);
  init_ivector nlen(1,nsp);        //Number of length bins for commercial and survey catch length frequencies
  ivector fage(1,nsp);              //Vector of the first age class for each species; used in trophic calcs
    !!for (i=1;i<=nsp;i++) {fage(i)=sum(nage(1,i))-(nage(i)-1);}
  ivector lage(1,nsp);              //Vector of the last age class for each species; used in trophic calcs
    !!for (i=1;i<=nsp;i++) {lage(i)=sum(nage(1,i));}

  init_int nMseg;                   //Number of segments of different natural mortality (M1) rates
  init_ivector nseg(1,nsp);      //Number of segments for the FIC data; for each segment, a separate q is estimated
  vector segtmp(1,nsp);
  !!segtmp = nseg;  for (i=1; i<=nsp; i++) {if (segtmp(i) > 1) {segtmp(i)--;}}
         //#breaks = FIC segments -1 unless there is only one segment, in which case #breaks= 1

  ivector nFt(1,nsp);                       //Number of fishing mortality parameters to be estimated
  ivector Fct(1,nsp);                       //Counter to assign initial F parameters, iFt, to the correct years
  init_ivector agePR(1,nsp);           //First age when partially recruited to the fishery; used to avoid hitting parameter bounds
  init_ivector ageFR(1,nsp);           //Age of full recruitment to the fishery
  init_ivector ficFR(1,nsp);            //Age of full recruitment to the fishery-independent survey, FIC
                                                     //ficFR == 0 for species that never become fully recruited to the survey
  init_vector FICs_lage(1,nsp);      //Survey selectivity coefficient for the last age class; used to anchor the curve
                                                     //FICs_lage == 1 for species that become fully recruited to the survey
  ivector nsel(1,nsp);                     //Number of fishery selectivity parameters estimated per species
  ivector FICnsel(1,nsp);               //Number of survey selectivity parameters estimated per species

  vector CPideal(1,nsp);                //Best possible value for the multinomial residuals of CAA proportions

  init_ivector M1yr(1,nMseg);      //The first year in each M1 segment, including the first segment

  //Data weightings
  init_vector TCwt(1,nsp);            //Total annual commercial catch in weight
  init_vector CPwt(1,nsp);            //Commercial catch proportions at age
  init_vector Bwt(1,nsp);             //Weight for Biomass Penalty Term, Bpen
  init_vector Ywt(1,nsp);             //Weight for Yr1 Penalty Term, Ypen)
  init_vector Rwt(1,nsp);             //Weight for Recruitment Penalty Term, Rpen
  init_vector Rwt2(1,nsp);           //Weight for Recruitment Penalty Term, Rpen2

  init_vector Bthres(1,nsp); //Biomass threshold used in the penalty function to avoid B == 0, which would cause M2 calc to crash
  ivector nBpen(1,nsp);      //Number of iterations where Bpen (within the calc_biomass_penalty function) was > 0
  ivector lBpen(1,nsp);       //Last iteration where the Bpen was > 0
  init_vector Rthres(1,nsp);//Threshold for the coefficient of variation of recruitment
  init_matrix TSwt(1,nsp,1,nFIC); //Total annual survey catch in number/tow
  init_matrix SPwt(1,nsp,1,nFIC); //Survey catch proportions-at-age

  init_matrix iM2(1,nsp,1,nage);              //M2 estimates (from previous model run) to use in non-FH years
  matrix SPideal(1,nsp,1,nFIC);              //Best possible value for the multiomial residuals of Survey, FIC, proportions
  
  init_imatrix FICfage(1,nsp,1,nFIC);      //First age captured for each trawl survey
  init_imatrix FIClage(1,nsp,1,nFIC);      //Last age captured for each trawl survey

  //Time series data
  init_matrix FICmon(1,nsp,1,nFIC);               //Month correpsonding to each trawl survey
  init_matrix FICyr(1,nsp,1,segtmp);               //Beginning year of FIC segments for species where nseg > 1; if nseg == 1, FICyr = 0
  init_matrix TotC(1,nsp,fyr,lyr);                     //Total commercial catch; metric tons
  matrix sumC(1,nsp,fyr,lyr);                           //Sum of obs CAA; Used to determine if age samples were taken in a particular yr
  init_3darray Wt(1,nsp,fyr,lyr,1,nage);           //Average individual weight-at-age; kg
  init_3darray ALKey(1,nsp,1,nlen,1,nage);    //Age-Length Key to convert age- to length- propz
  init_3darray Clen(1,nsp,fyr,lyr,1,nlen);         //Commercial Catch-at-Length; 10^6 number

  3darray CPlen(1,nsp,fyr,lyr,1,nlen);             //Proportion-at-length of the commercical catch, CAA
  init_3darray TotFIC(1,nsp,1,nFIC,fyr,lyr);    //Total annual survey catch; Necessary for years where age samples were not taken
  3darray sumFIC(1,nsp,1,nFIC,fyr,lyr);         //Sum of observed FIC; Used to determine if age samples were taken in a particular yr  

  init_3darray M1seg(1,nsp,1,nMseg,1,nage);  //Natural mortality rates or each segment, Mseg
  3darray M1(1,nsp,fyr,lyr,1,nage);                 //Natural mortality rate expanded to full 3darray

  init_4darray Slen(1,nsp2,1,nFIC2,fyr,lyr,1,nlen);   //Fishery-Independent trawl survey Catch (FIC) at length; number/tow
  4darray SPlen(1,nsp2,1,nFIC2,fyr,lyr,1,nlen);       //Proportion-at-length of the survey catch, FIC
  
  init_int eof;

	LOCAL_CALCS
    nFt.initialize();
    CPlen.initialize();  SPlen.initialize();  SPideal.initialize();  CPideal.initialize();
    nBpen.initialize();  lBpen.initialize();

    //Number of fishery selectivity parameters for each species:
    nsel =  ageFR-agePR;

    //Number of survey selectivity parameters for each species
    for (i = 1; i<=nsp; i++)  
      {
      if (ficFR(i) == 0)  {FICnsel(i) = nage(i) - 1;}
      else {FICnsel(i) = ficFR(i) - 1;}
      }

    for (i = 1; i<= nsp; i++)
      {
      sumC(i) = rowsum(Clen(i));        //Observed total catch in each year
      int M1ct = 1;                         //Count to correctly fill M1 array

      for (t = fyr(i); t<=lyr(i); t++)
        {
        //Fill M1 3darray
		    if (M1ct < nMseg && t == M1yr(M1ct + 1)) {M1ct ++;}
   		  M1(i,t) = M1seg(i,M1ct);

        //Determine number of F parameters to initialize; only initialize if fish were caught 
        if (TotC(i,t) > 0.) {nFt(i)++;}

        //Calculate the CLen proportions if the total catch in that year is not zero
        if (sumC(i,t) != 0) {CPlen(i)(t) = Clen(i)(t)/sumC(i)(t);}

        //Multinomial residuals for the CAA data for a perfect fit
        dvector CPobs = CPwt(i)*elem_prod( CPlen(i,t)+p , log(CPlen(i,t)+p) );
        CPideal(i) -= sum(CPobs);
        }  //end of year loop

      //Survey datasets
      for (j = 1; j<=nFIC; j++)
        {
        sumFIC(i,j) = rowsum( Slen(i,j) );  //Observed total survey catch in each year
        for (t = fyr(i); t<=lyr(i); t++)
        {
        if (sumFIC(i,j,t) != 0)
          {
          //Calculating the Survey proportions
          SPlen(i)(j)(t) = Slen(i)(j)(t)/sumFIC(i)(j)(t);
          //Multinomial residuals for the CAA data for a perfect fit
          dvector SPobs =  SPwt(i)(j)*elem_prod( SPlen(i,j,t)+p , log(SPlen(i,j,t)+p) );
          SPideal(i)(j) -= sum(SPobs);
          }  //end of sumFIC If Statement
        }  //end of year loop
      }  //end of nFIC loop

    }  //end of species loop

  if (debug == 1)
    {
    cout<<"nage\n"<<nage<<endl;
    cout<<"ageFR\n"<<ageFR<<endl;
    cout<<"nsel\n"<<nsel<<endl;
    cout<<"nFIC\n"<<nFIC<<endl;
    cout<<"FICfage\n"<<FICfage<<endl;
    cout<<"FIClage\n"<<FIClage<<endl;
    cout<<"FICmonth\n"<<FICmon<<endl;
    cout<<"ficFR\n"<<ficFR<<endl;
    cout<<"FICs_lage\n"<<FICs_lage<<endl;
    cout<<"FICnsel\n"<<FICnsel<<endl;
    cout<<"M1seg\n"<<M1seg<<endl;
    cout<<"nseg\n"<<nseg<<endl;
    cout<<"segtmp\n"<<segtmp<<endl;
    cout<<"FICyr\n"<<FICyr<<endl;
    cout<<"Clen\n"<<Clen<<endl;
    cout<<"Slen\n"<<Slen<<endl;
    cout<<"TotC\n"<<TotC<<endl;
    cout<<"sumFIC\n"<<sumFIC<<endl;
    cout<<"TotFIC\n"<<TotFIC<<endl;
    cout<<"CPlen\n"<<CPlen<<endl;
    cout<<"SPlen\n"<<SPlen<<endl;
    cout<<"nFt\n"<<nFt<<endl;
    cout<<"TSwt\n"<<TSwt<<endl;
    cout<<"SPwt\n"<<SPwt<<endl;
    cout<<"SPideal\n"<<SPideal<<endl;
    cout<<"eof\n"<<eof<<endl;
    }

  if(eof != 54321) {cout<<"Stop, data not inputted correctly"<<endl<<"eof: "<<eof<<endl; exit(1);}

  if (debug == 1) {cout<<"\nManually exiting at end of data section..."<<endl;  exit(-1);}

	END_CALCS


PARAMETER_SECTION
  objective_function_value ofv;
  number ofv_ideal;         //Total objective function value, subtracting out the ideal multinomial values

  vector Bpen(1,nsp);      //Penalty function for biomass and calculation of M2
  vector tBpen(1,nsp);     //Total Biomass penalty for each species
  vector Ypen(1,nsp);      //Penalty function for Year 1 abundances
  vector tYpen(1,nsp);     //Total Yr1 penalty for each species
  vector Rpen(1,nsp);      //Penalty function for the CV of recruitment
  vector Rpen2(1,nsp);      //Penalty function for recruitment deviations in the last three years
  vector tRpen(1,nsp);     //Total Recruitment penalty for each species
  vector tRpen2(1,nsp);     //Total Recruitment penalty for each species
  vector Devs(1,nsp);      //Square of summed deviations for each deviation initial parameter
  vector TCres(1,nsp);     //Objective function component: Total commercial catch
  vector CPmulti(1,nsp);  //Objective function component: Catch-at-age proportions
  vector ofvsp(1,nsp);          //Total objective function value for each species
  vector ofvsp_ideal(1,nsp); //Total objective function value for each species, subtracting out the ideal multinomial values

  init_vector aAge1(1,nsp,1);                                     //Average annual recruits, log space
  init_vector aFt(1,nsp,1);                                          //Average annual fishing mortality rates; averaged over years where TotC > 0, log space
  init_vector aYr1(1,nsp,-2);                                       //Scalar for Yr-1 abundance
  init_bounded_matrix dAge1(1,nsp,fyr+1,lyr,-10,5,3);   //Annual deviations in recruits, log space
  init_bounded_matrix dFt(1,nsp,1,nFt,-5,5,3);             //Annual deviations in species-specific Fs in yrs where TotC > 0, log space

  matrix Ft(1,nsp,fyr,lyr);                                          //Annual species-specific fishing mortality rates
  init_bounded_matrix FICsel(1,nsp,1,FICnsel,0,1,4)   //Age-specific survey selectivity parameters
  init_bounded_matrix agesel(1,nsp,1,nsel,0,1,-4);        //Age-specific fishery selectivity/vulnerability parameters

  init_matrix Yr1(1,nsp,1,nage,2);                              //Year 1 age-specific N's, 'means + deviations', log space
  matrix Age1(1,nsp,fyr+1,lyr);                                 //Annual recruits, 'means + deviations', log space
  matrix iFt(1,nsp,1,nFt);                                          //Fishing mortality rates, 'means + deviations',  in yrs where TotC > 0, log space

  matrix s(1,nsp,1,nage);                            //Fishery selectivity/vulnerability at age
  //Set s = 1 so that oldest ages will be fully recruited when the initial parameters are added to s
  !!s = 1;
  matrix FICs(1,nsp,1,nage);                      //Survey selectivity at age
  //Set selectivity coefficients to the inputted values for the final ages, FICs_lage, to anchor the curve
  //For species that become fully recruited to the survey, FICs_lage = 1; for others it represents the selectivity of the last age class 
  !!for (i=1; i<=nsp; i++)  {FICs(i)(FICnsel(i)+1,nage(i)) = FICs_lage(i);}

  matrix TotC_hat(1,nsp,fyr,lyr);     //Total commerical catch *in weight*; summed over ages
  matrix sumC_hat(1,nsp,fyr,lyr);    //Total commercial catch *in numbers*; summed over ages; for calc of CAA proportions
  matrix TCresid(1,nsp,fyr,lyr);       //Residuals of total commercial catch *in weight*; summed over ages
  matrix TSresA(1,nsp,1,nFIC);      //Objective function component: Total survey catch, method A
  matrix TSresB(1,nsp,1,nFIC);      //Objective function component: Total survey catch, method B
  matrix SPmulti(1,nsp,1,nFIC);      //Objective function component: Survey catch-at-age proportions

  3darray q(1,nsp,1,nFIC,1,nseg);                 //Survey catchability; age-invariant
  3darray toteps(1,nsp,1,nFIC,fyr,lyr);          //Residuals of total survey catch, TotFIC, Calculated by Method A
  3darray TotFIC_hat(1,nsp,1,nFIC,fyr,lyr); //Total survey (FIC) catch in number/tow; summed over ages
  3darray TSresid(1,nsp,1,nFIC,fyr,lyr);       //Residuals of total survey catch, TotFIC, Calculated by Method B

  3darray N(1,nsp,fyr,lyr,1,nage);                //Abundance-at-age
  3darray B(1,nsp,fyr,lyr,1,nage);                //Biomass-at-age = N * Wt
  3darray F(1,nsp,fyr,lyr,1,nage);                //Predicted realized F-at-age = Ft * s
  3darray Z(1,nsp,fyr,lyr,1,nage);                //Total mortality-at-age
  3darray M2(1,nsp,fyr,lyr,1,nage);             //Predation mortality-at-age
  3darray C_hat(1,nsp,fyr,lyr,1,nage);         //Predicted commercial catch-at-age
  3darray Cprop_hat(1,nsp,fyr,lyr,1,nage);  //Commercial catch proportions-at-age
  3darray CPlen_hat(1,nsp,fyr,lyr,1,nlen);    //Predicted commercial catch proportions-at-length

  4darray FIC_hat(1,nsp2,1,nFIC2,fyr,lyr,1,nage);     //Predicted survey catch, FIC, at age
  4darray Sprop_hat(1,nsp2,1,nFIC2,fyr,lyr,1,nage);  //Survey catch proportions-at-age
  4darray SPlen_hat(1,nsp2,1,nFIC2,fyr,lyr,1,nlen);    //Survey catch proportions-at-length

	LOCAL_CALCS
    if (debug == 2)
      {
      cout<<"Yr1\n"<<Yr1<<endl;
      cout<<"aAge1\n"<<aAge1<<endl<<"aFt\n"<<aFt<<endl;
      cout<<"dAge1\n"<<dAge1<<endl<<"dFt\n"<<dFt<<endl;
      cout<<"agesel\n"<<agesel<<endl<<"FICsel\n"<<FICsel<<endl;
      cout<<"\nManually exiting at the end of the parameter section...\n"<<endl;
      exit(-1);
      }
	END_CALCS


PROCEDURE_SECTION

  //Winter skate addition
  agesel = FICsel;
  
  calc_initial_states();  if (debug == 4) {cout<<"completed Initial States"<<endl;}

  for (t = minfyr; t<=maxlyr; t++) 
    {
    if (debug == 5) {cout<<"Year: "<<t<<endl;}
    calc_biomass_penalty();  if (debug == 5) {cout<<"completed Biomass Penalty"<<endl;}
    pop_dynamics();  if (debug == 5) {cout<<"completed Pop Dynamics"<<endl;}
    }  //end of year loop
  if (debug == 4) {cout<<"completed Year loop"<<endl;}

  calc_survey_abundance();  if (debug == 4) {cout<<"completed Survey Abundance"<<endl;}
  calc_year1_constraint();  if (debug == 4) {cout<<"completed Yr1 Constraint"<<endl;}
  calc_recruitment_penalty();  if (debug == 4) {cout<<"completed Recruitment Penalty"<<endl;}
  calc_length_props();  if(debug == 4) {cout<<"completed Length Props"<<endl;}
  calc_negative_loglikelihood();  if (debug == 4) {cout<<"completed Log Likelihood"<<endl;}

  if (debug == 6) {cout<<"Parameter estimates:\n"<<"Yr1\n"<<Yr1<<"\nAge1\n"<<Age1<<"\nagesel\n"<<agesel<<"\nFt\n"<<Ft
                                 <<"\nFICsel\n"<<FICsel<<endl<<endl;}
  if (debug == 7) {cout<<"N\n"<<N<<"\nC_hat\n"<<C_hat<<"\nZ\n"<<Z<<"\nM2\n"<<M2<<"\nF\n"<<F<<endl<<endl;}
  if (debug == 9) {cout<<"TotC_hat\n"<<TotC_hat<<"\nTotC\n"<<TotC<<"\nTCresid\n"<<TCresid
                                 <<"\nTotFIC_hat\n"<<TotFIC_hat<<"\nTotFIC\n"<<TotFIC<<"\ntoteps\n"<<toteps
                                 <<"\nCprop_hat\n"<<Cprop_hat<<"\nCPlen_hat\n"<<CPlen_hat<<"\nCPlen\n"<<CPlen
                                 <<"\nSprop_hat\n"<<Sprop_hat<<"\nSPlen_hat\n"<<SPlen_hat<<"\nSPlen\n"<<SPlen<<endl;  }

  nf++;  //Tallys the # of function evaluations
  if (debug == 66) {cout<<"\nmoving onto next iteration "<<"nf: "<<nf<<endl<<endl;}

  if (debug == 3) 
    {
    cout<<endl<<"ofv: "<<ofv<<endl<<"ofv_ideal: "<<ofv_ideal<<endl;
    cout<<endl<<"manually exiting at end of procedure section...\n"<<endl;
    exit(-1);
    }


GLOBALS_SECTION
  //Including C++ libraries
  #include <stats.cxx>


RUNTIME_SECTION
  convergence_criteria 1.e-3 ,  1.e-4
  maximum_function_evaluations 80000


TOP_OF_MAIN_SECTION
  arrmblsize = 8000000;  //Increase amount of available dvar memory


FUNCTION calc_initial_states
  /*
  Function 1) Initializes the differentiable arrays, 2) Calculates fishery recruitment at age (s),
    3)Initializes age-specific abundance and biomass in the first year, and 4) Expands the food-selection 
    initial parameter vectors into matrices [predator x prey species]
  */
  Fct.initialize();  //Counter to determine number of estimable fishing mortality rates
  N.initialize();  B.initialize();  
  Bpen.initialize();  Ypen.initialize();  Rpen.initialize();
  M2.initialize();  Z.initialize();  C_hat.initialize();  FIC_hat.initialize();  TotFIC_hat.initialize();  TotC_hat.initialize();
  Cprop_hat.initialize();  Sprop_hat.initialize();  CPlen_hat.initialize();  SPlen_hat.initialize();
  q.initialize();  toteps.initialize();  TSresid.initialize();
  Age1.initialize();  iFt.initialize();

  for (i=1; i<=nsp; i++)
    {
    //Expand species-specific averages and annual/age-specific deviations into full arrays
    Age1(i) = aAge1(i) + dAge1(i);
    iFt(i) = aFt(i) + dFt(i);

    //Add the fishery selectivity initial parameters to the selectivity matrix; method is a function of agePR
    if (agePR(i) == 1) {s(i)(1,nsel(i)) = agesel(i);}
    else {
      s(i)(1,agePR(i)-1) = 0;
      agesel(i).shift(agePR(i));  //Shift the index to ensure compatible array bounds
      s(i)(agePR(i),agePR(i)+nsel(i)-1) = agesel(i);
           }
    //Add the survey selectivity initial parameters to the survey selectivity matrix
    FICs(i)(1,FICnsel(i)) = FICsel(i);

    //Initializes the populations in the first year
    N(i)(fyr(i)) = mfexp(aYr1(i))*mfexp(Yr1(i));            //Abundance in first year
    B(i)(fyr(i)) = elem_prod(N(i)(fyr(i)),Wt(i)(fyr(i)));     //Biomass-at-age = N * Wt
    if (debug == 13) {cout<<"sp: "<<i<<endl<<"N(i,fyr)\n"<<N(i)(fyr(i))<<endl<<"Wt(i,fyr)\n"<<Wt(i)(fyr(i))<<endl
                                      <<"B(i,fyr)\n"<<B(i)(fyr(i))<<endl<<endl;}

    }  //end of species loop

  if (debug == 12) {cout<<"nage\n"<<nage<<endl<<"agePR\n"<<agePR<<endl<<"ageFR\n"<<ageFR<<endl
                                   <<"nsel\n"<<nsel<<endl<<"agesel\n"<<agesel<<endl<<"s\n"<<s<<endl
                                   <<"ficFR\n"<<ficFR<<endl<<"FICs_lage\n"<<FICs_lage<<endl<<"FICnsel\n"<<FICnsel<<endl
                                   <<"FICsel\n"<<FICsel<<endl<<"FICs\n"<<FICs<<endl; exit(-1);}
  if (debug == 13) {cout<<"N\n"<<N<<endl<<"B\n"<<B<<endl;  exit(-1);}
  if (debug == 14) {cout<<"Yr1\n"<<Yr1<<endl<<"aAge1\n"<<aAge1<<endl<<"dAge1\n"<<dAge1<<endl<<"Age1\n"<<Age1<<endl
                                   <<"aFt\n"<<aFt<<endl<<"dFt\n"<<dFt<<endl<<"iFt\n"<<iFt<<endl ;  exit(-1);}


FUNCTION pop_dynamics
   /*
  Function calculates 1) age-specific fishing mortality, F, 2) total mortality, Z, 3) state dynamics, N & B, and 4) commercial ,
  catch, including A) catch-at-age, C_hat, B) total catch in weight, TotC_hat, C) total catch in number, sumC_hat, and
  D) catch proportions at age, Cprop_hat.
  Function assumes N and C are in millions, TotC is in thousands of metric tons  
  */
  for (i=1; i<=nsp; i++)
    {
    //Only calculate state dynamics if the year is between the first and last yrs for that particular species
    if (t>=fyr(i) && t<=lyr(i))
      {
      //Expand initial F parameters, iFt, into a species x year matrix, Ft
      if (TotC(i,t) != 0) {
        Fct(i)++;
        Ft(i,t) = mfexp(iFt(i,Fct(i))); }
      else {Ft(i,t) = 0.;}

      //Age-specific fishing mortality
      F(i,t) = Ft(i,t)*s(i);

      //Total mortality
        //When either 1) trophic interactions are not turned on, or 2) data for all species are not available for year t, M2 is set
        //at a fixed/inputted species- and age- specific value.  Otherwise, M2 is set to the value calculated within the model
      Z(i)(t)=M1(i)(t)+F(i)(t)+iM2(i);
      /*  //SSP MODIFIED//
      if (trophic) {
        if (t>=FHfyr && t<=FHlyr) {Z(i)(t)=M1(i)(t)+F(i)(t)+M2(i)(t);}
        }  //end of Trophic If Statement
      */  //SSP MODIFIED//
      
      if (debug == 42) {cout<<"nf: "<<nf<<" sp: "<<i<<" yr: "<<t<<endl
                                       <<"M1\n"<<M1(i)(t)<<endl<<"iM2\n"<<iM2(i)<<endl<<"M2\n"<<M2(i)(t)<<endl
                                       <<"F\n"<<F(i)(t)<<endl<<"Z\n"<<Z(i)(t)<<endl<<endl;}
      //State dynamics
      if (t != lyr(i))
        {
        N(i)(t+1,1)=mfexp(Age1(i,t+1));
          //Annual recruits
        N(i)(t+1)(2,nage(i))=++ elem_prod(N(i)(t)(1,nage(i)-1),mfexp(-Z(i)(t)(1,nage(i)-1)));
          //Subsequent age classes
        N(i)(t+1)(nage(i)) += N(i)(t,nage(i))*mfexp(-Z(i)(t,nage(i)));
          //Plus group
        B(i)(t+1) = elem_prod(N(i)(t+1),Wt(i)(t+1)); 
        }  //end of lyr If statement

      //Commercial catch, Baranov catch equation
      dvar_vector Ntmp = elem_prod( N(i)(t),1. - exp(-Z(i)(t)) );
      C_hat(i)(t) = elem_prod( elem_div(F(i)(t),Z(i)(t)),Ntmp );                  //millions of fish
      //TotC_hat(i)(t) = sum(1000*elem_prod( C_hat(i)(t),Wt(i)(t) ));         //metric tons
      //TotC_hat(i)(t) = 1.e-6*sum(1000*elem_prod( C_hat(i)(t),Wt(i)(t) ));  //millions of metric tons
      TotC_hat(i)(t) = 1.e-3*sum(1000*elem_prod( C_hat(i)(t),Wt(i)(t) ));  //thousands of metric tons
      sumC_hat(i)(t) = sum(C_hat(i)(t));
      if (sumC_hat(i)(t) != 0) {Cprop_hat(i)(t) = C_hat(i)(t)/sumC_hat(i)(t);}

      } //end of fyr&&lyr If statement
    }  //end of species loop
  if (debug == 40) {cout<<"Z\n"<<Z<<endl<<"N\n"<<N<<endl<<"C_hat\n"<<C_hat<<endl<<"Cprop_hat\n"<<Cprop_hat<<endl;}
  if (debug == 41) {cout<<"Z\n"<<Z<<"N\n"<<N<<endl<<"C_hat\n"<<C_hat<<endl<<"Cprop_hat\n"<<Cprop_hat<<endl;  exit(-1);}
  if (debug == 43) {cout<<"M1yr\n"<<M1yr<<endl<<endl
                                   <<"Z\n"<<Z<<endl<<"M1\n"<<M1<<endl<<"iM2\n"<<iM2<<endl<<"M2\n"<<M2<<endl
                                   <<"F\n"<<F<<endl<<"Ft\n"<<Ft<<endl<<"s\n"<<s<<endl;  exit(-1);}


FUNCTION calc_survey_abundance
    /*
    Function calculates the residuals between the fishery-independent survey catch-at-age (FIC) data
    and the vulnerbale proportions-at-age to the survey
    */
    for (i=1; i <= nsp; i++)
      {
    for (j=1; j<=nFIC; j++)
      {
      int bage = FICfage(i,j);
      int eage = FIClage(i,j);
            
      //Species for which FIC data are not split into multiple segments
      if (nseg(i) == 1)
        {
        //The following code assumes one calculated q for every species (age-invariant) but calculates this q using total FIC and total 
          //N instead of age-specific indices.  Assumes time-invariant, but age-specific, survey selectivity (Method 2)

        //Incorporating time-invariant survey selectivity and the timing of the trawl survey;
            //selN represents the portion of N that has survived to the time of the survey and the proportion of that portion selected by the survey
 				    //selN = FICs * N * exp(-Z*FICmon/12):  Incorporates amount of Z the sp has incurred between Jan1 and FICmon
        dvar_matrix selN(fyr(i),lyr(i),bage,eage);
        for (t = fyr(i); t<=lyr(i); t++) {
          selN(t) = elem_prod(  elem_prod(FICs(i),N(i,t)) , mfexp(-Z(i,t)*(FICmon(i,j)/12))  ).sub(bage,eage);  }

       //Calculate q from deviations between FIC and N
       dvar_vector FICtmp = TotFIC(i)(j);
       dvar_vector Ntmp = rowsum(selN);
       if (debug == 50) {cout<<"\nLoop for species without any breaks in the survey data\nSpecies\n"<<i<<endl<<"Survey\n"<<j<<endl
                                        <<"q\n"<<q(i)(j)<<endl<<"FICtmp\n"<<FICtmp<<endl
                                        <<"FICs\n"<<FICs(i)<<endl<<"N\n"<<N(i)<<endl<<"selN\n"<<selN<<endl<<"Ntmp\n"<<Ntmp<<endl;}
                                        
       dvar_vector zt = log(FICtmp+p) - log(Ntmp+p);
       q(i,j,1) = mean(zt);
       if (debug == 50) {cout<<"zt\n"<<zt<<endl<<"q(i,j) = mean(zt)\n"<<q(i,j,1)<<endl;}
       toteps(i,j) = zt - q(i,j,1);
       if (debug == 50) {cout<<"toteps(i,j)\n"<<toteps(i,j)<<endl;}

       //Calculate predicted FIC
 			   //FIC_hat = q * FICs * N * exp(-Z*FICmon/12):  Incorporates amount of Z the sp has incurred between Jan1 and FICmon
       for (t = fyr(i); t<=lyr(i); t++) {FIC_hat(i,j,t) = mfexp(q(i,j,1))*elem_prod(  elem_prod(FICs(i),N(i,t)) , mfexp(-Z(i,t)*(FICmon(i,j)/12))); }

       if (debug == 52) {cout<<"sp: "<<i<<endl<<"survey: "<<j<<endl<<"mfexp(q)\n"<<mfexp(q(i,j,1))<<endl<<"FICs\n"<<FICs(i)<<endl
                                        <<"FIC_hat\n"<<FIC_hat(i)(j)<<endl<<"N\n"<<N(i)<<endl;}
      if (debug == 50) {cout<<"End of nseg=1 If statement\n"<<endl;}

        }  //end of nseg=1 If statement


      if (nseg(i) >1)
        {
        //The following code assumes multiple calculated q's for every species (age-invariant, but splitting the FIC time series) but
          //calculates this q using total FIC and total N instead of age-specific indices (Method 2A)

        //Incorporating time-invariant survey selectivity and the timing of the trawl survey;
            //selN represents the portion of N that has survived to the time of the survey and the proportion of that portion selected by the survey
 				    //selN = FICs * N * exp(-Z*FICmon/12):  Incorporates amount of Z the sp has incurred between Jan1 and FICmon

        dvar_matrix selN(fyr(i),lyr(i),bage,eage);
        for (t = fyr(i); t<=lyr(i); t++) {
          selN(t) = elem_prod(  elem_prod(FICs(i),N(i,t)) , mfexp(-Z(i,t)*(FICmon(i,j)/12))  ).sub(bage,eage);  }    //Incorp. survey timing

        //First and last years of each FIC segment
        if (debug == 51) {cout<<"\nLoop for species with breaks in the survey data\nSpecies'n"<<i<<endl<<"Survey\n"<<j<<endl
                                         <<"q\n"<<q(i)(j)<<endl<<"FICs\n"<<FICs(i)<<endl<<"N\n"<<N(i)<<endl
                                         <<"selN\n"<<selN<<endl<<"FICyr\n"<<FICyr(i)<<endl;}
        dvector FICfyr(1,nseg(i));  FICfyr.initialize();
        dvector FIClyr(1,nseg(i));  FIClyr.initialize();
        FICfyr(1) = fyr(i);
        FIClyr(nseg(i)) = lyr(i);

        for (seg=1; seg<nseg(i); seg++)
          {
          FICfyr(seg+1) = FICyr(i,seg);
          FIClyr(seg) = FICyr(i,seg)-1;
          }
        if (debug == 51) {cout<<"FICfyr\n"<<FICfyr<<endl<<"FIClyr\n"<<FIClyr<<endl;}

       //For each segment, calculate q from deviations between FIC and N
        for (seg=1; seg<=nseg(i); seg++)
          {
          int byr = FICfyr(seg);
          int eyr = FIClyr(seg);
          dvar_vector FICtmp = TotFIC(i)(j).sub(byr,eyr);
          dvar_vector Ntmp = rowsum(selN.sub(byr,eyr));
          if (debug == 51) {cout<<"\nseg: "<<seg<<endl<<"FICfyr: "<<byr<<"  FIClyr: "<<eyr<<endl<<"TotFIC\n"<<TotFIC(i,j)
                                          <<endl<<"FICtmp\n"<<FICtmp<<endl<<"N\n"<<N(i)<<endl<<"Ntmp\n"<<Ntmp<<endl;}
          dvar_vector zt = log(FICtmp+p) - log(Ntmp+p);
          q(i,j,seg) = mean(zt);
          if (debug == 51) {cout<<"zt\n"<<zt<<endl<<"q(i,j,seg) = mean(zt)\n"<<q(i,j,seg)<<endl;}
          toteps(i)(j).sub(byr,eyr) = zt - q(i,j,seg);
          if (debug == 51) {cout<<"toteps\n"<<toteps(i)(j).sub(byr,eyr)<<endl;}

         //Calculate predicted FIC
 				   //FIC_hat = q * FICs * N * exp(-Z*FICmon/12):  Incorporates amount of Z the sp has incurred between Jan1 and FICmon
         for (t = byr; t<=eyr; t++) {FIC_hat(i,j,t) = mfexp(q(i,j,seg))*elem_prod(  elem_prod(FICs(i),N(i,t)) , mfexp(-Z(i,t)*(FICmon(i,j)/12))); }
 
          } //end of seg loop
 
         if (debug == 52) {cout<<"sp: "<<i<<endl<<"survey: "<<j<<endl<<"mfexp(q)\n"<<mfexp(q(i,j,1))<<endl<<"FICs\n"<<FICs(i)<<endl
                                        <<"FIC_hat\n"<<FIC_hat(i)(j)<<endl<<"N\n"<<N(i)<<endl;}

        }  //end of nseg  > 1 If statement
      if (debug == 51) {cout<<"End of nseg>1 If statement\n"<<endl;}

      //Calculating Survey Catch, FIC, proportions-at-ge
      TotFIC_hat(i)(j) = colsum(trans(FIC_hat(i)(j)).sub(bage,eage));
      for (t = fyr(i); t<=lyr(i); t++)
        {
        if (TotFIC_hat(i,j,t) != 0) { Sprop_hat(i)(j)(t).sub(bage,eage) = FIC_hat(i)(j)(t).sub(bage,eage)/TotFIC_hat(i)(j)(t); } 
        }
      if (debug == 52)  {cout<<"TotFIC_hat\n"<<TotFIC_hat(i)(j)<<endl<<"Sprop_hat\n"<<Sprop_hat(i)(j)<<endl<<endl;}

      }  //end of nFIC loop

      if (debug == 50 | debug == 51) {cout<<"Exit at end of nFIC loop\n"<<endl; exit(-1);}

      }  //end of species loop

      if (debug == 52) {cout<<"Exit at end of species loop\n"<<endl;  exit(-1);}
      if (debug == 53) {cout<<"mfexp.q\n"<<mfexp(q)<<endl;}


FUNCTION calc_biomass_penalty
  /*
  Function calculates a penalty, Bpen, that is added to the objective function in the calc_likelihood function if the biomass for
  any species falls below a threshold value, Bthres
  */
  dvariable tmpBpen;
  for (i = 1; i<=nsp; i++)
    {
    if (t >= fyr(i) && t<=lyr(i))
      {
      if (debug == 80) {cout<<"t: "<<t<<"  sp: "<<i<<endl<<"B(i,t)\n"<<B(i,t)<<endl;}
      if (debug == 81 && min(B(i,t)) <= Bthres(i)) {cout<<"t: "<<t<<"  sp: "<<i<<endl<<"B(i,t)\n"<<B(i,t)<<endl;}

      for (a=1; a<=nage(i); a++)
        {
        tmpBpen.initialize();
        if (debug == 81 && B(i,t,a) <= Bthres(i)) {cout<<"a: "<<a<<"  pre Bpen: "<<Bpen<<endl;}  //09-14-09
        B(i,t,a) = posfun(B(i,t,a),Bthres(i),tmpBpen);
        Bpen(i) += tmpBpen;
        if (debug == 81 && B(i,t,a) <= Bthres(i)) {cout<<"a: "<<a<<"      Bpen: "<<Bpen<<"     tmpBpen: "<<tmpBpen<<endl;}
        }  //end of prey age loop

      if (debug == 80) {cout<<"Post B(i,t)\n"<<B(i,t)<<endl<<endl; }
      if (debug == 81 && min(B(i,t)) <= Bthres(i)) {cout<<"Post B(i,t)\n"<<B(i,t)<<endl<<endl; }
      }  //end of year If Statement
    }  //end of species loop


FUNCTION calc_year1_constraint
  /*
  Function calculates a penalty, Ypen, that is added to the objective function in the calc_likelihood function.  This penalty ensures
  that species-specific Yr1 abundances approximate an exponential decline with increasing age.  Specifically, the function
  1) calculates the average species- and age-specific total mortality, avgZ, observed over the time series, 2) uses this avgZ and
  the predicted age-1 abundance in year-1, Yr1(i,1), to create a synthetic, theoretical cohort in year 1, thYr1, and 3) calculates the 
  penalty, Ypen, as the sum of the squared deviations between the theoretical, thYr1, and predicted, Yr1, Year-1 abundances.
  */

  dvar_matrix avgZ(1,nsp,1,nage);   //Average age-specific total mortality; averaged across years
  dvar_matrix thYr1(1,nsp,1,nage);  //Theoretical year-1 abundances
  thYr1.initialize();
  for (i = 1; i<=nsp; i++)
    {
    avgZ(i) = colMeans(Z(i));
    if (debug == 82) {cout<<"sp: "<<i<<endl<<"Z\n"<<Z(i)<<endl<<"avgZ\n"<<avgZ(i)<<endl;}
    //Set Age-1 abundance equal to predicted Yr1 abundance
    thYr1(i,1) = mfexp(Yr1(i,1));
    //Loop over ages to calculate remaining N's instead of using elem_prod because age, in this case, is recursive
    for (a =2; a<=nage(i); a++)
      {thYr1(i,a) = thYr1(i,a-1)*mfexp(-avgZ(i,a-1));}
    if (debug == 82) {cout<<"mfexp(Yr1(i,1)): "<<mfexp(Yr1(i,1))<<endl<<"thYr1\n"<<thYr1(i)<<endl;}

    //Calculate penalty, Ypen, as the sum of squared deviations between theoretical and predicted Year-1 abundance
    Ypen(i) = norm2(mfexp(Yr1(i)) - thYr1(i));
    if (debug == 82) {cout<<"Ypen(i)\n"<<Ypen(i)<<endl<<endl;}
    }  //end of species loop

  if (debug == 82) {cout<<"Exiting at end of Yr1_constraint function..."<<endl;  exit(-1);}


FUNCTION calc_length_props
  /*
  Function calculates predicted length proportions for commercial catches (CPlen_hat) and survey catches (SPlen_hat) from:
    1) Predicted age-proportions, and 2) Age-length key
  */
  for (i=1; i<=nsp; i++)
    {
    CPlen_hat(i) = Cprop_hat(i) * trans(ALKey(i));
    for (j = 1; j<=nFIC; j++)
      {
      SPlen_hat(i)(j) = Sprop_hat(i)(j) * trans(ALKey(i));
      }  // end of survey loop
    }  // end of species loop


FUNCTION calc_recruitment_penalty
  /*
  Function calculates a penalty, Rpen, that is added to the objective function in the calc_likelihood function if the coefficient of variation for the recruitment of any species becomes greater than a threshold value, Rthres
  */
  int devyr = 3;

  for (i = 1; i<=nsp; i++)
    {
    if (debug == 83) {cout<<"sp: "<<i<<endl<<"Age1(i)\n"<<Age1(i)<<endl;
                               cout<<"std_dev: "<<std_dev(Age1(i))<<endl;
                               cout<<"avg: "<< mean(Age1(i))<<endl;  }
    dvariable cvAge1 = std_dev(Age1(i))/mean(Age1(i));
    if (debug == 83) {cout<<"cvAge1:  "<<cvAge1<<endl;}    
    if  (cvAge1>Rthres(i))  
      { Rpen(i) = .01*square(cvAge1-Rthres(i));  }

    int devfyr = lyr(i)-devyr+1;
    dvar_vector dev_tmp = dAge1(i)(devfyr,lyr(i));
    Rpen2(i) = norm2(dev_tmp);
    //cout<<"dev_tmp\n"<<dev_tmp<<endl;
    //cout<<"devfyr\n"<<devfyr<<endl;
    //exit(-1);

    if (debug == 83) {cout<<"Rpen(i):  "<<Rpen(i)<<endl<<endl;}
    }  //end of species loop

    if (debug == 83) {cout<<"Rpen\n"<<Rpen<<endl<<endl;}    
    


FUNCTION calc_negative_loglikelihood
  /*
  Function calculates the likelihood components for each data source:
    1.  Total commercial catch, in weight; TotC and TotC_hat
    2.  Commercial catch proportions-at-length; CPlen and CPlen_hat
    3.  Total fishery-independent survey catch in number/tow; toteps
    4.  Survey catch proportions-at-length; SPlen and SPlen_hat
    5.  Predator food-habits in proportion, by weight, at-age; FH and FH_hat

  For data sources 1 and 3, a lognormal distribution is assumed.  For data sources 2, 4, and 5, two different types of residuals were
  calculated - A) assuming a multinomial distribution, and B) assuming a multinomial distribution but first subtracting the best possible value of the residuals assuming a perfect fit (so that now the objective function goes to zero)
  In prevoius versions, I also tried an arcsine squareroot transformation, but I was unable to take the arcsine of the value 1
  when the element is part of a differentiable array.  Accordingly, I had to subtract a tiny number (5.e-8) from the observed 
  and predicted FH so that I could still take the arcsine of the element even if the predicted FH value == 1, but I was not
  confortable with this modification.

  The function also calculates the likelihood components for each of the initial parameters that represent deviations from means (dYr1,
  dAge1, and dFt).  The likelihood components for each of these parameters ensure that all species-specific deviations sum to zero.
  */

   //Total commercial catch component
       //init_matrix TotC(1,nsp,fyr,lyr);                //Total commercial catch; metric tons
       //matrix TotC_hat(1,nsp,fyr,lyr);                //Total commerical catch *in weight*; summed over ages
       //dvar_vector TCres(1,nsp);
  TCres.initialize();                                          //Lognormal distribution

   //Catch-at-age proportions component
       //3darray CPlen(1,nsp,fyr,lyr,1,nlen);         //Proportion-at-length of the commercical catch
       //3darray CPlen_hat(1,nsp,fyr,lyr,1,nlen);  //Commercial catch proportions-at-age
       //dvar_vector CPmulti(1,nsp);
  CPmulti.initialize();                                        //Multinomial distribution

  //Total survey catch component
      //3darray TotFIC(1,nsp,1,nFIC,fyr,lyr);        //Sum of observed FIC  
      //3darray TotFIC_hat(1,nsp,1,nFIC,fyr,lyr); //Total survey (FIC) catch in number/tow; summed over ages
      //matrix TSresA(1,nsp,1,nFIC);      //Objective function component: Total survey catch, method A
      //matrix TSresB(1,nsp,1,nFIC);      //Objective function component: Total survey catch, method B
  TSresA.initialize();  TSresB.initialize();           //Lognormal distribution

  //FIC, Survey catch-at-length proportions component
      //init_3darray sumFIC(1,nsp,1,nFIC,fyr,lyr);        //Sum of observed FIC  
      //3darray SPlen(1,nsp,1,nFIC,fyr,lyr,1,nlen);        //Proportion-at-length of the survey catch, FIC
      //3darray SPlen_hat(1,nsp,1,nFIC,fyr,lyr,1,nlen);  //Survey catch proportions-at-length
  SPmulti.initialize();                                                    //Multinomial distribution

  //Penalties
  tBpen.initialize();                                           //Total Biomass penalty for each species
  tYpen.initialize();                                          //Total Yr1 penalty for each species
  tRpen.initialize();                                           //Total Recruitment penalty for each species
  tRpen2.initialize();                                           //Total Recruitment penalty for each species

  //Deviations component
      //init_bounded_matrix dYr1(1,nsp,1,nage,-5,5,2);        //Age-specific deviations in Year 1 N's, log space
      //init_bounded_matrix dAge1(1,nsp,fyr+1,lyr,-5,5,2);   //Annual deviations in recruits, log space
      //init_bounded_matrix dFt(1,nsp,1,nFt,-5,5,2);             //Annual deviations in species-specific fishing mortality rates in yrs where TotC > 0
  Devs.initialize();

  //Total objective function value, some species-specific;
  ofvsp.initialize();  ofvsp_ideal.initialize();  ofv_ideal.initialize();

  //Species loop for each likelihood component
  for (i = 1; i<=nsp; i++)
    {

    //Total commercial catch component
    TCresid(i) = log(TotC(i)+o) - log(TotC_hat(i)+o);
    TCres(i) = TCwt(i) *norm2(TCresid(i));
      /* ******In future datasets, double check to make sure TotC and TotC_hat units are the same! --
            Best way to check is to make sure the numbers are of the same order of magnitude!******* */
    if (debug == 62) {cout<<"TotC_hat\n"<<TotC_hat(i)+o<<endl<<"log TotC_hat\n"<<log(TotC_hat(i)+o)<<endl;}
    if (debug == 60) {cout<<"Completed lognormal TotC,  sp="<<i<<endl;}

    //Total survey catch component; two possible ways to calculate TotFIC residuals
    //B
    TSresid(i) = log(TotFIC(i)+o) - log(TotFIC_hat(i)+o);
    for (j=1; j<=nFIC; j++)
      {
      //A
      TSresA(i)(j) = TSwt(i)(j) *norm2(toteps(i)(j));
      //B
      TSresB(i)(j) = TSwt(i)(j) *norm2(TSresid(i)(j));
      }  //end of nFIC loop
    if (debug == 62) {cout<<"TotFIC_hat\n"<<TotFIC_hat(i)+o<<endl<<"log TotFIC_hat\n"<<log(TotFIC_hat(i)+o)<<endl;}
    if (debug == 60) {cout<<"Completed lognormal TotFIC,  sp="<<i<<endl;}
    
    for (t = fyr(i); t<=lyr(i); t++)
      {

      //Catch-at-age proportions component
      if (sumC(i,t) != 0)  //only calculate CPmulti if sum of obs CAA > 0, i.e. if age samples were taken in a particular year
        {
        //Multinomial distribution
        dvar_vector CPvec = CPwt(i) *elem_prod( CPlen(i,t)+p, log(CPlen_hat(i,t)+p) );
          //It does not matter if CPwt is in- or outisde of the summation
        CPmulti(i) -= sum(CPvec);
        if (debug == 63) {cout<<"t: "<<t<<endl
                                         <<"CPlen_hat\n"<<CPlen_hat(i,t)+p<<endl<<"log CPlen_hat\n"<<log(CPlen_hat(i,t)+p)<<endl;}
        }  //end of CAA proportions component
      if (debug == 61) {cout<<"Completed multinomial CPlen,  sp="<<i<<"  yr="<<t<<endl;}

      //Survey catch-at-age proportions component
      for (j=1; j<=nFIC; j++)
        {
        if (sumFIC(i,j,t) != 0)
          {
          //Multinomial distribution
          dvar_vector SPvec = SPwt(i)(j) *elem_prod( SPlen(i,j,t)+p,log(SPlen_hat(i,j,t)+p) );
          SPmulti(i,j) -= sum(SPvec);
          if (debug == 64) {cout<<"t: "<<j<<"  j: "<<t<<endl
                                           <<"SPlen_hat\n"<<SPlen_hat(i,j,t)+p<<endl<<"log SPlen_hat\n"<<log(SPlen_hat(i,j,t)+p)<<endl;}
          }  //end of FIC proportions component
        }  //end of nFIC loop
      if (debug == 61) {cout<<"Completed multinomial SPlen,  sp="<<i<<"  yr="<<t<<endl;}

      }  //end of year loop
    if (debug == 60) {cout<<"Completed multinomial CPlen and SPlen,  sp="<<i<<endl;}

    //Putting all of the likelihood components together
    //Commercial and survey catch data
    ofvsp_ideal(i) = TCres(i) + (CPmulti(i) - CPideal(i)) + sum(TSresA(i)) + sum( (SPmulti(i) - SPideal(i)) );
    ofvsp(i) = TCres(i) + CPmulti(i) + sum(TSresA(i)) + sum(SPmulti(i));

    }  //end of species loop

  //Penalty terms
  tBpen = elem_prod(Bwt,Bpen);  //Total biomass penalty for each species
  tYpen = elem_prod(Ywt,Ypen);  //Total Yr1 penalty for each species
  tRpen = elem_prod(Rwt,Rpen);  //Total recruitment penalty for each species
  tRpen2 = elem_prod(Rwt2,Rpen2);  //Total recruitment penalty for each species
  //Diagnostics for biomass and recruitment penalty terms
  for (i=1; i<=nsp; i++)
    {
    if(Bpen(i)> 0 ) {
    nBpen(i)++;  lBpen(i) = nf;
    cout<<"\nSpecies: "<<i<<"  B penalty = "<<Bpen(i)<<"  tBpen =  "<<tBpen(i)<<"  nBpen = "<<nBpen(i)<<endl;}
    }  //end of species loop
 
  //Deviation terms
  Devs = 1.e6*( square(rowsum(dAge1)) + square(rowsum(dFt)) );
  if (debug == 66 && sum(Devs) > 0.001) {
    cout<<"Age1 devs: "<<square(rowsum(dAge1))<<endl<<"Ft devs: "<<square(rowsum(dFt))<<endl; }

  //Total objective function value
  ofv = sum(ofvsp) + sum(tBpen) + sum(tYpen) + sum(tRpen) + sum(tRpen2) + sum(Devs);
  ofv_ideal = sum(ofvsp_ideal) + sum(tBpen) + sum(tYpen) + sum(tRpen) + sum(tRpen2) + sum(Devs);

  if (debug == 66) {
    cout<<"ofv: "<<ofv<<"   phase: "<<current_phase()<<"   nf: "<<nf<<endl;
    cout<<setprecision(10)<<"ofvsp: "<<ofvsp<<"   ofv_ideal: "<<ofv_ideal<<"   ofvsp_ideal: "<<ofvsp_ideal<<endl
      <<"TotC.Res: "<<TCres<<endl
      <<"TotFIC.ResA\n"<<TSresA<<endl<<"TotFIC.ResB\n"<<TSresB<<endl
      <<"CAA.proportions: "<<CPmulti<<endl<<"FIC.proportions\n"<<SPmulti<<endl
      <<"CAA.proportions.revised: "<<CPmulti - CPideal<<endl<<"FIC.prop.revised\n"<<SPmulti - SPideal<<endl;
    cout<<"tYpen: "<<tYpen<<"  tRpen: "<<tRpen<<"  tRpen2: "<<tRpen2<<"   Devs: "<<Devs<<endl;
    cout<<"Last nf with Bpen: "<<lBpen<<"  # of nf's with Bpen: "<<nBpen<<endl<<endl;  }


BETWEEN_PHASES_SECTION
  cout<<"\nbetween_phases_section, moving to phase: "<<current_phase()<<endl<<endl;


FINAL_SECTION
  cout<<"\nofv\n"<<ofv<<endl<<"nf\n"<<nf<<endl;
  cout<<setprecision(10)<<"ofvsp\n"<<ofvsp<<endl<<"ofv_ideal\n"<<ofv_ideal<<endl<<"ofvsp_ideal\n"<<ofvsp_ideal<<endl
    <<"TotC.Res\n"<<TCres<<endl
    <<"TotFIC.ResA\n"<<TSresA<<endl<<"TotFIC.ResB\n"<<TSresB<<endl
    <<"CAA.proportions\n"<<CPmulti<<endl<<"FIC.proportions\n"<<SPmulti<<endl
    <<"CAA.proportions.revised\n"<<CPmulti - CPideal<<endl<<"FIC.prop.revised\n"<<SPmulti - SPideal<<endl;
    cout<<"tYpen: "<<tYpen<<"  tRpen: "<<tRpen<<"  tRpen2: "<<tRpen2<<"   Devs: "<<Devs<<endl;
  cout<<"Last nf with Bpen\n"<<lBpen<<endl<<"# of nf's with Bpen: "<<nBpen<<endl<<endl;


REPORT_SECTION
  report<<"mgc\n"<<objective_function_value::pobjfun->gmax<<endl;
  report<<"ofv\n"<<ofv<<endl;
  report<<"ofv_ideal\n"<<ofv_ideal<<endl;
  report<<"nsp\n"<<nsp<<endl;
  report<<"nf\n"<<nf<<endl;
  report<<"nFIC\n"<<nFIC<<endl;

  imatrix yrs(1,nsp,fyr,lyr);
  for (i = 1; i<=nsp; i++)  { yrs(i).fill_seqadd(fyr(i),1); }
  report<<"yrs\n"<<yrs<<endl;
  report<<"nyr\n"<<nyr<<endl;
  report<<"nage\n"<<nage<<endl;
  report<<"FICnseg\n"<<nseg<<endl;
  report<<"M1nseg\n"<<nMseg<<endl;
  report<<"M1yr\n"<<M1yr<<endl;
  report<<"agePR\n"<<agePR<<endl;
  report<<"ageFR\n"<<ageFR<<endl;
  report<<"ficFR\n"<<ficFR<<endl;
  report<<"FICfage\n"<<FICfage<<endl;
  report<<"FIClage\n"<<FIClage<<endl;
  report<<"FICs.lage\n"<<FICs_lage<<endl;

  report<<"Ft\n"<<Ft<<endl;
  report<<"nFt\n"<<nFt<<endl;
  report<<"nsel\n"<<nsel<<endl;
  report<<"FICnsel\n"<<FICnsel<<endl;

  report<<"obs.sumC\n"<<sumC<<endl;
  report<<"obs.totc\n"<<TotC<<endl;
  report<<"pred.totc\n"<<TotC_hat<<endl;
  report<<"resid.totc\n"<<TCresid<<endl;

  report<<"Age1\n"<<Age1<<endl;

  for (i=1; i<=nsp; i++)
    {
    report<<"FICmonth."<<i<<endl<<FICmon(i)<<endl;
    report<<"sel."<<i<<endl<<s(i)<<endl;
    report<<"FICsel."<<i<<endl<<FICs(i)<<endl;
    report<<"FIC.q."<<i<<endl<<mfexp(q(i))<<endl;
    report<<"FIC.yr."<<i<<endl<<FICyr(i)<<endl;
    report<<"M1."<<i<<endl<<M1seg(i)<<endl;

    report<<"obs.cplen."<<i<<endl<<CPlen(i)<<endl;
    report<<"pred.CAA."<<i<<endl<<C_hat(i)<<endl;
    report<<"pred.cprop."<<i<<endl<<Cprop_hat(i)<<endl;
    report<<"pred.cplen."<<i<<endl<<CPlen_hat(i)<<endl;

    report<<"N."<<i<<endl<<N(i)<<endl;
    report<<"Wt."<<i<<endl<<Wt(i)<<endl;
    report<<"iM2."<<i<<endl<<iM2(i)<<endl;
    report<<"M2."<<i<<endl<<M2(i)<<endl;

    report<<"obs.totFIC."<<i<<endl<<TotFIC(i)<<endl;
    report<<"pred.totFIC."<<i<<endl<<TotFIC_hat(i)<<endl;
    report<<"toteps."<<i<<endl<<toteps(i)<<endl;
    report<<"resid.totFIC."<<i<<endl<<TSresid(i)<<endl;

    report<<"Age1CV."<<i<<endl<<(std_dev(Age1(i))/mean(Age1(i)))<<endl;

    for (j = 1; j<= nFIC; j++)
      {
      report<<"obs.splen."<<j<<".sp."<<i<<endl<<SPlen(i)(j)<<endl;
      report<<"pred.FIC."<<j<<".sp."<<i<<endl<<FIC_hat(i)(j)<<endl;
      report<<"pred.sprop."<<j<<".sp."<<i<<endl<<Sprop_hat(i)(j)<<endl;
      report<<"pred.splen."<<j<<".sp."<<i<<endl<<SPlen_hat(i)(j)<<endl;
      }  //end of nFIC loop

    }  //end of species loop

  report<<"TCwt\n"<<TCwt<<endl;
  report<<"TSwt\n"<<TSwt<<endl;
  report<<"CPwt\n"<<CPwt<<endl;
  report<<"SPwt\n"<<SPwt<<endl;
  report<<"nBpen\n"<<nBpen<<endl;
  report<<"lBpen\n"<<lBpen<<endl;
  report<<"Bthres\n"<<Bthres<<endl;
  report<<"Bwt\n"<<Bwt<<endl;
  report<<"Ywt\n"<<Ywt<<endl;
  report<<"Rwt\n"<<Rwt<<endl;
  report<<"Rthres\n"<<Rthres<<endl;

  report<<"TCres\n"<<TCres<<endl;
  report<<"CPmulti\n"<<CPmulti<<endl;
  report<<"CPideal\n"<<CPideal<<endl;
  report<<"TSresA\n"<<TSresA<<endl;
  report<<"TSresB\n"<<TSresB<<endl;
  report<<"SPmulti\n"<<SPmulti<<endl;
  report<<"SPideal\n"<<SPideal<<endl;
  report<<"tYpen\n"<<tYpen<<endl;
  report<<"tRpen\n"<<tRpen<<endl;
  report<<"tRpen2\n"<<tRpen2<<endl;
  report<<"Devs\n"<<Devs<<endl;
  report<<"ofvsp\n"<<ofvsp<<endl;
  report<<"ofvsp_ideal\n"<<ofvsp_ideal<<endl;



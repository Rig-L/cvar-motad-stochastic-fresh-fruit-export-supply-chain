param zp>=0; #Number of potential management zones
param ns>=0; #Number of sample points in the field	
#param tl:=21;
param tl:=14; #Number of days to harvest
param nw:=27; #Number of scenarios

#########################################
#########################################
#Sets and indices
set Z:=1..zp; #Set of potential management zones
set Ce:=1..ns; #Sample points in the field
set T:=0..tl; #Planning horizon
set D:=1..tl; #Planning horizon
set S:=1..tl; #Planning horizon

set C:= {'BR','CA','CHI','COL','COS','ECU','US','EU','JA'}; #Set of customers
set P:= {'SA','VA','SV'}; #Set of ports
set W:= 1..nw; #Set of scenarios

#########################################
#########################################
#Parameters

#Waste rates
param delta:=0.02; #Waste rate of perishable product in orchards
param beta:=0.05; #Waste rate of unprocessed product harvested
param gamma:=0.02; #Waste rate of processed perishable product harvested

#Workers
param twmax:=66; #Maximum number of temporary workers
param pp:=80; #Productivity of permanent workers
param pt :=40; #Productivity of temporary workers
param sp :=500000; #Cost of permanent workers (monthly salary)
param st :=35000; #Cost of temporary workers (dayly salary)

#Transport
param vcmax :=3320; #Maximum transport capacity
param tcz{Z} >=0; #Transportation cost to zone z
param ttz{Z} >=0; #Transportation time to zone z
param th{T} >=0; #Available travel time in period t

#Inventory
param mscmax :=570000; #Maximum inventory capacity
param ic:=58; #Unit inventory cost

#Processing
param cpmax; #Maximum processing capacity [kg] per day
param pc :=316; #Unit cost of processing fruits

#Zone parameters
param alfa:=0.8; #Required homogeneity level
param totalvariance>=0; #Variance of the index in the field
param vcoef{Z}>=0; #Variance of the index for the potential zone z
param a{Z,Ce}>=0; #Binary parameter that takes a value of 1 if sample point m belongs to management zone z and 0 otherwise
param pz{Z} >=0; #Production of potential management zone z

#Suitability binary variable, that it takes a value 1 if management zone z can be harvested in period t and 0 otherwise
param bz{z in Z, t in T} := 
    if (t >= thsz[z] and t <= tf[z]) then 1
    else 0;

param tf{Z} >=0; #Final harvest period for the potential zone z
param thsz {Z} >=0; #Initial harvest period for the potential zone z
param c{Z}>=0; #Cost of choosing the potential zone z

#Demand
param dem{C}; #Demand of customer c
param ps{C, W}; #Unit profit of exported fruit sold at the optimal price to customer c
param pm{cl in C, w in W}:= ps[cl,w]*(10/10);

#Ports
param av{P,C}; #Binary parameter that takes a value of 1 if customer c can be supplied from port p and 0 otherwise
param cman:=200000; #Port selection cost for export operations

#Quality parameters
param qu:=0.0526; #Perishability rate
param sa:= 1.57894; #Increase in fruit shelf–life integrity due to processing
param days{C}; #Shelf–life days to meet the demand of customer c with the necessary quality
param qsu{t in T,d in D: t>=d}:= 1 - qu*(t-d); #Fruit shelf–life level in period t, given that it was harvested in period d
param qsp{t in T,d in D,s in S: t>=s}:= (1 - qu*(s-d))+sa-qu*(t-s); #Fruit shelf–life level in period t, given that it was harvested in period d and processed in period s
param qc{C}>=0; #Quality required by customer c
param prob{W}>=0; #Probability of scenario w
param rq{P,W}; #Logistic delay factor in port p in scenario w

#MOTAD Parameters
param aversion:= 0.90; #Risk aversion


#########################################
#########################################
#FIST STAGE VARIABLES
  
#HARVEST RELATED VARIABLES
var Q{Z} binary >=0;
var XU{Z,T,C} >= 0;
var XP{Z,T,C} >= 0;
var Y{Z,T} >= 0;

#OPERATIONAL RELATED VARIABLES
var NT{Z,T} integer >= 0;
var PW integer >= 0;
var TW{T} integer>= 0;
var CTW{T} integer >= 0;
var HTW{T} integer >= 0;
var PP{T,D,C} >= 0;
var IUP{T,D,C} >= 0;
var ITP{T,D,C} >= 0;
var IPP{T,D,S,C} >= 0;
var SP{P,D,S,C} >= 0;
var SU{P,D,C} >= 0;

#PORT RELATED VARIABLES
var CR{P,C} >= 0;
var CT{P,C} binary;

#DEMAND INDICATORS
var UD{C,W} >= 0;
var MD{C,W} >= 0;

#########################################
#SECOND STAGE VARIABLES

#Differentiation of shippers by price
var SHIPU{P,D,C,W} >= 0;
var SHIPP{P,D,S,C,W} >= 0;
var SMINU{P,D,C,W} >= 0;
var SMINP{P,D,S,C,W} >= 0;
var FQU{P,D,C,W} binary;
var FQP{P,D,S,C,W} binary;
var FAU{P,D,C,W} >= 0;
var SOU{P,D,C,W} >= 0;
var FAP{P,D,S,C,W} >= 0;
var SOP{P,D,S,C,W} >= 0;

#DEVIATIONS
var DEV{W} >= 0;

#OBJECTIVE FUNCTION
maximize FO:
(1-aversion)*(-sum{t in T,d in D, cl in C}(pc*PP[t,d,cl]) - 
sum{t in T,z in Z}(tcz[z]*NT[z,t])-                      
(sp*PW)-                                                 
sum{t in T}(st*TW[t])-                                  
sum{t in T, d in D, cl in C}(IUP[t,d,cl]+ITP[t,d,cl])*ic-  
sum{t in T, d in D, s in S, cl in C}(IPP[t,d,s,cl])*ic- 
sum {p in P, cl in C }(cman*CT[p,cl])-
sum{z in Z} (c[z])*Q[z]+                            
sum{re in W}(prob[re]*(
	sum{p in P, cl in C, d in D}(SHIPU[p,d,cl,re])*ps[cl,re]+
	sum{p in P, cl in C, d in D, s in S}(SHIPP[p,d,s,cl,re])*ps[cl,re]+
	sum{p in P, cl in C, d in D}(SMINU[p,d,cl,re]*pm[cl,re]*(qsu[tl,d]-rq[p,re]*days[cl]*qu)/qc[cl])+
	sum{p in P, cl in C, d in D, s in S}(SMINP[p,d,s,cl,re]*pm[cl,re]*(qsp[tl,d,s]-rq[p,re]*days[cl]*qu)/qc[cl])
)))
-aversion*sum{w in W}(prob[w]*DEV[w]);
	
#MOTAD CONSTRAINTS
s.t. desv1{w in W}:
	DEV[w]	>= (-sum{t in T,d in D, cl in C}(pc*PP[t,d,cl]) - 
			sum{t in T,z in Z}(tcz[z]*NT[z,t])-                      
			(sp*PW)-                                                 
			sum{t in T}(st*TW[t])-                                  
			sum{t in T, d in D, cl in C}(IUP[t,d,cl]+ITP[t,d,cl])*ic-  
			sum{t in T, d in D, s in S, cl in C}(IPP[t,d,s,cl])*ic- 
			sum {p in P, cl in C }(cman*CT[p,cl])-
			sum{z in Z} (c[z])*Q[z]+                            
			sum{re in W}(prob[re]*(
				sum{p in P, cl in C, d in D}(SHIPU[p,d,cl,re])*ps[cl,re]+
				sum{p in P, cl in C, d in D, s in S}(SHIPP[p,d,s,cl,re])*ps[cl,re]+
				sum{p in P, cl in C, d in D}(SMINU[p,d,cl,re]*pm[cl,re]*(qsu[tl,d]-rq[p,re]*days[cl]*qu)/qc[cl])+
				sum{p in P, cl in C, d in D, s in S}(SMINP[p,d,s,cl,re]*pm[cl,re]*(qsp[tl,d,s]-rq[p,re]*days[cl]*qu)/qc[cl])
			)))
			-
			(-sum{t in T,d in D, cl in C}(pc*PP[t,d,cl]) - 
			sum{t in T,z in Z}(tcz[z]*NT[z,t])-                      
			(sp*PW)-                                                 
			sum{t in T}(st*TW[t])-                                  
			sum{t in T, d in D, cl in C}(IUP[t,d,cl]+ITP[t,d,cl])*ic-  
			sum{t in T, d in D, s in S, cl in C}(IPP[t,d,s,cl])*ic- 
			sum {p in P, cl in C }(cman*CT[p,cl])-
			sum{z in Z} (c[z])*Q[z]+                            
			
			sum{p in P, cl in C, d in D}(SHIPU[p,d,cl,w])*ps[cl,w]+
			sum{p in P, cl in C, d in D, s in S}(SHIPP[p,d,s,cl,w])*ps[cl,w]+
			sum{p in P, cl in C, d in D}(SMINU[p,d,cl,w]*pm[cl,w]*(qsu[tl,d]-rq[p,w]*days[cl]*qu)/qc[cl])+
			sum{p in P, cl in C, d in D, s in S}(SMINP[p,d,s,cl,w]*pm[cl,w]*(qsp[tl,d,s]-rq[p,w]*days[cl]*qu)/qc[cl])
			);

#########################################
#FIRST STAGE CONSTRAINTS

###
s.t. RR1z{ce in Ce}:
	sum{z in Z}(a[z,ce]*Q[z])=1;

s.t. RR2z:
	sum{z in Z}(((1-alfa)*totalvariance + vcoef[z])*Q[z])<=(1-alfa)*totalvariance*ns;

s.t. RR2 {z in Z, t in T: t <thsz[z]}:
	Y[z,t] = 0;

s.t. RR4 {z in Z, t in T: t > tf[z]}:
	Y[z,t] = 0;
	
s.t. RR5a {z in Z,t in thsz[z]+1..tf[z]}:
     Y[z,t]=Y[z,t-1]*(1-delta)-sum{cl in C}(XU[z,t,cl]+XP[z,t,cl]);
     
s.t. RR5b {z in Z}:
     Y[z,thsz[z]]=pz[z]*Q[z]-sum{cl in C}(XU[z,thsz[z],cl]+XP[z,thsz[z],cl]);

s.t. RR6 {z in Z, t in 1..tl}:
   	sum{cl in C}(XU[z,t,cl]+XP[z,t,cl]) <= pz[z]*Q[z]*bz[z,t];
   
###
s.t. RR8 {t in T}:
    sum{z in Z,cl in C}(XU[z,t,cl]+XP[z,t,cl])<= pp*PW + pt*TW[t];

s.t. RR9 {t in T}:
    TW[t] <=twmax;

s.t. RR10 {t in 1..tl}:
    TW[t] = TW[t-1]-HTW[t-1]+CTW[t];

s.t. RR11:
	TW[0]=0;

s.t. RR12:
	HTW[0]=0;

###
s.t. RR19a {t in 1..tl, cl in C}:
   IUP[t,t,cl]=sum{z in Z}(XU[z,t,cl]);
   
s.t. RR19b {t in 1..tl, d in 1..tl, cl in C: t>d}:
   IUP[t,d,cl]=IUP[t-1,d,cl]*(1-beta);
   
s.t. RR19c {t in 1..tl, d in 1..tl, cl in C: t<d}:
   IUP[t,d,cl]=0;
   
s.t. RR19d {d in 1..tl, cl in C}:
   sum{p in P}(SU[p,d,cl])=IUP[tl,d,cl];
   
s.t. RR19e {d in 1..tl, cl in C}:
   qsu[tl,d]*IUP[tl,d,cl]>=qc[cl]*IUP[tl,d,cl];   

###
s.t. RR21a {t in 1..tl, cl in C}:
    ITP[t,t,cl]=sum{z in Z}(XP[z,t,cl])-PP[t,t,cl];
    
s.t. RR21b {t in 1..tl, d in 1..tl, cl in C: t>d}:
    ITP[t,d,cl]=ITP[t-1,d,cl]*(1-beta)-PP[t,d,cl];

s.t. RR21c {t in 1..tl, d in 1..tl, cl in C: t<d}:
    ITP[t,d,cl]=0;           
    
###
s.t. RR23a {t in 1..tl, d in 1..tl, cl in C}:
	IPP[t,d,t,cl]=PP[t,d,cl];

s.t. RR23b {t in 1..tl, d in 1..tl, s in 1..tl, cl in C: t>s}:
	IPP[t,d,s,cl]=IPP[t-1,d,s,cl]*(1-gamma);
	
s.t. RR23c {t in 1..tl, d in 1..tl, s in 1..tl, cl in C: t<s}:
	IPP[t,d,s,cl]=0;
	
s.t. RR23d {d in D, s in S, cl in C}:
	sum{p in P}SP[p,d,s,cl] = IPP[tl,d,s,cl];
	
s.t. RR23e {d in D, s in S, cl in C}:
	qsp[tl,d,s]*IPP[tl,d,s,cl]>=qc[cl]*IPP[tl,d,s,cl];

###
s.t. RR16 {d in 1..tl, cl in C}: 
       IUP[0,d,cl]=0;

s.t. RR17 {d in 1..tl, cl in C}: 
      ITP[0,d,cl]=0;
      
s.t. RR18 {d in 1..tl,s in 1..tl, cl in C}: 
      IPP[0,d,s,cl]=0;

###
s.t. RR25 {t in T}: 
     sum{d in D, cl in C}(IUP[t,d,cl]+ ITP[t,d,cl]+sum{s in S}(IPP[t,d,s,cl])) <= mscmax;

s.t. RR26 {t in T}: 
     sum{d in D, cl in C}(PP[t,d,cl])<= cpmax;

s.t. RR27 {t in 1..tl, d in 1..tl, cl in C: t<d}:
     PP[t,d,cl]=0;

s.t. RR28 {z in Z, t in T}:
  sum{cl in C}(XU[z,t,cl]+XP[z,t,cl])<= NT[z,t]*vcmax;

s.t. RR29 {t in 1..tl}:
	sum{z in Z}(NT[z,t]*ttz[z])<= th[t];

###
s.t. RR30 {cl in C}:
	sum{p in P}(CT[p,cl]) = 1; 
   
s.t. RRde {p in P, cl in C}: 
	sum{d in D}(SU[p,d,cl]+sum{s in S}(SP[p,d,s,cl]))=dem[cl]*CT[p,cl];
   
s.t. RR31b {p in P, cl in C}: 
   	CT[p,cl] <= av[p,cl];  

s.t. RR31a {p in P, cl in C}:
	CR[p,cl] = sum{d in D}(SU[p,d,cl]+sum{s in S}(SP[p,d,s,cl])); 	

#########################################
#SECOND STAGE CONSTRAINTS

###
s.t. RR33a {p in P, d in D, cl in C, w in W}: 
   SU[p,d,cl]*(qsu[tl,d]-rq[p,w]*days[cl]*qu)+FAU[p,d,cl,w]-SOU[p,d,cl,w] = qc[cl]*SU[p,d,cl]; 

s.t. RR33b {p in P, d in D, cl in C, w in W}: 
   FAU[p,d,cl,w]<=FQU[p,d,cl,w]*dem[cl]*qc[cl]; 
   
s.t. RR33b2 {p in P, d in D, cl in C, w in W}: 
   SOU[p,d,cl,w]<=(1-FQU[p,d,cl,w])*dem[cl]*qc[cl];   

s.t. RR33c {p in P, d in D, cl in C, w in W}: 
   SHIPU[p,d,cl,w]+SMINU[p,d,cl,w] = SU[p,d,cl]; 

s.t. RR33d {p in P, d in D, cl in C, w in W}: 
   SHIPU[p,d,cl,w] <= (1-FQU[p,d,cl,w])*dem[cl];

s.t. RR33e {p in P, d in D, cl in C, w in W}: 
   SMINU[p,d,cl,w] <= (FQU[p,d,cl,w])*dem[cl];
   
###
s.t. RR34a {p in P, d in D, s in S, cl in C, w in W}: 
   SP[p,d,s,cl]*(qsp[tl,d,s]-rq[p,w]*days[cl]*qu)+FAP[p,d,s,cl,w]-SOP[p,d,s,cl,w]=qc[cl]*SP[p,d,s,cl]; 

s.t. RR34b {p in P, d in D, s in S, cl in C, w in W}: 
   FAP[p,d,s,cl,w]<=FQP[p,d,s,cl,w]*dem[cl]*qc[cl];

s.t. RR34b2 {p in P, d in D, s in S, cl in C, w in W}: 
   SOP[p,d,s,cl,w]<=(1-FQP[p,d,s,cl,w])*dem[cl]*qc[cl];
  
s.t. RR34c {p in P, d in D, s in S, cl in C, w in W}: 
   SHIPP[p,d,s,cl,w]+SMINP[p,d,s,cl,w] = SP[p,d,s,cl]; 

s.t. RR34d {p in P, d in D, s in S, cl in C, w in W}: 
   SHIPP[p,d,s,cl,w] <= (1-FQP[p,d,s,cl,w])*dem[cl]; 

s.t. RR34e {p in P, d in D, s in S, cl in C, w in W}: 
   SMINP[p,d,s,cl,w] <= FQP[p,d,s,cl,w]*dem[cl];
   
###
s.t. RR35 {cl in C, w in W}: 
   MD[cl,w]=sum{p in P, d in D, s in S}(SHIPP[p,d,s,cl,w])+sum{p in P, d in D}(SHIPU[p,d,cl,w]); 
   
s.t. RR36 {cl in C, w in W}: 
   MD[cl,w]+UD[cl,w]=dem[cl]; 
   




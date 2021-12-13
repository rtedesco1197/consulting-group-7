/*STAT 488-002 Consulting
Mixed Models to investigate the effects of aromatase inhibitors on fish fertility
Joelle Strom
Updated: Dec. 13, 2021*/


%MACRO SortAndPlot(DSName, Var);
	Proc Sort data=&DSName;
	   By descending TRT TNK d;
	Run;
	 
	Proc Sgpanel data=&DSName;
	   Panelby TRT / columns=3 onepanel;
	   Series x=d y=Pred / group=TNK groupLC=TRT break lineattrs=(pattern=solid);
	   Keylegend / type=linecolor title="";
	Run;
	
	Proc Sgpanel data=&DSName;
		panelby TRT / columns=3 onepanel;
		vline d / response=&Var group=TRT stat=mean limitstat=stderr;
	Run;
%MEND;

Proc Import Datafile='/home/u49579191/Consulting/LetrozolelDat.csv' DBMS=CSV OUT=letrozolel;
 GETNAMES=YES;
 RUN;

data letrozolel;
set letrozolel;
t = d;                       /* discrete copy of time */
T1 = ifn(d<=10, 0, d - 10);  /* knot at day 10 for PWL analysis */
T2 = ifn(d<=14, 0, d - 14);  /* knot at day 14 for PWL analysis */
sqrtegg = sqrt(Eggs_gFem);
PctFert = PctFert/100;
PctVBL = PctVBL/100;
PctFert = (PctFert * 432 + 0.5) / 433;
PctVBL = (PctVBL * 432 + 0.5) / 433;
run;

proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model Eggs_gFem = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;

proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model Eggs_gFem = TRT d T1 TRT*d TRT*T1/ s chisq distribution=gaussian;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;

proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model sqrtegg = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;
/*Transforming the response variable decreases AIC over a model with untransformed variable and
also over the piece-wise regression (with knot at day 10)*/

proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model sqrtegg = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual type=CS subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;

proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model sqrtegg = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual type=AR(1) subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;
/*Log L decreased by changing variance matrix structure to AR(1)*/

proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model Eggs_gFem = TRT d TRT*d/ s chisq distribution=poisson;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;
/*Attempted poisson link but did not converge*/

/*FINAL MODEL*/
proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model sqrtegg = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual type=AR(1) subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
   output out=MixedOut1 pred=Pred;
run;

%SortAndPlot(MixedOut1, sqrtegg);

proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctFert = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;

proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctFert = TRT d T2 TRT*d TRT*T2/ s chisq distribution=beta;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;
/*Standard regression performs better than piece-wise with knot at day 14*/

proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctFert = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual type=CS subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;

proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctFert = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual type=AR(1) subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;
/*Changing variance matrix structure increases AIC*/

/*FINAL MODEL*/
proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctFert = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
   output out=MixedOut2 pred=Pred;
run;

%SortAndPlot(MixedOut2, PctFert);

proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctVBL = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
run;

proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctVBL = TRT d T2 TRT*d TRT*T2/ s chisq distribution=beta;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
run;
/*Piece-wise regression with knot at day 14 improves upon standard regression*/

proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctVBL = TRT d T2 TRT*d TRT*T2/ s chisq distribution=beta;
   random t / residual type=CS subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
run;

proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctVBL = TRT d T2 TRT*d TRT*T2/ s chisq distribution=beta;
   random t / residual type=AR(1) subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
run;
/*Log L reduced by changing variance matrix structure to AR(1)*/

/*FINAL MODEL*/
proc glimmix data=letrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctVBL = TRT d T2 TRT*d TRT*T2/ s chisq distribution=beta;
   random t / residual type=AR(1) subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   output out=MixedOut3 pred=Pred;
run;

%SortAndPlot(MixedOut3, PctVBL);

Proc Import Datafile='/home/u49579191/Consulting/AnastrozolelDat.csv' DBMS=CSV OUT=anastrozolel;
 GETNAMES=YES;
 RUN;

data anastrozolel;
set anastrozolel;
t = d;                       /* discrete copy of time */
T1 = ifn(d<=10, 0, d - 10);  /* knot at day 10 for PWL analysis */
sqrtegg = sqrt(Eggs_gFem);
PctFert = PctFert/100;
PctVBL = PctVBL/100;
PctFert = (PctFert * 432 + 0.5) / 433;
PctVBL = (PctVBL * 432 + 0.5) / 433;
run;

proc glimmix data=anastrozolel plots=all;
   class TNK t TRT(ref='control');
   model Eggs_gFem = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;

proc glimmix data=anastrozolel plots=all;
   class TNK t TRT(ref='control');
   model Eggs_gFem = TRT d T1 TRT*d TRT*T1/ s chisq distribution=gaussian;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;

proc glimmix data=anastrozolel plots=all;
   class TNK t TRT(ref='control');
   model sqrtegg = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;
/*Transforming the response variable decreases AIC over a model with untransformed variable and
also over the piece-wise regression (with knot at day 10)*/

proc glimmix data=anastrozolel plots=all;
   class TNK t TRT(ref='control');
   model sqrtegg = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual type=CS subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;

proc glimmix data=anastrozolel plots=all;
   class TNK t TRT(ref='control');
   model sqrtegg = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual type=AR(1) subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;
/*Log L decreased by changing variance matrix structure to AR(1)*/

/*FINAL MODEL*/
proc glimmix data=anastrozolel plots=all;
   class TNK t TRT(ref='control');
   model sqrtegg = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
   output out=MixedOut4 pred=Pred;
run;

%SortAndPlot(MixedOut4, sqrtegg);

proc glimmix data=anastrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctFert = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;

/*No obvious knot*/

proc glimmix data=anastrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctFert = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual type=CS subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;

proc glimmix data=anastrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctFert = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual type=AR(1) subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;
/*Changing variance matrix structure increases AIC*/

/*FINAL MODEL*/
proc glimmix data=anastrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctFert = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
   output out=MixedOut5 pred=Pred;
run;

%SortAndPlot(MixedOut5, PctFert);

proc glimmix data=anastrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctVBL = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;

/*No obvious knot*/

proc glimmix data=anastrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctVBL = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual type=CS subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;

proc glimmix data=anastrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctVBL = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual type=AR(1) subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;
/*No parameter estimates with different variance structures*/

/*FINAL MODEL*/
proc glimmix data=anastrozolel plots=all;
   class TNK t TRT(ref='control');
   model PctVBL = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
   output out=MixedOut6 pred=Pred;
run;

%SortAndPlot(MixedOut6, PctVBL);

Proc Import Datafile='/home/u49579191/Consulting/ExemestaneDat.csv' DBMS=CSV OUT=exemestane;
 GETNAMES=YES;
 RUN;

data exemestane;
set exemestane;
t = d;                       /* discrete copy of time */
sqrtegg = sqrt(Eggs_gFem);
PctFert = PctFert/100;
PctVBL = PctVBL/100;
PctFert = (PctFert * 432 + 0.5) / 433;
PctVBL = (PctVBL * 432 + 0.5) / 433;
run;

proc glimmix data=exemestane plots=all;
   class TNK t TRT(ref='control');
   model Eggs_gFem = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;

/*No obvious knot*/

proc glimmix data=exemestane plots=all;
   class TNK t TRT(ref='control');
   model sqrtegg = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;
/*Transforming the response variable decreases AIC over a model with untransformed variable*/

proc glimmix data=exemestane plots=all;
   class TNK t TRT(ref='control');
   model sqrtegg = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual type=CS subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;

proc glimmix data=exemestane plots=all;
   class TNK t TRT(ref='control');
   model sqrtegg = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual type=AR(1) subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
run;
/*Log L decreased by changing variance matrix structure to AR(1)*/

/*FINAL MODEL*/
proc glimmix data=exemestane plots=all;
   class TNK t TRT(ref='control');
   model sqrtegg = TRT d TRT*d/ s chisq distribution=gaussian;
   random t / residual type=AR(1) subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;    /* each subject gets its own intercept */
   output out=MixedOut7 pred=Pred;
run;

%SortAndPlot(MixedOut7, Eggs_gFem);

proc glimmix data=exemestane plots=all;
   class TNK t TRT(ref='control');
   model PctFert = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;

/*No obvious knot*/

proc glimmix data=exemestane plots=all;
   class TNK t TRT(ref='control');
   model PctFert = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual type=CS subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;

proc glimmix data=exemestane plots=all;
   class TNK t TRT(ref='control');
   model PctFert = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual type=AR(1) subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;
/*Changing variance matrix structure to compound symmetry decreases AIC*/

/*FINAL MODEL*/
proc glimmix data=exemestane plots=all;
   class TNK t TRT(ref='control');
   model PctFert = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual type=CS subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
   output out=MixedOut8 pred=Pred;
run;

%SortAndPlot(MixedOut8, PctFert);

proc glimmix data=exemestane plots=all;
   class TNK t TRT(ref='control');
   model PctVBL = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;

/*No obvious knot*/

proc glimmix data=exemestane plots=all;
   class TNK t TRT(ref='control');
   model PctVBL = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual type=CS subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;

proc glimmix data=exemestane plots=all;
   class TNK t TRT(ref='control');
   model PctVBL = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual type=AR(1) subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
run;
/*Log L reduced by changing variance matrix structure to compound symmetry*/

/*FINAL MODEL*/
proc glimmix data=exemestane plots=all;
   class TNK t TRT(ref='control');
   model PctVBL = TRT d TRT*d/ s chisq distribution=beta;
   random t / residual type=CS subject=TNK;   /* measurements are repeated for subjects */
   random intercept / subject=TNK;          /* each subject gets its own intercept */
   Nloptions maxiter=100 tech=nrridg;
   output out=MixedOut9 pred=Pred;
run;

%SortAndPlot(MixedOut9, PctVBL);



















**End;
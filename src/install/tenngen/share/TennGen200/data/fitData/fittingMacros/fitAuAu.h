#ifndef FITAUAU_H
#define FITAUAU_H


#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TROOT.h"
#include <iostream>
#include "TGraphErrors.h"
#include "TFile.h"
#include "Riostream.h"
#include <fstream>
#include <string>
#include <stdio.h>
#include <TSystem.h>
#include <sstream>


Double_t getRatio(Double_t* pT, Double_t* params){
	Double_t pt		= pT[0];
	Double_t A 	= params[0];
	Double_t B 	= params[1];
	Double_t C	= params[2];
	//Double_t D		= params[3];
	return A*TMath::Power(pt,B)*TMath::Exp(C*pt);
	//return A*TMath::Log(B*pt)+C;
	//return TMath::Log(C*pt)+ A*TMath::Power(pt,B);
}
Double_t getdNdptOverptIntegrand(Double_t* rad, Double_t* par){
	// (dN/dpt)/pt= r*dr*mt*I0((pt*sinh(rho))/T)*K1((mt*cosh(rho))/T)
	// rho=arctanh(beta); beta=betaMax*(r/R)^n
	//Double_t pT		= pt[0];
	Double_t r 			= rad[0];
	// parameters used to fit: mass, beta, temp, n, norm
	Double_t mass 		= par[0];// not mT
	Double_t pt			= par[1];
	Double_t betaMax 	= par[2];
	Double_t temp 		= par[3];
	Double_t n 			= par[4];
	

	Double_t beta = betaMax*TMath::Power(r,n);
	if(beta > 0.99999999999999999999) beta = 0.99999999999999999999;

	Double_t mT 	= TMath::Sqrt(mass*mass + pt*pt);
	Double_t rho0 	= TMath::ATanH(beta);
	Double_t avoidFPE = pt*TMath::SinH(rho0)/temp;
	if(avoidFPE > 700.) avoidFPE = 700.;
	Double_t bk1arg = mT*TMath::CosH(rho0)/temp;
	Double_t integrand = r*mT*TMath::BesselI0(avoidFPE)*TMath::BesselK1(bk1arg);

	return integrand;
}



Double_t getdNdpt(Double_t* pT, Double_t* params){

	TF1* dNdptOverptIntegrandFunc = new TF1("integrandFunc", getdNdptOverptIntegrand, 0, 1, 5 );
	
	Double_t pt		= pT[0];
	Double_t mass 	= params[0];
	Double_t beta 	= params[1];
	Double_t temp 	= params[2];
	Double_t n 		= params[3];
	Double_t norm	= params[4];
	
	dNdptOverptIntegrandFunc->SetParameters(mass,pt,beta,temp,n);

	Double_t dNdptOverpt 	= dNdptOverptIntegrandFunc->Integral(0,1);
	
	Double_t dNdpt_normalized			= pt* norm * dNdptOverpt;
	gSystem->ProcessEvents();
	gROOT->Reset();
	return dNdpt_normalized;
}


Double_t getdNdyIntegrand(Double_t* myPt, Double_t* par){
	Double_t pt   = myPt[0]; // x-axis of integration
	Double_t mass = par[0];
	
	return getdNdpt(myPt,par);
}

Double_t Multiplicity(Double_t* X, Double_t* par){
    Double_t x      = X[0];
    return 0.0569104-0.00129485*x+1.45861e-05*x*x-8.83828e-08*TMath::Power(x,3)+3.10488e-10*TMath::Power(x,4)-6.54566e-13*TMath::Power(x,5)+8.25612e-16*TMath::Power(x,6)-5.9263e-19*TMath::Power(x,7)+2.10337e-22*TMath::Power(x,8)-2.40376e-26*TMath::Power(x,9);
}

string concatenateHistoname(string centStr,string pName,string colSp,string colEn){
	string initText = "cent";
	string undScr = "_";//underscore
	//string enUnits = "GeV";
	string addedString = initText+centStr+undScr+pName+undScr+colSp+undScr+colEn;//+enUnits;
	return addedString; //type: const char*: to be done later
}

#endif

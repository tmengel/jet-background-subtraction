#include <iostream>
#include <string>
#include "TKey.h"
#include <sstream>
#include <fstream>
#include "fitAuAu.h"
#include "TClass.h"


using namespace std;

Double_t getdNdptOverptIntegrand(Double_t* rad, Double_t* par);// not used
Double_t getRatio(Double_t* pT, Double_t* params);



// main function:
int FitData(){
	
	TFile* myFile = new TFile("rhicAuAu200.root");
	TIter next(myFile->GetListOfKeys());
	TKey* mikey;
	TH1D* h;
	TCanvas* c1;
	TClass* class1;
	TF1* funcBGBW;
	TF1* ratioFit;
	
	int breakOutForTesting =0;
	int stop =70; // breakOut after this many iterations (if achieved); default: 140
	cout << "Flag" << endl;
	Double_t A_array[5][7];
	Double_t B_array[5][7];
	Double_t C_array[5][7];


	Double_t mass_array[5][7];
	Double_t beta_array[5][7];
	Double_t temp_array[5][7];
	Double_t n_array[5][7];
	Double_t norm_array[5][7];
	


	Double_t P0 [3][6][3];
	Double_t P1 [3][6][3];
	Double_t P2 [3][6][3];
	Double_t P3 [3][6][3];
	Double_t P4 [3][6][3];
	int ratioNum =0;
	while((mikey=(TKey*)next())){
	  
		///cout << "Histo iter: " << breakOutForTesting+1 << endl;
		class1 = gROOT->GetClass(mikey->GetClassName());
		if(!class1->InheritsFrom("TH1")){
			delete class1;
			mikey->DeleteBuffer();
			continue;
		}

		c1 = new TCanvas(); // a la Rademakers

		
		funcBGBW = new TF1("getdNdpt",getdNdpt,0.00000000000001,10.,5); // actually has 5 parameters
		ratioFit = new TF1("getRatio",getRatio, 0.000000000001, 10, 3);

		gPad->SetLogy();
	    gStyle->SetOptFit(1111);// display fit parameters; customizable
		//	gStyle->SetOptDate();// display date (at bottom left)
		gROOT-> SetBatch(kTRUE);// save canvases without displaying them
		c1->Update();

		// read histogram object for current iteration of key:
		h = (TH1D*)mikey->ReadObj();
		string histoName = h->GetName();
		Double_t collEn = 200.;// initialize
		string fitType = histoName.substr(0,5);
		cout << fitType <<endl;
		string particleID = histoName.substr(6,3);// starting position in array:6, 3 chars total
		string centrality = histoName.substr(histoName.length()-1,1);// starting position in array:4, 1 char total
		
		if(fitType == "Spect"){
			//------------ FIT PT SPECTRA -----------------//
			Double_t mass; // in GeV
			if		(particleID=="pi-"||particleID=="pi+")
					{mass = 0.13957;}
			else if	(particleID=="ka-"||particleID=="ka+")
					{mass = 0.49368;}
			else if	(particleID=="pro" ||particleID=="pba")
					{mass = 0.93827; }
			else if(particleID == "ALL"){mass = 0.49368; }
			int partIndx; 
			if		(particleID=="pi-"){partIndx =0;}
			else if	(particleID=="pi+"){partIndx =1;}
			else if	(particleID=="ka-"){partIndx =2;}
			else if (particleID=="ka+"){partIndx =3;}
			else if	(particleID=="pba"){partIndx =4;}
			else if (particleID=="pro"){partIndx =5;}		
			else if(particleID == "ALL"){partIndx =6; }
			

			// //------------- Begin BGBW fit --------------------------//
			if(histoName == "Spect_proton_cent0"||histoName=="Spect_proton_cent1"){
				funcBGBW->SetParameters(mass,0.5,0.5,1.,100.);
				cout << "TEST CHECK ALT PARMS";
			}
			else if(histoName == "Spect_pbar_cent1"){
				funcBGBW->SetParameters(mass,0.5,0.2,1.,10000.);
			}	
			else{funcBGBW->SetParameters(mass,0.95,0.05,0.1,1000000.);}

			funcBGBW->SetParNames("mass","beta (c)","temp","n","norm");
			funcBGBW->SetParLimits(1,0.35,0.999999999999999999999);//param 1
			funcBGBW->FixParameter(0,mass);// mass in GeV
					
			ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(20000);
			TFitResultPtr r = h->Fit("getdNdpt","S","",0.00000000000001,10.);
			
			
			Double_t chi2Prob = r->Prob();	
			cout << "chi-sq prob: " << chi2Prob <<endl;

			h->SetMaximum(5*(h->GetMaximum()));
			h-> GetXaxis()->SetRangeUser(0.,10.);
			TString xlabel = "p_{T}";
			TString ylabel = "#frac{d^{2}N}{dydp_{T}}";
			h-> SetXTitle(xlabel);
			h-> SetYTitle(ylabel);
					
			Double_t beta 			= funcBGBW->GetParameter(1);
			Double_t temp 			= funcBGBW->GetParameter(2);
			Double_t n	  			= funcBGBW->GetParameter(3);
			Double_t norm 			= funcBGBW->GetParameter(4);
			Double_t betaErr 		= funcBGBW->GetParError(1);
			Double_t tempErr 		= funcBGBW->GetParError(2);
			Double_t nErr 			= funcBGBW->GetParError(3);
			Double_t normErr 		= funcBGBW->GetParError(4);

			int centInd = std::atoi(centrality.c_str());
			mass_array[centInd][partIndx] = mass;
	        beta_array[centInd][partIndx] = beta;
	        temp_array[centInd][partIndx]= temp;
	        n_array[centInd][partIndx] = n;
	        norm_array[centInd][partIndx] = norm;
			//------------- end BGBW fit ----------------------------
		
			//-------- Find integrals left and right of data points -------//
			funcBGBW			 	-> SetParameters(mass,beta,temp,n,norm);
			//-- end- output results to file------------------------
			c1->Update();
			Double_t chi2BGBW = funcBGBW->GetChisquare();
			Double_t nDFBGBW = funcBGBW->GetNDF();
			Double_t p2 = funcBGBW->GetParameter(2);
			Double_t e2 = funcBGBW->GetParError(2);
			c1->Update();
			string imgPathAndName = "./fittedPlots/"+histoName+".png";
					//c1 -> SaveAs("./fittedPlots/trial1.png");
			TImage *png = TImage::Create();// FIXME try to use canvas method instead of png object
			png->FromPad(c1);
			const char* imgPathAndNameConstCharPtr = imgPathAndName.c_str();
			png->WriteImage(imgPathAndNameConstCharPtr);
			
			//mikey->DeleteBuffer();// works!

		}
		else if(fitType=="ratio"){
			ratioFit->SetParameters(10.,1.0,0.1);
			ratioFit->SetParNames("A","B","C");
			
			string ratioType = histoName.substr(6,3)+histoName.substr(10,3);
			if(ratioType == "pba/pi") {
				ratioType = "pbar"+histoName.substr(11,3);
				ratioNum = 5;
				}
			else if(ratioType == "pi-pi+"){
				ratioNum = 0;
				}
			else if(ratioType == "pba/pr") {
				ratioType = "pbar"+histoName.substr(11,3);
				ratioNum = 2;
			}
			else if(ratioType == "ka-ka+"){
				ratioNum=1;
			}
			else if(ratioType == "ka-pi-"){
				ratioNum=3;
			}
			else if(ratioType == "ka+pi+"){
				ratioNum=4;
			}
			else if(ratioType == "propi+"){
				ratioNum=6;
			}
			//A_array[std::atoi(centrality.c_str())][ratioNum] = 
			cout << ratioType <<endl;
			ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(20000);
			TFitResultPtr r = h->Fit("getRatio","S","",0.00000000000001,10.);
			Double_t chi2Prob = r->Prob();	
			cout << "chi-sq prob: " << chi2Prob <<endl;
			//h->Fit("pol6");
			h->SetMaximum(5*(h->GetMaximum()));
			h-> GetXaxis()->SetRangeUser(0.,10.);
			TString xlabel = "p_{T}";
			//TString ylabel = ratioType;
			h-> SetXTitle(xlabel);
			h-> SetYTitle(ratioType.c_str());
			c1->Update();
			histoName = "ratio_"+ratioType+"_cent"+centrality;
			
			Double_t A 			= ratioFit->GetParameter(0);
			Double_t B 			= ratioFit->GetParameter(1);
			Double_t C		= ratioFit->GetParameter(2);
			Double_t Aerr 			= ratioFit->GetParameter(0);
			Double_t Berr 		= ratioFit->GetParError(1);
			Double_t Cerr 		= ratioFit->GetParError(2);
			A_array[std::atoi(centrality.c_str())][ratioNum] = A;
			B_array[std::atoi(centrality.c_str())][ratioNum] = B;
			C_array[std::atoi(centrality.c_str())][ratioNum] = C;
			
			//------------- end BGBW fit ----------------------------
		
			//-------- Find integrals left and right of data points -------//
			ratioFit			 	-> SetParameters(A,B,C);
			//-- end- output results to file------------------------
			c1->Update();
			Double_t chi2BGBW = ratioFit->GetChisquare();
			Double_t nDFBGBW = ratioFit->GetNDF();
			Double_t p2 = ratioFit->GetParameter(2);
			Double_t e2 = ratioFit->GetParError(2);
			c1->Update();
			string imgPathAndName = "./fittedPlots/"+histoName+".png";
					//c1 -> SaveAs("./fittedPlots/trial1.png");
			TImage *png = TImage::Create();// FIXME try to use canvas method instead of png object
			png->FromPad(c1);
			const char* imgPathAndNameConstCharPtr = imgPathAndName.c_str();
			png->WriteImage(imgPathAndNameConstCharPtr);
			
			//mikey->DeleteBuffer();// works!

		}
				
		 mikey->DeleteBuffer();// works!
		
		

		gSystem->ProcessEvents();
		delete h;
		delete funcBGBW;
		delete ratioFit;
		delete c1;	// Rademakers
		//delete mikey; // FIXME 9 segmentation violation
		//delete class1; // segmentation violation
	}

	ifstream 			in;
	string harmonicName;
	in.open(Form("./datafiles-formatted/phenixAuAu200_partVn.dat"));
	string 				myString;// temporary string to hold ifstream instance
	string				collidingSpecies;/////////////
	string 				collisionEnergy;
	string 				collisionEnergyStr;////////////////////
	string 				particleName;// eg. pi+
	string ratioCent;
	Float_t 			myFloat; //temporary float to hold energy value
	Double_t			myDouble; // temporary double to hold values from histo row
	for(int i=0;i<4;i++){in>>myString;}
	for(int graphNum=0;graphNum<53;graphNum++){
		in>>myString;
		in>>particleName;
		in>>harmonicName;
		in>>ratioCent;
		std::vector<Double_t> ptVec;
		std::vector<Double_t> vnVec;
		std::vector<Double_t> ptErr;
		std::vector<Double_t> vnErr;
		while(in>>myDouble){
			ptVec.push_back(myDouble);
			in>>myDouble;
			vnVec.push_back(myDouble);
			in>>myDouble;
			ptErr.push_back(myDouble);
			in>>myDouble;
			vnErr.push_back(myDouble);
		}
		Double_t PT[ptVec.size()];
		Double_t VN[ptVec.size()];
		Double_t PTERR[ptVec.size()];
		Double_t VNERR[ptVec.size()];
		int numpt = 0;
		for(int j =0;j<ptVec.size();j++){
		PT[j] = ptVec[j];
		VN[j]= vnVec[j];
		PTERR[j]=ptErr[j];
		VNERR[j]=vnErr[j];
		numpt++;
		}
		//in>>myString;
		TGraphErrors *g = new TGraphErrors(numpt,PT,VN,PTERR,VNERR);
		string graphName = particleName+"_"+harmonicName+"_"+ratioCent;

		c1 = new TCanvas(); // a la Rademakers
		g->Draw();
		g->Fit("pol4");
		g->SetTitle(graphName.c_str());
		TF1 *fit = g->GetFunction("pol4");
		//g->GetXAxis()->SetXTitle("p_{T}");
		//g->SetYTitle("v_n");
		int partind;
		int centind;
		int harmind;
		if(particleName == "pion"){partind = 0;}
		if(particleName == "kaon"){partind = 1;}
		if(particleName == "proton"){partind = 2;}

		if(ratioCent == "cent0"){centind = 0;}
		if(ratioCent == "cent1"){centind = 1;}
		if(ratioCent == "cent2"){centind = 2;}
		if(ratioCent == "cent3"){centind = 3;}
		if(ratioCent == "cent4"){centind = 4;}
		if(ratioCent == "cent5"){centind = 5;}

		if(harmonicName == "v2"){harmind = 0;}
		if(harmonicName == "v3"){harmind = 1;}
		if(harmonicName == "v4"){harmind = 2;}
		
		P0 [partind][centind][harmind] = fit->GetParameter(0);
		P1 [partind][centind][harmind] = fit->GetParameter(1);
		P2 [partind][centind][harmind] = fit->GetParameter(2);
		P3 [partind][centind][harmind] = fit->GetParameter(3);
		P4 [partind][centind][harmind] = fit->GetParameter(4);
		// P1 [3][6][3];
		// P2 [3][6][3];
		// P3 [3][6][3];
		// P4 [3][6][3];
		c1->Update();
		string imgPathAndName = "./fittedGraphs/"+graphName+".png";
				//c1 -> SaveAs("./fittedPlots/trial1.png");
		TImage *png = TImage::Create();// FIXME try to use canvas method instead of png object
		png->FromPad(c1);
		const char* imgPathAndNameConstCharPtr = imgPathAndName.c_str();
		png->WriteImage(imgPathAndNameConstCharPtr);
		
		//mikey->DeleteBuffer();// works!
		cout << "Made " << graphName << " Graph number: " << graphNum <<endl;
		in.clear();
		delete c1;
		delete g;
		//g->Write();
	}
	Double_t Multip[10];
	in.close();
	in.open(Form("./datafiles-formatted/starAuAu200_multiplicity.dat"));
	for(int i=0;i<4;i++){in>>myString;}
	for(int graphNum=0;graphNum<1;graphNum++){
		in>>myString;
		std::vector<Double_t> ptVec;
		std::vector<Double_t> vnVec;
		std::vector<Double_t> ptErr;
		std::vector<Double_t> vnErr;
		while(in>>myDouble){
			ptVec.push_back(myDouble);
			in>>myDouble;
			vnVec.push_back(myDouble);
			in>>myDouble;
			ptErr.push_back(myDouble);
			in>>myDouble;
			vnErr.push_back(myDouble);
		}
		Double_t PT[ptVec.size()];
		Double_t VN[ptVec.size()];
		Double_t PTERR[ptVec.size()];
		Double_t VNERR[ptVec.size()];
		int numpt = 0;
		for(int j =0;j<ptVec.size();j++){
		PT[j] = ptVec[j];
		VN[j]= vnVec[j];
		VNERR[j]=TMath::Sqrt(ptErr[j]*ptErr[j]+vnErr[j]*vnErr[j]);
		PTERR[j]=0;
		numpt++;
		}
		//in>>myString;
		TGraphErrors *g = new TGraphErrors(numpt,PT,VN,PTERR,VNERR);
		string graphName = "multiplicity";
		c1 = new TCanvas(); // a la Rademakers
		g->Draw();
		c1->SetLogy();
		g->Fit("pol9");
		
		TF1 *fit = g->GetFunction("pol9");
		for(int i=0;i<10;i++){
			Multip[i]=fit->GetParameter(i);
		}
		TH1F* h1 = new TH1F(graphName.c_str(),graphName.c_str(),100, 0.,750);
		h1->FillRandom("multiplicity",10000);
		g->SetTitle(graphName.c_str());
		//g->GetXAxis()->SetXTitle("p_{T}");
		//g->SetYTitle("v_n");
		c1->Update();
		string imgPathAndName = "./fittedGraphs/"+graphName+".png";
				//c1 -> SaveAs("./fittedPlots/trial1.png");
		TImage *png = TImage::Create();// FIXME try to use canvas method instead of png object
		png->FromPad(c1);
		const char* imgPathAndNameConstCharPtr = imgPathAndName.c_str();
		png->WriteImage(imgPathAndNameConstCharPtr);
		
		//mikey->DeleteBuffer();// works!
		cout << "Made " << graphName << " Graph number: " << graphNum <<endl;
		in.clear();
		delete c1;
		delete g;
		//g->Write();
	}

	in.close();

	ofstream fout ("PrintedFunctions.h");
      if(!fout.is_open()){
        cout << "ERROR OPENING OUTPUT FILE" <<endl;
        return 0;
      }

	fout << "#ifndef PRINTEDFUNCTIONS_H\n#define PRINTEDFUNCTIONS_H\n\n";
	fout << "float RatioFunction(float pT, int Cent, int RatioType){\n";
	fout << "\tfloat pt\t= pT;\n";
	fout << "\tint cent\t= Cent;\n\tint ratio\t= RatioType;\n\n";
	for(int i=0;i<4;i++){
		if(i==0){fout<<"\tif(cent=="<<i<<"){\n";}
		else{fout<<"\telse if(cent=="<<i<<"){\n";}
		for(int j=0;j<7;j++){
			if(j==0){fout<< "\t\tif(ratio=="<<j<<"){ return "<< A_array[i][j] <<"*TMath::Power(pt,"<< B_array[i][j] << ")*TMath::Exp("<< C_array[i][j] <<"*pt); }\n";}
			else{fout<< "\t\telse if(ratio=="<<j<<"){ return "<< A_array[i][j] <<"*TMath::Power(pt,"<< B_array[i][j] << ")*TMath::Exp("<< C_array[i][j] <<"*pt); }\n";}			
		}
		fout << "\t}\n";
	}
	fout <<"}\n\n // End of RatioFunction\n\n";


	fout <<"float OmegaFunction(float pT, int Cent){\n";
	fout << "\tfloat pt\t= pT;\n";
	fout << "\tfloat cent\t= Cent;\n";
	for(int i=0;i<4;i++){
		if(i==0){fout<<"\tif(cent=="<<i<<"){\n";}
		else{fout<<"\telse if(cent=="<<i<<"){\n";}
		fout <<"\t return (1+"<< A_array[i][3] <<"*TMath::Power(pt,"<< B_array[i][3] << ")*TMath::Exp("<< C_array[i][3] <<"*pt)+" << A_array[i][5] <<"*TMath::Power(pt,"<< B_array[i][5] << ")*TMath::Exp("<< C_array[i][5] <<"*pt)+("<< A_array[i][0] <<"*TMath::Power(pt,"<< B_array[i][0] << ")*TMath::Exp("<< C_array[i][0] <<"*pt))*(1+"<< A_array[i][4] <<"*TMath::Power(pt,"<< B_array[i][4] << ")*TMath::Exp("<< C_array[i][4] <<"*pt)+" << A_array[i][6] <<"*TMath::Power(pt,"<< B_array[i][6] << ")*TMath::Exp("<< C_array[i][6] <<"*pt)));\n";
		fout << "\t}\n";
	}
	fout <<"}\n\n // End of OmegaFunction\n\n";


	fout << "Double_t Multiplicity(Double_t* X, Double_t* par){"<<endl;
			fout << "\tDouble_t x      = X[0];" << endl;
			fout<< "\treturn " << Multip[0];
			for(int i=1;i<10;i++){
					fout << "+" <<Multip[i];
					for(int j =0;j<i;j++){
						fout << "*x";
					}
			}
			fout << ";\n}\n // End of multiplicity\n\n";
		

	fout << "float HarmonicFunction(float pT, int Cent, int Harmonic, int KF){"<<endl;
	fout << "\tfloat pt;\n\tint part;\n\tif(pT > 4.0){  pt = 4.0;  }\n\telse if(pt < 0.5){ return 0.0; }\n\telse{ pt = pT; }\n\n\tint cent\t= Cent;\n\tif(KF ==211||KF==-211){ part =0; }\n\telse if(KF ==321||KF==-321){ part =1; }\n\telse if(KF ==2212||KF==-2212){ part =2; }\n\telse{ part =0; }\n\n\tint Vn\t=Harmonic;\n" << endl;
	for(int i=0;i<6;i++){
		if(i==0){fout<<"\tif(cent=="<<i<<"){\n";}
		else{fout<<"\telse if(cent=="<<i<<"){\n";}

		for(int j=0;j<3;j++){
			if(j==0){fout<<"\t\tif(part=="<<j<<"){\n";}
			else{fout<<"\t\telse if(part=="<<j<<"){\n";}
			for(int k=0;k<3;k++){
				if(k==0){fout<<"\t\t\tif(Vn=="<<k<<"){";}
				else{fout<<"\t\t\telse if(Vn=="<<k<<"){";}
				fout <<" return "<< P0[j][i][k] <<"+"<<P1[j][i][k]<<"*pt+"<< P2[j][i][k]<<"*pt*pt+"<<P3[j][i][k]<<"*pt*pt*pt+" <<P4[j][i][k] <<"*pt*pt*pt*pt;}"<<endl;
			}
			fout <<"\t\t}\n";

			
		}
		fout <<"\t}\n";
	
		
	}
	fout << "}"<<endl;

	fout << "Double_t HarmonicFunction(Double_t *pT, Int_t* par){"<<endl;
	fout << "\tDouble_t pt;\n\tInt_t part = par[0];\n\tif(pT[0] > 4.0){  pt = 4.0;  }\n\telse if(pT[0] < 0.5){ return 0.0; }\n\telse{ pt = pT[0]; }\n\n\tInt_t cent\t= par[1];\n\t\n\n\tint Vn\t=par[2];\n" << endl;
	for(int i=0;i<6;i++){
		if(i==0){fout<<"\tif(cent=="<<i<<"){\n";}
		else{fout<<"\telse if(cent=="<<i<<"){\n";}

		for(int j=0;j<3;j++){
			if(j==0){fout<<"\t\tif(part=="<<j<<"){\n";}
			else{fout<<"\t\telse if(part=="<<j<<"){\n";}
			for(int k=0;k<3;k++){
				if(k==0){fout<<"\t\t\tif(Vn=="<<k<<"){";}
				else{fout<<"\t\t\telse if(Vn=="<<k<<"){";}
				fout <<" return "<< P0[j][i][k] <<"+"<<P1[j][i][k]<<"*pt+"<< P2[j][i][k]<<"*pt*pt+"<<P3[j][i][k]<<"*pt*pt*pt+" <<P4[j][i][k] <<"*pt*pt*pt*pt;}"<<endl;
			}
			fout <<"\t\t}\n";

			
		}
		fout <<"\t}\n";
	
		
	}
	fout << "}"<<endl;
	fout <<"\n\n float MASS[4][6] ={";
	for(int i =0; i<4;i++){
		for(int j =0;j<6;j++){
			if(j==0){fout<<"{ ";}
			fout  << mass_array[i][j];
			if(j!=5) {fout<< " , ";}
			if(j==5){ fout<<" } ";}
		}
		if(i!=3){fout << " , ";}
		if(i==3){fout << " };\n";}
	}


	fout <<"\n\n float BETA[4][6] ={";
	for(int i =0; i<4;i++){
		for(int j =0;j<6;j++){
			if(j==0){fout<<"{ ";}
			fout  << beta_array[i][j];
			if(j!=5) {fout<< " , ";}
			if(j==5){ fout<<" } ";}
		}
		if(i!=3){fout << " , ";}
		if(i==3){fout << " };\n";}
	}

	fout <<"\n\n float TEMP[4][6] ={";
	for(int i =0; i<4;i++){
		for(int j =0;j<6;j++){
			if(j==0){fout<<"{ ";}
			fout  << temp_array[i][j];
			if(j!=5) {fout<< " , ";}
			if(j==5){ fout<<" } ";}
		}
		if(i!=3){fout << " , ";}
		if(i==3){fout << " };\n";}
	}

	fout <<"\n\n float NPAR[4][6] ={";
	for(int i =0; i<4;i++){
		for(int j =0;j<6;j++){
			if(j==0){fout<<"{ ";}
			fout  << n_array[i][j];
			if(j!=5) {fout<< " , ";}
			if(j==5){ fout<<" } ";}
		}
		if(i!=3){fout << " , ";}
		if(i==3){fout << " };\n";}
	}

	fout <<"\n\n float NORM[4][6] ={";
	for(int i =0; i<4;i++){
		for(int j =0;j<6;j++){
			if(j==0){fout<<"{ ";}
			fout  << norm_array[i][j];
			if(j!=5) {fout<< " , ";}
			if(j==5){ fout<<" } ";}
		}
		if(i!=3){fout << " , ";}
		if(i==3){fout << " };\n";}
	}
	fout << "#endif\n";
	fout.close();
	//	cout << ";\n\t}"<<endl;
	//}// end of while loop to iterate through every key
	gObjectTable->Print();
	delete myFile;
	return 0;
}

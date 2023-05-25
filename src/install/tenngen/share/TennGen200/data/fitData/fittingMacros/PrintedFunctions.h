#ifndef PRINTEDFUNCTIONS_H
#define PRINTEDFUNCTIONS_H

float RatioFunction(float pT, int Cent, int RatioType){
	float pt	= pT;
	int cent	= Cent;
	int ratio	= RatioType;

	if(cent==0){
		if(ratio==0){ return 1.05232*TMath::Power(pt,0.0444398)*TMath::Exp(-0.0391436*pt); }
		else if(ratio==1){ return 0.973497*TMath::Power(pt,0.00746373)*TMath::Exp(-0.0231327*pt); }
		else if(ratio==2){ return 0.714843*TMath::Power(pt,-0.0765232)*TMath::Exp(0.0302362*pt); }
		else if(ratio==3){ return 0.368629*TMath::Power(pt,0.943152)*TMath::Exp(-0.172881*pt); }
		else if(ratio==4){ return 0.378045*TMath::Power(pt,0.902791)*TMath::Exp(-0.139721*pt); }
		else if(ratio==5){ return 0.430421*TMath::Power(pt,3.51419)*TMath::Exp(-1.18616*pt); }
		else if(ratio==6){ return 0.623285*TMath::Power(pt,3.61225)*TMath::Exp(-1.24021*pt); }
	}
	else if(cent==1){
		if(ratio==0){ return 0.999801*TMath::Power(pt,-0.0249828)*TMath::Exp(0.0081137*pt); }
		else if(ratio==1){ return 0.963395*TMath::Power(pt,-0.0157204)*TMath::Exp(-0.0115677*pt); }
		else if(ratio==2){ return 0.746027*TMath::Power(pt,-0.0585174)*TMath::Exp(0.0057852*pt); }
		else if(ratio==3){ return 0.357072*TMath::Power(pt,0.933124)*TMath::Exp(-0.164785*pt); }
		else if(ratio==4){ return 0.369955*TMath::Power(pt,0.918328)*TMath::Exp(-0.143286*pt); }
		else if(ratio==5){ return 0.418842*TMath::Power(pt,3.43009)*TMath::Exp(-1.16914*pt); }
		else if(ratio==6){ return 0.571988*TMath::Power(pt,3.49185)*TMath::Exp(-1.18507*pt); }
	}
	else if(cent==2){
		if(ratio==0){ return 1.03665*TMath::Power(pt,0.0206189)*TMath::Exp(-0.0233864*pt); }
		else if(ratio==1){ return 0.979407*TMath::Power(pt,0.00336944)*TMath::Exp(-0.0266132*pt); }
		else if(ratio==2){ return 0.766358*TMath::Power(pt,-0.0469923)*TMath::Exp(-0.00814672*pt); }
		else if(ratio==3){ return 0.341224*TMath::Power(pt,0.885281)*TMath::Exp(-0.144538*pt); }
		else if(ratio==4){ return 0.346016*TMath::Power(pt,0.842213)*TMath::Exp(-0.101821*pt); }
		else if(ratio==5){ return 0.392903*TMath::Power(pt,3.21868)*TMath::Exp(-1.10371*pt); }
		else if(ratio==6){ return 0.536303*TMath::Power(pt,3.29701)*TMath::Exp(-1.12682*pt); }
	}
	else if(cent==3){
		if(ratio==0){ return 1.00893*TMath::Power(pt,-0.0147369)*TMath::Exp(0.00145578*pt); }
		else if(ratio==1){ return 0.961856*TMath::Power(pt,-0.0420376)*TMath::Exp(-0.00230289*pt); }
		else if(ratio==2){ return 0.777819*TMath::Power(pt,-0.0822537)*TMath::Exp(0.000508055*pt); }
		else if(ratio==3){ return 0.309516*TMath::Power(pt,0.775202)*TMath::Exp(-0.0930582*pt); }
		else if(ratio==4){ return 0.319566*TMath::Power(pt,0.776558)*TMath::Exp(-0.0740102*pt); }
		else if(ratio==5){ return 0.395983*TMath::Power(pt,2.91526)*TMath::Exp(-1.06911*pt); }
		else if(ratio==6){ return 0.517985*TMath::Power(pt,2.98984)*TMath::Exp(-1.07494*pt); }
	}
}

 // End of RatioFunction

float OmegaFunction(float pT, int Cent){
	float pt	= pT;
	float cent	= Cent;
	if(cent==0){
	 return (1+0.368629*TMath::Power(pt,0.943152)*TMath::Exp(-0.172881*pt)+0.430421*TMath::Power(pt,3.51419)*TMath::Exp(-1.18616*pt)+(1.05232*TMath::Power(pt,0.0444398)*TMath::Exp(-0.0391436*pt))*(1+0.378045*TMath::Power(pt,0.902791)*TMath::Exp(-0.139721*pt)+0.623285*TMath::Power(pt,3.61225)*TMath::Exp(-1.24021*pt)));
	}
	else if(cent==1){
	 return (1+0.357072*TMath::Power(pt,0.933124)*TMath::Exp(-0.164785*pt)+0.418842*TMath::Power(pt,3.43009)*TMath::Exp(-1.16914*pt)+(0.999801*TMath::Power(pt,-0.0249828)*TMath::Exp(0.0081137*pt))*(1+0.369955*TMath::Power(pt,0.918328)*TMath::Exp(-0.143286*pt)+0.571988*TMath::Power(pt,3.49185)*TMath::Exp(-1.18507*pt)));
	}
	else if(cent==2){
	 return (1+0.341224*TMath::Power(pt,0.885281)*TMath::Exp(-0.144538*pt)+0.392903*TMath::Power(pt,3.21868)*TMath::Exp(-1.10371*pt)+(1.03665*TMath::Power(pt,0.0206189)*TMath::Exp(-0.0233864*pt))*(1+0.346016*TMath::Power(pt,0.842213)*TMath::Exp(-0.101821*pt)+0.536303*TMath::Power(pt,3.29701)*TMath::Exp(-1.12682*pt)));
	}
	else if(cent==3){
	 return (1+0.309516*TMath::Power(pt,0.775202)*TMath::Exp(-0.0930582*pt)+0.395983*TMath::Power(pt,2.91526)*TMath::Exp(-1.06911*pt)+(1.00893*TMath::Power(pt,-0.0147369)*TMath::Exp(0.00145578*pt))*(1+0.319566*TMath::Power(pt,0.776558)*TMath::Exp(-0.0740102*pt)+0.517985*TMath::Power(pt,2.98984)*TMath::Exp(-1.07494*pt)));
	}
}

 // End of OmegaFunction

Double_t Multiplicity(Double_t* X, Double_t* par){
	Double_t x      = X[0];
	return 0.0569104+-0.00129485*x+1.45861e-05*x*x+-8.83828e-08*x*x*x+3.10488e-10*x*x*x*x+-6.54566e-13*x*x*x*x*x+8.25612e-16*x*x*x*x*x*x+-5.9263e-19*x*x*x*x*x*x*x+2.10337e-22*x*x*x*x*x*x*x*x+-2.40376e-26*x*x*x*x*x*x*x*x*x;
}
 // End of multiplicity

float HarmonicFunction(float pT, int Cent, int Harmonic, int KF){
	float pt;
	int part;
	if(pT > 4.0){  pt = 4.0;  }
	else if(pt < 0.5){ return 0.0; }
	else{ pt = pT; }

	int cent	= Cent;
	if(KF ==211||KF==-211){ part =0; }
	else if(KF ==321||KF==-321){ part =1; }
	else if(KF ==2212||KF==-2212){ part =2; }
	else{ part =0; }

	int Vn	=Harmonic;

	if(cent==0){
		if(part==0){
			if(Vn==0){ return -0.0143392+0.0759773*pt+-0.0154024*pt*pt+-0.000915513*pt*pt*pt+0.000264875*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.00558048+0.0136122*pt+0.0356842*pt*pt+-0.0164884*pt*pt*pt+0.00195171*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.00256332+-0.00483424*pt+0.0390473*pt*pt+-0.016438*pt*pt*pt+0.00205353*pt*pt*pt*pt;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0287711+0.0777838*pt+-0.0110908*pt*pt+-0.00251265*pt*pt*pt+0.000435634*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.00843333+0.00413273*pt+0.03708*pt*pt+-0.0140296*pt*pt*pt+0.00138853*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.0147847+0.0238129*pt+0.00149144*pt*pt+0.00122299*pt*pt*pt+-0.000583605*pt*pt*pt*pt;}
		}
		else if(part==2){
			if(Vn==0){ return 0.0139405+-0.0778337*pt+0.122313*pt*pt+-0.0441081*pt*pt*pt+0.00502149*pt*pt*pt*pt;}
			else if(Vn==1){ return 0.0204243+-0.0952528*pt+0.118049*pt*pt+-0.0375978*pt*pt*pt+0.00388916*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.0368641+-0.132059*pt+0.140577*pt*pt+-0.0465704*pt*pt*pt+0.00527664*pt*pt*pt*pt;}
		}
	}
	else if(cent==1){
		if(part==0){
			if(Vn==0){ return -0.0164604+0.119236*pt+-0.0140501*pt*pt+-0.00715798*pt*pt*pt+0.00137041*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.00615225+0.0205719*pt+0.0373663*pt*pt+-0.0176272*pt*pt*pt+0.00201875*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.00343526+-0.0156758*pt+0.058631*pt*pt+-0.0234438*pt*pt*pt+0.00258536*pt*pt*pt*pt;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0335949+0.108473*pt+0.00189956*pt*pt+-0.0119077*pt*pt*pt+0.00177362*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.0141866+0.0239085*pt+0.0233996*pt*pt+-0.0081926*pt*pt*pt+0.000430152*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.0178734+-0.0725294*pt+0.105302*pt*pt+-0.0391775*pt*pt*pt+0.00466704*pt*pt*pt*pt;}
		}
		else if(part==2){
			if(Vn==0){ return 0.0147481+-0.0885341*pt+0.16892*pt*pt+-0.0604128*pt*pt*pt+0.00661366*pt*pt*pt*pt;}
			else if(Vn==1){ return 0.020801+-0.0910493*pt+0.118184*pt*pt+-0.0365487*pt*pt*pt+0.0035618*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.0300511+-0.108966*pt+0.122315*pt*pt+-0.0365423*pt*pt*pt+0.00350489*pt*pt*pt*pt;}
		}
	}
	else if(cent==2){
		if(part==0){
			if(Vn==0){ return -0.0220529+0.172125*pt+-0.0353618*pt*pt+-0.003559*pt*pt*pt+0.00113968*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.0066372+0.0262161*pt+0.0372216*pt*pt+-0.0187145*pt*pt*pt+0.00228567*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.000152642+0.00135534*pt+0.0523496*pt*pt+-0.0225954*pt*pt*pt+0.0025451*pt*pt*pt*pt;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0424241+0.152629*pt+-0.00506494*pt*pt+-0.0151633*pt*pt*pt+0.00254353*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.0149325+0.0253627*pt+0.0329371*pt*pt+-0.0153877*pt*pt*pt+0.00170996*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.0171898+0.0261749*pt+0.032913*pt*pt+-0.0180592*pt*pt*pt+0.00240376*pt*pt*pt*pt;}
		}
		else if(part==2){
			if(Vn==0){ return 0.0128407+-0.0812974*pt+0.196424*pt*pt+-0.0729275*pt*pt*pt+0.0081403*pt*pt*pt*pt;}
			else if(Vn==1){ return 0.0216277+-0.0905268*pt+0.125852*pt*pt+-0.0410326*pt*pt*pt+0.00433817*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.0296393+-0.113592*pt+0.137947*pt*pt+-0.0424535*pt*pt*pt+0.00422479*pt*pt*pt*pt;}
		}
	}
	else if(cent==3){
		if(part==0){
			if(Vn==0){ return -0.0273469+0.215291*pt+-0.0580156*pt*pt+0.0015503*pt*pt*pt+0.00068957*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.00634738+0.0244379*pt+0.0472794*pt*pt+-0.0265474*pt*pt*pt+0.00383202*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.00529299+-0.0155944*pt+0.0851034*pt*pt+-0.0399046*pt*pt*pt+0.00537977*pt*pt*pt*pt;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0457415+0.184931*pt+-0.0184578*pt*pt+-0.0130774*pt*pt*pt+0.00241422*pt*pt*pt*pt;}
			else if(Vn==1){ return 0.00059818+-0.0174573*pt+0.0752039*pt*pt+-0.0318181*pt*pt*pt+0.00386052*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.00319935+-0.0357498*pt+0.0956003*pt*pt+-0.0389201*pt*pt*pt+0.0046787*pt*pt*pt*pt;}
		}
		else if(part==2){
			if(Vn==0){ return 0.00914554+-0.0597874*pt+0.203465*pt*pt+-0.0797661*pt*pt*pt+0.00929514*pt*pt*pt*pt;}
			else if(Vn==1){ return 0.0304227+-0.111558*pt+0.150866*pt*pt+-0.0511995*pt*pt*pt+0.00556649*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.0025491+-0.0227755*pt+0.0628781*pt*pt+-0.0165041*pt*pt*pt+0.00111185*pt*pt*pt*pt;}
		}
	}
	else if(cent==4){
		if(part==0){
			if(Vn==0){ return -0.0300557+0.23959*pt+-0.0712208*pt*pt+0.004233*pt*pt*pt+0.000504197*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.0047109+0.0195728*pt+0.0522525*pt*pt+-0.0282469*pt*pt*pt+0.00377098*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.0132215+0.0468*pt+0.0341852*pt*pt+-0.0206421*pt*pt*pt+0.00294137*pt*pt*pt*pt;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0425067+0.191418*pt+-0.0147714*pt*pt+-0.0177701*pt*pt*pt+0.00341417*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.000136675+-0.0175618*pt+0.0863983*pt*pt+-0.0430817*pt*pt*pt+0.00620464*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.0570229+0.17675*pt+-0.123802*pt*pt+0.0478088*pt*pt*pt+-0.00638515*pt*pt*pt*pt;}
		}
		else if(part==2){
			if(Vn==0){ return 0.0054852+-0.0327023*pt+0.19693*pt*pt+-0.0815048*pt*pt*pt+0.0098101*pt*pt*pt*pt;}
			else if(Vn==1){ return 0.0109575+-0.0600514*pt+0.115052*pt*pt+-0.0418587*pt*pt*pt+0.00470501*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.04566+-0.148386*pt+0.193706*pt*pt+-0.0675996*pt*pt*pt+0.00792379*pt*pt*pt*pt;}
		}
	}
	else if(cent==5){
		if(part==0){
			if(Vn==0){ return -0.0327924+0.250176*pt+-0.0765101*pt*pt+0.00390845*pt*pt*pt+0.000819225*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.003819+0.0112537*pt+0.0588299*pt*pt+-0.0333377*pt*pt*pt+0.0046983*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.00266408+0.00640717*pt+0.023585*pt*pt+-0.0121294*pt*pt*pt+0.00172664*pt*pt*pt*pt;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0131631+0.0325158*pt+-0.00707803*pt*pt+0.000728541*pt*pt*pt+-8.91768e-05*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.0571195+0.249406*pt+-0.074045*pt*pt+0.00463722*pt*pt*pt+0.000540562*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.0143528+0.0402737*pt+0.0106742*pt*pt+-0.00873702*pt*pt*pt+0.000978765*pt*pt*pt*pt;}
		}
		else if(part==2){
			if(Vn==0){ return -0.00897794+0.0202506*pt+0.159824*pt*pt+-0.0719297*pt*pt*pt+0.00894275*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.00854808+0.0237419*pt+0.00737678*pt*pt+0.00711372*pt*pt*pt+-0.00275382*pt*pt*pt*pt;}
			else if(Vn==2){ return 6.95282e-310+4.68399e-310*pt+4.684e-310*pt*pt+6.95282e-310*pt*pt*pt+5.05923e-321*pt*pt*pt*pt;}
		}
	}
}
Double_t HarmonicFunction(Double_t *pT, Int_t* par){
	Double_t pt;
	Int_t part = par[0];
	if(pT[0] > 4.0){  pt = 4.0;  }
	else if(pT[0] < 0.5){ return 0.0; }
	else{ pt = pT[0]; }

	Int_t cent	= par[1];
	

	Int_t Vn	=par[2];

	if(cent==0){
		if(part==0){
			if(Vn==0){ return -0.0143392+0.0759773*pt+-0.0154024*pt*pt+-0.000915513*pt*pt*pt+0.000264875*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.00558048+0.0136122*pt+0.0356842*pt*pt+-0.0164884*pt*pt*pt+0.00195171*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.00256332+-0.00483424*pt+0.0390473*pt*pt+-0.016438*pt*pt*pt+0.00205353*pt*pt*pt*pt;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0287711+0.0777838*pt+-0.0110908*pt*pt+-0.00251265*pt*pt*pt+0.000435634*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.00843333+0.00413273*pt+0.03708*pt*pt+-0.0140296*pt*pt*pt+0.00138853*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.0147847+0.0238129*pt+0.00149144*pt*pt+0.00122299*pt*pt*pt+-0.000583605*pt*pt*pt*pt;}
		}
		else if(part==2){
			if(Vn==0){ return 0.0139405+-0.0778337*pt+0.122313*pt*pt+-0.0441081*pt*pt*pt+0.00502149*pt*pt*pt*pt;}
			else if(Vn==1){ return 0.0204243+-0.0952528*pt+0.118049*pt*pt+-0.0375978*pt*pt*pt+0.00388916*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.0368641+-0.132059*pt+0.140577*pt*pt+-0.0465704*pt*pt*pt+0.00527664*pt*pt*pt*pt;}
		}
	}
	else if(cent==1){
		if(part==0){
			if(Vn==0){ return -0.0164604+0.119236*pt+-0.0140501*pt*pt+-0.00715798*pt*pt*pt+0.00137041*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.00615225+0.0205719*pt+0.0373663*pt*pt+-0.0176272*pt*pt*pt+0.00201875*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.00343526+-0.0156758*pt+0.058631*pt*pt+-0.0234438*pt*pt*pt+0.00258536*pt*pt*pt*pt;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0335949+0.108473*pt+0.00189956*pt*pt+-0.0119077*pt*pt*pt+0.00177362*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.0141866+0.0239085*pt+0.0233996*pt*pt+-0.0081926*pt*pt*pt+0.000430152*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.0178734+-0.0725294*pt+0.105302*pt*pt+-0.0391775*pt*pt*pt+0.00466704*pt*pt*pt*pt;}
		}
		else if(part==2){
			if(Vn==0){ return 0.0147481+-0.0885341*pt+0.16892*pt*pt+-0.0604128*pt*pt*pt+0.00661366*pt*pt*pt*pt;}
			else if(Vn==1){ return 0.020801+-0.0910493*pt+0.118184*pt*pt+-0.0365487*pt*pt*pt+0.0035618*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.0300511+-0.108966*pt+0.122315*pt*pt+-0.0365423*pt*pt*pt+0.00350489*pt*pt*pt*pt;}
		}
	}
	else if(cent==2){
		if(part==0){
			if(Vn==0){ return -0.0220529+0.172125*pt+-0.0353618*pt*pt+-0.003559*pt*pt*pt+0.00113968*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.0066372+0.0262161*pt+0.0372216*pt*pt+-0.0187145*pt*pt*pt+0.00228567*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.000152642+0.00135534*pt+0.0523496*pt*pt+-0.0225954*pt*pt*pt+0.0025451*pt*pt*pt*pt;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0424241+0.152629*pt+-0.00506494*pt*pt+-0.0151633*pt*pt*pt+0.00254353*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.0149325+0.0253627*pt+0.0329371*pt*pt+-0.0153877*pt*pt*pt+0.00170996*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.0171898+0.0261749*pt+0.032913*pt*pt+-0.0180592*pt*pt*pt+0.00240376*pt*pt*pt*pt;}
		}
		else if(part==2){
			if(Vn==0){ return 0.0128407+-0.0812974*pt+0.196424*pt*pt+-0.0729275*pt*pt*pt+0.0081403*pt*pt*pt*pt;}
			else if(Vn==1){ return 0.0216277+-0.0905268*pt+0.125852*pt*pt+-0.0410326*pt*pt*pt+0.00433817*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.0296393+-0.113592*pt+0.137947*pt*pt+-0.0424535*pt*pt*pt+0.00422479*pt*pt*pt*pt;}
		}
	}
	else if(cent==3){
		if(part==0){
			if(Vn==0){ return -0.0273469+0.215291*pt+-0.0580156*pt*pt+0.0015503*pt*pt*pt+0.00068957*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.00634738+0.0244379*pt+0.0472794*pt*pt+-0.0265474*pt*pt*pt+0.00383202*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.00529299+-0.0155944*pt+0.0851034*pt*pt+-0.0399046*pt*pt*pt+0.00537977*pt*pt*pt*pt;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0457415+0.184931*pt+-0.0184578*pt*pt+-0.0130774*pt*pt*pt+0.00241422*pt*pt*pt*pt;}
			else if(Vn==1){ return 0.00059818+-0.0174573*pt+0.0752039*pt*pt+-0.0318181*pt*pt*pt+0.00386052*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.00319935+-0.0357498*pt+0.0956003*pt*pt+-0.0389201*pt*pt*pt+0.0046787*pt*pt*pt*pt;}
		}
		else if(part==2){
			if(Vn==0){ return 0.00914554+-0.0597874*pt+0.203465*pt*pt+-0.0797661*pt*pt*pt+0.00929514*pt*pt*pt*pt;}
			else if(Vn==1){ return 0.0304227+-0.111558*pt+0.150866*pt*pt+-0.0511995*pt*pt*pt+0.00556649*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.0025491+-0.0227755*pt+0.0628781*pt*pt+-0.0165041*pt*pt*pt+0.00111185*pt*pt*pt*pt;}
		}
	}
	else if(cent==4){
		if(part==0){
			if(Vn==0){ return -0.0300557+0.23959*pt+-0.0712208*pt*pt+0.004233*pt*pt*pt+0.000504197*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.0047109+0.0195728*pt+0.0522525*pt*pt+-0.0282469*pt*pt*pt+0.00377098*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.0132215+0.0468*pt+0.0341852*pt*pt+-0.0206421*pt*pt*pt+0.00294137*pt*pt*pt*pt;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0425067+0.191418*pt+-0.0147714*pt*pt+-0.0177701*pt*pt*pt+0.00341417*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.000136675+-0.0175618*pt+0.0863983*pt*pt+-0.0430817*pt*pt*pt+0.00620464*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.0570229+0.17675*pt+-0.123802*pt*pt+0.0478088*pt*pt*pt+-0.00638515*pt*pt*pt*pt;}
		}
		else if(part==2){
			if(Vn==0){ return 0.0054852+-0.0327023*pt+0.19693*pt*pt+-0.0815048*pt*pt*pt+0.0098101*pt*pt*pt*pt;}
			else if(Vn==1){ return 0.0109575+-0.0600514*pt+0.115052*pt*pt+-0.0418587*pt*pt*pt+0.00470501*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.04566+-0.148386*pt+0.193706*pt*pt+-0.0675996*pt*pt*pt+0.00792379*pt*pt*pt*pt;}
		}
	}
	else if(cent==5){
		if(part==0){
			if(Vn==0){ return -0.0327924+0.250176*pt+-0.0765101*pt*pt+0.00390845*pt*pt*pt+0.000819225*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.003819+0.0112537*pt+0.0588299*pt*pt+-0.0333377*pt*pt*pt+0.0046983*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.00266408+0.00640717*pt+0.023585*pt*pt+-0.0121294*pt*pt*pt+0.00172664*pt*pt*pt*pt;}
		}
		else if(part==1){
			if(Vn==0){ return -0.0131631+0.0325158*pt+-0.00707803*pt*pt+0.000728541*pt*pt*pt+-8.91768e-05*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.0571195+0.249406*pt+-0.074045*pt*pt+0.00463722*pt*pt*pt+0.000540562*pt*pt*pt*pt;}
			else if(Vn==2){ return -0.0143528+0.0402737*pt+0.0106742*pt*pt+-0.00873702*pt*pt*pt+0.000978765*pt*pt*pt*pt;}
		}
		else if(part==2){
			if(Vn==0){ return -0.00897794+0.0202506*pt+0.159824*pt*pt+-0.0719297*pt*pt*pt+0.00894275*pt*pt*pt*pt;}
			else if(Vn==1){ return -0.00854808+0.0237419*pt+0.00737678*pt*pt+0.00711372*pt*pt*pt+-0.00275382*pt*pt*pt*pt;}
			else if(Vn==2){ return 0.0;}
		}
	}
	return 0.0;
}


 float MASS[4][6] ={{ 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  , { 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  , { 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  , { 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  };


 float BETA[4][6] ={{ 0.759913 , 0.760621 , 0.7319 , 0.732047 , 0.648344 , 0.680417 }  , { 0.765472 , 0.763144 , 0.739863 , 0.739893 , 0.641706 , 0.684773 }  , { 0.772956 , 0.773497 , 0.754914 , 0.756607 , 0.661272 , 0.692553 }  , { 0.913187 , 0.788666 , 0.780603 , 0.780076 , 0.681797 , 0.696818 }  };


 float TEMP[4][6] ={{ 0.175254 , 0.17647 , 0.163823 , 0.167908 , 0.25008 , 0.256333 }  , { 0.175317 , 0.17444 , 0.166239 , 0.170082 , 0.247217 , 0.256199 }  , { 0.173122 , 0.17437 , 0.164637 , 0.168712 , 0.240833 , 0.249546 }  , { 0.0835104 , 0.169124 , 0.16093 , 0.163796 , 0.218115 , 0.22822 }  };


 float NPAR[4][6] ={{ 32.371 , 33.1419 , 11.3702 , 12.1617 , 44.6847 , 74.0065 }  , { 30.7 , 28.7665 , 11.8083 , 12.3836 , 35.0342 , 65.2486 }  , { 28.408 , 29.0054 , 12.4268 , 13.1964 , 34.7224 , 54.8703 }  , { 27.5124 , 27.7284 , 14.806 , 14.6286 , 27.4966 , 36.6379 }  };


 float NORM[4][6] ={{ 10211.8 , 9733.9 , 7669.71 , 6958.68 , 714.188 , 833.707 }  , { 6923.66 , 7024.51 , 4679.75 , 4275.71 , 521.372 , 559.871 }  , { 4030.29 , 3828.84 , 2667.37 , 2413.47 , 341.159 , 360.574 }  , { 109206 , 1663.87 , 1117.2 , 1030.07 , 254.529 , 241.505 }  };
#endif

#ifndef TENNGEN_H
#define TENNGEN_H

#include "TGParticle.h"
#include "TGEvent.h"
#include "tenngenAuAu.h"
#include "TLegend.h"



TENNGEN_BEGIN_NAMESPACE

class TGSettings;
class TG200;

class TGSettings{
private:
    int energy;
    int nEvent, centbin;
    int dovn[5];
    float setpsi[5];
    float etarange;
    bool Histos;
    std::string outDirectory;
    bool dottree;
    bool QA;
    bool Batch;
    bool Stream;

public:

    ~TGSettings() {};
    TGSettings() {
        energy = 200;
        nEvent =1000;
        centbin = 0;
        etarange = 1.1;

        dovn[0] = 1;
        dovn[1] = 1;
        dovn[2] = 1;
        dovn[3] = 1;
        dovn[4] = 1;

        setpsi[0] = 0.;
        setpsi[1] = -1.;
        setpsi[2] = 0.;
        setpsi[3] = -1.;
        setpsi[4] = -1;

        Histos = true;
        dottree= true;
        Batch = true;
        QA = false;
        Stream = false;
        

        outDirectory = "./";
    }

    int getCollEn(){return energy;}
    int getCentBin(){ return centbin; }
    int getNevents(){ return nEvent; }
    float getEtaRange(){ return etarange;}
    int getVN(int n){ return dovn[n-1]; }
    float getPsiN(int n){ return setpsi[n-1]; }
    std::string getOutputDir(){return outDirectory;}
    bool getHistos(){return Histos;}
    bool getTTree(){return dottree; }
    bool qamode(){return QA;}
    bool batchmode(){return Batch;}
    bool streammode(){return Stream;}



    void setCollEn(int Energy){
        if(Energy == 200) energy=200;
        if(Energy == 5020) energy =5020;
    }
    void setCentBin(int cent){
        if(cent >=0 && cent<=5)   centbin = cent;
    }
    void setNevents(int numevents){nEvent = numevents;}
    void setEtaRange(float absEta){
        if(absEta >=0.0 && absEta<=1.1 )etarange = absEta;
    }
    void setVN(int n, bool On){
        if (On ) dovn[n-1] = 1;
        else if(!On) dovn[n-1] = 0; 
    }
    void setPsiN(int n, float userPsi){
        if((userPsi >=0 && userPsi <= 2.0*M_PI)||(userPsi == -1.0)){
            setpsi[n-1] = userPsi;
        }
    }
    void setOutputDir(std::string outDir){outDirectory  = outDir;}
    void doHistos(bool userInput){Histos = userInput;}
    void doTTree(bool userInput){dottree = userInput;}
    void setQA(bool userInput){
        QA = userInput;
        if(userInput){
        Histos = false;
        dottree = false;
        Batch = false; 
        Stream = false;
        }
    }
     void setBatch(bool userInput){
        Batch = userInput;
        if(userInput){
        QA = false; 
        Stream = false;
        }
     }
     void setStream(bool userInput){
        Stream = userInput;
        if(userInput){
        Histos = false;
        dottree = false;
        QA = false; 
        Batch = false;
        }
    }

};
class TG200{
private:
    const static int raw_high[4];
    const static int raw_low[4];
    const static int raw_bin[4];

    const static int correction_low[4];
    const static int correction_high[4];

    const static float piplusRatios[4];
    const static float piminusRatios[4];
    const static float kaplusRatios[4];
    const static float kaminusRatios[4];
    const static float proRatios[4];
    const static float pbarRatios[4];

    TGSettings settings;
    TGEventList tgEvents;

    TGPartList myparticles;
    TGEvent tmpEvent;

    int nEvent,setCent,doV1,doV2,doV3,doV4;
    float doPsi1,doPsi2,doPsi3,doPsi4;
    float etaRange;

    float ptPart, etaPart, phiPart;

    float v1, v2, v3, v4, maxPhi, dndphi;
    float psi1,psi2,psi3,psi4, tmpPhi;

    int multiplicity, spectraCent, harmonicCent, KF, CHECK;
    
    int partNumbers[6];

    TFile* myFile = new TFile(tenngen::AUAUFILEDIR);

    TH1F* ptDistroPip;
    TH1F* ptDistroPim;
    TH1F* ptDistroKap;
    TH1F* ptDistroKam;
    TH1F* ptDistroPro;
    TH1F* ptDistroPbar;
    TH1F* MultiDistro;

    string ptHistoPip; 
    string ptHistoPim; 
    string ptHistoKap; 
    string ptHistoKam; 
    string ptHistoPro; 
    string ptHistoPbar;
    string tFileHistos;
    string tFileTree;


    TRandom3* fRandom;

public:
    ~TG200() {};
    TG200() {     
        tgEvents.clear();
        clearEventBuffer();
        doPsi1 = 0;
        doPsi2 = 0;
        doPsi3 = 0;
        doPsi4 = 0;
    }
    TG200(TGSettings inputSettings, TRandom3* FRandom) : settings(inputSettings) , fRandom(FRandom) {
        

        tgEvents.clear();
        doPsi1 = 0;
        doPsi2 = 0;
        doPsi3 = 0;
        doPsi4 = 0;

        nEvent = settings.getNevents();
        etaRange = settings.getEtaRange();
        setCent = settings.getCentBin();
        doV1 = settings.getVN(1);
        doV2 = settings.getVN(2);
        doV3 = settings.getVN(3);
        doV4 = settings.getVN(4);

        psi1 = settings.getPsiN(1);
        if(psi1 ==-1.0 ) doPsi1 = -1.0;
        psi2 = settings.getPsiN(2);
        if(psi2 ==-1.0 ) doPsi2 = -1.0;
        psi3 = settings.getPsiN(3);
        if(psi3 ==-1.0 ) doPsi3 = -1.0;
        psi4 = settings.getPsiN(4);
        if(psi4 ==-1.0 ) doPsi4 = -1.0;
        
        getRootDistros();
        
        if(!settings.qamode()){genEvents();}
        if(settings.qamode()){genEventsQA();}

        tFileHistos=settings.getOutputDir()+"TG_"+std::to_string(nEvent)+"_eta"+std::to_string(etaRange)+"_"+"Cent"+std::to_string(setCent)+"_"+std::to_string(doV1)+std::to_string(doV2)+std::to_string(doV3)+std::to_string(doV4)+"_Histos.root";
        tFileTree=settings.getOutputDir()+"TG_"+std::to_string(nEvent)+"_eta"+std::to_string(etaRange)+"_"+"Cent"+std::to_string(setCent)+"_"+std::to_string(doV1)+std::to_string(doV2)+std::to_string(doV3)+std::to_string(doV4)+"_TTree.root";
        if(settings.getHistos()) tgEvents.writeHistos(tFileHistos);
        if(settings.getTTree()) tgEvents.writeTTree(tFileTree);

        
    

    }
   
    TGEventList events(){return tgEvents;}
    TGEvent& operator[](int i) {return tgEvents[i];}
    const TGEvent& operator[](int i) const {return tgEvents[i];}
    TGEvent& at(int i) {return tgEvents[i];}
    int size(){return tgEvents.size();}

    void config(TGSettings inputSettings);
    void seed(TRandom3* FRandom);
    void genEvents();
    void genEventsQA();
    void getRootDistros();
    void clearDistroBuffer();
    void clearEventBuffer();

    float dNdPhi();
    int SpectraCentBin();
    int HarmonicCentBin();
    float HarmonicFunction(const int harmonicN);
    

};

class TennGen{
private:

    friend class TGSettings;
    friend class TG200;
    TGSettings settings;
    TGEventList events;
    TRandom3* fRandom;
    TG200 tg;
    bool STREAMING;
    int StreamInt;

    string tFileHistos;
    string tFileTree;
    string dir;
  


public:

    ~TennGen() {};
    TennGen(const int collEn = 200, const int setEvents = 1000, const int setCent = 0, const float setEta =1.1, const int doV1=1, const int doV2=1,const int doV3=1,const int doV4=1,const int doV5 =1) {
        
        events.clear();
        settings.setCollEn(collEn);
        settings.setNevents(setEvents);
        settings.setCentBin(setCent);
        settings.setEtaRange(setEta);

        if(doV1 ==0) settings.setVN(1,false);
        if(doV2 ==0) settings.setVN(2,false);
        if(doV3 ==0) settings.setVN(3,false);
        if(doV4 ==0) settings.setVN(4,false);
        if(collEn == 5020){
            if(doV5 ==0) settings.setVN(5,false);
        }

        fRandom = new TRandom3(); 
        settings.setBatch(true);
       
    }
    

    void setcollen(int collen){ settings.setCollEn(collen);}
    void setnevent(int numevents){ settings.setNevents(numevents);}
    void setcent(int cent){settings.setCentBin(cent);}
    void seteta(float absEta){settings.setEtaRange(absEta);}
    
    
    void setvN(int n, bool On){ settings.setVN(n,On); }
    void setpsiN(int n, float userPsi){ settings.setPsiN(n,userPsi); }
    
    void do_Histos(bool userIn){settings.doHistos(userIn);}
    void do_TTree(bool userIn){settings.doTTree(userIn);}
    void set_QA(bool userIn){settings.setQA(userIn);}
    void set_Batch(bool userIn){settings.setBatch(userIn);}
    void set_Stream(bool userIn){settings.setStream(userIn);}
    
    void setOutputDir(std::string outDir){ settings.setOutputDir(outDir); }

    void defaultSettings200(bool userIn){

        if(userIn){
        settings.setCollEn(200);
        settings.setNevents(1000);
        settings.setCentBin(0);
        settings.setEtaRange(1.1);

        settings.setVN(1,true);
        settings.setVN(2,true);
        settings.setVN(3,true);
        settings.setVN(4,true);
        settings.setVN(5,false);

        settings.setPsiN(1,-1.0);
        settings.setPsiN(2,0.0);
        settings.setPsiN(3,-1.0);
        settings.setPsiN(4,0.0);


        settings.setOutputDir("./");

        settings.doHistos(true);
        settings.doTTree(true);
        settings.setQA(false);
        settings.setBatch(true);
        }

    }
    void defaultSettings5020(bool userIn){

        if(userIn){
        settings.setCollEn(5020);
        settings.setNevents(1000);
        settings.setCentBin(0);
        settings.setEtaRange(0.9);

        settings.setVN(1,true);
        settings.setVN(2,true);
        settings.setVN(3,true);
        settings.setVN(4,true);
        settings.setVN(5,true);

        settings.setPsiN(1,-1.0);
        settings.setPsiN(2,0.0);
        settings.setPsiN(3,-1.0);
        settings.setPsiN(4,0.0);
        settings.setPsiN(5,-1.0);

        settings.setOutputDir("./");

        settings.doHistos(true);
        settings.doTTree(true);
        settings.setQA(false);
        settings.setBatch(true);
        }

    }


    void init();
    void runBatch();
    void runQA();
    void runStream();

    TGEvent& next();


    TGEvent& operator[](int i) {return events[i];}
    const TGEvent& operator[](int i) const {return events[i];}
    TGEvent& at(int i) {return events[i];}
    int size(){return events.size();}

 
};






TENNGEN_END_NAMESPACE

#endif
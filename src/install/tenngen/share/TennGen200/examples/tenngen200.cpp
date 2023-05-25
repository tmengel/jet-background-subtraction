#include <iostream>
#include <string>

#include "TennGen.h"

using namespace tenngen;


const int MAXTRACKS = 2000;
float TrackPt[MAXTRACKS];
float TrackEta[MAXTRACKS];
float TrackPhi[MAXTRACKS];
int TrackKF[MAXTRACKS];

int main(){

    TennGen tg;
    TGEvent tmpEvent;


    tg.setcollen(200);
    tg.setnevent(10000);
    tg.setcent(0); // 0-10%
    tg.seteta(0.5); // abs value between 0 and 1.1

    tg.setvN(1,true);
    tg.setvN(2,true);
    tg.setvN(3,true);
    tg.setvN(4,true);
    
    tg.setpsiN(1,-1.0); // -1.0 means it will be randomly choosen between 0-2 pi
    tg.setpsiN(2,0.0);
    tg.setpsiN(3,-1.0);
    tg.setpsiN(4,0.0);
        
    tg.do_Histos(true); // will make debug root plots
    tg.do_TTree(true); // will make TTree 
    tg.set_QA(true);

    tg.init();

    
 

    return 0;

}
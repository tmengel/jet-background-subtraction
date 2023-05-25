#ifndef TENNGEN_BASE_H
#define TENNGEN_BASE_H

#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>

#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TRandom3.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "Riostream.h"
#include <TSystem.h>

#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>


#ifndef M_PI
#define M_PI 3.1415926535897932385
#endif

#ifndef M_C
#define M_C 2.99792458e8
#endif


/// \namespace tenngen
/// the TennGen namespace

// define this for easier readability (and obfuscation?) in
// a range of places
#define TENNGEN_BEGIN_NAMESPACE namespace tenngen {
#define TENNGEN_END_NAMESPACE   }

TENNGEN_BEGIN_NAMESPACE

using std::abs;
using std::string;
using std::to_string;
using std::vector;

using std::cout;
using std::istream;
using std::ostream;
using std::fstream;
using std::ofstream;

using std::endl;
using std::fixed;
using std::scientific;
using std::setprecision;
using std::setw;


TENNGEN_END_NAMESPACE

#endif
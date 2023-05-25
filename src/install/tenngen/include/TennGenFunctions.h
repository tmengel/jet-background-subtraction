#ifndef TENNGEN200AUAU_H
#define TENNGEN200AUAU_H

#include "tenngenbase.h"

TENNGEN_BEGIN_NAMESPACE


Double_t Multiplicity(Double_t* X, Double_t* par);

Double_t getdNdptOverptIntegrand(Double_t* rad, Double_t* par);

Double_t getdNdpt(Double_t* pT, Double_t* params);

std::string partcode(const int j);

void makeDistros200(std::string rootfilename);

void createAuAuHeader(std::string rootfileLocation);

TENNGEN_END_NAMESPACE




#endif
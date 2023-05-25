#ifndef TGPARTICLE_H
#define TGPARTICLE_H

#include "tenngenbase.h"

TENNGEN_BEGIN_NAMESPACE

class TGParticle {
private:

    double pt,eta,phi;
    int kf;
    int userindex;

public:
    
    ~TGParticle() {};

    TGParticle(const double PT =0, const double ETA=0, const double PHI=0, const int KF=0, const int u = -1) : pt(PT), eta(ETA), phi(PHI), kf(KF) , userindex(u){}
    TGParticle(const TGParticle& v) : pt(v.pt), eta(v.eta), phi(v.phi), kf(v.kf) , userindex(v.userindex){ }
    TGParticle& operator=(const TGParticle& v) { if (this != &v) { pt = v.pt; eta = v.eta;
    phi = v.phi; kf = v.kf, userindex = v.userindex; } return *this; }

 
    
    double Pt() { return pt; }
    double Eta(){ return eta; }
    double Phi(){ return phi; }
    int Kf(){ return kf; }

    void clear() {pt = 0.; eta = 0.; phi = 0.; kf = 0; setUserIndex(-1);}


    double M();
    double Px();
    double Py();
    double Pz();
    double E();
    double Rap();

    double getMass();

  
    void setUserIndex(const int u){ userindex = u;}
    int UserIndex(){return userindex;}

    friend ostream& operator<<(ostream&, const TGParticle& v);

 
    friend bool operator==(const TGParticle& a, const TGParticle& b);
    //inline bool operator!=(const TGParticle & a, const TGParticle & b) {return !(a==b);}

    friend bool operator==(const TGParticle & a, const int type);

};

ostream& operator<<(ostream&, const TGParticle& v) ;


bool operator==(const TGParticle& a, const TGParticle& b);
inline bool operator!=(const TGParticle & a, const TGParticle & b) {return !(a==b);}

bool operator==(const TGParticle& a, const int type);

//-------------------------------------------------------------------------------------------------

class TGPartList{
private:

    friend class TGParticle;

    std::vector<TGParticle> particle;

public:
    ~TGPartList() {};
    TGPartList() {}
    TGPartList& operator=(const TGPartList& List);
    TGPartList(const TGPartList& List) {*this = List;}


    void add(TGParticle newPart){particle.push_back(newPart);}
    int size() const {return particle.size();}
    void clear() { particle.clear();}

    void setListUserIndex(int u);

    TGParticle& operator[](int i) {return particle.at(i);}
    const TGParticle& operator[](int i) const {return particle.at(i);}

    TGParticle& at(int i) {return particle.at(i);}

    TGPartList& operator+=(const TGPartList& addList);
    TGPartList& operator+=(const TGParticle& addPart);

    friend TGPartList operator+(const TGPartList &a, const TGPartList &b);

};

  TGPartList operator+(const TGPartList &a, const TGPartList &b);

TENNGEN_END_NAMESPACE

#endif
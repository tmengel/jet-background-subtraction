#ifndef TGEVENT_H
#define TGEVENT_H

#include "TGParticle.h"

TENNGEN_BEGIN_NAMESPACE

class TGEventList;

class TGEvent {
private:
    friend class TGEventList;
    friend class TGParticle;
    friend class TGPartList;
    TGPartList particles;
    int centrality, multiplicity;
    float eta;

public:
    
    
    ~TGEvent() {};
    TGEvent() {}
    TGEvent(const TGPartList& Particles,const int Cent,const int Mult) : particles(Particles), centrality(Cent), multiplicity(Mult) {}
    TGEvent& operator=(const TGEvent& event);
    TGEvent(const TGEvent& event) {*this = event;}


    void add(TGParticle newPart){
    particles.add(newPart);
    multiplicity++;
    }
    void Eta(float in){eta = in;}
    void Centrality(int in){centrality =in;}
    void Multiplicity(int in){multiplicity = in;}

    int size() const {return particles.size();}
    void clear() { 
        particles.clear();
        centrality = -1;
        multiplicity = -1;
        eta = -1;
    }

    TGParticle& operator[](int i) {return particles[i];}
    const TGParticle& operator[](int i) const {return particles[i];}
    TGParticle& at(int i) {return particles[i];}

    int Centrality(){ return centrality; }
    int Multiplicity(){ return multiplicity; }
    float Eta() {return eta;}

    void setEventUserIndex(int in);


    TGEvent& operator+=(const TGEvent& addEvent);
    TGEvent& operator+=(const TGPartList& addList);
    TGEvent& operator+=(const TGParticle& addParticle);

    friend ostream& operator<<(ostream&, const TGEvent& v);
    friend TGEvent operator+(const TGEvent &a, const TGEvent &b);

    float avgPt();
    float recoPsiN(const int n);

    friend TGEventList operator+(const TGEventList &a, const TGEventList &b);
   


};

    ostream& operator<<(ostream&, const TGEvent& v) ;
    TGEvent operator+(const TGEvent &a, const TGEvent &b);

//-------------------------------------------------------------------------------------------------

class TGEventList{
private:

    friend class TGEvent;

    std::vector<TGEvent> events;

public:

    ~TGEventList() {};
    TGEventList() {}
    TGEventList& operator=(const TGEventList& eventList);
    TGEventList(const TGEventList& eventList) {*this = eventList;}


    void add(TGEvent newE){events.push_back(newE);}
    int size() const {return events.size();}
    void clear() { events.clear();}

    TGEvent& operator[](int i) {return events.at(i);}
    const TGEvent& operator[](int i) const {return events.at(i);}
    TGEvent& at(int i) {return events.at(i);}

    int Centrality(int i){ return events[i].Centrality(); }
    void Centrality(int i, int in){events[i].Centrality(in);}
    int Multiplicity(int i){ return events[i].Multiplicity(); }
    void Multiplicity(int i, int in){ events[i].Multiplicity(in);}
    float Eta(int i) {return  events[i].Eta();}
    void Eta(int i, float in){ events[i].Eta(in);}
    

    void writeTTree(std::string& outdir);
    void writeHistos(std::string& outdir);


    void setEventListUserIndex(int in);
    
    TGEventList& operator+=(const TGEventList& addList);
    TGEventList& operator+=(const TGEvent& addEvent);

    friend TGEventList operator+(const TGEventList &a, const TGEventList &b);
   

};

   TGEventList operator+(const TGEventList &a, const TGEventList &b);
   




TENNGEN_END_NAMESPACE

#endif
#include "TROOT.h"
#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH3D.h>
#include <TH1D.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TRandom.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <numeric>
//#include <array>        
//#include <random>       
//#include <chrono> 

#include <time.h>

#include <stdio.h>
#include <stdlib.h>

#include "Riostream.h"
#include "TString.h"
#include "TSystem.h"

#include "TMath.h"

using namespace std; 


// basic data structure to hold 3 masses
struct event_t {

    Double_t m1;
    Double_t m2; 
    Double_t m3;  

};


// distance-finding function
Double_t get_distance(event_t event1, event_t event2) {

   Double_t a = event1.m1 - event2.m1; 
   Double_t b = event1.m2 - event2.m2; 
   Double_t c = event1.m3 - event2.m3; 

   return sqrt(a*a + b*b + c*c);
}


Int_t ksbinned_tests(string infile1name, string infile2name, Int_t nbins, Int_t npermutations) {


// Opening data files and writing to a tree 
//-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0

std::ifstream infile1(infile1name.c_str(),std::ifstream::in);
std::ifstream infile2(infile2name.c_str(),std::ifstream::in);


int nlines1 = 0; 
int nlines2 = 0; 
std::string line; 
while(std::getline(infile1, line)) ++nlines1;
while(std::getline(infile2, line)) ++nlines2; 
const int nlinesc = nlines1 + nlines2; 

std::cout << "Read " << nlines1 << " events from file 1" << std::endl; 
std::cout << "Read " << nlines2 << " events from file 2" << std::endl; 

std::vector<event_t> events (nlinesc);

std::ifstream xinfile1(infile1name.c_str(),std::ifstream::in);
std::ifstream xinfile2(infile2name.c_str(),std::ifstream::in);


// file reading is fast so i'm not worried about having two loops.
// the files are read into the same array, one after the other, because 
// this way it is easier to randomly sample events when i permute it. 
int i=0;

Double_t a, b, c, d, e; 

while(xinfile1 >> a >> b >> c >> d >> e) {

    events[i].m1 = a;  
    events[i].m2 = b; 
    events[i].m3 = c;
    i++;
}

assert(i == nlines1); 
while(xinfile2 >> a >> b >> c >> d >> e) {
    events[i].m1 = a; 
    events[i].m2 = b; 
    events[i].m3 = c; 
    i++;
}


infile1.close(); 
infile2.close(); 
xinfile1.close(); 
xinfile2.close(); 


//-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0

vector<event_t>::const_iterator first = events.begin(); 
vector<event_t>::const_iterator middle = events.begin() + nlines1; 
vector<event_t>::const_iterator last = events.end(); 

vector<event_t> s1vec(first, middle); 
vector<event_t> s2vec(middle, last); 

assert(s1vec.size() == s2vec.size());


// do the nominal test here 

TH1D * s1hist = new TH1D("s1","",nbins,0,sqrt(3));
TH1D * s2hist = new TH1D("s2","",nbins,0,sqrt(3));

vector<Double_t> ksvector (npermutations+1);  

clock_t start = clock(); 

for(unsigned j=0; j<s1vec.size(); j++) {

    for(unsigned k=j+1; k<s1vec.size(); k++) {

        Double_t d1, d2;
        
        d1 = get_distance(s1vec[j], s1vec[k]); 
        d2 = get_distance(s2vec[j], s2vec[k]); 

        s1hist->Fill(d1);
        s2hist->Fill(d2); 
    }
}

clock_t end = clock(); 

std::cout << "that took " <<  (end - start)/CLOCKS_PER_SEC << std::endl; 

Double_t ksvaluenominal = s1hist->KolmogorovTest(s2hist, "M");
std::cout << "the nominal KS value is " << ksvaluenominal << std::endl; 

clock_t end2 = clock(); 

//std::cout << "the test then took " <<  (end2 - end)/CLOCKS_PER_SEC << std::endl;

std::cout << ksvaluenominal << std::endl; 

ksvector[0] = ksvaluenominal; 

//-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0

TH1D * kshist = new TH1D("ks","ks statistic distribution",80,0,0.1);

Int_t pcounter = 0;  

//main perumutation loop
for (Int_t p = 0; p<npermutations; p++) {

    if (p%20 == 0) {
        std::cout << "Reached permutation # " << p << std::endl; 
    }

    s1hist->Reset(); 
    s2hist->Reset(); 

    std::random_shuffle(events.begin(), events.end()); 

    //vector<event_t>::const_iterator first = events.begin(); 
    //vector<event_t>::const_iterator middle = events.begin() + nlines1; 
    //vector<event_t>::const_iterator last = events.end(); 
    
    //vector<event_t> s1vec(0, nlines1); 
    //vector<event_t> s2vec(nlines1, nlines1+nlines2); 

    vector<event_t> s1vec(first, middle); 
    vector<event_t> s2vec(middle, last);

    for(unsigned j=0; j<s1vec.size(); j++) {

        for(unsigned k=j+1; k<s1vec.size(); k++) {

            Double_t d1, d2;
        
            d1 = get_distance(s1vec[j], s1vec[k]); 
            d2 = get_distance(s2vec[j], s2vec[k]); 

            s1hist->Fill(d1);
            s2hist->Fill(d2); 
        }
    }

    Double_t ksvalue = s1hist->KolmogorovTest(s2hist, "M");
    
    ksvector[p+1] = ksvalue; 

    kshist->Fill(ksvalue);

    if (ksvalue > ksvaluenominal) {
        ++pcounter;
    }

    //std::cout << pcounter << std::endl; 

}
std::cout << pcounter << std::endl; 
//s1hist->Reset(); 
Double_t pvalue = (Double_t)pcounter/(Double_t)npermutations; 

std::cout << "The p-value is " << pvalue << std::endl; 

std::cout << ksvector.size() << std::endl; 
std::cout << ksvector[0] << std::endl; 

TCanvas * canvas = new TCanvas("canvas","canvas",1000,700); 
//canvas->Setlogy(); 


canvas->cd(); 

//gStyle->SetOptStat(""); 


kshist->GetYaxis()->SetTitle("");
kshist->GetXaxis()->SetTitle("");  

kshist->Draw("hist");


delete s1hist;
delete s2hist; 

return 69420;
}


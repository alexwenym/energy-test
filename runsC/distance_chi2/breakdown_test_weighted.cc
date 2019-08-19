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
#include <cmath>

#include <time.h>

#include <stdio.h>
#include <stdlib.h>

#include "Riostream.h"
#include "TString.h"
#include "TSystem.h"

#include "TMath.h"

using namespace std; 

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

// basic data structure to hold 3 masses. 
// Each event is defined as a collection of three masses. 
// If more masses are needed, you only need to modify this function, 
// the distance finding function, and the data-reading process. 
struct event_t {

    Double_t m1;
    Double_t m2; 
    Double_t m3;  

};

// distance-finding function - finds the distance between two events and uses the delta parameter. 
// gaussian metric 
Double_t get_distance(event_t event1, event_t event2, Double_t delta) {

   Double_t a = event1.m1 - event2.m1; 
   Double_t b = event1.m2 - event2.m2; 
   Double_t c = event1.m3 - event2.m3; 

   Double_t dij = sqrt(a*a + b*b + c*c);    
   
   //return sqrt(a*a + b*b + c*c);
   return exp(-(dij*dij)/(2*delta*delta));
}


void fill_distances (
    vector<event_t> s1vec, 
    vector<event_t> s2vec, 
    TH1D * shist, 
    TH1D * chist, 
    Double_t delta) {
/*
    for(unsigned j=0; j<s1vec.size(); j++) {
        for(unsigned k=0; k<s1vec.size(); k++) {
            if (k>j) { 
                Double_t d1, d2;
                d1 = get_distance(s1vec[j], s1vec[k], delta); 
                d2 = get_distance(s2vec[j], s2vec[k], delta); 
                shist->Fill(d1);
                shist->Fill(d2);
            }
            Double_t cross = get_distance(s1vec[j], s2vec[k], delta); 
            chist->Fill(cross); 
        }
    }
*/ 

    TH1D * shist1 = (TH1D*)shist->Clone(); 
    TH1D * shist2 = (TH1D*)shist->Clone(); 

    for(unsigned j=0; j<s1vec.size(); j++) {
        for(unsigned k=0; k<s1vec.size(); k++) {
            if (k>j) { 
                Double_t d1, d2;
                d1 = get_distance(s1vec[j], s1vec[k], delta); 
    
                shist1->Fill(d1);
    
            }
        }
    }

    Double_t si1 = (Double_t)s1vec.size(); 

    shist1->Scale(1.0/(si1*(si1-1)));

    for(unsigned j=0; j<s2vec.size(); j++) {
        for(unsigned k=0; k<s2vec.size(); k++) {
            if (k>j) { 
                Double_t d1, d2;
                d1 = get_distance(s2vec[j], s2vec[k], delta); 
    
                shist2->Fill(d1);
    
            }
        }
    }
        
    Double_t si2 = (Double_t)s2vec.size(); 

    shist2->Scale(1.0/(si2*(si2-1)));

    for(unsigned j=0; j<s1vec.size(); j++) {
        for(unsigned k=0; k<s2vec.size(); k++) {
            
            Double_t cross = get_distance(s1vec[j], s2vec[k], delta); 
            chist->Fill(cross); 

        }
    }
    
    *shist = (*shist1)+(*shist2); 

    chist->Scale(1.0/(si1*si2));


}



// main function
void breakdown_test_weighted(string infile1name, string infile2name, Int_t nbins, Int_t npermutations, Double_t delta) {
// two input files 
// number of bins in the binned KS test
// number of permutations to run when calculating p-value. Note, this needs to be high enough! 
// delta parameter for the gaussian distance metric

//-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0
// Opening data files 
//-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0

// define ifstream
std::ifstream infile1(infile1name.c_str(),std::ifstream::in);
std::ifstream infile2(infile2name.c_str(),std::ifstream::in);

// make histograms which will be filled with distances and which the ks-test will be applied to
TH1D * shist = new TH1D("pure","",nbins,0,1);
TH1D * chist = new TH1D("crossterms","",nbins,0,1);

// counting the number of lines in each file read in for verification and to get the number of events
Int_t nlines1 = 0; 
Int_t nlines2 = 0; 
std::string line; 
while(std::getline(infile1, line)) ++nlines1;
while(std::getline(infile2, line)) ++nlines2; 
const Int_t nlinesc = nlines1 + nlines2; 

std::cout << "Read " << nlines1 << " events from file 1" << std::endl; 
std::cout << "Read " << nlines2 << " events from file 2" << std::endl; 
std::cout << "The delta metric is " << delta << std::endl;
std::cout << "The number of bins is " << nbins << std::endl; 

// main vector holding all events from both samples. 
// the reason both samples are in one vector is that this is the most efficient way to shuffle the data
// to run permutations
std::vector<event_t> events (nlinesc);

// redefine ifstream after the last instance is used to count events
std::ifstream xinfile1(infile1name.c_str(),std::ifstream::in);
std::ifstream xinfile2(infile2name.c_str(),std::ifstream::in);

// two loops to read all the events in the file 
int i=0;
Double_t a, b, c, d, e; 
while(xinfile1 >> a >> b >> c >> d >> e) {
    events[i].m1 = a;  
    events[i].m2 = b; 
    events[i].m3 = c;
    i++;
}
// using one single iterator so that after s1 is filled to the event vector, s2 follows immediately
assert(i == nlines1); 
while(xinfile2 >> a >> b >> c >> d >> e) {
    events[i].m1 = a; 
    events[i].m2 = b; 
    events[i].m3 = c; 
    i++;
}

// clean up! 
infile1.close(); 
infile2.close(); 
xinfile1.close(); 
xinfile2.close(); 


// iterators to navigate to important locations in the main event vector
vector<event_t>::const_iterator first = events.begin(); 
vector<event_t>::const_iterator middle = events.begin() + nlines1; 
vector<event_t>::const_iterator last = events.end(); 


//-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0
// Finding nominal value
//-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0

// separate out samples from the main event vector
vector<event_t> s1vec(first, middle); 
vector<event_t> s2vec(middle, last); 
// currently, different size samples are not supported, but it should be easy to implement
//assert(s1vec.size() == s2vec.size());

// main loop that finds distances and fills the two histograms. 
fill_distances(s1vec, s2vec, shist, chist, delta);  

TH1F *disthist1 = (TH1F*)shist->Clone();
TH1F *disthist2 = (TH1F*)chist->Clone();

unsigned shist_xbins = shist->GetNbinsX(); 
unsigned chist_xbins = chist->GetNbinsX(); 
assert(shist_xbins == chist_xbins); 

Double_t increment = 1.0/(Double_t)shist_xbins; 
vector<Double_t> tvalues (shist_xbins); 

TH1F* tcontributions = new TH1F("tcontributions1","",shist_xbins,0,1);
TH1F* tcontributions_cu = new TH1F("tcontributions_cu", "", shist_xbins, 0,1);

for (unsigned i=1; i<=shist_xbins; i++) {
    Double_t psi1 = ((i-0.5)*increment)*(shist->GetBinContent(i));
    Double_t psi2 = ((i-0.5)*increment)*(chist->GetBinContent(i));

    Double_t t = psi1 - psi2; 

    tvalues[i-1] = t; 
    tcontributions->SetBinContent(i,t);
}

TH1F * tdistnominal = (TH1F*)tcontributions->Clone();
TH1F * tdistnominal_cu = (TH1F*)tcontributions->GetCumulative();

Double_t tareanominal = (tdistnominal_cu->Integral());
Double_t tvaluenominal = (tdistnominal_cu->GetBinContent(shist_xbins));
Double_t tprodnominal = tareanominal * tvaluenominal; 

cout << "The T test statistic is: " << tvaluenominal << endl; 
cout << "The nominal CDF area is: " << tareanominal << endl; 
cout << "The nominal area*T is: " << tprodnominal << endl; 


//-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0
// Starting permutations
//-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0

// make a histogram to contain the permutation values. Not essential, good for plotting 
//TH1D * thist = new TH1D("th","",1000,10e4,10e9);

TFile * histofile = new TFile("extra_histograms.root", "RECREATE");

//main perumutation loop
cout << "Starting " << npermutations << " permutations." << endl; 
Int_t pcounter = 0;  
for (Int_t p = 0; p<npermutations; p++) {

    if (p%1000 == 0) {
        std::cout << "Reached permutation # " << p << std::endl; 
    }
    
    shist->Reset(); 
    chist->Reset(); 
    
    std::random_shuffle(events.begin(), events.end()); 

    vector<event_t> s1vec(first, middle); 
    vector<event_t> s2vec(middle, last);
   
    fill_distances(s1vec, s2vec, shist, chist, delta); 
 
    unsigned shist_xbins = shist->GetNbinsX(); 
    unsigned chist_xbins = chist->GetNbinsX(); 
    assert(shist_xbins == chist_xbins); 

    Double_t increment = 1.0/(Double_t)shist_xbins; 
    vector<Double_t> tvalues (shist_xbins); 
    
    tcontributions->Reset(); 
    tcontributions_cu->Reset(); 

    for (unsigned i=1; i<=shist_xbins; i++) {
        Double_t psi1 = ((i-0.5)*increment)*(shist->GetBinContent(i));
        Double_t psi2 = ((i-0.5)*increment)*(chist->GetBinContent(i));
        Double_t t = psi1 - psi2; 

        tcontributions->SetBinContent(i,t); 
    }

    tcontributions_cu = (TH1F*)tcontributions->GetCumulative();

    Double_t tarea = tcontributions_cu->Integral();
    Double_t tvalue = tcontributions_cu->GetBinContent(shist_xbins);
    Double_t tprod = tarea*tvalue; 
 
    if (tvalue > tvaluenominal) {
        ++pcounter;
        TString ss = TString::Format("histo_%i", p);
        tcontributions_cu->Write(ss);
        cout << tvalue << endl; 
    }
}

// p value calculated here
Double_t pvalue = (Double_t)pcounter/(Double_t)npermutations; 

std::cout << "The p-value is " << pvalue << std::endl; 

std::cout << "Finished!" << std::endl; 


/*
Double_t s = 3.0;

TCanvas * c1 = new TCanvas("c1", "", 1000,700);
c1->Divide(2,2); 
c1->cd(1);
tdistnominal->SetTitle("nominal");
tdistnominal->GetYaxis()->SetRangeUser(-s*5e-6, s*5e-6); 
tdistnominal->Draw("hist");

//TH1F * tdistnominal_cu = (TH1F*)tdistnominal->GetCumulative();
c1->cd(3); 
tdistnominal_cu->SetLineColor(2);
tdistnominal_cu->GetYaxis()->SetRangeUser(-s*2e-5, s*8e-5); 
tdistnominal_cu->Draw("hist"); 
Double_t area1 = tdistnominal_cu->Integral(); 
//cout << area1 << endl; 
 

c1->cd(2); 
tcontributions->SetTitle("permutation");
tcontributions->GetYaxis()->SetRangeUser(-s*5e-6, s*5e-6); 
tcontributions->Draw("hist");

//TH1F * tcontributions_cu = (TH1F*)tcontributions->GetCumulative();
c1->cd(4); 
tcontributions_cu->SetLineColor(2);
tcontributions_cu->GetYaxis()->SetRangeUser(-s*2e-5, s*8e-5);
tcontributions_cu->Draw("hist"); 
Double_t area2 = tcontributions_cu->Integral(); 
//cout << area2 << endl; 

cout << "Finished again!" << endl; 

*/

// cleanup! (not really necessary if only run once)
histofile->Write(); 
//histofile->Close(); 
delete shist;
delete chist;
//delete tcontributions; 
//delete kshist; 


}


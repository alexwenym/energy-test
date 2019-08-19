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
//#include <bits/stdc++.h> 

#include <time.h>

#include <stdio.h>
#include <stdlib.h>

#include "Riostream.h"
#include "TString.h"
#include "TSystem.h"

#include "TMath.h"

using namespace std; 

// This script can take in two samples, make distributions of the distances, and run a binned KS-test on the
// distributions. This is a more fundamental way to probe the processes that underlie the energy test. 
// - Alex Wen, July 2019


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


// main function
Int_t breakdown_test(string infile1name, string infile2name, Int_t nbins, Int_t npermutations, Double_t delta) {
// two input files 
// number of bins in the binned KS test
// number of permutations to run when calculating p-value. Note, this needs to be high enough! 
// delta parameter for the gaussian distance metric


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
assert(s1vec.size() == s2vec.size());

// main loop that finds distances and fills the two histograms. 
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

// Run the KS test to get the nominal KS value. 
// if the M option is not specificed, the test returns an (inaccurate) p value instead of the test statistic

shist->Scale(1.0/3998000.0); 
chist->Scale(1.0/4000000.0); 

TH1F *disthist1 = (TH1F*)shist->Clone();
TH1F *disthist2 = (TH1F*)chist->Clone();

unsigned shist_xbins = shist->GetNbinsX(); 
unsigned chist_xbins = chist->GetNbinsX(); 
assert(shist_xbins == chist_xbins); 

Double_t increment = 1.0/(Double_t)shist_xbins; 
vector<Double_t> tvalues (shist_xbins); 

TH1F* tcontributions = new TH1F("tcontributions1","",shist_xbins,0,1);

for (unsigned i=1; i<=shist_xbins; i++) {
    Double_t psi1 = ((i-0.5)*increment)*(shist->GetBinContent(i));
    Double_t psi2 = ((i-0.5)*increment)*(chist->GetBinContent(i));
    Double_t sum = (shist->GetBinContent(i))+(chist->GetBinContent(i));
    Double_t t = psi1 - psi2; 

    tvalues[i-1] = t; 
    tcontributions->SetBinContent(i,t); 
}



Double_t tsum = std::accumulate(tvalues.begin(),tvalues.end(), 0.0);
std::cout << "The total T is: " << tsum << std::endl; 
std::cout << "The max T is: " << *max_element(tvalues.begin(), tvalues.end()) << std::endl; 
std::cout << "The min T is: " << *min_element(tvalues.begin(), tvalues.end()) << std::endl; 


/*
disthist1->Add(disthist2, -1); 
//*disthist1 = (*disthist1)*(*disthist1);
//assert(disthist1->GetEntries() > 1);

Double_t tvaluenominal = disthist1->Integral(); 

cout << "The nominal test statistic is " << tvaluenominal << endl; 
*/

//TCanvas * tcanvas = new TCanvas("tcanvas", "tcanvas", 600,400);
//tcanvas->cd(); 
//tcontributions->Draw("hist"); 

//for (unsigned i = 0; i < shist_xbins; ++i) {
//    std::cout << tvalues[i] << ' ';
//}

// create a vector that will contain all the ks values obtained from the permutation test. 
// The first entry is the nominal one, so there are npermutations+1 entries in the vector. 
// This is mainly to facilitate saving/writing the KS distribution to a file should it be needed. 
//vector<Double_t> tvector (npermutations+1);
//ksvector[0] = ksvaluenominal; 


//-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0
// Starting permutations
//-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0-----0



// make a histogram to contain the permutation values. Not essential, good for plotting 
//TH1D * thist = new TH1D("th","",1000,10e4,10e9);

//main perumutation loop
Int_t pcounter = 0;  
for (Int_t p = 0; p<npermutations; p++) {

    if (p%2000 == 0) {
        std::cout << "Reached permutation # " << p << std::endl; 
    }
    
    // histograms are reused, so they are cleared of the last loop's contents and re-filled
    shist->Reset(); 
    chist->Reset(); 
    
    // shuffling the events for permutation. Even C++ documentation doesn't say what the 
    // exact randomness mechanism is, but this works....
    std::random_shuffle(events.begin(), events.end()); 

    vector<event_t> s1vec(first, middle); 
    vector<event_t> s2vec(middle, last);
    
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

    shist->Scale(1.0/3998000.0); 
    chist->Scale(1.0/4000000.0); 
    
    unsigned shist_xbins = shist->GetNbinsX(); 
    unsigned chist_xbins = chist->GetNbinsX(); 
    assert(shist_xbins == chist_xbins); 

    Double_t increment = 1.0/(Double_t)shist_xbins; 
    vector<Double_t> tvalues (shist_xbins); 
    
    tcontributions->Reset(); 
    

    for (unsigned i=1; i<=shist_xbins; i++) {
        Double_t psi1 = ((i-0.5)*increment)*(shist->GetBinContent(i));
        Double_t psi2 = ((i-0.5)*increment)*(chist->GetBinContent(i));
        Double_t sum = (shist->GetBinContent(i))+(chist->GetBinContent(i));
        Double_t t = psi1 - psi2; 

        tvalues[i-1] = t; 
        tcontributions->SetBinContent(i,t); 
    }


    Double_t tvalue = std::accumulate(tvalues.begin(),tvalues.end(), 0.0);
    
    //ksvector[p+1] = ksvalue; 

    //thist->Fill(tvalue);
    
    // rudimentary way to count how many ks values are greater than the nominal one, to calculate the p-value
    if (tvalue > tsum) {
         ++pcounter;
    }
}


// p value calculated here
Double_t pvalue = (Double_t)pcounter/(Double_t)npermutations; 

std::cout << "The p-value is " << pvalue << std::endl; 

std::cout << "Finished!" << std::endl; 



/*
// plotting the ks distribution
TCanvas * canvas = new TCanvas("canvas","canvas",1000,700);  
canvas->cd(); 
gPad->SetLogx(); 
//gStyle->SetOptStat(""); 
thist->GetYaxis()->SetTitle("");
thist->GetXaxis()->SetTitle("");  
thist->Draw("hist");


// some plotting code, mainly for debugging. Shows the distance distribution
TCanvas * canvas0 = new TCanvas("canvas0n","canvas0",1000,700); 
canvas0->Divide(2,1); 
canvas0->cd(1); 
gPad->SetLogy(); 
disthist1->Draw("hist"); 
canvas0->cd(2); 
gPad->SetLogy(); 
disthist2->Draw("hist");
*/

//TCanvas * c1 = new TCanvas("c1", "", 1000,700); 
//c1->cd(); 
//tcontributions->Draw();

/*
std::cout << disthist1->GetNbinsX() << std::endl; 
std::cout << disthist2->GetNbinsX() << std::endl; 

TCanvas * diffcanvas = new TCanvas("diffcanvas","",700,500); 
diffcanvas->cd();
disthist1->Add(disthist2, -1);
//TH1F disthistsq = (*disthist1); 
//TH1F diffcuhist;
//diffcuhist.Add(disthist1->GetCumulative(), disthist2->GetCumulative(), 1, -1);
*disthist1 = (*disthist1)*(*disthist1); 
disthist1->Draw("hist");

std::cout << disthist1->GetEntries() << std::endl; 
*/

// cleanup! 
delete shist;
delete chist;
delete tcontributions; 
//delete kshist; 

 
return 1729;
}


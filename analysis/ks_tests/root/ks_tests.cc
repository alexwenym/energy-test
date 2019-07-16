#include "TROOT.h"
#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH3.h>
#include <TH3F.h>
#include <TH3D.h>
#include <TH1D.h>
#include <TH2D.h>
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

#include <stdio.h>
#include <stdlib.h>

#include "Riostream.h"
#include "TString.h"
#include "TSystem.h"

#include "TMath.h"



// helper functions here


int ks_tests() {

ifstream infile1("../../../data_generation/s1_low.txt",std::ifstream::in);
ifstream infile2("../../../data_generation/s2_low.txt",std::ifstream::in);

TH1D * histo = new TH1D("Energy","",105,0,11.5); 
TCanvas * canvas = new TCanvas("canvas","canvas",1000,700); 

double m11, m12, m13, m21, m22, m23

while(!infile1.eof()) {

    infile1 >> 
    m11 >> 
    m12 >>
    m13 ; 

    infile2 >> 
    m21 >>
    m22 >> 
    m23 ; 

    //histo->Fill(energy,freq); 

}

canvas->cd(); 

gStyle->SetOptStat(""); 

histo->GetYaxis()->SetTitle("Relative Counts");
histo->GetXaxis()->SetTitle("Neutron Energy (MeV)");  

histo->Draw("hist C");

return 0;
}

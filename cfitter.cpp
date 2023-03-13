#include <iostream>
#include <fstream>
#include <string>

void extractor(TH2 *hist, const char* name){

    TCanvas *c1 = new TCanvas("c1","c1",1200,1000);
    c1->cd();
    c1->Print((std::string(name)+ std::string("_plots.pdf[")).c_str());
    TH1D *iter;
    TFitResultPtr res;
    double p1[hist->GetNbinsX()];
    double error[hist->GetNbinsX()];
    TF1 *f1;
    for(int i = 0; i < hist->GetNbinsX(); i++){
        iter = hist->ProjectionY("",i,i);
        iter->Scale(1/iter->Integral());
        f1 = new TF1("f1", "[0] + [1]*cos(x)");
        iter->Fit(f1,"","",-3.14,3.14);
        p1[i] = f1->GetParameter(1);
        error[i] = f1->GetParError(1);
        iter->Draw();
        c1->Print((std::string(name)+ std::string("_plots.pdf")).c_str());
    }

    ofstream paramfile;
    paramfile.open((std::string(name)+std::string("_output.txt")).c_str());
    for(int i =0; i < hist->GetNbinsX(); i++){
        paramfile << p1[i] <<endl;
    }
    paramfile.close();
    ofstream errfile;
    errfile.open((std::string(name)+std::string("_err.txt")).c_str());
    for(int i =0; i < hist->GetNbinsX(); i++){
        errfile << error[i] <<endl;
    }
    errfile.close();
    c1->Print((std::string(name)+std::string("_plots.pdf]")).c_str());
}

void cfitter(){
    TFile *f = new TFile("RootFiles/cp_phase_pi_half.root");
    TH2F *bellhist = (TH2F*) f->Get("cp_pi_half_bellInequality;1");
    TH2F *mockbellhist = (TH2F*) f->Get("cp_pi_half_mockBellInequality;1");
    

    bellhist->SetDirectory(0);
    mockbellhist->SetDirectory(0);
    extractor(bellhist,"bell");
    extractor(mockbellhist,"mock");
    f->Close();
}

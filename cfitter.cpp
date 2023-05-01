#include <iostream>
#include <fstream>
#include <string>

void extractor(TH2 *hist, const char* name){

    TCanvas *c1 = new TCanvas("c1","c1",1200,1000);
    c1->cd();
    c1->Print((std::string("plots/")+ std::string(name)+ std::string("_plots.pdf[")).c_str());
    TH1D *iter;
    TFitResultPtr res;
    double p1[hist->GetNbinsX()];
    double p0[hist->GetNbinsX()];
    double error[hist->GetNbinsX()];
    TF1 *f1;
    for(int i = 0; i < hist->GetNbinsX(); i++){
        iter = hist->ProjectionY("",i,i);
        iter->Scale(1/iter->Integral());
        f1 = new TF1("f1", "[0] + [1]*cos(x)");
        iter->Fit(f1,"","",-3.14,3.14);
        p0[i] = f1->GetParameter(0);
        p1[i] = f1->GetParameter(1);
        error[i] = f1->GetParError(1);
        iter->Draw();
        c1->Print((std::string("plots/")+std::string(name)+ std::string("_plots.pdf")).c_str());
    }

    ofstream paramfile;
    paramfile.open((std::string("data/")+std::string(name)+std::string("_output.txt")).c_str());
    for(int i =0; i < hist->GetNbinsX(); i++){
        paramfile << p1[i] <<endl;
    }
    paramfile.close();
    ofstream param2file;
    param2file.open((std::string("data/")+std::string(name)+std::string("_offset.txt")).c_str());
    for(int i =0; i < hist->GetNbinsX(); i++){
        param2file << p0[i] <<endl;
    }
    param2file.close();
    ofstream errfile;
    errfile.open((std::string("data/")+std::string(name)+std::string("_err.txt")).c_str());
    for(int i =0; i < hist->GetNbinsX(); i++){
        errfile << error[i] <<endl;
    }
    errfile.close();
    c1->Print((std::string("plots/")+std::string(name)+std::string("_plots.pdf]")).c_str());
}

void cfitter(){
    TFile *f = new TFile("RootFiles/cp_phase_0.root");
    TH2F *bellhist = (TH2F*) f->Get("cp_0_speedBellInequality;1");
    TH2F *mockbellhist = (TH2F*) f->Get("cp_0_speedMockInequality;1");
    TH2F *mock2bellhist = (TH2F*) f->Get("cp_0_speedMock2Inequality;1");
    TH2F *smbellhist = (TH2F*) f->Get("cp_0_smspeedBellInequality;1");
    TH2F *smmockbellhist = (TH2F*) f->Get("cp_0_smspeedMockInequality;1");
    TH2F *smmock2bellhist = (TH2F*) f->Get("cp_0_smspeedMock2Inequality;1");
    

    bellhist->SetDirectory(0);
    mockbellhist->SetDirectory(0);
    mock2bellhist->SetDirectory(0);
    smbellhist->SetDirectory(0);
    smmockbellhist->SetDirectory(0);
    smmock2bellhist->SetDirectory(0);
    extractor(bellhist,"bell");
    extractor(mockbellhist,"mock");
    extractor(mock2bellhist,"mock2");
    extractor(smbellhist,"smbell");
    extractor(smmockbellhist,"smmock");
    extractor(smmock2bellhist,"smmock2");
    f->Close();
}

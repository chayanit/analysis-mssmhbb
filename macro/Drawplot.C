#include <memory>
#include <fstream>
#include <ostream>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

void drawcanvas(TH1D* h_novo, TH1D* h_syst)
{
  	//gStyle->SetOptTitle(0);
  	gStyle->SetOptStat("");
  	gStyle->SetLegendBorderSize(0);

	TH1D* h_central = (TH1D*)h_novo->Clone("h_central");
	//h_central->Draw("hist");
	TH1D* h_ratio = (TH1D*)h_novo->Clone("h_ratio");

	bool logy = 1;
	TCanvas canvas;
	canvas.SetCanvasSize(500,500);

	TPad* pad1;
   	pad1 = new TPad("pad1","",0,0.1,1,1);
    	pad1->SetBottomMargin(0.2);
    	pad1->SetRightMargin(0.05); // The ratio plot below inherits the right and left margins settings here!
    	pad1->SetLeftMargin(0.16); 
  	pad1->Draw();
  	pad1->cd();
  	if (logy) { 	pad1->SetLogy();
			h_central->SetMaximum(h_central->GetMaximum()*10);
			h_central->SetMinimum(h_central->GetMinimum()/10);
			h_central->SetMinimum(0.1);
	}
 	else {		h_central->SetMaximum(h_central->GetMaximum()+200.);
			h_central->SetMinimum(0.);
	}
	//h_central->SetTitle(h_syst->GetName());
      	h_central->GetXaxis()->SetTitleOffset(999); //Effectively turn off x axis title on main plot
      	h_central->GetXaxis()->SetLabelOffset(999); //Effectively turn off x axis label on main plot
      	h_central->GetYaxis()->SetTitleSize(0.041);
      	h_central->GetYaxis()->SetTitleOffset(1.20);
      	h_central->GetYaxis()->SetLabelSize(0.04);
 	h_central->Draw("hist");
	h_syst->SetLineColor(kBlack);
	h_syst->Draw("same hist");	

        TLegend *leg1 = new TLegend(0.60,0.60,0.86,0.87);
        leg1->SetFillColor(kWhite);
        leg1->SetLineColor(kWhite);
        leg1->AddEntry("h_central","Rereco C,D,E","L");
        leg1->AddEntry("h_syst","Rereco F","L");
        leg1->Draw("same");

	canvas.cd();
	//bottom ratio
	TPad *pad2 = new TPad("pad2","",0,0.0,1,0.25);
	pad2->SetTopMargin(1);
	pad2->SetBottomMargin(0.33);
	pad2->SetLeftMargin(pad1->GetLeftMargin());
	pad2->SetRightMargin(pad1->GetRightMargin());
     	pad2->SetGridy();
	pad2->Draw();
	pad2->cd();

	h_ratio->SetTitle("");
	h_ratio->Divide(h_syst);
    	h_ratio->SetMarkerStyle(8);
  	h_ratio->SetMarkerSize(0.5);	
	//h_ratio->SetMarkerColor(46);
	h_ratio->SetLineColor(kBlack);
    	h_ratio->GetXaxis()->SetTitleSize(0.15);
    	h_ratio->GetXaxis()->SetTitleOffset(0.85);
    	h_ratio->GetXaxis()->SetLabelSize(0.12);
    	h_ratio->GetXaxis()->SetLabelOffset(0.008);
    	//h_ratio->SetYTitle("Ratio");
	h_ratio->SetYTitle("#frac{CDE}{F}");
    	h_ratio->GetYaxis()->CenterTitle(kTRUE);
    	h_ratio->GetYaxis()->SetTitleSize(0.15);
    	h_ratio->GetYaxis()->SetTitleOffset(0.3);
    	h_ratio->GetYaxis()->SetNdivisions(3,5,0);
    	h_ratio->GetYaxis()->SetLabelSize(0.12);
    	h_ratio->GetYaxis()->SetLabelOffset(0.011);

    	h_ratio->Draw("p");
    	h_ratio->SetMaximum(5.0);
    	h_ratio->SetMinimum(-5.0);

	//TString pdfname;
	//pdfname = h_syst->GetName();
	//canvas.SaveAs(pdfname+"_withratio.pdf");
	canvas.SaveAs("M12_CDEF_CRwithratio.pdf");

}

void Drawplot() 
{
	TH1::SetDefaultSumw2();
 
	gStyle->SetOptStat("");
        gStyle->SetLegendBorderSize(0);

  	TFile *f1 = new TFile("/nfs/dust/cms/user/chayanit/MSSMHBB_Myfork2017/CMSSW_9_4_4/src/Analysis/Mssmhbb/test/examples/allhadronic/histogram_C_CR.root") ;
 	TFile *f2 = new TFile("/nfs/dust/cms/user/chayanit/MSSMHBB_Myfork2017/CMSSW_9_4_4/src/Analysis/Mssmhbb/test/examples/allhadronic/histogram_D_CR.root") ;
	TFile *f3 = new TFile("/nfs/dust/cms/user/chayanit/MSSMHBB_Myfork2017/CMSSW_9_4_4/src/Analysis/Mssmhbb/test/examples/allhadronic/histogram_E_CR.root") ;
	TFile *f4 = new TFile("/nfs/dust/cms/user/chayanit/MSSMHBB_Myfork2017/CMSSW_9_4_4/src/Analysis/Mssmhbb/test/examples/allhadronic/histogram_F_CR.root") ;
 
	TTree *t1 = (TTree*)f1->Get("MssmHbb_13TeV");
	TTree *t2 = (TTree*)f2->Get("MssmHbb_13TeV");
	TTree *t3 = (TTree*)f3->Get("MssmHbb_13TeV");
	TTree *t4 = (TTree*)f4->Get("MssmHbb_13TeV");

	TH1D *h1 = new TH1D("h1", ";mbb [GeV];Entries", 50, 0., 2500.);
	TH1D *h2 = new TH1D("h2", ";mbb [GeV];Entries", 50, 0., 2500.);
	TH1D *h3 = new TH1D("h3", ";mbb [GeV];Entries", 50, 0., 2500.);
	TH1D *h4 = new TH1D("h4", ";mbb [GeV];Entries", 50, 0., 2500.);
 	
	double mbb1,mbb2,mbb3,mbb4;
	t1->SetBranchAddress("mbb", &mbb1);
	t2->SetBranchAddress("mbb", &mbb2);
	t3->SetBranchAddress("mbb", &mbb3);
	t4->SetBranchAddress("mbb", &mbb4);
 
	for (int i = 0; i < t1->GetEntriesFast(); ++i) 
	{	t1->GetEntry(i);
		if(mbb1 > 0.) h1->Fill(mbb1);
	}
	//h1->Scale(1/h1->Integral());
	h1->SetLineWidth(2);
	h1->SetLineColor(kBlack);

        for (int i = 0; i < t2->GetEntriesFast(); ++i)
        {       t2->GetEntry(i);
                if(mbb2 > 0.) h2->Fill(mbb2);
        }
	//h2->Scale(1/h2->Integral());
	h2->SetLineWidth(2);
	h2->SetLineColor(kRed);
	
        for (int i = 0; i < t3->GetEntriesFast(); ++i)
        {       t3->GetEntry(i);
                if(mbb3 > 0.) h3->Fill(mbb3);
        }
	//h3->Scale(1/h3->Integral());
	h3->SetLineWidth(2);
        h3->SetLineColor(kGreen+2);

	TH1D *h_CDE = (TH1D*)h1->Clone("h_CDE");
	h_CDE->Add(h2);
	h_CDE->Add(h3);
	h_CDE->SetLineWidth(2);
        h_CDE->SetLineColor(46);
	//h_CDE->Draw("hist");	
	//std::cout << "h_CDE integral = " << h_CDE->Integral() << std::endl;

        for (int i = 0; i < t4->GetEntriesFast(); ++i)
        {       t4->GetEntry(i);
                if(mbb4 > 0.) h4->Fill(mbb4);
        }	
	//std::cout << "h_F integral = " << h4->Integral() << std::endl;
        h4->SetLineWidth(2);
        h4->SetLineColor(9);
 
	drawcanvas(h_CDE,h4);

        h1->Scale(1/h1->Integral());
        h2->Scale(1/h2->Integral());
        h3->Scale(1/h3->Integral());
        h_CDE->Scale(1/h_CDE->Integral());
	h4->Scale(1/h4->Integral());	
	//drawcanvas(h_CDE,h4);

  	TCanvas *canvas = new TCanvas("canvas","canvas",1000,600);
  	//canvas->Divide(3,2);
  	//canvas->cd(1);
	canvas->SetLogy();
  	h1->Draw("hist");
	h2->Draw("hist same");
	h3->Draw("hist same");
	h4->Draw("hist same");
	
	TLegend *leg1 = new TLegend(0.60,0.60,0.86,0.87);
        leg1->SetFillColor(kWhite);
        leg1->SetLineColor(kWhite);
	leg1->AddEntry("h1","Rereco C","L");
	leg1->AddEntry("h2","Rereco D","L");
	leg1->AddEntry("h3","Rereco E","L");
	leg1->AddEntry("h4","Rereco F","L");
	leg1->Draw("same");	

  	//drawcanvas(h_CDE,h4);
  	canvas->SaveAs("M12_CDEF_CR.pdf");
  	cout << "End of code :)" << endl;

}



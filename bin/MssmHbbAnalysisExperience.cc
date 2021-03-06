#include "boost/program_options.hpp"
#include "boost/algorithm/string.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h" 

#include "Analysis/MssmHbb/interface/MssmHbbAnalyser.h"

using namespace std;
using namespace analysis;
using namespace analysis::tools;
using namespace analysis::mssmhbb;

using TH1s = std::map<std::string, TH1F*>;
using TH2s = std::map<std::string, TH2F*>;
     

// float GetBTag(const Jet & jet, const std::string & algo);
// void CreateHistograms(TH1s &, const int & n = 3);
// void CreateHistograms(TH2s &, const int & n = 3);
// template <typename T>
// void WriteHistograms(T & );


// =============================================================================================   
int main(int argc, char ** argv)
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   MssmHbbAnalyser mssmhbb(argc,argv);
   
   if ( mssmhbb.config()->override() ) 
      std::cout << "*w* The event selection was overridden" << std::endl;
      
   // Analysis of events
   std::cout << "This analysis has " << mssmhbb.analysis()->size() << " events" << std::endl;
   
// 
   int counter = 0;
      
   int maxevt = mssmhbb.config()->nEventsMax();
   if ( maxevt < 0 ) maxevt = mssmhbb.analysis()->size();
   for ( int i = 0 ; i < maxevt ; ++i )
   {
      if ( i > 0 && i%100000==0 ) std::cout << i << "  events processed! " << std::endl;
      bool goodEvent = mssmhbb.event(i);
      if ( ! goodEvent ) continue;
      
      if ( ! mssmhbb.config()->override() )
      {
         ++counter;  
         continue;
      }
   
      // histograms
//      auto h1 = mssmhbb.H1Fs();
      auto h_cut = mssmhbb.H1F("cutflow");
         
      if ( ! mssmhbb.selectionTrigger() ) continue;
      h_cut -> Fill(0);
         
      if ( mssmhbb.analysisWithJets() )
      {
         if ( ! mssmhbb.selectionJetId() ) continue;
         h_cut -> Fill(1);
         
         // you now have a list of selectedJets that you can play with
//         std::vector<Jet*> selectedJets = mssmhbb.selectedJets();  FIXME
               
         // standard jet selection
         if ( ! mssmhbb.Analyser::selectionJet(1) ) continue;
         if ( ! mssmhbb.Analyser::selectionJet(2) ) continue;
         h_cut -> Fill(2);
         
         // additional jet selection for MssmHbb
         if ( ! mssmhbb.selectionJetDeta(1,2,1.5) ) continue;
         h_cut -> Fill(3);
         
         // matching to online jets
         if ( ! mssmhbb.onlineJetMatching() ) continue;
         h_cut -> Fill(4);
         
         // btag of two leading jets
         if ( ! mssmhbb.selectionBJet(1) ) continue;
         if ( ! mssmhbb.selectionBJet(2) ) continue;
         h_cut -> Fill(5);
         
         // matching to online btag objects
         if ( ! mssmhbb.onlineBJetMatching() ) continue;
         h_cut -> Fill(6);
         
         // 3rd jet selection
         if ( ! mssmhbb.Analyser::selectionJet(3) ) continue;
         h_cut -> Fill(7);
         
         // delta R jet selection
         if ( ! mssmhbb.Analyser::selectionJetDr(1,2,-1.) ) continue;
         if ( ! mssmhbb.Analyser::selectionJetDr(1,3,-1.) ) continue;
         if ( ! mssmhbb.Analyser::selectionJetDr(2,3,-1.) ) continue;
         h_cut -> Fill(8);
         
         
         if ( mssmhbb.config()->signalRegion() )
         {
            if ( ! mssmhbb.selectionBJet(3) ) continue;
         }
         else
         {
            if ( ! mssmhbb.selectionNonBJet(3) ) continue;
         }
         h_cut -> Fill(9);
         
         mssmhbb.fillJetHistograms();
      }
      
   }
      
   mssmhbb.end();
   
   
   std::cout << counter << " events passed the selection" << std::endl;
   
} //end main


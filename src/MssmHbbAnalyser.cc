/**\class MssmHbb MssmHbbAnalyser.cc Analysis/Tools/src/MssmHbbAnalyser.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberval Walsh Bastos Rangel
//         Created:  Mon, 20 Oct 2014 14:24:08 GMT
//
//

// system include files
#include <iostream>
// 
// user include files
#include "Analysis/MssmHbb/interface/MssmHbbAnalyser.h"


//
// class declaration
//

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::mssmhbb;

//
// constructors and destructor
//
MssmHbbAnalyser::MssmHbbAnalyser()
{
}

MssmHbbAnalyser::MssmHbbAnalyser(int argc, char ** argv) : Analyser(argc,argv)
{
   histograms("jet",config_->nJetsMin());
   for ( int i = 0; i < 20; ++i ) cutflow_.push_back(0);
   
}

MssmHbbAnalyser::~MssmHbbAnalyser()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//
// ------------ method called for each event  ------------


bool MssmHbbAnalyser::event(const int & i)
{
   // parent function checks only json and run range validity
   if ( ! Analyser::event(i) ) return false;
   
   if ( config_->override() )  return true;
   
   if ( ! selectionTrigger() ) return false;
   h1_["cutflow"] -> Fill(0);
      
   if ( analysisWithJets() )
   {
      if ( ! selectionJetId() ) return false;
      h1_["cutflow"] -> Fill(1);
            
      // standard jet selection
      if ( ! Analyser::selectionJet() ) return false;
      h1_["cutflow"] -> Fill(2);
      
      // additional jet selection for MssmHbb
      if ( ! selectionJet() ) return false;
      h1_["cutflow"] -> Fill(3);
      
      // matching to online jets
      if ( ! onlineJetMatching() ) return false;
      h1_["cutflow"] -> Fill(4);
      
      // btag of two leading jets
      if ( ! selectionBJet(1) ) return false;
      if ( ! selectionBJet(2) ) return false;
      h1_["cutflow"] -> Fill(5);
      
      // matching to online btag objects
      if ( ! onlineBJetMatching() ) return false;
      h1_["cutflow"] -> Fill(6);
      
      if ( config_->signalRegion() )
      {
         if ( ! selectionBJet(3) ) return false;
      }
      else
      {
         if ( ! selectionNonBJet(3) ) return false;
      }
      h1_["cutflow"] -> Fill(7);
      
      fillJetHistograms();
      
   }
      
   return true;
}

bool MssmHbbAnalyser::selectionJet()
{
   bool isgood = true;
   
   // jet kinematics and btag
   std::map<std::string,bool> isOk;
   for ( int j = 0; j < config_->nJetsMin() ; ++j )
   {
      for ( int k = j+1; k < config_->nJetsMin() && j < config_->nJetsMin(); ++k )
      {
         isOk[Form("dr%d%d",j,k)]   = true;
         isOk[Form("deta%d%d",j,k)] = true;
      }
   }
   // kinematic selection
   for ( int j = 0 ; j < config_->nJetsMin() ; ++j )
   {
      // delta R between jets
      for ( int k = j+1; k < config_->nJetsMin() && j < config_->nJetsMin(); ++k )
         if ( selectedJets_[j]->deltaR(*selectedJets_[k]) < config_->drmin_ )                                            isOk[Form("dr%d%d",j,k)]   = false;
   }
   // delta eta 2 leading jets
   if ( config_->nJetsMin() > 1 )
      if ( fabs(selectedJets_[0]->eta() - selectedJets_[1]->eta()) > config_->detamax_ && !(config_->detamax_ < 0) )     isOk[Form("deta%d%d",0,1)] = false;
   
   for ( auto & ok : isOk )
      isgood = ( isgood && ok.second );
   
   return isgood;
}

void MssmHbbAnalyser::histograms(const std::string & obj, const int & n)
{
   Analyser::histograms(obj,n);
   
}


void MssmHbbAnalyser::end()
{
   Analyser::end();
   
}


#ifndef Analysis_MssmHbb_MssmHbbAnalyser_h
#define Analysis_MssmHbb_MssmHbbAnalyser_h 1

// -*- C++ -*-
//
// Package:    Analysis/MssmHbb
// Class:      Analysis
// 
/**\class Analysis MssmHbbAnalyser.cc Analysis/MssmHbb/src/MssmHbbAnalyser.cc

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
#include <memory>
#include <vector>
#include <string>
// 
// user include files

#include "Analysis/Tools/interface/Analyser.h"

//
// class declaration
//

namespace analysis {
   namespace mssmhbb {

      class MssmHbbAnalyser : public analysis::tools::Analyser {
         
         public:
            MssmHbbAnalyser();
            MssmHbbAnalyser(int argc, char ** argv);
           ~MssmHbbAnalyser();
           
            virtual bool event(const int &i);
            bool preselection();
            virtual bool muonJet(const bool & swap = false);
            void fillMssmHbbTree();
            void mssmHbbTree();
            void fillMssmHbbHistograms();
            bool jetCorrections();
            bool triggerSelection();
            bool jetSelection();
            bool muonSelection();
            bool btagSelection();
            bool endSelection();
            // ----------member data ---------------------------
         protected:
            bool do_tree_;
            std::shared_ptr<TTree> mssmhbb_tree_;
            double mbb_;
            double mbbw_;
               
         private:
               

      };
   }
}

#endif  // Analysis_MssmHbb_MssmHbbAnalyser_h

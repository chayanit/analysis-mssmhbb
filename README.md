# `Analysis/MssmHbb`

**Extension of the Analysis Framework for the MSSM Hbb all-hadronic and semileptonic analyses in CMS**

See also the code documentation [Doxygen](https://www.desy.de/~walsh/analysis-framework/doxygen/pro) page


**This package inherits from the [Analysis/Tools](https://github.com/desy-cms/analysis-tools) package**

* [Installation](#Installation)
* [Running the analysis](#Running-the-analysis)
  * [Semileptonic channel](#Semileptonic-channel)
  * [All-hadronic channel](#All-hadronic-channel) (to do)
  * [Main macro](#Main-macro)
  * [Configuration file](#Configuration-file)
  * [Workflow](#Workflow)
  * [Output](#Output)
* [MssmHbbAnalyser class](#MssmHbbAnalyser-class)



## Installation

The codes as well as the ntuples are independent of CMSSW. However, in order to compile it uses `scram`. So the latest version in any architecture should be fine.


```bash
cmsrel CMSSW_X_Y_Z
cd CMSSW_X_Y_Z/src
cmsenv

git clone https://github.com/desy-cms/analysis-tools.git Analysis/Tools
git clone https://github.com/desy-cms/analysis-mssmhbb.git Analysis/MssmHbb
git clone https://github.com/desy-cms/analysis-calibrations.git Analysis/Tools/data/calibrations
git clone https://github.com/desy-cms/analysis-ntuples.git Analysis/MssmHbb/test/ntuples

scram b -j4 USER_CXXFLAGS="-Wno-misleading-indentation"
hash -r
```

:zzz: The `USER_CXXFLAGS="-Wno-misleading-indentation"` prevents a large number of warnings from misleading identation in modules of the boost library.

:zzz: More information about the additional packages, `calibration`, `ntuples` etc can be found in the `Analysis/Tools` [README.md](https://github.com/desy-cms/analysis-tools/blob/develop/README.md).


## Running the analysis

A macro with the full analysis workflow for both semileptonic and all-hadronic channels can be found in [Analysis/MssmHbb/bin/MssmHbbAnalysis.cc](bin/MssmHbbAnalysis.cc)

### Semileptonic channel

The configuration file for the semileptonic analysis with 2017 data is available in the file [Analysis/MssmHbb/test/mssmhbb_semilep_2017.cfg](test/mssmhbb_semilep_2017.cfg).

#### The semileptonic standard selection

:warning: Unfortunately Github markdown does not render Latex equations :rage:

* Trigger
  * HLT: HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33_v
  * L1: L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6)
* Jet
  * Jet energy corrections applied
  * b-jet energy regression applied
  * Jet ID: tight
  * Jet PU ID: loose
  * Jet 1: $p_{T} > 60$ GeV; $|\eta| < 2.2$; BTag: DeepFlavour Medium WP; matched to trigger objects
  * Jet 2: $p_{T} > 50$ GeV; $|\eta| < 2.2$; BTag: DeepFlavour Medium WP; matched to trigger objects
  * Jet 3: $p_{T} > 30$ GeV; $|\eta| < 2.2$; BTag: DeepFlavour Medium WP
    * Reverse btag Loose WP in control region
  * $\Delta\eta(jet_1,jet_2) < 1.5$
  * $\Delta R (jet_{i},jet_{j}) > 1$ ; $i,j=1,2,3$ ; $i\neq j$
* Muon
  * Muon ID: tight
  * Muon: $p_{T} > 13$ GeV; $|\eta| < 2.2$; matched to trigger objects
* Muon-Jet
  * At least one selected muon must be within one of the two leading jet: $\Delta R(\mu,jet)<0.4$
  

### All-hadronic channel

The configuration file for the semileptonic analysis is available in the file [Analysis/MssmHbb/test/mssmhbb_semilep_2017.cfg](test/mssmhbb_semilep_2017.cfg).

#### The all-hadronic standard selection

**TO DO** 

:warning: Unfortunately Github markdown does not render Latex equations :rage:


### Main macro
The structure of the macro is as follows:

```cpp
#include "Analysis/MssmHbb/interface/MssmHbbAnalyser.h"

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::mssmhbb;

int main(int argc, char ** argv)
{
    TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

    MssmHbbAnalyser mssmhbb(argc,argv);

    for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
    {
        if ( ! mssmhbb.event(i)             )  continue;    // read event, run selection/json
        if ( ! mssmhbb.triggerSelection ()  )  continue;    // trigger
        if ( ! mssmhbb.preselection()       )  continue;    // preselection
        if ( ! mssmhbb.jetCorrections()     )  continue;    // jet corrections
        if ( ! mssmhbb.muonSelection()      )  continue;    // muon
        if ( ! mssmhbb.jetSelection()       )  continue;    // jets
        if ( ! mssmhbb.muonJet()            )  continue;    // muon-jet association
        if ( ! mssmhbb.btagSelection()      )  continue;    // btagging
        mssmhbb.endSelection();
    }
} // end main
```

Compared to the Analysis/Tools [example](https://github.com/desy-cms/analysis-tools/blob/develop/README.md#example) this is much simplified, for the individual selections were grouped and included in methods of the `MssmHbbAnalyser` class. However some flexibility to re-order the workflow is still allowed. Of course, the user can also replace any of the `if` statements above and enable more detailed steps of the selections, or add other steps to the workflow.

#### Event selection detailed description

1. `MSSMHBBAnalyser::event(i)`

When calling the method `event(i)`, the entry `i` is read from the `TTree` of the ntuple. If running on data, the event is checked whether it is valid from a list of certified data given by the provided JSON file.

When running on Monte Carlo, generated event weights are applied. Pileup re-weighing is also applied if re-weigh file is provided in the configuration file in the `[Corrections]` section, e.g.:

```ini
[Corrections]
Pileup.reweight  = calibrations/2017/pileup/PileupWeight_Run2017_Mix_2017.root
```


2. `MSSMHBBAnalyser::triggerSelection()`

It returns `true` if the event passes both the HLT and L1 triggers, which are given in the configuration file section and whose results are in the `TriggerResults` collection, thus

```ini
[Trigger]
results   = TriggerResults
hltPath   = HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33_v 
l1Seed    = L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6
```

3. `MSSMHBBAnalyser::preselection()`

This method will prepare a vector of muons and a vector of jets that pass the identification criteria, namely muon id for muons, and jet id and jet pileup id for jets. It will return `false` if any of those vectors has zero entries.

The identification criteria is defined in the configuration file as, for example

```ini
[Jets]
id = tight
puId = loose
```
```ini
[Muons]
id = tight
```

4. `MSSMHBBAnalyser::jetCorrections()`

This method will apply b-jet energy regression corrections to all the jets in the event, if requested in the configuration file:
```ini
[Corrections]
Jets.bRegression = true
```
It will also smear the jet $p_{t}$ according to the jet energy resolution scale factors, if the relevant calibration files are provided in the configuration, e.g.:
```ini
[Corrections]
Jets.bRegression = true
Jets.jerPtRes    = calibrations/2017/jer/Fall17_V3_MC_PtResolution_AK4PFchs.txt
Jets.jerSF       = calibrations/2017/jer/Fall17_V3_MC_SF_AK4PFchs.txt
```

In the `jetCorrections` method the returned value is always `true`.

5. `MSSMHBBAnalyser::muonSelection()`

The `muonSelection` method selects events according to the minimum number of muons required. It also skims the vector of muons to contain only muons with the required kinematics and matching to muon trigger objects.

Example of the requirements give in the configuration:
```ini
[Muons]
ptMin = 13.
etaMax = 2.2
```
```ini
[Trigger.Objects]
Muons.L1  = hltL1Mu12er2p3Jet40er2p3dRMax0p4DoubleJet40er2p3dEtaMax1p6_l1muon
Muons.L3  = hltL3fL1sMu12Dijet40L1f0L2f8L3Filtered12
```

6. `MSSMHBBAnalyser::jetSelection()`

The `jetSelection` method selects events according to the minimum number of jets required. It also checks whether the event fulfils the required kinematic selection and matching of the two leading jets to jet trigger objects. It returns `false` if any of the requirements is not fulfilled.

For a selection, in which:
  * Jet 1: $p_{t} > 60$ GeV; $|\eta| < 2.2$
  * Jet 2: $p_{t} > 50$ GeV; $|\eta| < 2.2$
  * Jet 3: $p_{t} > 30$ GeV; $|\eta| < 2.2$
  * $\Delta\eta(jet_1,jet_2) < 1.5$
  * $\Delta R (jet_{i},jet_{j}) > 1$ ; $i,j=1,2,3$ ; $i\neq j$
  

the configuration reads

```ini
[Jets]
ptMin = 60.
ptMin = 50.
ptMin = 30.
etaMax = 2.2
etaMax = 2.2
etaMax = 2.2
...
```
:warning: Notice that the order matters! First `ptMin` entry is for the leading jet, second `ptMin` is for the second leading jet, ...

```ini
[Trigger.Objects]
Jets.L1   = hltL1Mu12er2p3Jet40er2p3dRMax0p4DoubleJet40er2p3dEtaMax1p6_l1jet
Jets.Calo = hltDoubleCaloBJets30eta2p3
Jets.PF   = hltDoublePFBJets40Eta2p3 
```

7. `MSSMHBBAnalyser::muonJet()`

The `muonJet` method returns true if either of the two leading jets contains at least one of the selected muons.


8. `MSSMHBBAnalyser::btagSelection()`

If the leading jets fulfil the required b-tag working point selection and matching of the two leading jets to the corresponding trigger objects, the `btagSelection` method returns `true`. The requirements are controlled in the configuration file. For example, select events with three leading jets passing the `medium` working point of the deepflavour algorithm.

```ini
[BTag]
nMin  = 3
wp = medium
wp = medium
wp = medium
algorithm = deepflavour
```
where `wp` if the working point of each required b-tagged jet.

:warning: Notice that the order matters! First `wp` entry is for the leading jet, second `wp` is for the second leading jet, ...


The values of the WP's must be given in the configuration file:
```ini
[BTag]
loose  = 0.0521
medium = 0.3033
tight  = 0.7489
```

Because in the MSSM Hbb analysis, the background modelling is performed in a control region, in which on the b-tagged jets has its b-tagging requirement reversed, one can apply this possibility in the configuration files by setting the flag

```ini
[Info]
signalRegion = false
```
and defininig reverse b tag WP and the rank of the jet to have the reverse b-tag applied, e.g. reverse b tag of the 3rd leading jet with `loose` WP:
```ini
[BTag]
revWP = loose
revBJet = 3
```


The matching to trigger object is done given the appropriate trigger object collection, e.g.
```ini
[Trigger.Objects]
BTag.Calo = hltBTagCalo30x8CSVp0p92DoubleWithMatching
```

### Configuration file

The structure of the configuration file is defined by the `Analysis/Tools` package. More details can be found in the `Analysis/Tools` [README.md](https://github.com/desy-cms/analysis-tools/blob/develop/README.md#detailed-description).


### Workflow

The workflow of the main macro for DATA is

```
+----------------------------------------------------------------------------------------------+
| Workflow DATA                                                                                |
+----------------------------------------------------------------------------------------------+
|  0 - Total events read                                                                       |
|  1 - Certified data: Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt            |
|  2 - Using Jet collection: MssmHbb/Events/updatedPatJets                                     |
|  3 - Using Muon collection: MssmHbb/Events/slimmedMuons                                      |
|  4 - HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33_v                               |
|  5 - L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6                          |
|  6 - MuonId: tight                                                                           |
|  7 - JetId: tight                                                                            |
|  8 - JetPileupId: loose                                                                      |
|  9 - b jet energy regression                                                                 |
| 10 - NMuons >= 1                                                                             |
| 11 - Muons selected: pt >  13.0 GeV and |eta| < 2.2                                          |
| 12 - Muons selected: Online muon matching (L1 and L3)                                        |
| 13 - NJets >= 3                                                                              |
| 14 - Jet 1: pt >  60.0 GeV and |eta| < 2.2                                                   |
| 15 - Jet 1: online jet match                                                                 |
| 16 - Jet 2: pt >  50.0 GeV and |eta| < 2.2                                                   |
| 17 - Jet 2: online jet match                                                                 |
| 18 - Jet 3: pt >  30.0 GeV and |eta| < 2.2                                                   |
| 19 - MSSMHbb Semileptonic: Jet-muon association                                              |
| 20 - Jet 1: deepflavour btag > 0.3033 (medium)                                               |
| 21 - Jet 1: online b jet match                                                               |
| 22 - Jet 2: deepflavour btag > 0.3033 (medium)                                               |
| 23 - Jet 2: online b jet match                                                               |
| 24 - Jet 3: deepflavour btag < 0.0521 (loose) [reverse btag]                                 |
| 25 - *** Filling jets histograms - final_selection                                           |
| 26 - Fill MssmHbb Histograms                                                                 |
| 27 - Fill MssmHbb tree                                                                       |
+----------------------------------------------------------------------------------------------+

```

where as for the Monte Carlo it is

```
+----------------------------------------------------------------------------------------------+
| Workflow Monte Carlo                                                                         |
+----------------------------------------------------------------------------------------------+
|  0 - Total events read                                                                       |
|  1 - Generated weighted events (sign of weights)                                             |
|  2 - Using GenParticles collection: MssmHbb/Events/prunedGenParticles                        |
|  3 - Using GenJets collection: MssmHbb/Events/slimmedGenJets                                 |
|  4 - Using Jet collection: MssmHbb/Events/updatedPatJets                                     |
|  5 - Using Muon collection: MssmHbb/Events/slimmedMuons                                      |
|  6 - Pileup weight (PileupWeight_Run2017_Mix_2017.root)                                      |
|  7 - HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33_v                               |
|  8 - L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6                          |
|  9 - MuonId: tight                                                                           |
| 10 - JetId: tight                                                                            |
| 11 - JetPileupId: loose                                                                      |
| 12 - b jet energy regression                                                                 |
| 13 - JER smearing (Fall17_V3_MC_PtResolution_AK4PFchs.txt,Fall17_V3_MC_SF_AK4PFchs.txt)      |
| 14 - NMuons >= 1                                                                             |
| 15 - Muons selected: pt >  13.0 GeV and |eta| < 2.2                                          |
| 16 - Muons selected: Online muon matching (L1 and L3)                                        |
| 17 - NJets >= 3                                                                              |
| 18 - Jet 1: pt >  60.0 GeV and |eta| < 2.2                                                   |
| 19 - Jet 1: online jet match                                                                 |
| 20 - Jet 2: pt >  50.0 GeV and |eta| < 2.2                                                   |
| 21 - Jet 2: online jet match                                                                 |
| 22 - Jet 3: pt >  30.0 GeV and |eta| < 2.2                                                   |
| 23 - MSSMHbb Semileptonic: Jet-muon association                                              |
| 24 - Jet 1: deepflavour btag > 0.3033 (medium)                                               |
| 25 - Jet 1: btag SF applied (deepflavour medium WP)                                          |
| 26 - Jet 1: online b jet match                                                               |
| 27 - Jet 2: deepflavour btag > 0.3033 (medium)                                               |
| 28 - Jet 2: btag SF applied (deepflavour medium WP)                                          |
| 29 - Jet 2: online b jet match                                                               |
+----------------------------------------------------------------------------------------------+

```

### Output

The macro produces a root file containing a directory with jet-related histograms and a TTree containing the dijet mass $m_{bb}$  of the two leading b-tagged jets and the event weights after the full selection.

## MssmHbbAnalyser class

:warning: To do





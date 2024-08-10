/******************************************************************************
 *  [File] 
 *      MultiLepPAT.cc 
 *  [Class]      
 *      MultiLepPAT 
 *  [Directory] 
 *      myAnalyzers/MultiLepPAT/src/MultiLepPAT.cc
 *  [Description]
 *      Make rootTuple for JPsiKKK reconstruction
 *  [Implementation]
 *     <Notes on implementation>
 *  [Original Author]
 *  
 *  [Update Log]
 *      20240626 [Eric Wang]
 *          Added necessary annotation for the code
 *          Modified for TPS analysis on pp > Jpsi + Jpsi + Ups
 *              (Decay channels: Jpsi > mu- + mu+, Ups > mu- + mu+)
 *  [Note]
 *      20240626 [Eric Wang]
 *          Upsilon is abbreviated as "Ups"
******************************************************************************/



#ifndef _MultiLepPAT_h
#define _MultiLepPAT_h

// system include files
#include <memory>
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // xining MINIAODtest

// user include files
#include "../interface/VertexReProducer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
//#include "RecoVertex/V0Producer/interface/V0Producer.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "PhysicsTools/RecoUtils/interface/CheckHitPattern.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>
#include <utility>

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

//
// class decleration
//

using std::vector;
using namespace edm;
using namespace reco;
using namespace std;

class MultiLepPAT : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit MultiLepPAT(const ParameterSet&);
    ~MultiLepPAT();

  
private:
    virtual void beginJob() ;
    virtual void beginRun(Run const & iRun, EventSetup const& iSetup);
    virtual void analyze(const Event&, const EventSetup&);
    virtual void endJob() ;

 
    //add token here
    edm::EDGetTokenT< L1GlobalTriggerReadoutRecord    > gtRecordToken_;
    edm::EDGetTokenT< BeamSpot                        > gtbeamspotToken_;
    edm::EDGetTokenT< VertexCollection                > gtprimaryVtxToken_;
    edm::EDGetTokenT< edm::View<pat::Muon>            > gtpatmuonToken_; // MINIAOD
    edm::EDGetTokenT< edm::TriggerResults             > gttriggerToken_;
    edm::EDGetTokenT< edm::View<pat::PackedCandidate> > trackToken_; // MINIAOD
    edm::EDGetTokenT< reco::GenParticleCollection     > genParticlesToken_;

    // ----------member data ---------------------------

    virtual int GetHitsBefore(const EventSetup& setup, const vector<reco::Track>& trks, RefCountedKinematicTree& fitTr){
        return -1;
    }

    virtual int GetMissesAfter(const EventSetup& setup, const vector<reco::Track>& trks, RefCountedKinematicTree& VrtxTree){
        return -1;
    }

    //get ctau from beamspot
    virtual double GetcTau( RefCountedKinematicVertex&   decayVrtx, 
                            RefCountedKinematicParticle& kinePart, 
                            Vertex&                             bs ){	
        TVector3 vtx;
        TVector3 pvtx;
        vtx.SetXYZ((*decayVrtx).position().x(), (*decayVrtx).position().y(), 0);
        pvtx.SetXYZ(bs.position().x(), bs.position().y(), 0);
        VertexDistanceXY vdistXY;
        TVector3 pperp(kinePart->currentState().globalMomentum().x(),
	    	   kinePart->currentState().globalMomentum().y(), 0);

        TVector3 vdiff = vtx - pvtx;
        double cosAlpha = vdiff.Dot(pperp) / (vdiff.Perp() * pperp.Perp());
        Measurement1D distXY = vdistXY.distance(Vertex(*decayVrtx), Vertex(bs));
        double ctauPV = distXY.value() * cosAlpha * kinePart->currentState().mass() / pperp.Perp();
        return ctauPV;    
    }

    virtual double GetcTauErr(  RefCountedKinematicVertex& decayVrtx, 
                                RefCountedKinematicParticle& kinePart, 
                                Vertex& bs                              ){       
        TVector3 pperp(kinePart->currentState().globalMomentum().x(),
		               kinePart->currentState().globalMomentum().y(), 
                       0                                              );
        AlgebraicVector vpperp(3);
        vpperp[0] = pperp.x();
        vpperp[1] = pperp.y();
        vpperp[2] = 0.;

        GlobalError v1e = (Vertex(*decayVrtx)).error();
        GlobalError v2e = bs.error();
        AlgebraicSymMatrix vXYe = asHepMatrix(v1e.matrix()) + asHepMatrix(v2e.matrix());
        double ctauErrPV = sqrt(vXYe.similarity(vpperp)) * kinePart->currentState().mass() / (pperp.Perp2());

        return ctauErrPV;    
    }

    virtual bool tracksToMuonPair(vector<RefCountedKinematicParticle>&        arg_MuonResults,
                                  KinematicParticleFactoryFromTransientTrack& arg_MuFactory,
                                  const MagneticField&                        arg_bField,
                                  const TrackRef& arg_Trk1,   const TrackRef& arg_Trk2       );

    double deltaR(double eta1, double phi1, double eta2, double phi2) {
        double deta = eta1 - eta2;
        double dphi = phi1 - phi2;
        while (dphi >   M_PI) dphi -= 2*M_PI;
        while (dphi <= -M_PI) dphi += 2*M_PI;
        return sqrt(deta*deta + dphi*dphi);
    }
    
    // Member data

    // Essentials [Annotation by Eric Wang, 20240626]
    string proccessName_;
    HLTConfigProvider hltConfig_;
    
    InputTag hlTriggerResults_;
    InputTag inputGEN_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord>  magneticFieldToken_;
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord>  theTTBuilderToken_;
    string  vtxSample;
    bool    doMC;
    int     MCParticle;
    bool    doJPsiMassCost;
    int     MuPixHits_c;
    int     MuSiHits_c;
    double  MuNormChi_c;
    double  MuD0_c;

    // Limits for secondary particles [Annotation by Eric Wang, 20240626]
    
    double JMaxM_c;
    double JMinM_c;
    int PiSiHits_c;
    double MuPt_c;
    double JPiPiDR_c;
    double XPiPiDR_c;
    bool UseXDr_c;	
    double JPiPiMax_c;
    double JPiPiMin_c;

    // Trigger info [Annotation by Eric Wang, 20240626]
    
    bool resolveAmbiguity_; 
    bool addXlessPrimaryVertex_;
    vector<string>      TriggersForJpsi_;
    vector<string>      FiltersForJpsi_;
    vector<string>      TriggersForUpsilon_;
    vector<string>      FiltersForUpsilon_;
    
    int JpsiMatchTrig[50], UpsilonMatchTrig[50];
    
    vector<string>      TriggersForMatching_;
    vector<string>      FiltersForMatching_;
    int  MatchingTriggerResult[50];
    bool Debug_;
    double Chi_Track_;

    // Constructing TTree object [Annotation by Eric Wang, 20240626]
    
    TTree* X_One_Tree_;
    
    unsigned int        runNum, evtNum, lumiNum;
    unsigned int        nGoodPrimVtx;
    
    vector<unsigned int>*   trigRes;
    vector<std::string>*    trigNames;
    vector<unsigned int>*   L1TT;
    vector<std::string>*    MatchTriggerNames;

    // primary vertices [Annotation by Eric Wang, 20240626]
    float               priVtxX,    priVtxY,    priVtxZ, 
                        priVtxXE,   priVtxYE,   priVtxZE, 
                        priVtxChiNorm, priVtxChi, priVtxCL;
    vector<float>       *PriVtxXCorrX, *PriVtxXCorrY, *PriVtxXCorrZ;
    vector<double>      *PriVtxXCorrEX, *PriVtxXCorrEY, *PriVtxXCorrEZ;
    vector<float>	    *PriVtxXCorrC2, *PriVtxXCorrCL;
    
    // all muons: kinematic and other information [Annotation by Eric Wang, 20240626]
    unsigned int        nMu;
    vector<float>       *muPx, *muPy, *muPz, 
                        *muD0, *muD0E, *muDz, 
                        *muChi2, *muGlChi2,  *mufHits;   
    vector<bool>        *muFirstBarrel, *muFirstEndCap;
    vector<float>       *muDzVtx, *muDxyVtx;   
    vector<int>         *muNDF, *muGlNDF, *muPhits, *muShits, *muGlMuHits, *muType, *muQual;
    vector<int>         *muTrack;
    vector<float>       *muCharge;
    vector<float>       *muIsoratio;

    // all muons: selection result [Annotation by Eric Wang, 20240626]
    vector<int>         *muIsGoodLooseMuon,         *muIsGoodLooseMuonNew, 
                        *muIsGoodSoftMuonNewIlse,   *muIsGoodSoftMuonNewIlseMod, 
                        *muIsGoodTightMuon,         *muIsJpsiTrigMatch,         
                        *muIsUpsTrigMatch,          *munMatchedSeg;
    vector<int>         *muIsPatLooseMuon, *muIsPatTightMuon, *muIsPatSoftMuon, *muIsPatMediumMuon;
    
    //for Maksat trigger match [Annotation by Eric Wang, 20240626]
    vector<int> *muUpsVrtxMatch, *muL3TriggerMatch;
    
    //added by zhenhu for MuonID
    vector<float>  *muMVAMuonID, *musegmentCompatibility; 
    //for Stoyan slope pull
    vector<float>  *mupulldXdZ_pos_noArb, *mupulldYdZ_pos_noArb;
    //addition 3,4,5
    vector<float>  *mupulldXdZ_pos_ArbDef, *mupulldYdZ_pos_ArbDef;
    vector<float>  *mupulldXdZ_pos_ArbST, *mupulldYdZ_pos_ArbST;
    vector<float>  *mupulldXdZ_pos_noArb_any, *mupulldYdZ_pos_noArb_any;


    // Index of selected muons in one event [Added by Eric Wang, 20240626]
    // Reference: code from Wang Xining
    vector<float> *Jpsi1_mu1_Idx,   *Jpsi1_mu2_Idx, *Jpsi2_mu1_Idx, *Jpsi2_mu2_Idx;
    vector<float> *Ups1_mu1_Idx,    *Ups1_mu2_Idx;

    // Vertex fitting [Modified by Eric Wang, 20240626]
    vector<float> *VtxProb,      *Chi2,           *ndof,        *VtxPt2;

    // Secondary Jpsi's and Upsilon's kinematics [Modified by Eric Wang, 20240626]
    vector<float> *Jpsi1_mass,   *Jpsi1_VtxProb,  *Jpsi1_Chi2,  *Jpsi1_ndof,
                  *Jpsi2_mass,   *Jpsi2_VtxProb,  *Jpsi2_Chi2,  *Jpsi2_ndof,
                  *Jpsi1_px,     *Jpsi1_py,       *Jpsi1_pz,    *Jpsi1_massErr, 
                  *Jpsi2_px,     *Jpsi2_py,       *Jpsi2_pz,    *Jpsi2_massErr;
    vector<float> *Ups_mass,     *Ups_VtxProb,    *Ups_Chi2,    *Ups_ndof,
                  *Ups_px,       *Ups_py,         *Ups_pz,      *Ups_massErr;

    // Secondary Jpsi's and Upsilon's kinematics fitted under constraints 
    //      [Modified by Eric Wang, 20240626]
    vector<float> *CS_Jpsi1_mass,   *CS_Jpsi1_VtxProb,  *CS_Jpsi1_Chi2,  *CS_Jpsi1_ndof,
                  *CS_Jpsi2_mass,   *CS_Jpsi2_VtxProb,  *CS_Jpsi2_Chi2,  *CS_Jpsi2_ndof,
                  *CS_Jpsi1_px,     *CS_Jpsi1_py,       *CS_Jpsi1_pz,    *CS_Jpsi1_massErr, 
                  *CS_Jpsi2_px,     *CS_Jpsi2_py,       *CS_Jpsi2_pz,    *CS_Jpsi2_massErr;
    vector<float> *CS_Ups_mass,     *CS_Ups_VtxProb,    *CS_Ups_Chi2,    *CS_Ups_ndof,
                  *CS_Ups_px,       *CS_Ups_py,         *CS_Ups_pz,      *CS_Ups_massErr;

    // MC results [Modified by Eric Wang, 20240626]
    // Kinematics
    vector<float> 
        *MC_G1_Jpsi1_px, *MC_G1_Jpsi1_py, *MC_G1_Jpsi1_pz, *MC_G1_Jpsi1_mass,
        *MC_G1_Jpsi2_px, *MC_G1_Jpsi2_py, *MC_G1_Jpsi2_pz, *MC_G1_Jpsi2_mass,
        *MC_G1_Ups1_px,  *MC_G1_Ups1_py,  *MC_G1_Ups1_pz,  *MC_G1_Ups1_mass;
    vector<float> 
        *MC_G2_Mu1_px, *MC_G2_Mu1_py, *MC_G2_Mu1_pz, *MC_G2_Mu1_mass,
        *MC_G2_Mu2_px, *MC_G2_Mu2_py, *MC_G2_Mu2_pz, *MC_G2_Mu2_mass,
        *MC_G2_Mu3_px, *MC_G2_Mu3_py, *MC_G2_Mu3_pz, *MC_G2_Mu3_mass,
        *MC_G2_Mu4_px, *MC_G2_Mu4_py, *MC_G2_Mu4_pz, *MC_G2_Mu4_mass,
        *MC_G2_Mu5_px, *MC_G2_Mu5_py, *MC_G2_Mu5_pz, *MC_G2_Mu5_mass,
        *MC_G2_Mu6_px, *MC_G2_Mu6_py, *MC_G2_Mu6_pz, *MC_G2_Mu6_mass ;

    // Particle ID

    vector<float> 
        *MC_G1_Jpsi1_PDG_ID,
        *MC_G1_Jpsi2_PDG_ID,
        *MC_G1_Ups1_PDG_ID;
    vector<float> 
        *MC_G2_Mu1_PDG_ID, 
        *MC_G2_Mu2_PDG_ID, 
        *MC_G2_Mu3_PDG_ID, 
        *MC_G2_Mu4_PDG_ID, 
        *MC_G2_Mu5_PDG_ID, 
        *MC_G2_Mu6_PDG_ID;
    vector<float> 
        *Match_mu1px, *Match_mu1py, *Match_mu1pz,
        *Match_mu2px, *Match_mu2py, *Match_mu2pz,
        *Match_mu3px, *Match_mu3py, *Match_mu3pz,
        *Match_mu4px, *Match_mu4py, *Match_mu4pz,
        *Match_mu5px, *Match_mu5py, *Match_mu5pz,
        *Match_mu6px, *Match_mu6py, *Match_mu6pz ;

    

// below are removed from the main code [Modified by Eric Wang, 20240626]
/*
    //xining
    // 
    vector<float> 
        *X_mu1Idx,       *X_mu2Idx,         *X_mu3Idx,      *X_mu4Idx,
        *X_mass,         *X_VtxProb,        *X_Chi2,        *X_ndof, 
        *X_px,           *X_py,             *X_pz,          *X_massErr, 
        *X_JPiPi_mass,   *X_JPiPi_VtxProb,  *X_JPiPi_Chi2,  *X_JPiPi_ndof,
        *X_Jpsi1_mass,   *X_Jpsi1_VtxProb,  *X_Jpsi1_Chi2,  *X_Jpsi1_ndof,
        *X_Jpsi2_mass,   *X_Jpsi2_VtxProb,  *X_Jpsi2_Chi2,  *X_Jpsi2_ndof,
        *X_JPiPi_px,     *X_JPiPi_py,       *X_JPiPi_pz,    *X_JPiPi_massErr, 
        *X_Jpsi1_px,     *X_Jpsi1_py,       *X_Jpsi1_pz,    *X_Jpsi1_massErr, 
        *X_Jpsi2_px,     *X_Jpsi2_py,       *X_Jpsi2_pz,    *X_Jpsi2_massErr, 
        *X_JPiPi_Pi1Idx, *X_JPiPi_Pi2Idx,
        *X_JPiPi_Pi1px,  *X_JPiPi_Pi1py,    *X_JPiPi_Pi1pz,
        *X_JPiPi_Pi2px,  *X_JPiPi_Pi2py,    *X_JPiPi_Pi2pz; 

    //added for mass constrain on 1208, cs stands for some kind of mass constraint
    // Variables for Jpsi and Pion
    vector <float>  *cs_X_Jpsi1_mass, *cs_X_Jpsi1_VtxProb,  *cs_X_Jpsi1_Chi2,   *cs_X_Jpsi1_ndof, 
                    *cs_X_Jpsi2_mass, *cs_X_Jpsi2_VtxProb,  *cs_X_Jpsi2_Chi2,   *cs_X_Jpsi2_ndof, 
                    *cs_X_Jpsi1_px,   *cs_X_Jpsi1_py,       *cs_X_Jpsi1_pz,     *cs_X_Jpsi1_massErr,
                    *cs_X_Jpsi2_px,   *cs_X_Jpsi2_py,       *cs_X_Jpsi2_pz,     *cs_X_Jpsi2_massErr        ;
    // Notice that Jpsi and Pion information is actually the same as before, 
    // because JPiPi MC fit doesn't change them by default.  
    // So Jpsi information retrived from Jpsi MC fit but not JPiPi MC fit. 
    // And Pion information is actually the same as that in 4mu2pi X fit so there are no variables for them
    // If you do what to change them, please check JPiPi MC fit and information stored carefully.

    // Variables for JPiPi fitted with MassConstraint Jpsi1
    vector<float> *cs_X_JPiPi_mass, *cs_X_JPiPi_VtxProb,    *cs_X_JPiPi_Chi2,   *cs_X_JPiPi_ndof, 
                  *cs_X_JPiPi_px ,  *cs_X_JPiPi_py,         *cs_X_JPiPi_pz,     *cs_X_JPiPi_massErr ;

    // Variables for X and JPiPi; MassConstraint to Psi2S and X3782 
    vector<float> *cs_X_mass_Psi2S,         *cs_X_VtxProb_Psi2S,        *cs_X_Chi2_Psi2S,       *cs_X_ndof_Psi2S, 
                  *cs_X_px_Psi2S,           *cs_X_py_Psi2S,             *cs_X_pz_Psi2S,         *cs_X_massErr_Psi2S,
                  *cs_X_JPiPi_mass_Psi2S,   *cs_X_JPiPi_VtxProb_Psi2S,  *cs_X_JPiPi_Chi2_Psi2S, *cs_X_JPiPi_ndof_Psi2S, 
                  *cs_X_JPiPi_px_Psi2S,     *cs_X_JPiPi_py_Psi2S,       *cs_X_JPiPi_pz_Psi2S,   *cs_X_JPiPi_massErr_Psi2S;

    vector<float> *cs_X_mass_X3872,         *cs_X_VtxProb_X3872,        *cs_X_Chi2_X3872,       *cs_X_ndof_X3872, 
                  *cs_X_px_X3872,           *cs_X_py_X3872,             *cs_X_pz_X3872,         *cs_X_massErr_X3872,
                  *cs_X_JPiPi_mass_X3872,   *cs_X_JPiPi_VtxProb_X3872,  *cs_X_JPiPi_Chi2_X3872, *cs_X_JPiPi_ndof_X3872, 
                  *cs_X_JPiPi_px_X3872,     *cs_X_JPiPi_py_X3872,       *cs_X_JPiPi_pz_X3872,   *cs_X_JPiPi_massErr_X3872;


//doMC
vector<float> 
*MC_X_px,
*MC_X_py,
*MC_X_pz,
*MC_X_mass,
*MC_Dau_Jpsipx,
*MC_Dau_Jpsipy,
*MC_Dau_Jpsipz,
*MC_Dau_Jpsimass,
*MC_Dau_psi2spx,
*MC_Dau_psi2spy,
*MC_Dau_psi2spz,
*MC_Dau_psi2smass,
*MC_Granddau_mu1px,
*MC_Granddau_mu1py,
*MC_Granddau_mu1pz,
*MC_Granddau_mu2px,
*MC_Granddau_mu2py,
*MC_Granddau_mu2pz,
*MC_Granddau_Jpsipx,
*MC_Granddau_Jpsipy,
*MC_Granddau_Jpsipz,
*MC_Granddau_Jpsimass,
*MC_Granddau_pi1px,
*MC_Granddau_pi1py,
*MC_Granddau_pi1pz,
*MC_Granddau_pi2px,
*MC_Granddau_pi2py,
*MC_Granddau_pi2pz,
*MC_Grandgranddau_mu3px,
*MC_Grandgranddau_mu3py,
*MC_Grandgranddau_mu3pz,
*MC_Grandgranddau_mu4px,
*MC_Grandgranddau_mu4py,
*MC_Grandgranddau_mu4pz;

vector<int>
*MC_X_chg,
*MC_Dau_JpsipdgId,
*MC_Dau_psi2spdgId,
*MC_Granddau_mu1pdgId,
*MC_Granddau_mu2pdgId,
*MC_Granddau_JpsipdgId,
*MC_Granddau_pi1pdgId,
*MC_Granddau_pi2pdgId,
*MC_Grandgranddau_mu3pdgId,
*MC_Grandgranddau_mu4pdgId;

vector<float> 
*Match_mu1px,
*Match_mu1py,
*Match_mu1pz,
*Match_mu2px,
*Match_mu2py,
*Match_mu2pz,
*Match_mu3px,
*Match_mu3py,
*Match_mu3pz,
*Match_mu4px,
*Match_mu4py,
*Match_mu4pz,

*Match_pi1px,
*Match_pi1py,
*Match_pi1pz,
*Match_pi2px,
*Match_pi2py,
*Match_pi2pz; 

  // float mybxlumicorr,myrawbxlumi;
  ////////////////////////  

  


*/
  
};

#endif

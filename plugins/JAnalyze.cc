// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
//
// class declaration
//

class MiniAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MiniAnalyzer(const edm::ParameterSet&);
      ~MiniAnalyzer();

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;


      // ----------member function-----------------------
      bool JElectronID(const edm::Ptr<pat::Electron> & electron, const int & elecIDLevel, const std::vector<edm::Handle<edm::ValueMap<unsigned int>>> & electron_bitmaps, const std::vector<edm::Handle<edm::ValueMap<vid::CutFlowResult>>> & electron_cutflows); 
      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Electron> > electronViewToken_;
      std::vector<edm::EDGetTokenT<edm::ValueMap<unsigned int>>> electron_bitmaps_;
      std::vector<edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>>> electron_cutflows_;
      unsigned electron_nWP_;
};

MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig):
    electronViewToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons")))
{
    // Strings for working points to use
    std::vector<std::string> vwp = {
      "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto",
      "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose",
      "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium",
      "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight",
    };
    // There will be 4 working points
    electron_nWP_ = vwp.size();
    // Add ValueMaps from the VIDElectronIDProducer
    for (auto const& wp : vwp) {
      electron_bitmaps_.push_back(consumes<edm::ValueMap<unsigned int>>(edm::InputTag(wp + std::string("Bitmap"))));
      electron_cutflows_.push_back(consumes<edm::ValueMap<vid::CutFlowResult>>(edm::InputTag(wp)));
    }
}

MiniAnalyzer::~MiniAnalyzer()
{
}


void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::vector<edm::Handle<edm::ValueMap<unsigned int>>> electron_bitmaps(electron_nWP_);
    for (unsigned iWP = 0; iWP < electron_nWP_; ++iWP) iEvent.getByToken(electron_bitmaps_[iWP], electron_bitmaps[iWP]);
    std::vector<edm::Handle<edm::ValueMap<vid::CutFlowResult>>> electron_cutflows(electron_nWP_);
    for (unsigned iWP = 0; iWP < electron_nWP_; ++iWP) iEvent.getByToken(electron_cutflows_[iWP], electron_cutflows[iWP]);

    edm::Handle<edm::View<pat::Electron> > electronsView;
    iEvent.getByToken(electronViewToken_, electronsView);
    for (auto const& obj : electronsView->ptrs()) {
      // elecIDLevel: 0: veto, 1: loose, 2: medium, 3: tight
      int t_elecIDLevel = 0;
      auto bitmap = (*(electron_bitmaps[t_elecIDLevel]))[obj];
      auto cutflow = (*(electron_cutflows[t_elecIDLevel]))[obj];
      std::cout<<"pt: "<<obj->pt()<<" eta: "<<obj->eta()<<" bitmap: "<<bitmap;
      for (unsigned iCut = 0; iCut < cutflow.cutFlowSize(); ++iCut) {
        std::cout<<" "<<cutflow.getNameAtIndex(iCut)<<": "<<(bitmap>>iCut&1);
      }

      // elecIDLevel: 0: veto, 1: loose, 2: medium, 3: tight
      int elecIDLevel = 0;
      std::cout<<" ID pass: "<<JElectronID(obj, elecIDLevel, electron_bitmaps, electron_cutflows);

      std::cout<<std::endl;
    }


    printf("\n");
}

bool MiniAnalyzer::JElectronID(const edm::Ptr<pat::Electron> & electron, const int & elecIDLevel, const std::vector<edm::Handle<edm::ValueMap<unsigned int>>> & electron_bitmaps, const std::vector<edm::Handle<edm::ValueMap<vid::CutFlowResult>>> & electron_cutflows) 
{
  auto electron_bitmap = (*(electron_bitmaps[elecIDLevel]))[electron];
  auto electron_cutflow = (*(electron_cutflows[elecIDLevel]))[electron];
  bool pass = true;
  for (unsigned iCut = 0; iCut < electron_cutflow.cutFlowSize(); ++iCut) {
    if (iCut == 7) continue; // Ignore isolation cut
    bool pass_cut = electron_bitmap>>iCut&1;
    if (!pass_cut) pass = false;
  }
  return pass;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);

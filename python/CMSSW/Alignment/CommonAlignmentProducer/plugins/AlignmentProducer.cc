/// \file AlignmentProducer.cc
///
///  \author    : Frederic Ronga
///  Revision   : $Revision: 1.43 $
///  last update: $Date: 2010/09/10 11:46:17 $
///  by         : $Author: mussgill $

#include "AlignmentProducer.h"
#include "FWCore/Framework/interface/LooperFactory.h" 
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterBuilder.h" 
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterStore.h" 
#include "Alignment/CommonAlignment/interface/Alignable.h" 

#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

// System include files
#include <memory>
#include <sstream>

// Framework
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/Run.h"

// Conditions database
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

// Geometry
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeomBuilderFromGeometricDet.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/MuonNumberingRecord.h"
#include "Geometry/DTGeometryBuilder/src/DTGeometryBuilderFromDDD.h"
#include "Geometry/CSCGeometryBuilder/src/CSCGeometryBuilderFromDDD.h"
#include "Geometry/TrackingGeometryAligner/interface/GeometryAligner.h"
#include "CondFormats/AlignmentRecord/interface/TrackerAlignmentRcd.h"
#include "CondFormats/AlignmentRecord/interface/TrackerAlignmentErrorRcd.h"
#include "CondFormats/AlignmentRecord/interface/DTAlignmentRcd.h"
#include "CondFormats/AlignmentRecord/interface/DTAlignmentErrorRcd.h"
#include "CondFormats/AlignmentRecord/interface/CSCAlignmentRcd.h"
#include "CondFormats/AlignmentRecord/interface/CSCAlignmentErrorRcd.h"
#include "CondFormats/AlignmentRecord/interface/TrackerSurveyRcd.h"
#include "CondFormats/AlignmentRecord/interface/TrackerSurveyErrorRcd.h"
#include "CondFormats/AlignmentRecord/interface/DTSurveyRcd.h"
#include "CondFormats/AlignmentRecord/interface/DTSurveyErrorRcd.h"
#include "CondFormats/AlignmentRecord/interface/CSCSurveyRcd.h"
#include "CondFormats/AlignmentRecord/interface/CSCSurveyErrorRcd.h"
#include "CondFormats/AlignmentRecord/interface/GlobalPositionRcd.h"
#include "CondFormats/Alignment/interface/DetectorGlobalPosition.h"

// Tracking, LAS and cluster flag map (fwd is enough!) 
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Alignment/interface/AliClusterValueMapFwd.h"
#include "DataFormats/Alignment/interface/TkFittedLasBeamCollectionFwd.h"
#include "Alignment/LaserAlignment/interface/TsosVectorCollection.h"

// Alignment
#include "CondFormats/Alignment/interface/SurveyErrors.h"
#include "Alignment/TrackerAlignment/interface/TrackerScenarioBuilder.h"
#include "Alignment/MuonAlignment/interface/MuonScenarioBuilder.h"
#include "Alignment/CommonAlignment/interface/SurveyDet.h"
#include "Alignment/CommonAlignmentParametrization/interface/RigidBodyAlignmentParameters.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentAlgorithmPluginFactory.h"
#include "Alignment/CommonAlignmentMonitor/interface/AlignmentMonitorPluginFactory.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterSelector.h"

//_____________________________________________________________________________
AlignmentProducer::AlignmentProducer(const edm::ParameterSet& iConfig) :
  theAlignmentAlgo(0), theAlignmentParameterStore(0),
  theAlignableExtras(0), theAlignableTracker(0), theAlignableMuon(0), 
  globalPositions_(0),
  nevent_(0), theParameterSet(iConfig),
  theMaxLoops( iConfig.getUntrackedParameter<unsigned int>("maxLoops") ),
  stNFixAlignables_(iConfig.getParameter<int>("nFixAlignables") ),
  stRandomShift_(iConfig.getParameter<double>("randomShift")),
  stRandomRotation_(iConfig.getParameter<double>("randomRotation")),
  applyDbAlignment_( iConfig.getUntrackedParameter<bool>("applyDbAlignment")),
  doMisalignmentScenario_(iConfig.getParameter<bool>("doMisalignmentScenario")),
  saveToDB_(iConfig.getParameter<bool>("saveToDB")),
  saveApeToDB_(iConfig.getParameter<bool>("saveApeToDB")),
  doTracker_( iConfig.getUntrackedParameter<bool>("doTracker") ),
  doMuon_( iConfig.getUntrackedParameter<bool>("doMuon") ),
  useExtras_( iConfig.getUntrackedParameter<bool>("useExtras") ),
  useSurvey_( iConfig.getParameter<bool>("useSurvey") ),
  tjTkAssociationMapTag_(iConfig.getParameter<edm::InputTag>("tjTkAssociationMapTag")),
  beamSpotTag_(iConfig.getParameter<edm::InputTag>("beamSpotTag")),
  tkLasBeamTag_(iConfig.getParameter<edm::InputTag>("tkLasBeamTag")),
  clusterValueMapTag_(iConfig.getParameter<edm::InputTag>("hitPrescaleMapTag"))
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::AlignmentProducer";

  // Tell the framework what data is being produced
  if (doTracker_) {
     setWhatProduced(this, &AlignmentProducer::produceTracker);
  }
  if (doMuon_) {
     setWhatProduced(this, &AlignmentProducer::produceDT);
     setWhatProduced(this, &AlignmentProducer::produceCSC);
  }

  // Create the alignment algorithm
  edm::ParameterSet algoConfig = iConfig.getParameter<edm::ParameterSet>( "algoConfig" );
  std::string algoName = algoConfig.getParameter<std::string>("algoName");
  theAlignmentAlgo = AlignmentAlgorithmPluginFactory::get( )->create( algoName, algoConfig  );

  // Check if found
  if ( !theAlignmentAlgo )
	throw cms::Exception("BadConfig") << "Couldn't find algorithm called " << algoName;

  edm::ParameterSet monitorConfig = iConfig.getParameter<edm::ParameterSet>( "monitorConfig" );
  std::vector<std::string> monitors = monitorConfig.getUntrackedParameter<std::vector<std::string> >( "monitors" );

  for (std::vector<std::string>::const_iterator miter = monitors.begin();  miter != monitors.end();  ++miter) {
    AlignmentMonitorBase* newMonitor = AlignmentMonitorPluginFactory::get()->create(*miter, monitorConfig.getUntrackedParameter<edm::ParameterSet>(*miter));

    if (!newMonitor) throw cms::Exception("BadConfig") << "Couldn't find monitor named " << *miter;

    theMonitors.push_back(newMonitor);
  }
}


//_____________________________________________________________________________
// Delete new objects
AlignmentProducer::~AlignmentProducer()
{
  delete theAlignmentParameterStore;
  delete theAlignableExtras;
  delete theAlignableTracker;
  delete theAlignableMuon;

  delete globalPositions_;
}


//_____________________________________________________________________________
// Produce tracker geometry
boost::shared_ptr<TrackerGeometry>
AlignmentProducer::produceTracker( const TrackerDigiGeometryRecord& iRecord )
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::produceTracker";
  return theTracker;
}

//_____________________________________________________________________________
// Produce muonDT geometry
boost::shared_ptr<DTGeometry>
AlignmentProducer::produceDT( const MuonGeometryRecord& iRecord )
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::produceDT";
  return theMuonDT;
}

//_____________________________________________________________________________
// Produce muonCSC geometry
boost::shared_ptr<CSCGeometry>
AlignmentProducer::produceCSC( const MuonGeometryRecord& iRecord )
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::produceCSC";
  return theMuonCSC;  
}


//_____________________________________________________________________________
// Initialize algorithm
void AlignmentProducer::beginOfJob( const edm::EventSetup& iSetup )
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::beginOfJob";

  nevent_ = 0;

  // Create the geometries from the ideal geometries (first time only)
  this->createGeometries_( iSetup );
  
  // Retrieve and apply alignments, if requested (requires DB setup)
  if ( applyDbAlignment_ ) {
    // we need GlobalPositionRcd - and have to keep track for later removal
    // before writing again to DB...
    edm::ESHandle<Alignments> globalPositionRcd;
    iSetup.get<GlobalPositionRcd>().get(globalPositionRcd);
    globalPositions_ = new Alignments(*globalPositionRcd);

    if ( doTracker_ ) {     // apply to tracker
      this->applyDB<TrackerGeometry,TrackerAlignmentRcd,TrackerAlignmentErrorRcd>
	(&(*theTracker), iSetup,  
	 align::DetectorGlobalPosition(*globalPositions_, DetId(DetId::Tracker)));
    }

    if ( doMuon_ ) { // apply to tracker
      this->applyDB<DTGeometry,DTAlignmentRcd,DTAlignmentErrorRcd>
	(&(*theMuonDT), iSetup,
	 align::DetectorGlobalPosition(*globalPositions_, DetId(DetId::Muon)));
      this->applyDB<CSCGeometry,CSCAlignmentRcd,CSCAlignmentErrorRcd>
	(&(*theMuonCSC), iSetup,
	 align::DetectorGlobalPosition(*globalPositions_, DetId(DetId::Muon)));
    }
  }

  // Create alignable tracker and muon 
  if (doTracker_) {
    theAlignableTracker = new AlignableTracker( &(*theTracker) );
  }

  if (doMuon_) {
     theAlignableMuon = new AlignableMuon( &(*theMuonDT), &(*theMuonCSC) );
  }

  if (useExtras_) {
    theAlignableExtras = new AlignableExtras();
  }

  // Create alignment parameter builder
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::beginOfJob" 
                            << "Creating AlignmentParameterBuilder";
  edm::ParameterSet aliParamBuildCfg = 
    theParameterSet.getParameter<edm::ParameterSet>("ParameterBuilder");
  AlignmentParameterBuilder alignmentParameterBuilder(theAlignableTracker,
                                                      theAlignableMuon,
                                                      theAlignableExtras,
						      aliParamBuildCfg );
  // Fix alignables if requested
  if (stNFixAlignables_>0) alignmentParameterBuilder.fixAlignables(stNFixAlignables_);

  // Get list of alignables
  Alignables theAlignables = alignmentParameterBuilder.alignables();
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::beginOfJob" 
                            << "got " << theAlignables.size() << " alignables";

  // Create AlignmentParameterStore 
  edm::ParameterSet aliParamStoreCfg = 
    theParameterSet.getParameter<edm::ParameterSet>("ParameterStore");
  theAlignmentParameterStore = new AlignmentParameterStore(theAlignables, aliParamStoreCfg);
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::beginOfJob" 
                            << "AlignmentParameterStore created!";

  // Apply misalignment scenario to alignable tracker and muon if requested
  // WARNING: this assumes scenarioConfig can be passed to both muon and tracker
  if (doMisalignmentScenario_ && (doTracker_ || doMuon_)) {
    edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::beginOfJob" 
                              << "Applying misalignment scenario to "
                              << (doTracker_ ? "tracker" : "")
                              << (doMuon_    ? (doTracker_ ? " and muon" : "muon") : ".");
    edm::ParameterSet scenarioConfig 
      = theParameterSet.getParameter<edm::ParameterSet>( "MisalignmentScenario" );
    if (doTracker_) {
      TrackerScenarioBuilder scenarioBuilder( theAlignableTracker );
      scenarioBuilder.applyScenario( scenarioConfig );
    }
    if (doMuon_) {
      MuonScenarioBuilder muonScenarioBuilder( theAlignableMuon );
      muonScenarioBuilder.applyScenario( scenarioConfig );
    }
  } else {
    edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::beginOfJob" 
                              << "NOT applying misalignment scenario!";
  }

  // Apply simple misalignment
  const std::string sParSel(theParameterSet.getParameter<std::string>("parameterSelectorSimple"));
  this->simpleMisalignment_(theAlignables, sParSel, stRandomShift_, stRandomRotation_, true);

  // Initialize alignment algorithm
  theAlignmentAlgo->initialize( iSetup, 
				theAlignableTracker, theAlignableMuon, theAlignableExtras,
				theAlignmentParameterStore );

  for (std::vector<AlignmentMonitorBase*>::const_iterator monitor = theMonitors.begin();
       monitor != theMonitors.end();  ++monitor) {
     (*monitor)->beginOfJob(theAlignableTracker, theAlignableMuon, theAlignmentParameterStore);
  }
}

//_____________________________________________________________________________
// Terminate algorithm
void AlignmentProducer::endOfJob()
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::endOfJob";

  for (std::vector<AlignmentMonitorBase*>::const_iterator monitor = theMonitors.begin();  monitor != theMonitors.end();  ++monitor) {
     (*monitor)->endOfJob();
  }

  // Save alignments to database
  if (saveToDB_ || saveApeToDB_) {
    if ( doTracker_ ) { // first tracker
      const AlignTransform *trackerGlobal = 0; // will be 'removed' from constants 
      if (globalPositions_) { // i.e. applied before in applyDB
	trackerGlobal = &align::DetectorGlobalPosition(*globalPositions_,
						       DetId(DetId::Tracker));
      }
      // Get alignments+errors - ownership taken over by writeDB(..), so no delete
      Alignments *alignments = theAlignableTracker->alignments();
      AlignmentErrors *alignmentErrors = theAlignableTracker->alignmentErrors();
      this->writeDB(alignments, "TrackerAlignmentRcd",
		    alignmentErrors, "TrackerAlignmentErrorRcd", trackerGlobal);
    }
    
    if ( doMuon_ ) { // now muon
      const AlignTransform *muonGlobal = 0; // will be 'removed' from constants 
      if (globalPositions_) { // i.e. applied before in applyDB
	muonGlobal = &align::DetectorGlobalPosition(*globalPositions_,
						    DetId(DetId::Muon));
      }
      // Get alignments+errors, first DT - ownership taken over by writeDB(..), so no delete
      Alignments      *alignments       = theAlignableMuon->dtAlignments();
      AlignmentErrors *alignmentErrors  = theAlignableMuon->dtAlignmentErrors();
      this->writeDB(alignments, "DTAlignmentRcd",
		    alignmentErrors, "DTAlignmentErrorRcd", muonGlobal);
      
      // Get alignments+errors, now CSC - ownership taken over by writeDB(..), so no delete
      alignments       = theAlignableMuon->cscAlignments();
      alignmentErrors  = theAlignableMuon->cscAlignmentErrors();
      this->writeDB(alignments, "CSCAlignmentRcd",
		    alignmentErrors, "CSCAlignmentErrorRcd", muonGlobal);
    }
  }

  if (theAlignableExtras) theAlignableExtras->dump();
}

//_____________________________________________________________________________
// Called at beginning of loop
void AlignmentProducer::startingNewLoop(unsigned int iLoop )
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::startingNewLoop" 
                            << "Starting loop number " << iLoop;

  theAlignmentAlgo->startNewLoop();

  for (std::vector<AlignmentMonitorBase*>::const_iterator monitor = theMonitors.begin();  monitor != theMonitors.end();  ++monitor) {
     (*monitor)->startingNewLoop();
  }

  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::startingNewLoop" 
                            << "Now physically apply alignments to  geometry...";


  // Propagate changes to reconstruction geometry (from initialisation or iteration)
  GeometryAligner aligner;
  if ( doTracker_ ) {
    std::auto_ptr<Alignments> alignments(theAlignableTracker->alignments());
    std::auto_ptr<AlignmentErrors> alignmentErrors(theAlignableTracker->alignmentErrors());
    aligner.applyAlignments<TrackerGeometry>( &(*theTracker),&(*alignments),&(*alignmentErrors), AlignTransform() ); // don't apply global a second time!
  }
  if ( doMuon_ ) {
    std::auto_ptr<Alignments>      dtAlignments(       theAlignableMuon->dtAlignments());
    std::auto_ptr<AlignmentErrors> dtAlignmentErrors(  theAlignableMuon->dtAlignmentErrors());
    std::auto_ptr<Alignments>      cscAlignments(      theAlignableMuon->cscAlignments());
    std::auto_ptr<AlignmentErrors> cscAlignmentErrors( theAlignableMuon->cscAlignmentErrors());

    aligner.applyAlignments<DTGeometry>( &(*theMuonDT), &(*dtAlignments), &(*dtAlignmentErrors), AlignTransform() ); // don't apply global a second time!
    aligner.applyAlignments<CSCGeometry>( &(*theMuonCSC), &(*cscAlignments), &(*cscAlignmentErrors), AlignTransform() ); // nope!
  }
}


//_____________________________________________________________________________
// Called at end of loop
edm::EDLooper::Status 
AlignmentProducer::endOfLoop(const edm::EventSetup& iSetup, unsigned int iLoop)
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::endOfLoop" 
                            << "Ending loop " << iLoop;

  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::endOfLoop" 
                            << "Terminating algorithm.";
  theAlignmentAlgo->terminate();

  for (std::vector<AlignmentMonitorBase*>::const_iterator monitor = theMonitors.begin();  monitor != theMonitors.end();  ++monitor) {
     (*monitor)->endOfLoop(iSetup);
  }

  if ( iLoop == theMaxLoops-1 || iLoop >= theMaxLoops ) return kStop;
  else return kContinue;
}

//_____________________________________________________________________________
// Called at each event
edm::EDLooper::Status 
AlignmentProducer::duringLoop( const edm::Event& event, 
  const edm::EventSetup& setup )
{
  nevent_++;

  // reading in survey records
  this->readInSurveyRcds(setup);
	
  // Printout event number
  for ( int i=10; i<10000000; i*=10 )
    if ( nevent_<10*i && (nevent_%i)==0 )
      edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::duringLoop" 
                                << "Events processed: " << nevent_;
  
  // Retrieve trajectories and tracks from the event
  // -> merely skip if collection is empty
  edm::Handle<TrajTrackAssociationCollection> m_TrajTracksMap;
  if (event.getByLabel(tjTkAssociationMapTag_, m_TrajTracksMap)) {
    
    // Form pairs of trajectories and tracks
    ConstTrajTrackPairCollection trajTracks;
    for ( TrajTrackAssociationCollection::const_iterator iPair = m_TrajTracksMap->begin();
          iPair != m_TrajTracksMap->end(); ++iPair) {
      trajTracks.push_back( ConstTrajTrackPair( &(*(*iPair).key), &(*(*iPair).val) ) );
    }
    edm::Handle<reco::BeamSpot> beamSpot;
    event.getByLabel(beamSpotTag_, beamSpot);

    if (nevent_==1 && theAlignableExtras) {
      edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::duringLoop"
				<< "initializing AlignableBeamSpot" << std::endl;
      theAlignableExtras->initializeBeamSpot(beamSpot->x0(), beamSpot->y0(), beamSpot->z0(),
					     beamSpot->dxdz(), beamSpot->dydz());
    }

    // Run the alignment algorithm with its input
    const AliClusterValueMap *clusterValueMapPtr = 0;
    if(clusterValueMapTag_.encode().size()){//check that the input tag is not empty
      edm::Handle<AliClusterValueMap> clusterValueMap;
      event.getByLabel(clusterValueMapTag_, clusterValueMap);
      clusterValueMapPtr = &(*clusterValueMap);
    }

    const AlignmentAlgorithmBase::EventInfo eventInfo(event.id(), trajTracks, *beamSpot,
						      clusterValueMapPtr);
    theAlignmentAlgo->run(setup, eventInfo);


    for (std::vector<AlignmentMonitorBase*>::const_iterator monitor = theMonitors.begin();
	 monitor != theMonitors.end();  ++monitor) {
      (*monitor)->duringLoop(event, setup, trajTracks); // forward eventInfo?
    }
  } else {
    edm::LogError("Alignment") << "@SUB=AlignmentProducer::duringLoop" 
			       << "No track collection found: skipping event";
  }
  

  return kContinue;
}

// ----------------------------------------------------------------------------
void AlignmentProducer::beginRun(const edm::Run &run, const edm::EventSetup &setup)
{
  theAlignmentAlgo->beginRun(setup); // do not forward edm::Run...
}

// ----------------------------------------------------------------------------
void AlignmentProducer::endRun(const edm::Run &run, const edm::EventSetup &setup)
{
  // call with or without las beam info...
  typedef AlignmentAlgorithmBase::EndRunInfo EndRunInfo;
  if (tkLasBeamTag_.encode().size()) { // non-empty InputTag
    edm::Handle<TkFittedLasBeamCollection> lasBeams;
    edm::Handle<TsosVectorCollection> tsoses;
    run.getByLabel(tkLasBeamTag_, lasBeams);
    run.getByLabel(tkLasBeamTag_, tsoses);
    
    theAlignmentAlgo->endRun(EndRunInfo(run.id(), &(*lasBeams), &(*tsoses)), setup);
  } else {
    edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::endRun"
			      << "No Tk LAS beams to forward to algorithm.";
    theAlignmentAlgo->endRun(EndRunInfo(run.id(), 0, 0), setup);
  }
}

// ----------------------------------------------------------------------------
void AlignmentProducer::beginLuminosityBlock(const edm::LuminosityBlock &lumiBlock,
				    const edm::EventSetup &setup)
{
  theAlignmentAlgo->beginLuminosityBlock(setup); // do not forward edm::LuminosityBlock
}

// ----------------------------------------------------------------------------
void AlignmentProducer::endLuminosityBlock(const edm::LuminosityBlock &lumiBlock,
				  const edm::EventSetup &setup)
{
  theAlignmentAlgo->endLuminosityBlock(setup); // do not forward edm::LuminosityBlock
}

// ----------------------------------------------------------------------------

void AlignmentProducer::simpleMisalignment_(const Alignables &alivec, const std::string &selection, 
                                            float shift, float rot, bool local)
{

  std::ostringstream output; // collecting output

  if (shift > 0. || rot > 0.) {
    output << "Adding random flat shift of max size " << shift
           << " and adding random flat rotation of max size " << rot <<" to ";

    std::vector<bool> commSel(0);
    if (selection != "-1") {
      AlignmentParameterSelector aSelector(0,0); // no alignable needed here...
      const std::vector<char> cSel(aSelector.convertParamSel(selection));
      if (cSel.size() < RigidBodyAlignmentParameters::N_PARAM) {
	throw cms::Exception("BadConfig") 
	  << "[AlignmentProducer::simpleMisalignment_]\n"
	  << "Expect selection string '" << selection << "' to be at least of length " 
	  << RigidBodyAlignmentParameters::N_PARAM << " or to be '-1'.\n"
	  << "(Most probably you have to adjust the parameter 'parameterSelectorSimple'.)";
      }
      for (std::vector<char>::const_iterator cIter = cSel.begin(); cIter != cSel.end(); ++cIter) {
        commSel.push_back(*cIter == '0' ? false : true);
      }
      output << "parameters defined by (" << selection 
             << "), representing (x,y,z,alpha,beta,gamma),";
    } else {
      output << "the active parameters of each alignable,";
    }
    output << " in " << (local ? "local" : "global") << " frame.";

    for (std::vector<Alignable*>::const_iterator it = alivec.begin(); it != alivec.end(); ++it) {
      Alignable* ali=(*it);
      std::vector<bool> mysel(commSel.empty() ? ali->alignmentParameters()->selector() : commSel);
      
      if (std::abs(shift)>0.00001) {
	double s0 = 0., s1 = 0., s2 = 0.;
        if (mysel[RigidBodyAlignmentParameters::dx]) s0 = shift * double(random()%1000-500)/500.;
        if (mysel[RigidBodyAlignmentParameters::dy]) s1 = shift * double(random()%1000-500)/500.;
        if (mysel[RigidBodyAlignmentParameters::dz]) s2 = shift * double(random()%1000-500)/500.;
        
        if (local) ali->move( ali->surface().toGlobal(align::LocalVector(s0,s1,s2)) );
	else       ali->move( align::GlobalVector(s0,s1,s2) );

      //AlignmentPositionError ape(dx,dy,dz);
      //ali->addAlignmentPositionError(ape);
      }

      if (std::abs(rot)>0.00001) {
	align::EulerAngles r(3);
        if (mysel[RigidBodyAlignmentParameters::dalpha]) r(1)=rot*double(random()%1000-500)/500.;
        if (mysel[RigidBodyAlignmentParameters::dbeta])  r(2)=rot*double(random()%1000-500)/500.;
        if (mysel[RigidBodyAlignmentParameters::dgamma]) r(3)=rot*double(random()%1000-500)/500.;

        const align::RotationType mrot = align::toMatrix(r);
        if (local) ali->rotateInLocalFrame(mrot);
        else       ali->rotateInGlobalFrame(mrot);
        
      //ali->addAlignmentPositionErrorFromRotation(mrot);
      }
    } // end loop on alignables
  } else {
    output << "No simple misalignment added!";
  }
  edm::LogInfo("Alignment")  << "@SUB=AlignmentProducer::simpleMisalignment_" << output.str();
}


//__________________________________________________________________________________________________
void AlignmentProducer::createGeometries_( const edm::EventSetup& iSetup )
{
   edm::ESTransientHandle<DDCompactView> cpv;
   iSetup.get<IdealGeometryRecord>().get( cpv );

   if (doTracker_) {
     edm::ESHandle<GeometricDet> geometricDet;
     iSetup.get<IdealGeometryRecord>().get( geometricDet );
     TrackerGeomBuilderFromGeometricDet trackerBuilder;
     theTracker = boost::shared_ptr<TrackerGeometry>( trackerBuilder.build(&(*geometricDet)) );
   }

   if (doMuon_) {
     edm::ESHandle<MuonDDDConstants> mdc;
     iSetup.get<MuonNumberingRecord>().get(mdc);
     DTGeometryBuilderFromDDD DTGeometryBuilder;
     CSCGeometryBuilderFromDDD CSCGeometryBuilder;
     theMuonDT = boost::shared_ptr<DTGeometry>(new DTGeometry );
     DTGeometryBuilder.build( theMuonDT, &(*cpv), *mdc);
     theMuonCSC = boost::shared_ptr<CSCGeometry>( new CSCGeometry );
     CSCGeometryBuilder.build( theMuonCSC, &(*cpv), *mdc );
   }
}

void AlignmentProducer::addSurveyInfo_(Alignable* ali)
{
  const std::vector<Alignable*>& comp = ali->components();

  unsigned int nComp = comp.size();

  for (unsigned int i = 0; i < nComp; ++i) addSurveyInfo_(comp[i]);

  const SurveyError& error = theSurveyErrors->m_surveyErrors[theSurveyIndex];

  if ( ali->id() != error.rawId() ||
       ali->alignableObjectId() != error.structureType() )
  {
    throw cms::Exception("DatabaseError")
      << "Error reading survey info from DB. Mismatched id!";
  }

  const CLHEP::Hep3Vector&  pos = theSurveyValues->m_align[theSurveyIndex].translation();
  const CLHEP::HepRotation& rot = theSurveyValues->m_align[theSurveyIndex].rotation();

  AlignableSurface surf( align::PositionType( pos.x(), pos.y(), pos.z() ),
			 align::RotationType( rot.xx(), rot.xy(), rot.xz(),
					      rot.yx(), rot.yy(), rot.yz(),
					      rot.zx(), rot.zy(), rot.zz() ) );

  surf.setWidth( ali->surface().width() );
  surf.setLength( ali->surface().length() );

  ali->setSurvey( new SurveyDet( surf, error.matrix() ) );

  ++theSurveyIndex;
}

void AlignmentProducer::readInSurveyRcds( const edm::EventSetup& iSetup ){
	
  // Get Survey Rcds and add Survey Info
  if ( doTracker_ && useSurvey_ ){
    bool tkSurveyBool = watchTkSurveyRcd_.check(iSetup);
    bool tkSurveyErrBool = watchTkSurveyErrRcd_.check(iSetup);
    edm::LogInfo("Alignment") << "watcher tksurveyrcd: " << tkSurveyBool;
    edm::LogInfo("Alignment") << "watcher tksurveyerrrcd: " << tkSurveyErrBool;
    if ( tkSurveyBool || tkSurveyErrBool){
      
      edm::LogInfo("Alignment") << "ADDING THE SURVEY INFORMATION";
      edm::ESHandle<Alignments> surveys;
      edm::ESHandle<SurveyErrors> surveyErrors;
      
      iSetup.get<TrackerSurveyRcd>().get(surveys);
      iSetup.get<TrackerSurveyErrorRcd>().get(surveyErrors);
      
      theSurveyIndex  = 0;
      theSurveyValues = &*surveys;
      theSurveyErrors = &*surveyErrors;
      addSurveyInfo_(theAlignableTracker);
    }
  }
  
  if ( doMuon_ && useSurvey_) {
    bool DTSurveyBool = watchTkSurveyRcd_.check(iSetup);
    bool DTSurveyErrBool = watchTkSurveyErrRcd_.check(iSetup);
    bool CSCSurveyBool = watchTkSurveyRcd_.check(iSetup);
    bool CSCSurveyErrBool = watchTkSurveyErrRcd_.check(iSetup);
    
    if ( DTSurveyBool || DTSurveyErrBool || CSCSurveyBool || CSCSurveyErrBool ){
      edm::ESHandle<Alignments> dtSurveys;
      edm::ESHandle<SurveyErrors> dtSurveyErrors;
      edm::ESHandle<Alignments> cscSurveys;
      edm::ESHandle<SurveyErrors> cscSurveyErrors;
      
      iSetup.get<DTSurveyRcd>().get(dtSurveys);
      iSetup.get<DTSurveyErrorRcd>().get(dtSurveyErrors);
      iSetup.get<CSCSurveyRcd>().get(cscSurveys);
      iSetup.get<CSCSurveyErrorRcd>().get(cscSurveyErrors);
      
      theSurveyIndex  = 0;
      theSurveyValues = &*dtSurveys;
      theSurveyErrors = &*dtSurveyErrors;
      std::vector<Alignable*> barrels = theAlignableMuon->DTBarrel();
      for (std::vector<Alignable*>::const_iterator iter = barrels.begin();  iter != barrels.end();  ++iter) {
	addSurveyInfo_(*iter);
      }
      
      theSurveyIndex  = 0;
      theSurveyValues = &*cscSurveys;
      theSurveyErrors = &*cscSurveyErrors;
      std::vector<Alignable*> endcaps = theAlignableMuon->CSCEndcaps();
      for (std::vector<Alignable*>::const_iterator iter = endcaps.begin();  iter != endcaps.end();  ++iter) {
	addSurveyInfo_(*iter);
      }
    }
  }

}


//////////////////////////////////////////////////
// a templated method - but private, so not accessible from outside
// ==> does not have to be in header file
template<class G, class Rcd, class ErrRcd>
void AlignmentProducer::applyDB(G* geometry, const edm::EventSetup &iSetup,
				const AlignTransform &globalCoordinates) const
{
  // 'G' is the geometry class for that DB should be applied,
  // 'Rcd' is the record class for its Alignments 
  // 'ErrRcd' is the record class for its AlignmentErrors
  // 'globalCoordinates' are global transformation for this geometry
  edm::ESHandle<Alignments> alignments;
  iSetup.get<Rcd>().get(alignments);

  edm::ESHandle<AlignmentErrors> alignmentErrors;
  iSetup.get<ErrRcd>().get(alignmentErrors);

  GeometryAligner aligner;
  aligner.applyAlignments<G>(geometry, &(*alignments), &(*alignmentErrors),
			     globalCoordinates);
}


//////////////////////////////////////////////////
void AlignmentProducer::writeDB(Alignments *alignments,
				const std::string &alignRcd,
				AlignmentErrors *alignmentErrors,
				const std::string &errRcd,
				const AlignTransform *globalCoordinates) const
{
  Alignments * tempAlignments = alignments;
  AlignmentErrors * tempAlignmentErrors = alignmentErrors;

  // Call service
  edm::Service<cond::service::PoolDBOutputService> poolDb;
  if (!poolDb.isAvailable()) { // Die if not available
    delete tempAlignments;      // promised to take over ownership...
    delete tempAlignmentErrors; // dito
    throw cms::Exception("NotAvailable") << "PoolDBOutputService not available";
  }

  if (globalCoordinates  // happens only if (applyDbAlignment_ == true)
      && globalCoordinates->transform() != AlignTransform::Transform::Identity) {

    tempAlignments = new Alignments();            // temporary storage for
    tempAlignmentErrors = new AlignmentErrors();  // final alignments and errors

    GeometryAligner aligner;
    aligner.removeGlobalTransform(alignments, alignmentErrors,
                                  *globalCoordinates,
                                  tempAlignments, tempAlignmentErrors);
    
    delete alignments;       // have to delete original alignments
    delete alignmentErrors;  // same thing for the errors

    edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::writeDB"
			      << "globalCoordinates removed from alignments (" << alignRcd
			      << ") and errors (" << alignRcd << ").";
  }
  
  if (saveToDB_) {
    edm::LogInfo("Alignment") << "Writing Alignments to " << alignRcd << ".";
    poolDb->writeOne<Alignments>(tempAlignments, poolDb->beginOfTime(), alignRcd);
  } else { // poolDb->writeOne(..) takes over 'alignments' ownership,...
    delete tempAlignments; // ...otherwise we have to delete, as promised!
  }

  if (saveApeToDB_) {
    edm::LogInfo("Alignment") << "Writing AlignmentErrors to " << errRcd << ".";
    poolDb->writeOne<AlignmentErrors>(tempAlignmentErrors, poolDb->beginOfTime(), errRcd);
  } else { // poolDb->writeOne(..) takes over 'alignmentErrors' ownership,...
    delete tempAlignmentErrors; // ...otherwise we have to delete, as promised!
  }
}


DEFINE_FWK_LOOPER( AlignmentProducer );

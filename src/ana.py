"""
GaudiPython script for eta -> mu mu (gamma) analysis
Author: Michael Peters
Usage: lb-run DaVinci/v45r8 ipython src/ana.py

This analysis ntupling script is built for Run 2 MC and Turbo data, both for
local samples and analysis production with LHCb. To be even more specific, the
script is tailored towards 2018 files, which was the last year of Run 2 data.

Note that local samples require already being installed (in ntuples/).

You can find DaVinci configs using
`lb-dirac dirac-bookkeeping-production-information DATA_ID`, e.g., `00169948`
for eta -> mu mu gamma
"""

# DaVinci configuration.
from Configurables import DaVinci
from GaudiConf import IOHelper

# Analyisis Production configuration.
try: from ProdConf import ProdConf
except ImportError: pass
from pathlib import Path
import argparse
import sys
import importlib.machinery
import importlib.util
from datetime import datetime
import os

# =============================================================================

def parseArgs() -> bool:
    """
    Parser method required for GaudiPython to work with AnalysisProductions.

    Argument parser method required for GaudiPython to work with
    AnalysisProductions because the option files get unsorted.
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__,
    )
    parser.add_argument(
        "options_paths",
        nargs="+",
        type=Path,
        help="Additional options files to load before starting Gaudi Python",
    )
    args = parser.parse_args()

    # Options file loading logic
    opt = None
    for options_path in args.options_paths:
        print("Adding options file:", options_path)
        loader = importlib.machinery.SourceFileLoader(
            "option", str(options_path.resolve()))
        spec = importlib.util.spec_from_loader(loader.name, loader)
        mod = importlib.util.module_from_spec(spec)
        loader.exec_module(mod)
        try: opt = mod.backwards
        except AttributeError: pass

    if opt is None:
        opt = 'False'
        print("Warning: No options file set BDTTagger.Backwards." +
              "Defaulting to 'False'.")
    else: print(f"Running with BDTTagger.Backwards == '{opt}'.")

    if opt == 'True': backwards = True
    elif opt == 'False': backwards = False
    else: raise ValueError("Invalid input, must be 'True' or 'False'")

    return backwards

# =============================================================================

def parseArgsLocal():
    """
    Parser method for local samples.
    """
    parser = argparse.ArgumentParser(description='Local sample analysis script')
    parser.add_argument('evtmax', type=int, default=-1,
                        help='Maximum number of events to process. -1 for all events.')
    return parser.parse_args()

# =============================================================================

# Possible decay options
DECAYS = ['eta2mumu', 'eta2mumugamma', 'eta2mumumumu', 'eta2mumuee']

# Set flags
# True = MC | False = Turbo Run 2 data
IS_MC = False
# True = signal | False = minbias
# Only relevant if IS_MC is True
IS_SIGNAL = False
# True = local sample | False = analysis production
IS_SAMPLE = False
# Decay type
DECAY = 'eta2mumu'
if DECAY not in DECAYS: 
    raise ValueError(f"Invalid decay mode. Must be one of {DECAYS}.")

DaVinci().DataType = '2018'
DaVinci().Lumi = False  # Processing luminosity data
# Local sample
args = parseArgsLocal()
if IS_SAMPLE:
    # MC
    if IS_MC:
        DaVinci().Lumi = False  # No luminosity data for MC.
        DaVinci().Simulation = True  # MC simulation data.
        # Signal sample, i.e., at least one candidate guaranteed per event
        if IS_SIGNAL:
            if DECAY == 'eta2mumugamma':
                DaVinci().DDDBtag = 'dddb-20210528-8'
                DaVinci().CondDBtag = 'sim-20201113-8-vc-md100-Sim10'
                # event type 39112231 (MC 2018), 00169948 root files
                data_paths = [
                    'data/eta2mumugamma/00169948_00000003_7.AllStreams.dst',
                ]
            elif DECAY == 'eta2mumu':
                DaVinci().DDDBtag = '2018-v03.06'
                DaVinci().CondDBtag = 'sim-20201113-8-vc-md100-Sim10'
                # event type 39112031 (MC 2018), 00358503 root files
                data_paths = [
                    'data/eta2mumu/00358503_00000016_1.allstreams.dst'
                ]
            else:
                # TODO: add other decay modes for local sample tests
                raise ValueError("Invalid decay mode.")
        # Minbias
        else:  # 00090844 root files
            DaVinci().DDDBtag = 'dddb-20170721-3'
            DaVinci().CondDBtag = 'sim-20190128-vc-md100'
            data_paths = [
                'data/minbias/00090844_00000001_7.AllStreams.dst',
            ]
    # Sample data
    else: 
        data_paths = [
            'data/00080042_00003916_1.leptons.mdst'
        ]
    

    # Get input data.
    IOHelper('ROOT').inputFiles(data_paths, clear=True)
    # Use current date + '.root' as file suffix.
    extension = "_" + str(datetime.now().strftime("%Y%m%d")) + ".root"
# Analysis production
else:
    backwards = parseArgs()
    # Output file
    outfile = f"{ProdConf().OutputFilePrefix}.{ProdConf().OutputFileTypes[0]}"
    
# Reconstruction.
from Configurables import CombineParticles
from StandardParticles import StdAllLooseMuons as muons
from StandardParticles import StdLooseAllPhotons as photons
from StandardParticles import StdLooseElectrons as electrons
from PhysSelPython.Wrappers import Selection, SelectionSequence
# Data configuration
DaVinci().RootInTES = ''  # MC default; overridden below for Turbo data.
if not IS_MC:
    from PhysConf.Filters import LoKi_Filters
    DaVinci().Simulation = False
    DaVinci().InputType = 'MDST'
    DaVinci().RootInTES = '/Event/Leptons/Turbo'
    DaVinci().Turbo = True
    from Configurables import TurboConf
    TurboConf().RunPersistRecoUnpacking = True
    DaVinci().DDDBtag = 'dddb-20171030-3'
    DaVinci().CondDBtag = 'cond-20180202'
    hlt = LoKi_Filters(HLT2_Code =
        "HLT_PASS_RE('Hlt2ExoticaPrmptDiMuonTurboDecision') | "
        "HLT_PASS_RE('Hlt2ExoticaDisplDiMuonDecision') | "
        "HLT_PASS_RE('Hlt2ExoticaDiMuonNoIPTurboDecision')")
    DaVinci().EventPreFilters = hlt.filters('TriggerFilters')
    # For Turbo data, RebuildSelection is required so that standard particle
    # makers source their inputs from the Turbo persistent reco containers.
    from PhysConf.Selections import RebuildSelection
    muons   = RebuildSelection(muons)
    photons = RebuildSelection(photons)
    electrons = RebuildSelection(electrons)

# --- Decay mode config --------------------------------------------------------
# --- Combination cuts ---
combination_cuts = (
    "(AM > 398*MeV) & (AM < 1108*MeV) & "  # 398 = eta - 150, 1108 = eta_prime + 150
    "(AMAXDOCA('') < 0.4*mm) & "  # doca between children
    # possibly change TRCHI2DOF to 2.5
    "(AMAXCHILD('mu-' == ABSID, TRCHI2DOF) < 3) & "  # track
    "(AMINCHILD('mu-' == ABSID, PROBNNmu) > 0.4)"  # muon weights
)

# --- Daughter cuts ---
daughter_cuts = {}
required_selections = [muons]
if DECAY == 'eta2mumugamma':
    daughter_cuts["mu+"] = "(PT > 500*MeV) & (P > 3*GeV)"
    daughter_cuts["mu-"] = "(PT > 500*MeV) & (P > 3*GeV)"
    daughter_cuts["gamma"] = "(PT > 500*MeV) & (CL > 0.2)"
    required_selections.append(photons)
    if IS_SAMPLE: outfile = 'ntuples/eta2MuMuGamma' + ('_mc' if IS_MC else '') + extension
    decay_descriptor = "eta -> mu+ mu- gamma"
elif DECAY == 'eta2mumu':
    daughter_cuts["mu+"] = "(PT > 500*MeV) & (P > 3*GeV)"
    daughter_cuts["mu-"] = "(PT > 500*MeV) & (P > 3*GeV)"
    if IS_SAMPLE: outfile = 'ntuples/eta2MuMu' + ('_mc' if IS_MC else '') + extension
    decay_descriptor = "eta -> mu+ mu-"
elif DECAY == 'eta2mumumumu':
    # Half PT requirement for twice the number of final state particles
    daughter_cuts["mu+"] = "(PT > 250*MeV) & (P > 3*GeV)"
    daughter_cuts["mu-"] = "(PT > 250*MeV) & (P > 3*GeV)"
    if IS_SAMPLE: outfile = 'ntuples/eta2MuMuMuMu' + ('_mc' if IS_MC else '') + extension
    # Apply cuts to ensure at least 2 muons pass the trigger, then persist reco
    # will save the whole event and the other dimuon pair can be picked up in
    # reconstruction. If all 4 had a PT > 500 MeV cut, we would lose a lot of
    # signal unecessarily.
    combination_cuts += " & (ANUM((ABSID == 'mu+') & (PT > 500*MeV)) > 0)"
    combination_cuts += " & (ANUM((ABSID == 'mu-') & (PT > 500*MeV)) > 0)"
    decay_descriptor = "eta -> mu+ mu- mu+ mu-"
elif DECAY == 'eta2mumuee':
    # TODO: need to apply eta_4mu logic here as well
    # TODO: need special trigger lines for electron
    # Half PT requirement for twice the number of final state particles
    daughter_cuts["mu+"] = "(PT > 250*MeV) & (P > 3*GeV)"
    daughter_cuts["mu-"] = "(PT > 250*MeV) & (P > 3*GeV)"
    daughter_cuts["e+"] = "(PT > 250*MeV) & (CL > 0.2)"
    daughter_cuts["e-"] = "(PT > 250*MeV) & (CL > 0.2)"
    required_selections.append(electrons)
    if IS_SAMPLE: outfile = 'ntuples/eta2MuMuEE' + ('_mc' if IS_MC else '') + extension
    decay_descriptor = "eta -> mu+ mu- e+ e-"

print(f"Writing output to {outfile}")  # debug

# --- Apply cuts ---
comb = CombineParticles(
    'comb',
    DecayDescriptor=decay_descriptor,
    DaughtersCuts=daughter_cuts,
    CombinationCut=combination_cuts,
    # vertex fit can fail, must pass
    # vertex, hard to model this cut eff from data
    MotherCut="(HASVERTEX) & (VFASPF(VCHI2PDOF) < 10)")

# Selection
sel_comb = Selection(
    'sel',
    Algorithm=comb,
    RequiredSelections=required_selections)

# Final selection sequence
seq = SelectionSequence('seq', TopSelection=sel_comb)

# DaVinci algorithm sequence.
DaVinci().appendToMainSequence([seq])

# TisTos configuration.
from Configurables import ToolSvc, TriggerTisTos
for stage in ('Hlt1', 'Hlt2'):
    ToolSvc().addTool(TriggerTisTos, stage + "TriggerTisTos")
    tool = getattr(ToolSvc(), stage + "TriggerTisTos")
    tool.HltDecReportsLocation = os.path.join(DaVinci().RootInTES, stage, 'DecReports')
    tool.HltSelReportsLocation = os.path.join(DaVinci().RootInTES, stage, 'SelReports')

# GaudiPython configuration.
import GaudiPython
import ROOT
gaudi = GaudiPython.AppMgr()
tes = gaudi.evtsvc()  # Transient Event Storage

# Tools.
rndTool = ROOT.TRandom3(0)
genTool = gaudi.toolsvc().create(
    'DaVinciSmartAssociator',
    interface='IParticle2MCWeightedAssociator')
rftTool = gaudi.toolsvc().create(
    'PVOfflineTool',
    interface='IPVOfflineTool')
pvrTool = gaudi.toolsvc().create(
    'GenericParticle2PVRelator<_p2PVWithIPChi2, '
    'OfflineDistanceCalculatorName>/P2PVWithIPChi2',
    interface='IRelatedPVFinder')
dstTool = gaudi.toolsvc().create(
    'LoKi::TrgDistanceCalculator',
    interface='IDistanceCalculator')
trkTool = gaudi.toolsvc().create(
    'TrackMasterExtrapolator',
    interface='ITrackExtrapolator')
l0Tool = gaudi.toolsvc().create(
    'L0TriggerTisTos',
    interface='ITriggerTisTos')
hlt1Tool = gaudi.toolsvc().create(
    'TriggerTisTos/Hlt1TriggerTisTos',
    interface='ITriggerTisTos')
hlt2Tool = gaudi.toolsvc().create(
    'TriggerTisTos/Hlt2TriggerTisTos',
    interface='ITriggerTisTos')
physTool = gaudi.toolsvc().create(
    'TriggerTisTos/Strip/PhysTriggerTisTos',
    interface='ITriggerTisTos')
docaTool = GaudiPython.gbl.LoKi.Particles.DOCA(0, 0, dstTool)

# Initialize the tuple.
try: from ntuple import Ntuple
except: from src.ntuple import Ntuple
ntuple = Ntuple(outfile, IS_MC, DECAY, tes, genTool, rftTool, pvrTool,
                None, dstTool, None, trkTool, l0Tool, hlt1Tool, hlt2Tool)

# Run.
# Local sample
try: evtmax = args.evtmax if args.evtmax > 0 else float("inf")
# Analysis production
except: 
    try: evtmax = int(sys.argv[1])
    except: evtmax = float("inf")
evtnum = 0
while evtnum < evtmax:
    gaudi.run(1)  # Advance Gaudi by one event
    if not bool(tes['/Event']): break  # Exit if no data found in TES
    evtnum += 1
    ntuple.clear()

    # Fill event info.
    daq = tes['DAQ/ODIN']
    try:
        ntuple.ntuple['run_n'][0] = daq.runNumber()
        ntuple.ntuple['evt_n'][0] = daq.eventNumber()
        ntuple.ntuple['evt_tck'][0] = daq.triggerConfigurationKey()
    except: continue
    # Save number of primary vertices
    try: ntuple.ntuple['pvr_n'][0] = len(tes[os.path.join(DaVinci().RootInTES, 'Rec/Vertex/Primary')])
    except: pass
    # Scintilator pad multiplicity info from L0DUReport
    try: ntuple.ntuple['evt_spd'][0] = GaudiPython.gbl.LoKi.L0.DataValue('Spd(Mult)')(tes['Trig/L0/L0DUReport'])
    except: pass

    # Create tools.
    if not ntuple.detTool: ntuple.detTool = gaudi.detSvc()[
        '/dd/Structure/LHCb/BeforeMagnetRegion/Velo']

    # Fill candidates.
    fill = False

    # Fill MC.
    if IS_MC:
        mcps = tes['MC/Particles']
        try: len(mcps); run = True
        except: run = False
        if run:
            for mcp in mcps:
                if abs(mcp.particleID().pid()) == 221:
                    ntuple.fillMcp(mcp)
                    fill = True

    # Get particles and primary vertices
    # 20260407: removed trks
    prts = tes[os.path.join(DaVinci().RootInTES, seq.outputLocation())]
    pvrs = tes[os.path.join(DaVinci().RootInTES, 'Rec/Vertex/Primary')]

    # Fill tag and prt info.
    sigs = []
    try: len(prts); run = True
    except: run = False

    if run:
        for prt in prts:
            sigs += [prt]
            ntuple.fillPrt(prt, pvrs)
            fill = True

    # Fill ntuple if there is information for this event
    if fill: ntuple.fill()  # Debug

# Close and write the output.
ntuple.close()

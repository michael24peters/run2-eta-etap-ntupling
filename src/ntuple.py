################################################################################
# Ntuple class for eta analyses                                                #
# Author: Michael Peters                                                       #
################################################################################

import ROOT
import array
import os
import GaudiPython
from GaudiPython.Bindings import gbl
from collections import OrderedDict
import array
from LoKiArrayFunctors.decorators import AMAXDOCA  # For DOCA info
from Configurables import DaVinci  # For DaVinci().RootInTES

STD = gbl.std
LHCB = gbl.LHCb

# Tag configuration.
TrkCats = [('ve', 1), ('tt', 2), ('it', 3), ('ot', 4), ('mu', 7)]
l0Trgs = [
    'L0DiMuonDecision',
    'L0MuonDecision',
    'L0HadronDecision'  # TIS line
]
hlt1Trgs = [
    'Hlt1DiMuonNoIPDecision',
    'Hlt1DiMuonLowMassDecision',
    'Hlt1MultiDiMuonNoIPDecision',
    'Hlt1.*TrackMVA.*'  # TIS line
]
hlt2Trgs = [
    'Hlt2ExoticaPrmptDiMuonTurbo',
    'Hlt2ExoticaDiMuonNoIPTurbo',
]

# =============================================================================

class Ntuple:
    """
    Class to store an ntuple.
    """

    def __init__(self, name, is_mc, decay, tes, genTool, rftTool, pvrTool, velTool,
                 dstTool, detTool, trkTool, l0Tool, hlt1Tool, hlt2Tool):
        self.IS_MC = is_mc
        self.decay = decay
        self.tes = tes
        self.genTool = genTool
        self.rftTool = rftTool
        self.pvrTool = pvrTool
        self.velTool = velTool
        self.dstTool = dstTool
        self.detTool = detTool
        self.trkTool = trkTool
        self.l0Tool = l0Tool
        self.hlt1Tool = hlt1Tool
        self.hlt2Tool = hlt2Tool
        self.saved = {}
        self.ntuple = OrderedDict()
        self.tfile = ROOT.TFile(name, 'RECREATE')
        self.ttree = ROOT.TTree('tree', 'data')
        vrsVrt = ['x', 'y', 'z', 'dx', 'dy', 'dz']
        vrsMom = ['pid', 'px', 'py', 'pz', 'e']
        vrsPrt = ['iso', 'mu', 'loose_mu', 'tight_mu', 'pnn_mu', 'pnn_pi',
                  'pnn_k', 'pnn_p', 'pnn_ghost', 'prb_ghost', 'ip', 'ip_chi2',
                  'x0', 'y0', 'z0', 't0', 'p0', 'id0', 'z1', 'id1', 'id2',
                  'id3', 'xm2', 'ym2', 'zm2']
        # NOTE not used: prt_iso, tag_ve_ns, tag_tt_ns, tag_it_ns, tag_ot_ns,
        # tag_mu_ns
        # self.vrsInit('prt_iso')
        # self.vrsInit('tag', ['ve_ns', 'tt_ns', 'it_ns', 'ot_ns', 'mu_ns'])
        vrsTag = ['m', 'dtf_m', 'dtf_dm', 'ip', 'ip_chi2', 'fd', 'fd_chi2',
                  'doca', 'dtf_chi2', 'chi2']
        vrsTrg = (
            ['l0_tos%i' % i for i in range(len(l0Trgs))] +
            ['l0_tis%i' % i for i in range(len(l0Trgs))] +
            ['hlt1_tos%i' % i for i in range(len(hlt1Trgs))] +
            ['hlt1_tis%i' % i for i in range(len(hlt1Trgs))] +
            ['hlt2_tos%i' % i for i in range(len(hlt2Trgs))] +
            ['hlt2_tis%i' % i for i in range(len(hlt2Trgs))] +
            ['hlt2_tos_topo',  'hlt2_tis_topo']
        )
        self.vrsInit('pvr', vrsVrt)
        self.vrsInit('tag', ['idx_pvr'] + vrsMom + vrsVrt + vrsTag + vrsTrg)
        # NOTE not used: tag_ve_iso0, tag_ve_iso1, tag_ln_iso0, tag_ln_iso1,
        # self.vrsInit('tag', ['ve_iso0', 've_iso1', 'ln_iso0', 'ln_iso1'])
        self.vrsInit('prt', ['idx_pvr', 'deltar', 'idx_mom'] + vrsMom + vrsPrt)

        # MC data.
        if self.IS_MC:
            # Reco to gen-level mapping
            self.vrsInit('prt', ['idx_gen'] + vrsMom + vrsPrt)
            
            # Gen-level info
            vrsMcp = ['pid', 'q', 'px', 'py', 'pz', 'e', 'x', 'y', 'z']
            self.vrsInit('mcpvr', ['x', 'y', 'z'])
            self.vrsInit('mc', ['idx_pvr', 'idx_mom'] + vrsMcp)

        self.ntuple['pvr_n'] = array.array('d', [-1])
        self.ntuple['run_n'] = array.array('d', [-1])
        self.ntuple['evt_n'] = array.array('d', [-1])
        self.ntuple['evt_tck'] = array.array('d', [-1])
        self.ntuple['evt_spd'] = array.array('d', [-1])

        for key, val in self.ntuple.items():
            if type(val) is array.array:
                self.ttree.Branch(key, val, key + '/D')
            else: self.ttree.Branch(key, val)

    # ---------------------------------------------------------------------------

    def vrsInit(self, pre, vrs):
        for v in vrs: self.ntuple['%s_%s' % (pre, v)] = ROOT.vector('double')()

    # ---------------------------------------------------------------------------

    def key(self, obj):
        """
        Generate the key for an object.
        Return key for particle based on its physical properties, in order of
        momentum, track, calorimeter.
        obj : Object, usually of type Particle or MCParticle
        """

        key = None
        try:
            key = (obj.momentum().Px(), obj.momentum().Py(),
                   obj.momentum().Pz())
            try:
                trk = obj.proto().track()
                key = (trk.momentum().X(), trk.momentum().Y(),
                       trk.momentum().Z())
            except:
                try:
                    pos = obj.proto().calo()[0].position()
                    key = (pos.x(), pos.y(), pos.z(), pos.e())
                except: pass
        except:
            try: key = (obj.position().X(), obj.position().Y(),
                        obj.position().Z())
            except: pass
        return key

    # ---------------------------------------------------------------------------

    def close(self):
        """
        Close the ntuple.
        """

        self.tfile.Write('data', ROOT.TObject.kOverwrite)
        self.tfile.Close()

    # ---------------------------------------------------------------------------

    def clear(self):
        """
        Clear the ntuple.
        """

        self.saved.clear()
        for key, val in self.ntuple.items():
            if isinstance(val, array.array): val[0] = -1
            else: val.clear()

    # ---------------------------------------------------------------------------

    def is_event_empty(self):
        """
        Check if there is any data in the current event for the ntuple to
        fill. Return True if empty, False otherwise.
        """
        for key, val in self.ntuple.items():
            if hasattr(val, 'size') and val.size() > 0: return False
        return True

    # ---------------------------------------------------------------------------

    def fill(self, key=None, val=None, idx=None, vrs=None):
        """
        Fill the ntuple for either an event or an object.
        key : float, creates identifier for particle from its momentum
        val : value to key-value pair, usually the particle information
        idx : index of the particle (?)
        vrs : not in use
        """

        if key is None or val is None: 
            # Do not fill empty events
            if self.is_event_empty(): return
            self.tfile.Cd('')
            self.ttree.Fill()
        elif key in self.ntuple:
            if idx is None: self.ntuple[key].push_back(val)
            elif idx < len(self.ntuple[key]): self.ntuple[key][idx] = val

    # ---------------------------------------------------------------------------

    def fillVrt(self, pre, vrt, cov=None, pos=None):
        """
        Fill an object's vertex information (vrt : Vertex), including
        Covariant Matrix (cov) and Position (pos) information based on its
        ntuple prefix (pre).
        """

        from math import sqrt
        if cov is None: 
            try: cov = vrt.covMatrix()
            except: pass
        if pos is None: pos = vrt.position()
        self.fill('%s_x' % pre, pos.X())
        self.fill('%s_y' % pre, pos.Y())
        self.fill('%s_z' % pre, pos.Z())
        if cov:
            self.fill('%s_dx' % pre, sqrt(abs(cov[0][0])))
            self.fill('%s_dy' % pre, sqrt(abs(cov[1][1])))
            self.fill('%s_dz' % pre, sqrt(abs(cov[2][2])))

    # ---------------------------------------------------------------------------

    def fillMom(self, pre, mom):
        """
        Fill an object's Momentum (mom : LorentzVector) information based on
        its ntuple prefix (pre : String).
        """

        self.fill('%s_px' % pre, mom.Px())
        self.fill('%s_py' % pre, mom.Py())
        self.fill('%s_pz' % pre, mom.Pz())
        self.fill('%s_e' % pre, mom.E())
        try:
            self.fill('%s_dtf_m' % pre, mom.m().value())
            self.fill('%s_dtf_dm' % pre, mom.m().error())
        except:
            self.fill('%s_dtf_m' % pre, -1)
            self.fill('%s_dtf_dm' % pre, -1)

    # ---------------------------------------------------------------------------

    def turboTISTOS(self, obj, line, mode='tos'):
        """
        Compute TOS/TIS for Turbo HLT2 lines by comparing offline daughter track
        LHCb IDs against online candidate daughters in the Turbo container.
        For non-Turbo lines whose container is absent, falls back to the DecReport
        decision (TOS=False, TIS=True when fired).
        """
        if not obj: return False
        if   mode == 'tos': threshold = 0.7;  flag = False
        elif mode == 'tis': threshold = 0.01; flag = True
        else: return False
        try:
            reports = self.tes[os.path.join(DaVinci().RootInTES, 'Hlt2/DecReports')]
            if not reports.decReport(line + 'Decision').decision(): return False
            online_cands = self.tes[os.path.join(DaVinci().RootInTES, line, 'Particles')]
            trackPassed = [False] * len(obj.daughtersVector())
            for i, track in enumerate(obj.daughtersVector()):
                offlineIDs = set(id.lhcbID() for id in track.proto().track().lhcbIDs())
                for cand in (online_cands or []):
                    for dtr in cand.daughtersVector():
                        if dtr.particleID().pid() != track.particleID().pid(): continue
                        try: onlineIDs = set(id.lhcbID() for id in dtr.proto().track().lhcbIDs())
                        except: continue
                        if offlineIDs and len(offlineIDs & onlineIDs) / len(offlineIDs) > threshold:
                            trackPassed[i] = True
            return False if trackPassed.count(flag) else True
        except: pass

    # ---------------------------------------------------------------------------

    def fillPrt(self, prt, pvrs=None):
        """
        Fill all available particle (prt) information into ntuple, including
        its primary vertex information (pvrs).
        """

        pid = prt.particleID().pid()
        pro = prt.proto()  # ProtoParticle
        # Try to get raw Particle data
        try: prt = prt.data()
        except: pass
        # Try to get primary vertex associated with particle
        pvr_loc = os.path.join(DaVinci().RootInTES, 'Rec/Vertex/Primary')
        try: pvr = self.pvrTool.relatedPV(prt, pvr_loc)
        except: pvr = None
        # Recursive base case; check if a composite particle that decays.
        vrt = prt.endVertex()
        mom = prt.momentum()
        key = self.key(prt)  # Generate unique key for particle.
        pre = 'tag' if vrt else 'prt'

        # Save current index of TTree array for daughters to reference.
        idx = self.ntuple['%s_px' % pre].size()
        self.saved[key] = idx

        # Daughters.
        if pre == 'tag':  # tag => candidate => has daughters
            trks = []
            hits = [0] * 13  # Default (to prevent segfault)
            for dtr in prt.daughters():
                # Recursively loop thru daughters
                (dtrPre, dtrIdx) = self.fillPrt(dtr)

                # Mother (tag) indexing for daughters (prt)
                try: self.fill('%s_idx_mom' % dtrPre, idx
                               if idx is not None else -1)
                except: pass

                # Extract daughter track if able
                if dtr.proto() and dtr.proto().track():
                    trk, mu = dtr.proto().track(), None
                    # Store each hit in each sub-detector in a dictionary
                    if trk: hits = self.hits(trk)
                # 7 == CALO. If dtr has no hits in calo but does in muon
                # chambers, add the muon track. (Muons should skip the calo.)
                if hits[7] == 0 and dtr.proto().muonPID():
                    mu = dtr.proto().muonPID().muonTrack()
                    trks += [(trk, mu)]

        if pre == 'tag':
            # Find shared hits.
            # n = self.share(trks)
            # Maximum Distance of Closest Approach (DOCA) among daughters
            # from LoKiArrayFunctors.decorators import AMAXDOCA
            # Use PV as constraint
            if pvr:
                dtf = GaudiPython.gbl.DecayTreeFitter.Fitter(prt, pvr)
                dtf.fit()
            else: dtf = None
            if dtf and dtf.status() == 0:
                par = dtf.fitParams(prt)  # Vertex information
                mom = par.momentum()  # Fitted momentum
                # Save fitted values, including uncertainty
                self.fill('%s_dtf_chi2' % pre, dtf.chiSquare())
                self.fillVrt(pre, prt, par.posCovMatrix(), par.position())
            # Use original vertex info from Particle if no dtf.
            else:
                self.fill('%s_dtf_chi2' % pre, -1)
                self.fillVrt(pre, vrt)
            self.fill('%s_chi2' % pre, vrt.chi2())
            self.fill('%s_doca' % pre, AMAXDOCA('')(prt.daughters()))

        # Momentum and mass.
        self.fill('%s_m' % pre, prt.measuredMass())
        self.fillMom(pre, mom)

        # Trigger.
        if pre == 'tag':
            # TOS == Trigger On Signal, 
            # TIS == Trigger Independent of Signal,
            # TOB == Trigger On Both.
            # setOfflineInput(prt) sets candidate to analyze for the trigger tool.
            # First, clear the input for each line.
            self.l0Tool.setOfflineInput()
            self.hlt1Tool.setOfflineInput()
            self.hlt2Tool.setOfflineInput()
            # Now set trigger input for each line and fill TOS and TIS info for each line.
            self.l0Tool.setOfflineInput(prt)
            self.hlt1Tool.setOfflineInput(prt)
            self.hlt2Tool.setOfflineInput(prt)
            # Fill L0 TIS and TOS info.
            for i, name in enumerate(l0Trgs):
                self.l0Tool.setTriggerInput(name)
                self.fill('%s_l0_tos%i' % (pre, i), self.l0Tool.tisTosTobTrigger().tos())
                self.fill('%s_l0_tis%i' % (pre, i), self.l0Tool.tisTosTobTrigger().tis())
            # Fill HLT1 TIS and TOS info.
            # if pre == 'tag':
            for i, name in enumerate(hlt1Trgs):
                self.hlt1Tool.setTriggerInput(name)
                self.fill('%s_hlt1_tos%i' % (pre, i), self.hlt1Tool.tisTosTobTrigger().tos())
                self.fill('%s_hlt1_tis%i' % (pre, i), self.hlt1Tool.tisTosTobTrigger().tis())
            # Fill HLT2 TIS and TOS info.
            for i, name in enumerate(hlt2Trgs):
                if self.IS_MC:
                    self.hlt2Tool.setTriggerInput(name + 'Decision')
                    self.fill('%s_hlt2_tos%i' % (pre, i), self.hlt2Tool.tisTosTobTrigger().tos())
                    self.fill('%s_hlt2_tis%i' % (pre, i), self.hlt2Tool.tisTosTobTrigger().tis())
                else:
                    self.fill('%s_hlt2_tos%i' % (pre, i), self.turboTISTOS(prt, name))
                    self.fill('%s_hlt2_tis%i' % (pre, i), self.turboTISTOS(prt, name, mode='tis'))

            # Non-Turbo HLT2 lines persist SelReports; use TriggerTisTos.
            self.hlt2Tool.setTriggerInput('Hlt2ExoticaDisplDiMuon' + 'Decision')
            self.fill('%s_hlt2_tos_displ' % pre, self.hlt2Tool.tisTosTobTrigger().tos())
            self.fill('%s_hlt2_tis_displ' % pre, self.hlt2Tool.tisTosTobTrigger().tis())
            self.hlt2Tool.setTriggerInput('Hlt2Topo.*')
            self.fill('%s_hlt2_tos_topo' % pre, self.hlt2Tool.tisTosTobTrigger().tos())
            self.fill('%s_hlt2_tis_topo' % pre, self.hlt2Tool.tisTosTobTrigger().tis())

        # Particle ID.
        self.fill('%s_pid' % pre, pid)
        # If there is a proto track, try to fill what information is available.
        if pro:
            # Fill muon information
            try:
                # Muon hypothesis
                self.fill('%s_mu' % pre, pro.muonPID().IsMuon())
                # More inclusive muon hypothesis
                self.fill('%s_loose_mu' % pre, pro.muonPID().IsMuonLoose())
                # More exclusive muon hypothesis
                self.fill('%s_tight_mu' % pre, pro.muonPID().IsMuonTight())
            except: pass
            # Likelihood particle is muon, pion, kaon, proton, or ghost
            self.fill('%s_pnn_mu' % pre, pro.info(701, -100))
            self.fill('%s_pnn_pi' % pre, pro.info(702, -100))
            self.fill('%s_pnn_k' % pre, pro.info(703, -100))
            self.fill('%s_pnn_p' % pre, pro.info(704, -100))
            self.fill('%s_pnn_ghost' % pre, pro.info(705, -100))
            self.fill('%s_prb_ghost' % pre, pro.track().ghostProbability()
                      if pro.track() else -100)  # Info unavailable

            # Track
            trk = pro.track()
            ids = []
            # Fill VELO hits.
            if trk:
                for i in trk.lhcbIDs():
                    if i.isVelo():
                        # ?? What does detTool do? sensor()?
                        d = self.detTool.sensor(i.veloID())
                        ids += [(d.z(), d, i)]
                    ids.sort()
            # Up to 4 VELO hits allowed per track. If less found, fill with -1
            for hit in range(0, 4):
                if hit >= len(ids):
                    self.fill('%s_x%i' % (pre, hit), -1)
                    self.fill('%s_y%i' % (pre, hit), -1)
                    self.fill('%s_z%i' % (pre, hit), -1)
                    self.fill('%s_id%i' % (pre, hit), -1)
                    self.fill('%s_t%i' % (pre, hit), -1)
                    self.fill('%s_p%i' % (pre, hit), -1)
                    continue
                # ?? ADDITIONAL EXPLANATION NEEDED ??
                # Calculate track parameters based on z position
                z, d, i = ids[hit]
                s = i.veloID()
                v = GaudiPython.gbl.LHCb.StateVector()
                self.trkTool.propagate(trk, z, v, prt.particleID())
                # Fill track coordinate vertex info
                self.fill('%s_x%i' % (pre, hit), v.x())
                self.fill('%s_y%i' % (pre, hit), v.y())
                self.fill('%s_z%i' % (pre, hit), v.z())
                # Phi type strip == -1, r type strip == 1.
                self.fill('%s_id%i' % (pre, hit), i.lhcbID() *
                          (-1 if s.isPhiType() else 1))
                # Fill spatial coordinate and strip width
                if (s.isPhiType()):
                    self.fill('%s_t%i' % (pre, hit), d.globalPhi(s.strip(), 0))
                    self.fill('%s_p%i' % (pre, hit), d.phiPitch(s.strip()))
                if (s.isRType()):
                    self.fill('%s_t%i' % (pre, hit), d.globalR(s.strip(), 0))
                    self.fill('%s_p%i' % (pre, hit), d.rPitch(s.strip()))

            # Muon extrapolation.
            z = 15270.0  # Muon system region
            # Propagated state vector
            v = GaudiPython.gbl.LHCb.StateVector()
            # Save extrapolated position
            if trk:  # prevents segfault error
                self.trkTool.propagate(trk, z, v, prt.particleID())
                self.fill('%s_xm2' % pre, v.x())
                self.fill('%s_ym2' % pre, v.y())
                self.fill('%s_zm2' % pre, v.z())

        # Find linked MC particle matches only for daughters using 
        # DaVinciSmartAssociator
        if pre == 'prt' and self.IS_MC:
            try:
                # Initial dr value.
                deltar = -1.0
                # Relate reconstructed particle to generator-level particle.
                gen = None; wgt = 0; rels = self.genTool.relatedMCPs(prt)
                # Select match with heighest weight
                for rel in rels:
                    if rel.weight() > wgt:
                        gen = rel.to()
                        wgt = rel.weight()
                if gen: (genPre, genIdx) = self.fillMcp(gen)
                # If DaVinciSmartAssociator fails (no match), use delta r
                # matching
                else:
                    from math import sqrt
                    # store delta r minimum
                    mindr = float('inf'); relp = None
                    mcps = self.tes['MC/Particles']
                    for mcp in mcps:
                        dphi = mcp.momentum().phi() - prt.momentum().phi()
                        deta = mcp.momentum().eta() - prt.momentum().eta()
                        # Calculate delta r
                        dr = sqrt(dphi**2 + deta**2)
                        # Check if this is smaller than the current 
                        # minimum delta r. If so, update mindr and 
                        # relp. Add info to ntuple for this daughter's
                        # linked MCParticle
                        if dr < mindr: mindr = dr; relp = mcp
                    if relp: 
                        (genPre, genIdx) = self.fillMcp(relp)
                        deltar = mindr
                    else: genIdx = -1
                self.fill('%s_deltar' % pre, deltar)
                self.fill('%s_idx_gen' % pre, genIdx)
            except: pass  # Could fill deltar and idx_gen with -1, won't for now

        # IP.
        from ctypes import c_double
        ip, ipChi2 = c_double(-1.0), c_double(-1.0)
        # Compute Impact Parameter (IP)
        # Measures shortest distance from particle to PV. High IP means
        # particle likely isn't prompt (PV) eta and instead displaced (SV)
        # other particle (e.g. B or D meson) producing decay that fakes an eta.
        if pvr: self.dstTool.distance(prt, pvr, ip, ipChi2)
        self.fill('%s_ip' % pre, ip.value)  # need to convert c_double to py
        self.fill('%s_ip_chi2' % pre, ipChi2.value)  # same here
        fd, fdChi2 = c_double(-1.0), c_double(-1.0)
        # Compute Flight Distance (FD)
        # Measures how far reconstructed candidate traveled from PV to SV.
        # Should be approximately zero since prompt eta decays are basically at
        # production point. High FD means desired decay products likely coming
        # from LLP (e.g. B or D meson) and eta candidate is random combo of
        # displaced tracks that happen to look right. Or it could mean the
        # tracks were mis-constructed and really should have been identified as
        # an eta candidate.
        if pvr and vrt: self.dstTool.distance(vrt, pvr, fd, fdChi2)
        self.fill('%s_fd' % pre, fd.value)  # need to convert c_double to py
        self.fill('%s_fd_chi2' % pre, fdChi2.value)  # same here

        # Isolation.
        # !! Removed for now. See DiMuon.py for implementation.

        # Primary vertex.
        if pvr:
            key = self.key(pvr)
            if key not in self.saved:
                self.saved[key] = self.ntuple['pvr_x'].size()
                self.fillVrt('pvr', pvr)
            # Save index of primary vertex
            self.fill('%s_idx_pvr' % pre, self.saved[key])
        else: self.fill('%s_idx_pvr' % pre, -1)

        # Return prefix used for labeling and index where the particle's data
        # is stored.
        return (pre, idx)

    # ---------------------------------------------------------------------------

    def hits(self, trk1, trk2=None, mu1=None, mu2=None, dct=None):
        """
        hits() and share() combined find the tracks that have common hits. If
        they have many shared hits, they may be an identical track and thus are
        thus a fake particle. This is useful for 'clone rejection'.
        """

        if dct is None: dct = {n: 0 for n in range(0, 13)}
        if not trk1: return dct
        id1s = [id1 for id1 in trk1.lhcbIDs()]
        id2s = [id2 for id2 in trk2.lhcbIDs()] if trk2 else []
        if mu1: id1s += [id1 for id1 in mu1.lhcbIDs()]
        if mu2: id2s += [id2 for id2 in mu2.lhcbIDs()]
        for id1 in id1s:
            cat1 = id1.detectorType()
            chn1 = id1.channelID()
            if not trk2: dct[cat1] += 1; dct[0] += 1
            for id2 in id2s:
                cat2 = id2.detectorType()
                chn2 = id2.channelID()
                if chn1 == chn2:
                    dct[0] += 1
                    if cat1 == cat2: dct[cat1] += 1
        return dct

    # ---------------------------------------------------------------------------

    def share(self, trks):
        ns = {n: 0 for n in range(0, 13)}
        for trk1 in trks:
            for trk2 in trks:
                if trk1[0] == trk2[0]: continue
                hits = self.hits(trk1[0], trk2[0], trk1[1], trk2[1])
                for key, val in hits.items():
                    if val > ns[key]: ns[key] = val
        for cat in TrkCats:
            n = ns[cat[1]]
            self.fill('tag_%s_ns' % (cat[0]), n)
        return ns[TrkCats[0][1]]

    # ---------------------------------------------------------------------------

    def fillMcp(self, prt):
        """
        Fill MC truth entries:
        - If prt is an eta(221): only fill the eta and exactly its three
          daughters [-13, 13, 22]. Set mc_idx_mom of daughters to the eta
          index. mc_idx_rec = -1 for the eta.
        - If prt is matched to a rec particle: fill only if it hasn't been filled
          already. If it has been filled, point mc_idx_rec to the rec particle
          index but do not fill again.
        - For other gen-level particles: do not fill.
        """
        pid = prt.particleID().pid()
        mom = prt.momentum()
        pos = None
        key = self.key(prt)
        pre = 'mc'
        if key in self.saved: return (pre, self.saved[key])

        # If eta, fill eta and its daughters
        if pid == 221:
            dtrs = []
            for vrt in prt.endVertices():
                for dtr in vrt.products():
                    # eta -> mu+ mu- gamma
                    if self.decay == 'eta2mumugamma' and dtr.particleID().pid() in [-13, 13, 22]:
                        dtrs.append(dtr)
                    # eta -> mu+ mu-
                    elif self.decay == 'eta2mumu' and dtr.particleID().pid() in [-13, 13]:
                        dtrs.append(dtr)
                    # eta -> mu+ mu- mu+ mu-
                    elif self.decay == 'eta2mumumumu' and dtr.particleID().pid() in [-13, 13]:
                        dtrs.append(dtr)
                    # eta -> mu+ mu- e+ e-
                    elif self.decay == 'eta2mumuee' and dtr.particleID().pid() in [-13, 13, -11, 11]:
                        dtrs.append(dtr)
            # Sort daughter pids to create consistent output
            dtrs = sorted(dtrs, key=lambda d: d[0].particleID().pid())
            pids = [d[0].particleID().pid() for d in dtrs]
            # Check for target decay
            # eta -> mu+ mu- gamma
            if self.decay == 'eta2mumugamma' and pids != [-13, 13, 22]:
                return (None, None)
            # eta -> mu+ mu-
            elif self.decay == 'eta2mumu' and pids != [-13, 13]:
                return (None, None)
            # eta -> mu+ mu- mu+ mu-
            elif self.decay == 'eta2mumumumu' and pids != [-13, -13, 13, 13]:
                return (None, None)
            # eta -> mu+ mu- e+ e-
            elif self.decay == 'eta2mumuee' and pids != [-13, -11, 11, 13]:
                return (None, None)
            decay = [prt] + dtrs

            # Fill.
            idx_eta = -1
            # Loop through particles in decay
            for p in decay:
                pid = p.particleID().pid()
                mom = p.momentum()
                key = self.key(p)
                
                idx = self.ntuple['%s_px' % pre].size()
                self.saved[key] = idx

                # Save eta as mother index
                if pid == 221:
                    idx_eta = idx
                    self.fill('%s_idx_mom' % pre, -1)
                # Point daughters to mother index
                else:
                    self.fill('%s_idx_mom' % pre, idx_eta)
                # Momentum
                self.fillMom(pre, mom)
                self.fill('%s_q' % pre, float(p.particleID().threeCharge()) / 3.0)
                self.fill('%s_pid' % pre, pid)
                # Vertex.
                self.fillVrt(pre, p.originVertex())
                # Primary vertex.
                pvr = p.primaryVertex()
                if pvr:
                    key = self.key(pvr)
                    if key not in self.saved:
                        self.saved[key] = self.ntuple['mcpvr_x'].size()
                        self.fill('mcpvr_x', pvr.position().X())
                        self.fill('mcpvr_y', pvr.position().Y())
                        self.fill('mcpvr_z', pvr.position().Z())
                    self.fill('%s_idx_pvr' % pre, self.saved[key])
                else: self.fill('%s_idx_pvr' % pre, -1)
            return (pre, idx_eta) # Return just the eta.
        
        # If not a potential eta candidate, just fill normally.
        # Save particle index.
        idx = self.ntuple['%s_px' % pre].size()
        self.saved[key] = idx

        # Momentum.
        self.fillMom(pre, mom)
        # PID.
        self.fill('%s_q' % pre, float(prt.particleID().threeCharge()) / 3.0)
        self.fill('%s_pid' % pre, pid)
        # Vertex.
        self.fillVrt(pre, prt.originVertex())
        # Primary vertex.
        pvr = prt.primaryVertex()
        if pvr:
            key = self.key(pvr)
            if key not in self.saved:
                self.saved[key] = self.ntuple['mcpvr_x'].size()
                self.fill('mcpvr_x', pvr.position().X())
                self.fill('mcpvr_y', pvr.position().Y())
                self.fill('mcpvr_z', pvr.position().Z())
            self.fill('%s_idx_pvr' % pre, self.saved[key])
        else: self.fill('%s_idx_pvr' % pre, -1)

        # Mother.
        self.fill('%s_idx_mom' % pre, -1)

        return (pre, idx)

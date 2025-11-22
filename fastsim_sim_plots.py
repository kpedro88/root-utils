from glob import glob
import os
import math
import sys
sys.path.insert(0, os.path.expandvars('$PWD/root-utils'))
from utils import PlotFactory, HistoSample, TreeSample, Variable,  Ratio

fastsim_aliases = {
"EcalHitsEB": "PCaloHits_fastSimProducer_EcalHitsEB_SIM.obj",
"EcalHitsEE": "PCaloHits_fastSimProducer_EcalHitsEE_SIM.obj",
"EcalHitsES": "PCaloHits_fastSimProducer_EcalHitsES_SIM.obj",
"HcalHits":   "PCaloHits_fastSimProducer_HcalHits_SIM.obj",
"HcalIeta":   "(2*((HcalHits.detId & 0x80000)>0) - 1)*((HcalHits.detId/(2^10)) & 0x1FF)",
"HcalSubdet": "(HcalHits.detId/(2^25))&7",
"MuonCSCHits": "PSimHits_MuonSimHits_MuonCSCHits_SIM.obj.momentumAtEntry()",
"MuonDTHits": "PSimHits_MuonSimHits_MuonDTHits_SIM.obj.momentumAtEntry()",
"MuonRPCHits": "PSimHits_MuonSimHits_MuonRPCHits_SIM.obj.momentumAtEntry()",
"TrackerHits": "PSimHits_fastSimProducer_TrackerHits_SIM.obj.momentumAtEntry()",
}

simvars = {
"EcalHitsEB":  ["EcalHitsEB.energy()", "EcalHitsEB",  (0, 50   ), 1, 100],
"EcalHitsEE":  ["EcalHitsEE.energy()", "EcalHitsEE",  (0, 80   ), 1, 100],
"EcalHitsES":  ["EcalHitsES.energy()", "EcalHitsES",  (0, 0.008), 1, 100],
"HcalHits":    ["HcalHits.energy()",   "HcalHits",    (0, 50   ), 1, 100],
"HcalIeta":    ["HcalIeta",            "HcalIeta",    (-41.5,41.5),   1, 83],
"HcalIetaWt":  ["HcalIeta",            "HcalIetaWt",  (-41.5,41.5),   1, 83],
"HcalSubdet":  ["HcalSubdet",          "HcalSubdet",  (0, 8    ), 1, 8],
"HFHits":      ["HcalHits.energy()",   "HFHits",      (0, 50   ), 1, 50],
"MuonCSCHits": ["MuonCSCHits.perp()",  "MuonCSCHits", (0, 50   ), 1, 100],
"MuonDTHits":  ["MuonDTHits.perp()",   "MuonDTHits",  (0, 100  ), 1, 100],
"MuonRPCHits": ["MuonRPCHits.perp()",  "MuonRPCHits", (0, 100  ), 1, 100],
"TrackerHits": ["TrackerHits.perp()",  "TrackerHits", (0, 150  ), 1, 100],
}

simwts = {
"HcalIetaWt": "HcalHits.energy()",
"HFHits": "HcalSubdet==4",
}

import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--old", type=str, required=True, help="old file for comparison")
parser.add_argument("--new", type=str, required=True, help="new file for comparison")
parser.add_argument("--out", type=str, default=os.path.join(os.getcwd(),'plots_sim'), help="output dir for plots")
parser.add_argument("--var", type=str, nargs='*', help="variable(s) to compare (empty: make all plots)")
parser.add_argument("--mode", type=str, default="plot", choices=["dryrun","plot","valid"], help="mode of operation")
parser.add_argument("--verbose", default=False, action="store_true", help="verbose printouts")
args = parser.parse_args()

import ROOT
ROOT.gROOT.SetBatch(True)

ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.gSystem.Load("libDataFormatsPatCandidates.so")
ROOT.FWLiteEnabler.enable()

ntestfiles = 0

basepath = [args.old]
altpath = [args.new]

vectorselection='1'
eventselection=''

def _noop(self):
    pass

PlotFactory._complete_dfs_and_weights = _noop

def _get_histos_ttree(self):
    for s in self.stacksamples + self.markersamples + self.linesamples:
        for v in self.variables:
            hname = self.inputpattern.replace('VARIABLE', v.title).replace('SAMPLE', str(s))
            self.histos[v + s] = ROOT.TH1D(hname,"",v.nbins,v.axisrange[0],v.axisrange[1])
            s.chain.Draw(f"{v.vartoplot}>>{hname}", v.wt, "goff")
    for iv, v in enumerate(self.variables):
        if iv == 0:
            for s in self.markersamples[::-1]:
                if len(s.title) > 0 and s.group is None:
                    if self.legend: self.legend.AddEntry(self.histos[v + s], s.title, 'pe')

            for s in self.stacksamples[::-1]:
                if len(s.title) > 0 and s.group is None:
                    if self.legend: self.legend.AddEntry(self.histos[v + s], s.title, 'f')

            for s in self.linesamples:
                if len(s.title) > 0 and s.group is None:
                    if self.legend: self.legend.AddEntry(self.histos[v + s], s.title, 'l')

PlotFactory._get_histos = _get_histos_ttree

pf = PlotFactory(
    outputpath=args.out,
    outputpattern='VARIABLE',
    outputformat=['png'],
    axes='log',

    normalize=False,
    ylabel='Entries',
    ylabelratio='New/Old',

    yaxisrangeratio=(0.50001, 1.49999),
    uoflowbins=True,

    ncolumnslegend=1,
    linewidth=5,

    text='',
    extratext='#splitline{Work in progress}{Simulation}',
)

def set_aliases(sample, adict):
    for key,val in adict.items():
        sample.chain.SetAlias(key,val)

oldsample = TreeSample(ntestfiles=ntestfiles, category='line', name='old', title='Old',
               modifyvarname=lambda varname: varname,
               tree='Events', files=basepath,
               eventselection=eventselection,
               vectorselection=vectorselection,
               color=ROOT.kBlack)
set_aliases(oldsample, fastsim_aliases)

newsample = TreeSample(ntestfiles=ntestfiles, category='marker', name='new', title='New',
               tree='Events', files=altpath,
               eventselection=eventselection,
               vectorselection=vectorselection,
               color=ROOT.kBlue)
set_aliases(newsample, fastsim_aliases)

pf.add_samples([
    oldsample,
    newsample
])

ratio_cat = 'ratio'
if args.mode=='valid': ratio_cat = 'diff'
pf.add_ratios([
    Ratio(category=ratio_cat, name='new:old'),
])

available_vars = [v for v in simvars]
if args.var: available_vars = args.var
pf.add_variables([
    Variable.fromlist(simvars[v]) for v in simvars if v in available_vars
])

# manually append weights (for TTree::Draw(...,wt,...)
for var in pf.variables:
    if var.title in simwts:
        var.wt = simwts[var.title]
    else:
        var.wt = ""

pf.process(mode=args.mode)

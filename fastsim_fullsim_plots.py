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
}

fullsim_aliases = {
"EcalHitsEB": "PCaloHits_g4SimHits_EcalHitsEB_SIM.obj",
"EcalHitsEE": "PCaloHits_g4SimHits_EcalHitsEE_SIM.obj",
"EcalHitsES": "PCaloHits_g4SimHits_EcalHitsES_SIM.obj",
"HcalHits":   "PCaloHits_g4SimHits_HcalHits_SIM.obj",
"HcalIeta":   "(2*((HcalHits.detId/(2^20)&1)>0) - 1)*((HcalHits.detId/(2^10)) & 1023)",
"HcalSubdet": "(HcalHits.detId/(2^28))&15",
}

simvars = {
"EcalHitsEB":  ["EcalHitsEB.energy()", "EcalHitsEB",  (0, 50   ), 1, 100],
"EcalHitsEE":  ["EcalHitsEE.energy()", "EcalHitsEE",  (0, 80   ), 1, 100],
"EcalHitsES":  ["EcalHitsES.energy()", "EcalHitsES",  (0, 0.008), 1, 100],
"HcalHits":    ["HcalHits.energy()",   "HcalHits",    (0, 50   ), 1, 100],
"HcalIeta":    ["HcalIeta",            "HcalIeta",    (-41.5,41.5),   1, 83],
"HcalSubdet":  ["HcalSubdet",          "HcalSubdet",  (0, 8    ), 1, 8],
"HFHits":      ["HcalHits.energy()",   "HFHits",      (0, 50   ), 1, 50],
}

simwts = {
"HcalIeta": "HcalHits.energy()",
"HFHits": "HcalSubdet==4",
}

import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--full", type=str, required=True, help="fullsim file for comparison")
parser.add_argument("--fast", type=str, required=True, help="fastsim file for comparison")
parser.add_argument("--out", type=str, default=os.path.join(os.getcwd(),'plots_sim_fullfast'), help="output dir for plots")
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

fullpath = [args.full]
fastpath = [args.fast]

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
    ylabelratio='Fast/Full',

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

fullsample = TreeSample(ntestfiles=ntestfiles, category='line', name='fullsim', title='FullSim',
               modifyvarname=lambda varname: varname,
               tree='Events', files=fullpath,
               eventselection=eventselection,
               vectorselection=vectorselection,
               color=ROOT.kBlack)
set_aliases(fullsample, fullsim_aliases)

fastsample = TreeSample(ntestfiles=ntestfiles, category='marker', name='fastsim', title='FastSim',
               tree='Events', files=fastpath,
               eventselection=eventselection,
               vectorselection=vectorselection,
               color=ROOT.kBlue)
set_aliases(fastsample, fastsim_aliases)

pf.add_samples([
    fullsample,
    fastsample
])

ratio_cat = 'ratio'
if args.mode=='valid': ratio_cat = 'diff'
pf.add_ratios([
    Ratio(category=ratio_cat, name='fastsim:fullsim'),
])

available_vars = [v for v in simvars]
if args.var: available_vars = args.var
pf.add_variables([
    Variable.fromlist(simvars[v]) for v in simvars if v in available_vars
])

# manually append weights (for TTree::Draw(...,wt,...)
for var in pf.variables:
    if var.name in simwts:
        var.wt = simwts[var.name]
    else:
        var.wt = ""

pf.process(mode=args.mode)

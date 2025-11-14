from glob import glob
import os
import math
import sys
sys.path.insert(0, os.path.expandvars('$PWD/root-utils'))
from utils import PlotFactory, HistoSample, TreeSample, Variable,  Ratio
from nanovariables import nanovariables, nanovariables_CMSSW_10_6_22, nanovariables_CMSSW_13_0_13

import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--old", type=str, required=True, help="old file for comparison")
parser.add_argument("--new", type=str, required=True, help="new file for comparison")
parser.add_argument("--out", type=str, default=os.path.join(os.getcwd(),'plots'), help="output dir for plots")
parser.add_argument("--var", type=str, nargs='*', help="variable(s) to compare (empty: make all plots)")
parser.add_argument("--mode", type=str, default="plot", choices=["dryrun","plot","valid"], help="mode of operation")
parser.add_argument("--verbose", default=False, action="store_true", help="verbose printouts")
args = parser.parse_args()

import ROOT
ROOT.gROOT.SetBatch(True)

ntestfiles = 0

basepath = [args.old]
altpath = [args.new]

vectorselection='1'

pf = PlotFactory(
    outputpath=args.out,
    outputpattern='VARIABLE',
    outputformat=['png'],

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

pf.add_samples([
    TreeSample(ntestfiles=ntestfiles, category='line', name='old', title='Old',
               modifyvarname=lambda varname: varname,
               tree='Events', files=basepath,
               eventselection='',
               vectorselection=vectorselection,
               color=ROOT.kBlack),

    TreeSample(ntestfiles=ntestfiles, category='marker', name='new', title='New',
               tree='Events', files=altpath,
               eventselection='',
               vectorselection=vectorselection,
               color=ROOT.kBlue),

])

pf.add_ratios([
    Ratio(category='ratio', name='new:old'),
])

available_vars = list(ROOT.RDataFrame('Events', basepath[0]).GetColumnNames())
if args.verbose: print('available vars', available_vars)
if args.var: available_vars = args.var
pf.add_variables([
    Variable.fromlist(nanovariables_CMSSW_13_0_13[v]) for v in nanovariables_CMSSW_13_0_13 if v in available_vars
])

pf.process(mode=args.mode)

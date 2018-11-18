#!/usr/bin/env python

import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("--dhlT", "-td", type=str, required=True)
parser.add_argument("--dhlP", "-pd", type=str, required=True)
parser.add_argument("--ecoliT", "-te", type=str, required=True)

args = parser.parse_args()

Dhl_trans = pd.read_table(args.dhlT)
Dhl_prot = pd.read_table(args.dhlP)
Ecoli_Trans = pd.read_table(args.ecoliT)

{system("sed s/" $1 "/" $2 "/ test.txt")}
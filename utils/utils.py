import argparse
import uuid
import pandas as pd
import os

import warnings
warnings.filterwarnings("ignore")
def to_bed(df,out):
	df.to_csv(out,sep="\t",header=False,index=False)

def delete_file(f):
	os.system("rm %s"%(f))

def read_bed(f,d):
	df = pd.read_csv(f,sep="\t",header=None)
	df[1] = df[1].astype(int)
	df[2] = df[2].astype(int)
	df[1] = df[1]-d
	df[1][df[1] < 0] = 0
	df[2] = df[2]+d
	return df




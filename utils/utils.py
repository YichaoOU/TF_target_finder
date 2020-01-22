import argparse
import uuid
import pandas as pd
import os
import pandas as pd
import sys
import matplotlib
matplotlib.use('agg')
import os
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
from matplotlib.colors import ListedColormap
import argparse
import datetime
import getpass
import uuid
from scipy.interpolate import interp1d
import re, string
from matplotlib_venn import venn2
import subprocess
from joblib import Parallel, delayed
import warnings
import scipy
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


def wccount(filename):
	df = pd.read_csv(filename,sep="\t",header=None)
	df['name'] = df[0]+":"+df[1].astype(str)+"-"+df[2].astype(str)
	df = df.drop_duplicates("name")
	# out = subprocess.Popen(['wc', '-l', filename],stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0]
	return df.shape[0]


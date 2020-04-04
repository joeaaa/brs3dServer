'''
function : merge obabel transformed canonical SMILES file and BRS-3D 300 matrix file by row, using shared compound id.
author：Li-Zhaoyang
---
file format:

1. BRS-3D 300 matrix
ID	BRS1	BRS2	...	BRS300
AP-616/40879343	0.359	0.419	...	0.5	0.515

2. canonical SMILES file(generated by command obabel -imol2 demo.mol2 -ocan -O demo.mol2.can)
smiles	ID
Cc1ccc(cc1)c1coc2c1cc1c(cc(=O)oc1c2C)c1ccccc1	AP-616/40879343

3. after merged:
ID	BRS1	BRS2	BRS299	BRS300	smiles
AP-616/40879343	0.359	0.419	0.5	0.515	Cc1ccc(cc1)c1coc2c1cc1c(cc(=O)oc1c2C)c1ccccc1
---
'''

import pandas as pd
import sys
import os

try:
      brs3d_file = sys.argv[1]
      obabel_can_file = sys.argv[2]

      current_localtion = os.getcwd() + "\merge_result"
      if not os.path.exists(current_localtion):
            os.makedirs(current_localtion)

      (filepath, tempfilename) = os.path.split(brs3d_file)
      (filename, extension) = os.path.splitext(tempfilename)

      merged_file = current_localtion + "\\"+ filename +'_merged.csv'

      brs3d_file = pd.read_table(brs3d_file)
      obabel_can_file = pd.read_table(obabel_can_file)

      pd.merge(brs3d_file, obabel_can_file, on="cmpdid", how="left",sort=False).to_csv(merged_file, index=False) #how=left，就是让df1保留所有的行列数据，df2根据df1的行列进行补全。Sort：默认为True，将合并的数据进行排序，设置为False可以提高性能；
except ImportError:
      print("Make sure the pandas module in Python Path")
except KeyError:
      print("Make sure you have added the header on both of the two files.")
      print("BRS-3D file header:")
      print("cmpdid	brs1	brs2	brs299	brs300")
      print("obabel transformed result header:")
      print("smiles	cmpdid")
      print("Header label separated by table ('\t')")
except IndexError:
      print("Usage: merge_script.py brs3d_file_with_header obabel_result_with_header")
except IOError or FileNotFoundError:
      print("File not exist, check your input")
else:
      print("BRS-3D result and obabel transformed result merged on chemical ID successfully! ")
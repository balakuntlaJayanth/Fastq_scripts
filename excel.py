import sys
import pandas as pd
import numpy as np

f = sys.argv[1]
f1 = sys.argv[2]
output = sys.argv[3]


newdf1 = pd.read_csv(f , sep='\t')
newdf2 = pd.read_csv(f1 , sep='\t')


writer = pd.ExcelWriter(output)

newdf1.to_excel(writer,'Sheet1',startcol=1,index=False)
newdf2.to_excel(writer,'Sheet2',startcol=1,index=False)
writer.save()

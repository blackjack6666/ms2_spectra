# spec_lib.parquet_file_example created by bathy at 9/19/2022

import pyarrow.parquet as pq
import pandas as pd
import json


filename = r'F:\alanine_tailing\2022_10_28\raw_data/In_gel_lane_100ng_5.parquet'

# open parquet file
dt=pq.read_table(filename)
# open the spectra index file
spec_idx=json.load(open(filename.replace('.parquet','.idx.json')))
# inspect data
df_raw_data = dt.to_pandas()
# print(df_raw_data.iloc[:300,])
# get all the valid spectral number
# spec_no_list = dt.column('spec_no').unique().to_numpy()

# locate the 5000th spectrum from the file

spectra_num = 4492

# first get the start and end from spec_idx dictionary, spec_idx dictionary use string as key
spec_start, spec_end = spec_idx[str(spectra_num)]

# then let's extract the spec data
df_ms_unique=df_raw_data.iloc[spec_start:spec_end,]
mz, intensity = df_ms_unique['mz'], df_ms_unique['int']
for i, j in zip(mz,intensity):
    print (i,j)
# visualize the data
# df_ms_5000.plot(x='mz',y='int')
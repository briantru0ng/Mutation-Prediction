import pandas as pd

df1 = pd.read_csv('sift_predictions_vep.csv')
df2 = pd.read_csv('genedataset_vep.csv')

combined = pd.concat([df1, df2]).drop_duplicates()

combined.to_csv('whole_dataset.csv', index=False)
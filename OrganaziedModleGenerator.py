import pickle
import sys
import pandas as pd
from sklearn.feature_selection import VarianceThreshold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score

def remove_low_variance(input_data, threshold=0.1):
    selection = VarianceThreshold(threshold)
    selection.fit(input_data)
    return input_data[input_data.columns[selection.get_support(indices=True)]]

def run(path):
    dataset_url = path + '/06_bioactivity_data_3class_pIC50_pubchem_fp.csv'
    dataset = pd.read_csv(dataset_url)
    X = dataset.drop(['pIC50'], axis=1)
    Y = dataset.iloc[:,-1]

    X = remove_low_variance(X, threshold=0.1)
    X.to_csv(path+'/descriptor_list.csv', index=False)
    model = RandomForestRegressor(n_estimators=500, random_state=42)
    model.fit(X, Y)
    r2 = model.score(X, Y)
    Y_pred = model.predict(X)

    # Data visualizations were skipped
    pickle.dump(model, open(path+'/modle.pkl', 'wb'))

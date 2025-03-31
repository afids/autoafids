from sklearn.metrics import mean_squared_error 
import pandas as pd
import argparse

# compute the mean squared error between 2 fcsv
def compute_mse(t1w_fcsv, t2w_fcsv):
    
    # read the data, skipping the first 2 lines (headers and metadata)
    t1w_data = pd.read_csv(t1w_fcsv, comment='#', header=None)
    t2w_data = pd.read_csv(t2w_fcsv, comment='#', header=None)
    
    # extract the columns for coordinates (x, y, z) from both files (columns 1, 2, 3)
    t1w_coords = t1w_data.iloc[:, 1:4].values
    t2w_coords = t2w_data.iloc[:, 1:4].values
    
    # make sure the files are the same shape
    if t1w_coords.shape != t2w_coords.shape:
        raise ValueError("The files have a different number of points. MSE computation requires equal data points.")
    
    # compute total mse
    mse = mean_squared_error(t1w_coords, t2w_coords)
    
    return mse

parser = argparse.ArgumentParser(description="Set paths for the computing the mse for the fcsv files.")
parser.add_argument("--autoafids_fcsv", type=str, required=True, help="Path to autoafids generated fcsv")
parser.add_argument("--baseline_fcsv", type=str, required=True, help="Path to the baseline fcsv")
args = parser.parse_args()

mse = compute_mse(args.autoafids_fcsv, args.baseline_fcsv)
print(f"mse: {mse:.4f}")

if mse > 1.2:
    raise ValueError("MSE is too large")


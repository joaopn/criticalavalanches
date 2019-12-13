%% Computes avalanche statistics from thresholded data using the NCC Toolbox

function avalanche_statistics(filepath, datafolder)

data = hdf5read(filepath,datafolder);


DATADIR="./Analysis/Data"
if [ ! -d "$DATADIR" ]; then
    echo "Creating data folders..."
    mkdir ./Analysis/Data
    mkdir ./Analysis/Data/standard
    mkdir ./Analysis/Data/standard/onebody_density
    mkdir ./Analysis/Data/standard/singlerun
    mkdir ./Analysis/Data/standard/varying_dt
    mkdir ./Analysis/Data/standard/varying_N
    mkdir ./Analysis/Data/parallel
    mkdir ./Analysis/Data/parallel/onebody_density
    mkdir ./Analysis/Data/parallel/singlerun
    mkdir ./Analysis/Data/parallel/varying_dt
    mkdir ./Analysis/Data/parallel/varying_N
else
    echo "Data folders already exist"
fi
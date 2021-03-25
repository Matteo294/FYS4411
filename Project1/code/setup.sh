DATADIR="./Analysis/Data"

if [ ! -d "$DATADIR" ]; then
    echo "Creating data folders..."
    mkdir ./Analysis/Data
    mkdir ./Analysis/Data/standard
    mkdir ./Analysis/Data/standard/onebody_density
    mkdir ./Analysis/Data/standard/singlerun
    mkdir ./Analysis/Data/standard/varying_dt
    mkdir ./Analysis/Data/standard/varying_N
    mkdir ./Analysis/Data/standard/varying_alpha
    mkdir ./Analysis/Data/parallel
    mkdir ./Analysis/Data/parallel/onebody_density
    mkdir ./Analysis/Data/parallel/singlerun
    mkdir ./Analysis/Data/parallel/varying_dt
    mkdir ./Analysis/Data/parallel/varying_N
    mkdir ./Analysis/Data/parallel/varying_alpha
else
    echo "Folders are ready"
    if [ -d "$DATADIR/parallel" ]; then

        if [ ! -d "$DATADIR/parallel/onebody_density" ]; then
            mkdir ./Analysis/Data/parallel/onebody_density
        fi

        if [ ! -d "$DATADIR/parallel/singlerun" ]; then
            mkdir ./Analysis/Data/parallel/singlerun
        fi

        if [ ! -d "$DATADIR/parallel/varying_alpha" ]; then
            mkdir ./Analysis/Data/parallel/varying_alpha
        fi

        if [ ! -d "$DATADIR/parallel/varying_dt" ]; then
            mkdir ./Analysis/Data/parallel/varying_dt
        fi

        if [ ! -d "$DATADIR/parallel/varying_N" ]; then
            mkdir ./Analysis/Data/parallel/varying_N
        fi

    else
        mkdir ./Analysis/Data/parallel
        mkdir ./Analysis/Data/parallel/onebody_density
        mkdir ./Analysis/Data/parallel/singlerun
        mkdir ./Analysis/Data/parallel/varying_dt
        mkdir ./Analysis/Data/parallel/varying_N
        mkdir ./Analysis/Data/parallel/varying_alpha
    fi

    if [ -d "$DATADIR/standard" ]; then

        if [ ! -d "$DATADIR/standard/onebody_density" ]; then
            mkdir ./Analysis/Data/standard/onebody_density
        fi

        if [ ! -d "$DATADIR/standard/singlerun" ]; then
            mkdir ./Analysis/Data/standard/singlerun
        fi

        if [ ! -d "$DATADIR/standard/varying_alpha" ]; then
            mkdir ./Analysis/Data/standard/varying_alpha
        fi

        if [ ! -d "$DATADIR/standard/varying_dt" ]; then
            mkdir ./Analysis/Data/standard/varying_dt
        fi

        if [ ! -d "$DATADIR/standard/varying_N" ]; then
            mkdir ./Analysis/Data/standard/varying_N
        fi

    else 
        mkdir ./Analysis/Data/standard
        mkdir ./Analysis/Data/standard/onebody_density
        mkdir ./Analysis/Data/standard/singlerun
        mkdir ./Analysis/Data/standard/varying_dt
        mkdir ./Analysis/Data/standard/varying_N
        mkdir ./Analysis/Data/standard/varying_alpha
    fi
    
fi
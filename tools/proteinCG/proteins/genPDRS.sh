echo $1.pdb
pdb2pqr --ff=charmm --ffout=amber --chain $1.pdb $1.pqr
freesasa -n 500 --format=pdb --hydrogen $1.pqr > $1.pdrs


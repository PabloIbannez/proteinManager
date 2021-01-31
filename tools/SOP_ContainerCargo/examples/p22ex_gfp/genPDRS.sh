#!/bin/sh

if [ "$#" -ne 1 ]; then
    echo "Syntax: ./genPDRS pdbFile "
    exit 1
fi

if [ -f $1 ]; then
    
    in=`echo "$1" | cut -f 1 -d '.'`
    echo $in
    
    ./sepMod $in
    
    for f in *_model_*
    do
        echo $f 
        id=`echo $f | grep -P -o "[0-9]+\."` #get model id from filename
        id=${id%?} #remove last character from id (it's a dot)
        n=${f%?}"rs" #change pdb->pdrs and store it in n
        pdbfixer $f --keep-heterogens=none --ph=7.0 --output=$f"H" #add hydrogen
        freesasa --format=pdb -n 500 --hydrogen $f"H" > $n
        rm $f"H"
        sed -i '1,5d' $n #remove first and second line
        sed -i "1s/^/MODEL        $id\n/" $n #add model id to the beginning of the file
    done
    
    o=$in".pdrs"
    
    echo $o
    rm -f $o
    
    for f in $(ls -1v | grep ".pdrs");do cat $f >> $o ;done

    rm *_model_*

else
   echo "File $1 does not exist."
fi

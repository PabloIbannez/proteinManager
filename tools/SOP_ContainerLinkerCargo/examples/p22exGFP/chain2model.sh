file=$1

#echo $file

first=`awk '/^ATOM *1 / {print FNR}' $file`

#echo $first

sed -i "${first}iMODEL     1" $file 

modelCount=2

endChain=`awk '/CG2 THR A 230/ {print FNR}' $file`

#echo $endChain

endChainArray=($endChain)
unset 'endChainArray[${#endChainArray[@]} -1]'

#echo ${endChainArray[2]}

for line in ${endChainArray[*]}
do
    echo $modelCount
    let "line=line+1+(modelCount-2)*2"
    sed -i "${line}iENDMDL\nMODEL     ${modelCount}" $file
    let "modelCount++"
done

sed -i '$iENDMDL' $file 



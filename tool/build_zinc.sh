#!/bin/bash

# check if we need to get zinc
echo "Do you want/need to get zinc? (y/n)"
read GET_ZINC

if [ "$GET_ZINC" = "y" ]; then
    echo "All ZINC or instock only? (a/i)"
    read ZINC_TYPE

    if [ "$ZINC_TYPE" = "a" ]; then

        echo "Downloading ZINC20"
        for i in `seq 1 20`; do wget https://files.docking.org/zinc20-ML/smiles/ZINC20_smiles_chunk_$i.tar.gz; done;

        echo "Unpacking ZINC20"
        for i in `seq 1 20`; do tar -xvf ZINC20_smiles_chunk_$i.tar.gz; done;

        echo "Removing ZINC20 tarballs"
        for i in `seq 1 20`; do rm ZINC20_smiles_chunk_$i.tar.gz; done;
    else
        echo "Downloading ZINC20 instock"
        for i in `cat instock.txt`; do wget $(echo $i | tr -d '\r\n'); done;

        echo "Removing headers"
        for i in *smi; do sed -i '1d' $i; done;
    fi

else
    echo "What directory is ZINC in?"
    read ZINC_DIR
fi

echo "Resetting Filter"
rm -rf main.bloom

echo "Building Filter"
./molbloom-bloom 2000 zinc20 1006651037 `ls ${ZINC_DIR}/ZINC20_smiles_chunk_*/smiles_all_*.txt`

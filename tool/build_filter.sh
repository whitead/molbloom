#!/bin/bash

# check if we need to get zinc
echo "Do you want/need to get zinc? (y/n)"
read GET_ZINC

if [ "$GET_ZINC" = "y" ]; then
    echo "Getting zinc..."
    echo "Downloading ZINC20"
    for i in `seq 1 20`; do wget https://files.docking.org/zinc20-ML/smiles/ZINC20_smiles_chunk_$i.tar.gz; done;

    echo "Unpacking ZINC20"
    for i in `seq 1 20`; do tar -xvf ZINC20_smiles_chunk_$i.tar.gz; done;

    echo "Removing ZINC20 tarballs"
    for i in `seq 1 20`; do rm ZINC20_smiles_chunk_$i.tar.gz; done;
else
    echo "What directory is ZINC in?"
    read ZINC_DIR
fi


echo "Resetting Filter"
rm -rf main.bloom

echo "Building Filter"
./fastz-bloom $f 2000 zinc20 1006651037 `ls ${ZINC_DIR}/ZINC20_smiles_chunk_*/smiles_all_*.txt`
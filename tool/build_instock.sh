#!/bin/bash

# check if we need to get zinc
echo "Do you want/need to get zinc? (y/n)"
read GET_ZINC

if [ "$GET_ZINC" = "y" ]; then
    echo "What directory do you want to download to?"
    read ZINC_DIR
    mkdir -p ${ZINC_DIR}
    echo "Downloading ZINC20 instock"
    for i in `cat instock.txt`; do wget ${i::-1} -P "${ZINC_DIR}" ; done;

    echo "Removing headers"
    for i in "${ZINC_DIR}/*smi"; do sed -i '1d' $i; done;
else
    echo "What directory is ZINC instock in?"
    read ZINC_DIR
fi


echo "Resetting Filter"
rm -rf instock.bloom

echo "Building Filter"
./molbloom-bloom 100 instock 9227726 `ls ${ZINC_DIR}/*smi`

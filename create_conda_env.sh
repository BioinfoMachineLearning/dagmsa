#!/bin/bash

echo "Attempting to create conda environment 'dagmsa_env' ..."
conda create --name dagmsa_env python==3.8

if [ $? != 0 ]; then
   echo "Enable to create conda environment. Quitting!!!"
   exit 1
fi
echo "conda environment 'dagmsa_env' successfully created."

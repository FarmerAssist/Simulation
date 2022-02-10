#!/bin/bash
LAT=$1
LON=$2

echo "The LAT is: $1"
echo "The LON is: $2"

Date=$(date +"%Y_%m_%d")
DIRNAME="Results/$Date"

if [ ! -d $DIRNAME ]; then
  mkdir -p $DIRNAME;
fi


Date2=$(date +"%H_%M")
DIRNAME2="$DIRNAME/${LAT}_${LON}_$Date"
echo $DIRNAME2

mkdir "/pysims/data/${DIRNAME2}"
cd "/pysims/data/${DIRNAME2}"
python /pysims/data/pysims/pysims.py --param /pysims/data/Sims/Champaign/params.apsim.sample --campaign /pysims/data/Sims/Champaign/campaign_apsim --tlatidx $LAT --tlonidx $LON

cd ..

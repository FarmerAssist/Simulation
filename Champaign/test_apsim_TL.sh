LAT="0025"
LON="0046"
Date=$(date +"%d-%m-%Y")
DIRNAME="run_${LAT}_${LON}_$Date"
echo $DIRNAME

mkdir "/pysims/data/${DIRNAME}"
cd "/pysims/data/${DIRNAME}"
python /pysims/data/pysims/pysims.py --param ../params.apsim.sample --campaign ../campaign_apsim --tlatidx $LAT --tlonidx $LON

cd ..

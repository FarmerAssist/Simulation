LAT="0025"
LON="0046"
DIRNAME="run_${LAT}_${LON}"
echo $DIRNAME
rm -r "/pysims/data/${DIRNAME}"
mkdir "/pysims/data/${DIRNAME}"
cd "/pysims/data/${DIRNAME}"
/pysims/data/pysims/pysims.py --param ../params.apsim.sample --campaign ../campaign_apsim --tlatidx $LAT --tlonidx $LON

cd ..

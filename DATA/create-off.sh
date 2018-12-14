DATADIR="non-rigid-world-DATA/"
OFFDIR="non-rigid-world-OFF/"
VERTF=$DATADIR$1".vert"
TRIF=$DATADIR$1".tri"
OFFF=$OFFDIR$1".off"

touch $OFFF
echo "OFF" > $OFFF
echo $(wc -l < $VERTF)" "$(wc -l < $TRIF)" 0" >> $OFFF

cat $VERTF >> $OFFF
echo "" >> $OFFF

awk '{for(i=1;i<=NF;i++){$i=$i-1}}1' $TRIF | awk '$0="3 "$0' >> $OFFF

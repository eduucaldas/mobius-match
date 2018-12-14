DATADIR="non-rigid-world-DATA/"
OFFDIR="non-rigid-world-OFF/"
VERTF=$DATADIR$1".vert"
TRIF=$DATADIR$1".tri"
OFFF=$OFFDIR$1".off"

touch $OFFF  
nLines=$(wc -l < $VERTF)" "$(wc -l < $TRIF) 
echo $(wc -l < $VERTF)" "$(wc -l < $TRIF) > $OFFF

cat $VERTF $TRIF >> $OFFF

#!/bin/bash

WHAT="selectZmm"

LIST=$@
[ "X$LIST" == "X"] &&  LIST="selectZmm selectZmmGen"

echo "LIST is '$LIST'"

for WHAT in selectZmm selectZmmGen;
do
   while sleep 3s;
   do
   	ps aux | grep -v grep | grep "$USER" | grep $WHAT  || 
   		{ echo "Done $WHAT" ; echo "DONE" | mail -s "done $WHAT" $USER@cern.ch ;  break;   }
   done

done

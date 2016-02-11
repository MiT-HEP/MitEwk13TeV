#!/bin/bash

WHAT="selectZmm"

LIST="$@"
[ "X$LIST" == "X" ] &&  LIST="selectZmm selectZmmGen"

echo "LIST is '$LIST'"

for WHAT in $LIST ;
do
   while sleep 3s;
   do
   	ps aux | grep -v grep | grep "$USER" | grep $WHAT >& /dev/null || 
   		{ echo "Done $WHAT" ; echo "DONE" | mail -s "done $WHAT" $USER@cern.ch ;  break;   }
   done

done

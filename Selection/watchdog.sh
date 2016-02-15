#!/bin/bash

WHAT="selectZmm"

LIST="$@"
[ "X$LIST" == "X" ] &&  LIST="selectZmm selectZmmGen"

echo "LIST is '$LIST'"

for WHAT in $LIST ;
do
   echo "-> Checking for $WHAT"
   N=0
   while sleep 3s;
   do
   	ps aux | grep -v grep | grep "$USER" | grep $WHAT >& /dev/null && N=0 || N=$((N+1))
   	[ $N -gt 3 ] && { echo "Done $WHAT" ; echo "DONE" | mail -s "done $WHAT" $USER@cern.ch ;  break;   }
   done
   sleep 1m;
done

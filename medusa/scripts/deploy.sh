#!/bin/bash

# This script deploys docs to IJS server
TMPFILE=`mktemp`
cat > $TMPFILE << CODE
lcd docs
cd /var/www/html/ParallelAndDistributedSystems/medusa/docs/
#put latex/refman.pdf
put -r html/* html/
bye
CODE
sftp -b $TMPFILE e6.ijs.si
rm $TMPFILE

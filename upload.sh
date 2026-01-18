#!/bin/sh

dir="01_literature_protein_data/"

if (cd "$dir" && xsltproc apibaso.xml > apibaso.html); then
	svn commit
else
	echo "Please edit syntax error."
fi

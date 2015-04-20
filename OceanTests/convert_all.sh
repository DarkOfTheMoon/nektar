#!/bin/bash

FILE="test2.xml";
for f in Session*.chk ; do
	FldToVtk "${FILE}" "${f}"
done
rm -rf Save
mkdir Save
mv Session* Save/

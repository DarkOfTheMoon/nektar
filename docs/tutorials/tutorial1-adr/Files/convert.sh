#!/bin/bash/
mesh_file=ADR_mesh_aligned
condition_file=ADR_mesh
NekHome=/Users/ssherw/nektarpp/
path=${NekHome}/builds/utilities/FieldConvert/FieldConvert
for ((i=0; i<=10; i++)) do
    ${path} ${mesh_file}.xml ${condition_file}.xml ${mesh_file}_${i}.chk ${mesh_file}_${i}.vtu
done

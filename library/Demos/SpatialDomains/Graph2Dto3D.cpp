#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/MeshGraph3D.h>


using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);


    if(argc < 2)
    {
        fprintf(stderr,"Usage: Graph2Dto3D meshfile   \n");
        exit(1);
    }
    try
    {
        //----------------------------------------------
        // Read in mesh from input file
        SpatialDomains::MeshGraphSharedPtr graph2D = 
            SpatialDomains::MeshGraph::Read(vSession);
        //----------------------------------------------

        //----------------------------------------------
        // Create 3D Graph from 2D mesh
        SpatialDomains::MeshGraphSharedPtr graph3D  = 
            MemoryManager<SpatialDomains::MeshGraph3D>::
            AllocateSharedPtr(graph2D,1.0);
        //----------------------------------------------

        //----------------------------------------------
        std::string out = "output.xml";
        graph3D->WriteGeometry(out);
        //----------------------------------------------
    }
    catch (const std::runtime_error&)
    {
        cout << "Caught an error" << endl;
        return 1;
    }
}


#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Communication/CommMPI.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/MeshGraph3D.h>


using namespace Nektar;
using namespace LibUtilities;

int main(int argc, char *argv[])
{


    if((argc < 2)||(argc > 3))
    {
        fprintf(stderr,"Usage: Graph2Dto3D meshfile   \n");
        exit(1);
    }
    try
    {

        //----------------------------------------------
        // Read in mesh from input file

        LibUtilities::SessionReaderSharedPtr  vSession = 
            LibUtilities::SessionReader::CreateInstance(argc, argv);

        SpatialDomains::MeshGraphSharedPtr graph2D  = 
            SpatialDomains::MeshGraph::Read(vSession);
        //----------------------------------------------

        //----------------------------------------------
        // Create 3D Graph from 2D mesh (currently using 2D session)
        SpatialDomains::MeshGraphSharedPtr graph3D  = 
            MemoryManager<SpatialDomains::MeshGraph3D>::
            AllocateSharedPtr(vSession,graph2D,1.0);
        //----------------------------------------------

        //----------------------------------------------
        // Create empty TinyXML document.
        TiXmlDocument doc;
        TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "utf-8", "");
        doc.LinkEndChild(decl);

        // Write out geometry information.
        graph3D->WriteGeometry(doc);
        //----------------------------------------------

        //----------------------------------------------
        // Save file.
        std::string infile = argv[argc-1];
        std::string outfile;
        if(vSession->GetComm()->GetSize() == 1)
        {
            outfile = infile.substr(0,infile.find_last_of(".")) + "_extrude.xml";
        }
        else
        {
            outfile = infile.substr(0,infile.find_last_of(".")) + 
                "_P"+ boost::lexical_cast<std::string>(
                          vSession->GetComm()->GetRank()) + "_extrude.xml";
        }
        
        // Add Expansion details
        TiXmlElement *root = doc.FirstChildElement("NEKTAR");
        TiXmlElement *expTag = new TiXmlElement("EXPANSIONS");
        root->LinkEndChild(expTag);
        
        SpatialDomains::CompositeMap domain = graph3D->GetDomain()[0];
        SpatialDomains::CompositeMap::iterator dit;
        for(dit = domain.begin(); dit != domain.end(); ++dit)
        {
            std::string compstr  = "C[" + boost::lexical_cast<string>(dit->first) + "]";
            TiXmlElement * e = new TiXmlElement("E");
            e->SetAttribute("COMPOSITE",compstr.c_str());
            e->SetAttribute("NUMMODES",4);
            e->SetAttribute("TYPE","MODIFIED");
            e->SetAttribute("FIELDS","u");
            expTag->LinkEndChild(e);
        }
        
        cout << "Writing file " << outfile << endl;
        doc.SaveFile(outfile);
        //----------------------------------------------

        vSession->Finalise();

    }
    catch (const std::runtime_error&)
    {
        cout << "Caught an error" << endl;
        return 1;
    }
}


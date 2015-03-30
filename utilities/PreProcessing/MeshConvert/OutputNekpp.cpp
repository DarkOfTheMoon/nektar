////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputNekpp.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Nektar++ file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
namespace io = boost::iostreams;

#include <tinyxml.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SpatialDomains/MeshGraph.h>

#include "MeshElements.h"
#include "OutputNekpp.h"

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey OutputNekpp::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eOutputModule, "xml"), OutputNekpp::create,
                "Writes a Nektar++ xml file.");

        OutputNekpp::OutputNekpp(MeshSharedPtr m) : OutputXmlBase(m)
        {
            m_config["z"] = ConfigOption(true, "0",
                "Compress output file and append a .gz extension.");
            m_config["test"] = ConfigOption(true, "0",
                "Attempt to load resulting mesh and create meshgraph.");
        }

        OutputNekpp::~OutputNekpp()
        {

        }
        
        void OutputNekpp::Process()
        {
            if (m_mesh->m_verbose)
            {
                cout << "OutputNekpp: Writing file..." << endl;
            }

            TiXmlDocument doc;
            TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "utf-8", "");
            doc.LinkEndChild( decl );

            TiXmlElement * root = new TiXmlElement( "NEKTAR" );
            doc.LinkEndChild( root );

            // Begin <GEOMETRY> section
            TiXmlElement * geomTag = new TiXmlElement( "GEOMETRY" );
            geomTag->SetAttribute("DIM", m_mesh->m_expDim);
            geomTag->SetAttribute("SPACE", m_mesh->m_spaceDim);
            root->LinkEndChild( geomTag );

            WriteXmlNodes     (geomTag);
            WriteXmlEdges     (geomTag);
            WriteXmlFaces     (geomTag);
            WriteXmlElements  (geomTag);
            WriteXmlCurves    (geomTag);
            WriteXmlComposites(geomTag);
            WriteXmlDomain    (geomTag);
            WriteXmlExpansions(root);
            WriteXmlConditions(root);
            
            // Extract the output filename and extension
            string filename = m_config["outfile"].as<string>();

            // Compress output and append .gz extension
            if (m_config["z"].as<bool>())
            {
                filename += ".gz";
                ofstream fout(filename.c_str(),
                              std::ios_base::out | std::ios_base::binary);

                std::stringstream decompressed;
                decompressed << doc;
                io::filtering_streambuf<io::output> out;
                out.push(io::gzip_compressor());
                out.push(fout);
                io::copy(decompressed, out);

                fout.close();
            }
            else
            {
                doc.SaveFile(filename);
            }

            // Test the resulting XML file by loading it with the session reader
            // and generating the meshgraph.
            if (m_config["test"].beenSet)
            {
                vector<string> filenames(1);
                filenames[0] = filename;

                LibUtilities::SessionReaderSharedPtr vSession
                    = LibUtilities::SessionReader::CreateInstance(
                        0, NULL, filenames);
                SpatialDomains::MeshGraphSharedPtr graphShPt =
                    SpatialDomains::MeshGraph::Read(vSession);
            }
        }

        void OutputNekpp::WriteXmlNodes(TiXmlElement * pRoot)
        {
            TiXmlElement* verTag = new TiXmlElement( "VERTEX" );
            std::set<NodeSharedPtr>::iterator it;

            std::set<NodeSharedPtr> tmp(
                    m_mesh->m_vertexSet.begin(),
                    m_mesh->m_vertexSet.end());

            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                NodeSharedPtr n = *it;
                stringstream s;
                s << scientific << setprecision(8) 
                  << n->m_x << " " << n->m_y << " " << n->m_z;
                TiXmlElement * v = new TiXmlElement( "V" );
                v->SetAttribute("ID",n->m_id);
                v->LinkEndChild(new TiXmlText(s.str()));
                verTag->LinkEndChild(v);
            }
            pRoot->LinkEndChild(verTag);
        }

        void OutputNekpp::WriteXmlEdges(TiXmlElement * pRoot)
        {
            if (m_mesh->m_expDim >= 2)
            {
                TiXmlElement* verTag = new TiXmlElement( "EDGE" );
                std::set<EdgeSharedPtr>::iterator it;
                std::set<EdgeSharedPtr> tmp(m_mesh->m_edgeSet.begin(),
                                            m_mesh->m_edgeSet.end());
                for (it = tmp.begin(); it != tmp.end(); ++it)
                {
                    EdgeSharedPtr ed = *it;
                    stringstream s;

                    s << setw(5) << ed->m_n1->m_id << "  " << ed->m_n2->m_id << "   ";
                    TiXmlElement * e = new TiXmlElement( "E" );
                    e->SetAttribute("ID",ed->m_id);
                    e->LinkEndChild( new TiXmlText(s.str()) );
                    verTag->LinkEndChild(e);
                }
                pRoot->LinkEndChild( verTag );
            }
        }

        void OutputNekpp::WriteXmlFaces(TiXmlElement * pRoot)
        {
            if (m_mesh->m_expDim == 3)
            {
                TiXmlElement* verTag = new TiXmlElement( "FACE" );
                std::set<FaceSharedPtr>::iterator it;
                std::set<FaceSharedPtr> tmp(
                        m_mesh->m_faceSet.begin(),
                        m_mesh->m_faceSet.end());

                for (it = tmp.begin(); it != tmp.end(); ++it)
                {
                    stringstream s;
                    FaceSharedPtr fa = *it;

                    for (int j = 0; j < fa->m_edgeList.size(); ++j)
                    {
                        s << setw(10) << fa->m_edgeList[j]->m_id;
                    }
                    TiXmlElement * f;
                    switch(fa->m_vertexList.size())
                    {
                        case 3:
                            f = new TiXmlElement("T");
                            break;
                        case 4:
                            f = new TiXmlElement("Q");
                            break;
                        default:
                            abort();
                    }
                    f->SetAttribute("ID", fa->m_id);
                    f->LinkEndChild( new TiXmlText(s.str()));
                    verTag->LinkEndChild(f);
                }
                pRoot->LinkEndChild( verTag );
            }
        }

        void OutputNekpp::WriteXmlElements(TiXmlElement * pRoot)
        {
            TiXmlElement* verTag = new TiXmlElement( "ELEMENT" );
            vector<ElementSharedPtr> &elmt = m_mesh->m_element[m_mesh->m_expDim];

            for(int i = 0; i < elmt.size(); ++i)
            {
                TiXmlElement *elm_tag = new TiXmlElement(elmt[i]->GetTag());
                elm_tag->SetAttribute("ID", elmt[i]->GetId());
                elm_tag->LinkEndChild(new TiXmlText(elmt[i]->GetXmlString()));
                verTag->LinkEndChild(elm_tag);
            }
            pRoot->LinkEndChild(verTag);
        }

    }
}

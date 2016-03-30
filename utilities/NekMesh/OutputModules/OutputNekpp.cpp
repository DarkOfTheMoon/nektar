///////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
namespace io = boost::iostreams;

#include <tinyxml.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/MeshEntities.hpp>
#include <SpatialDomains/MeshGraph.h>

#include <NekMeshUtils/MeshElements/Element.h>
#include "OutputNekpp.h"

using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{
ModuleKey OutputNekpp::className =
    GetModuleFactory().RegisterCreatorFunction(ModuleKey(eOutputModule, "xml"),
                                               OutputNekpp::create,
                                               "Writes a Nektar++ xml file.");

OutputNekpp::OutputNekpp(MeshSharedPtr m) : OutputModule(m)
{
    m_config["z"] = ConfigOption(
        true, "0", "Compress output file and append a .gz extension.");
    m_config["test"] = ConfigOption(
        true, "0", "Attempt to load resulting mesh and create meshgraph.");
    m_config["uncompress"] = ConfigOption(true, "0", "Uncompress xml sections");
}

OutputNekpp::~OutputNekpp()
{
}

template <typename T> void TestElmts(SpatialDomains::MeshGraphSharedPtr &graph)
{
    const std::map<int, boost::shared_ptr<T> > &tmp =
        graph->GetAllElementsOfType<T>();
    typename std::map<int, boost::shared_ptr<T> >::const_iterator it1, it2;

    SpatialDomains::CurveMap &curvedEdges = graph->GetCurvedEdges();
    SpatialDomains::CurveMap &curvedFaces = graph->GetCurvedFaces();

    for (it1 = tmp.begin(), it2 = tmp.end(); it1 != it2; ++it1)
    {
        SpatialDomains::GeometrySharedPtr geom = it1->second;
        geom->FillGeom();
        geom->Reset(curvedEdges, curvedFaces);
    }
}

void OutputNekpp::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "OutputNekpp: Writing file..." << endl;
    }

    TiXmlDocument doc;
    TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "utf-8", "");
    doc.LinkEndChild(decl);

    TiXmlElement *root = new TiXmlElement("NEKTAR");
    doc.LinkEndChild(root);

    // Begin <GEOMETRY> section
    TiXmlElement *geomTag = new TiXmlElement("GEOMETRY");
    geomTag->SetAttribute("DIM", m_mesh->m_expDim);
    geomTag->SetAttribute("SPACE", m_mesh->m_spaceDim);
    root->LinkEndChild(geomTag);

    WriteXmlNodes(geomTag);
    WriteXmlEdges(geomTag);
    WriteXmlFaces(geomTag);
    WriteXmlElements(geomTag);
    WriteXmlCurves(geomTag);
    WriteXmlComposites(geomTag);
    WriteXmlDomain(geomTag);
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

    // Test the resulting XML file (with a basic test) by loading it
    // with the session reader, generating the MeshGraph and testing if
    // each element is valid.
    if (m_config["test"].beenSet)
    {
        vector<string> filenames(1);
        filenames[0] = filename;

        LibUtilities::SessionReaderSharedPtr vSession =
            LibUtilities::SessionReader::CreateInstance(0, NULL, filenames);
        SpatialDomains::MeshGraphSharedPtr graphShPt =
            SpatialDomains::MeshGraph::Read(vSession);

        TestElmts<SpatialDomains::SegGeom>(graphShPt);
        TestElmts<SpatialDomains::TriGeom>(graphShPt);
        TestElmts<SpatialDomains::QuadGeom>(graphShPt);
        TestElmts<SpatialDomains::TetGeom>(graphShPt);
        TestElmts<SpatialDomains::PrismGeom>(graphShPt);
        TestElmts<SpatialDomains::PyrGeom>(graphShPt);
        TestElmts<SpatialDomains::HexGeom>(graphShPt);
    }
}

void OutputNekpp::WriteXmlNodes(TiXmlElement *pRoot)
{
    bool UnCompressed = m_config["uncompress"].as<bool>();

    TiXmlElement *verTag = new TiXmlElement("VERTEX");
    std::set<NodeSharedPtr>::iterator it;

    std::set<NodeSharedPtr> tmp(m_mesh->m_vertexSet.begin(),
                                m_mesh->m_vertexSet.end());

    if (UnCompressed)
    {
        for (it = tmp.begin(); it != tmp.end(); ++it)
        {
            NodeSharedPtr n = *it;
            stringstream s;
            s << scientific << setprecision(8) << n->m_x << " " << n->m_y << " "
              << n->m_z;
            TiXmlElement *v = new TiXmlElement("V");
            v->SetAttribute("ID", n->m_id);
            v->LinkEndChild(new TiXmlText(s.str()));
            verTag->LinkEndChild(v);
        }
    }
    else
    {
        std::vector<LibUtilities::MeshVertex> vertInfo;
        for (it = tmp.begin(); it != tmp.end(); ++it)
        {
            LibUtilities::MeshVertex v;
            NodeSharedPtr n = *it;
            v.id            = n->m_id;
            v.x             = n->m_x;
            v.y             = n->m_y;
            v.z = n->m_z;
            vertInfo.push_back(v);
        }
        std::string vertStr;
        LibUtilities::CompressData::ZlibEncodeToBase64Str(vertInfo, vertStr);
        verTag->SetAttribute("COMPRESSED",
                             LibUtilities::CompressData::GetCompressString());
        verTag->SetAttribute("BITSIZE",
                             LibUtilities::CompressData::GetBitSizeStr());

        verTag->LinkEndChild(new TiXmlText(vertStr));
    }

    pRoot->LinkEndChild(verTag);
}

void OutputNekpp::WriteXmlEdges(TiXmlElement *pRoot)
{
    bool UnCompressed = m_config["uncompress"].as<bool>();

    if (m_mesh->m_expDim >= 2)
    {
        TiXmlElement *verTag = new TiXmlElement("EDGE");

        std::set<EdgeSharedPtr>::iterator it;
        std::set<EdgeSharedPtr> tmp(m_mesh->m_edgeSet.begin(),
                                    m_mesh->m_edgeSet.end());
        if (UnCompressed)
        {
            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                EdgeSharedPtr ed = *it;
                stringstream s;

                s << setw(5) << ed->m_n1->m_id << "  " << ed->m_n2->m_id
                  << "   ";
                TiXmlElement *e = new TiXmlElement("E");
                e->SetAttribute("ID", ed->m_id);
                e->LinkEndChild(new TiXmlText(s.str()));
                verTag->LinkEndChild(e);
            }
        }
        else
        {
            std::vector<LibUtilities::MeshEdge> edgeInfo;
            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                LibUtilities::MeshEdge e;
                EdgeSharedPtr ed = *it;

                e.id = ed->m_id;
                e.v0 = ed->m_n1->m_id;
                e.v1 = ed->m_n2->m_id;

                edgeInfo.push_back(e);
            }
            std::string edgeStr;
            LibUtilities::CompressData::ZlibEncodeToBase64Str(edgeInfo,
                                                              edgeStr);
            verTag->SetAttribute(
                "COMPRESSED", LibUtilities::CompressData::GetCompressString());
            verTag->SetAttribute("BITSIZE",
                                 LibUtilities::CompressData::GetBitSizeStr());
            verTag->LinkEndChild(new TiXmlText(edgeStr));
        }
        pRoot->LinkEndChild(verTag);
    }
}

void OutputNekpp::WriteXmlFaces(TiXmlElement *pRoot)
{
    bool UnCompressed = m_config["uncompress"].as<bool>();

    if (m_mesh->m_expDim == 3)
    {
        TiXmlElement *verTag = new TiXmlElement("FACE");
        std::set<FaceSharedPtr>::iterator it;
        std::set<FaceSharedPtr> tmp(m_mesh->m_faceSet.begin(),
                                    m_mesh->m_faceSet.end());

        if (UnCompressed)
        {
            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                stringstream s;
                FaceSharedPtr fa = *it;

                for (int j = 0; j < fa->m_edgeList.size(); ++j)
                {
                    s << setw(10) << fa->m_edgeList[j]->m_id;
                }
                TiXmlElement *f;
                switch (fa->m_vertexList.size())
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
                f->LinkEndChild(new TiXmlText(s.str()));
                verTag->LinkEndChild(f);
            }
        }
        else
        {
            std::vector<LibUtilities::MeshTri> TriFaceInfo;
            std::vector<LibUtilities::MeshQuad> QuadFaceInfo;

            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                FaceSharedPtr fa = *it;

                switch (fa->m_edgeList.size())
                {
                    case 3:
                    {
                        LibUtilities::MeshTri f;
                        f.id = fa->m_id;
                        for (int i = 0; i < 3; ++i)
                        {
                            f.e[i] = fa->m_edgeList[i]->m_id;
                        }
                        TriFaceInfo.push_back(f);
                    }
                    break;
                    case 4:
                    {
                        LibUtilities::MeshQuad f;
                        f.id = fa->m_id;
                        for (int i = 0; i < 4; ++i)
                        {
                            f.e[i] = fa->m_edgeList[i]->m_id;
                        }
                        QuadFaceInfo.push_back(f);
                    }
                    break;
                    default:
                        ASSERTL0(false, "Unkonwn face type");
                }
            }

            if (TriFaceInfo.size())
            {
                std::string vType("T");
                TiXmlElement *x = new TiXmlElement(vType);
                std::string faceStr;
                LibUtilities::CompressData::ZlibEncodeToBase64Str(TriFaceInfo,
                                                                  faceStr);
                x->SetAttribute(
                    "COMPRESSED",
                    LibUtilities::CompressData::GetCompressString());
                x->SetAttribute("BITSIZE",
                                LibUtilities::CompressData::GetBitSizeStr());
                x->LinkEndChild(new TiXmlText(faceStr));
                verTag->LinkEndChild(x);
            }

            if (QuadFaceInfo.size())
            {
                std::string vType("Q");
                TiXmlElement *x = new TiXmlElement(vType);
                std::string faceStr;
                LibUtilities::CompressData::ZlibEncodeToBase64Str(QuadFaceInfo,
                                                                  faceStr);
                x->SetAttribute(
                    "COMPRESSED",
                    LibUtilities::CompressData::GetCompressString());
                x->SetAttribute("BITSIZE",
                                LibUtilities::CompressData::GetBitSizeStr());
                x->LinkEndChild(new TiXmlText(faceStr));
                verTag->LinkEndChild(x);
            }
        }
        pRoot->LinkEndChild(verTag);
    }
}

void OutputNekpp::WriteXmlElements(TiXmlElement *pRoot)
{
    bool UnCompressed = m_config["uncompress"].as<bool>();

    TiXmlElement *verTag           = new TiXmlElement("ELEMENT");
    vector<ElementSharedPtr> &elmt = m_mesh->m_element[m_mesh->m_expDim];

    if (UnCompressed)
    {
        for (int i = 0; i < elmt.size(); ++i)
        {
            TiXmlElement *elm_tag = new TiXmlElement(elmt[i]->GetTag());
            elm_tag->SetAttribute("ID", elmt[i]->GetId());
            elm_tag->LinkEndChild(new TiXmlText(elmt[i]->GetXmlString()));
            verTag->LinkEndChild(elm_tag);
        }
    }
    else
    {
        std::vector<LibUtilities::MeshEdge> SegInfo;
        std::vector<LibUtilities::MeshTri> TriInfo;
        std::vector<LibUtilities::MeshQuad> QuadInfo;
        std::vector<LibUtilities::MeshTet> TetInfo;
        std::vector<LibUtilities::MeshPyr> PyrInfo;
        std::vector<LibUtilities::MeshPrism> PrismInfo;
        std::vector<LibUtilities::MeshHex> HexInfo;

        for (int i = 0; i < elmt.size(); ++i)
        {
            switch (elmt[i]->GetTag()[0])
            {
                case 'S':
                {
                    LibUtilities::MeshEdge e;
                    e.id = elmt[i]->GetId();
                    e.v0 = elmt[i]->GetVertex(0)->m_id;
                    e.v1 = elmt[i]->GetVertex(1)->m_id;
                    SegInfo.push_back(e);
                }
                break;
                case 'T':
                {
                    LibUtilities::MeshTri e;
                    e.id = elmt[i]->GetId();
                    for (int j = 0; j < 3; ++j)
                    {
                        e.e[j] = elmt[i]->GetEdge(j)->m_id;
                    }
                    TriInfo.push_back(e);
                }
                break;
                case 'Q':
                {
                    LibUtilities::MeshQuad e;
                    e.id = elmt[i]->GetId();
                    for (int j = 0; j < 4; ++j)
                    {
                        e.e[j] = elmt[i]->GetEdge(j)->m_id;
                    }
                    QuadInfo.push_back(e);
                }
                break;
                case 'A':
                {
                    LibUtilities::MeshTet e;
                    e.id = elmt[i]->GetId();
                    for (int j = 0; j < 4; ++j)
                    {
                        e.f[j] = elmt[i]->GetFace(j)->m_id;
                    }
                    TetInfo.push_back(e);
                }
                break;
                case 'P':
                {
                    LibUtilities::MeshPyr e;
                    e.id = elmt[i]->GetId();
                    for (int j = 0; j < 5; ++j)
                    {
                        e.f[j] = elmt[i]->GetFace(j)->m_id;
                    }
                    PyrInfo.push_back(e);
                }
                break;
                case 'R':
                {
                    LibUtilities::MeshPrism e;
                    e.id = elmt[i]->GetId();
                    for (int j = 0; j < 5; ++j)
                    {
                        e.f[j] = elmt[i]->GetFace(j)->m_id;
                    }
                    PrismInfo.push_back(e);
                }
                break;
                case 'H':
                {
                    LibUtilities::MeshHex e;
                    e.id = elmt[i]->GetId();
                    for (int j = 0; j < 6; ++j)
                    {
                        e.f[j] = elmt[i]->GetFace(j)->m_id;
                    }
                    HexInfo.push_back(e);
                }
                break;
                default:
                    ASSERTL0(false, "Unknown element type");
            }
        }

        if (SegInfo.size())
        {
            std::string vType("S");
            TiXmlElement *x = new TiXmlElement(vType);
            std::string Str;
            LibUtilities::CompressData::ZlibEncodeToBase64Str(SegInfo, Str);
            x->SetAttribute("COMPRESSED",
                            LibUtilities::CompressData::GetCompressString());
            x->SetAttribute("BITSIZE",
                            LibUtilities::CompressData::GetBitSizeStr());
            x->LinkEndChild(new TiXmlText(Str));
            verTag->LinkEndChild(x);
        }

        if (TriInfo.size())
        {
            std::string vType("T");
            TiXmlElement *x = new TiXmlElement(vType);
            std::string Str;
            LibUtilities::CompressData::ZlibEncodeToBase64Str(TriInfo, Str);
            x->SetAttribute("COMPRESSED",
                            LibUtilities::CompressData::GetCompressString());
            x->SetAttribute("BITSIZE",
                            LibUtilities::CompressData::GetBitSizeStr());
            x->LinkEndChild(new TiXmlText(Str));
            verTag->LinkEndChild(x);
        }

        if (QuadInfo.size())
        {
            std::string vType("Q");
            TiXmlElement *x = new TiXmlElement(vType);
            std::string Str;
            LibUtilities::CompressData::ZlibEncodeToBase64Str(QuadInfo, Str);
            x->SetAttribute("COMPRESSED",
                            LibUtilities::CompressData::GetCompressString());
            x->SetAttribute("BITSIZE",
                            LibUtilities::CompressData::GetBitSizeStr());
            x->LinkEndChild(new TiXmlText(Str));
            verTag->LinkEndChild(x);
        }

        if (TetInfo.size())
        {
            std::string vType("A");
            TiXmlElement *x = new TiXmlElement(vType);
            std::string Str;
            LibUtilities::CompressData::ZlibEncodeToBase64Str(TetInfo, Str);
            x->SetAttribute("COMPRESSED",
                            LibUtilities::CompressData::GetCompressString());
            x->SetAttribute("BITSIZE",
                            LibUtilities::CompressData::GetBitSizeStr());
            x->LinkEndChild(new TiXmlText(Str));
            verTag->LinkEndChild(x);
        }

        if (PyrInfo.size())
        {
            std::string vType("P");
            TiXmlElement *x = new TiXmlElement(vType);
            std::string Str;
            LibUtilities::CompressData::ZlibEncodeToBase64Str(PyrInfo, Str);
            x->SetAttribute("COMPRESSED",
                            LibUtilities::CompressData::GetCompressString());
            x->SetAttribute("BITSIZE",
                            LibUtilities::CompressData::GetBitSizeStr());
            x->LinkEndChild(new TiXmlText(Str));
            verTag->LinkEndChild(x);
        }

        if (PrismInfo.size())
        {
            std::string vType("R");
            TiXmlElement *x = new TiXmlElement(vType);
            std::string Str;
            LibUtilities::CompressData::ZlibEncodeToBase64Str(PrismInfo, Str);
            x->SetAttribute("COMPRESSED",
                            LibUtilities::CompressData::GetCompressString());
            x->SetAttribute("BITSIZE",
                            LibUtilities::CompressData::GetBitSizeStr());
            x->LinkEndChild(new TiXmlText(Str));
            verTag->LinkEndChild(x);
        }

        if (HexInfo.size())
        {
            std::string vType("H");
            TiXmlElement *x = new TiXmlElement(vType);
            std::string Str;
            LibUtilities::CompressData::ZlibEncodeToBase64Str(HexInfo, Str);
            x->SetAttribute("COMPRESSED",
                            LibUtilities::CompressData::GetCompressString());
            x->SetAttribute("BITSIZE",
                            LibUtilities::CompressData::GetBitSizeStr());
            x->LinkEndChild(new TiXmlText(Str));
            verTag->LinkEndChild(x);
        }
    }

    pRoot->LinkEndChild(verTag);
}

}
}

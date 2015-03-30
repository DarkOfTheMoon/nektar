////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputXmlBase.cpp
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
#include "OutputXmlBase.h"

namespace Nektar
{
    namespace Utilities
    {

        OutputXmlBase::OutputXmlBase(MeshSharedPtr m) : OutputModule(m)
        {
        }

        OutputXmlBase::~OutputXmlBase()
        {
        }
        

        void OutputXmlBase::WriteXmlCurves(TiXmlElement * pRoot)
        {
            int edgecnt = 0;

            bool curve = false;
            EdgeSet::iterator it;
            for (it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); ++it)
            {
                if ((*it)->m_edgeNodes.size() > 0) 
                {
                    curve = true;
                    break;
                }
            }
            if (!curve) return;

            TiXmlElement * curved = new TiXmlElement ("CURVED" );

            for (it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); ++it)
            {
                if ((*it)->m_edgeNodes.size() > 0)
                {
                    TiXmlElement * e = new TiXmlElement( "E" );
                    e->SetAttribute("ID",        edgecnt++);
                    e->SetAttribute("EDGEID",    (*it)->m_id);
                    e->SetAttribute("NUMPOINTS", (*it)->GetNodeCount());
                    e->SetAttribute("TYPE", 
                        LibUtilities::kPointsTypeStr[(*it)->m_curveType]);
                    TiXmlText * t0 = new TiXmlText((*it)->GetXmlCurveString());
                    e->LinkEndChild(t0);
                    curved->LinkEndChild(e);
                }
            }

            int facecnt = 0;

            // 2D elements in 3-space, output face curvature information
            if (m_mesh->m_expDim == 2 && m_mesh->m_spaceDim == 3)
            {
                vector<ElementSharedPtr>::iterator it;
                for (it  = m_mesh->m_element[m_mesh->m_expDim].begin();
                     it != m_mesh->m_element[m_mesh->m_expDim].end(); ++it)
                {
                    // Only generate face curve if there are volume nodes
                    if ((*it)->GetVolumeNodes().size() > 0)
                    {
                        TiXmlElement * e = new TiXmlElement( "F" );
                        e->SetAttribute("ID",        facecnt++);
                        e->SetAttribute("FACEID",    (*it)->GetId());
                        e->SetAttribute("NUMPOINTS", (*it)->GetNodeCount());
                        e->SetAttribute("TYPE",
                           LibUtilities::kPointsTypeStr[(*it)->GetCurveType()]);

                        TiXmlText * t0 = new TiXmlText((*it)->GetXmlCurveString());
                        e->LinkEndChild(t0);
                        curved->LinkEndChild(e);
                    }
                }
            }
            else if (m_mesh->m_expDim == 3)
            {
                FaceSet::iterator it2;
                for (it2 = m_mesh->m_faceSet.begin(); it2 != m_mesh->m_faceSet.end(); ++it2)
                {
                    if ((*it2)->m_faceNodes.size() > 0)
                    {
                        TiXmlElement * f = new TiXmlElement( "F" );
                        f->SetAttribute("ID",       facecnt++);
                        f->SetAttribute("FACEID",   (*it2)->m_id);
                        f->SetAttribute("NUMPOINTS",(*it2)->GetNodeCount());
                        f->SetAttribute("TYPE",
                                        LibUtilities::kPointsTypeStr[(*it2)->m_curveType]);
                        TiXmlText * t0 = new TiXmlText((*it2)->GetXmlCurveString());
                        f->LinkEndChild(t0);
                        curved->LinkEndChild(f);
                    }
                }
            }

            pRoot->LinkEndChild( curved );
        }

        void OutputXmlBase::WriteXmlComposites(TiXmlElement * pRoot)
        {
            TiXmlElement* verTag = new TiXmlElement("COMPOSITE");
            CompositeMap::iterator it;
            ConditionMap::iterator it2;
            int j = 0;

            for (it = m_mesh->m_composite.begin(); it != m_mesh->m_composite.end(); ++it, ++j)
            {
                if (it->second->m_items.size() > 0) 
                {
                    TiXmlElement *comp_tag = new TiXmlElement("C"); // Composite
                    bool doSort = true;
                    
                    // Ensure that this composite is not used for periodic BCs!
                    for (it2  = m_mesh->m_condition.begin(); 
                         it2 != m_mesh->m_condition.end(); ++it2)
                    {
                        ConditionSharedPtr c = it2->second;
                        
                        // Ignore non-periodic boundary conditions.
                        if (find(c->type.begin(), c->type.end(), ePeriodic) ==
                            c->type.end())
                        {
                            continue;
                        }

                        for (int i = 0; i < c->m_composite.size(); ++i)
                        {
                            if (c->m_composite[i] == j)
                            {
                                doSort = false;
                            }
                        }
                    }

                    doSort = doSort && it->second->m_reorder;
                    comp_tag->SetAttribute("ID", it->second->m_id);
                    comp_tag->LinkEndChild(
                        new TiXmlText(it->second->GetXmlString(doSort)));
                    verTag->LinkEndChild(comp_tag);
                }
                else
                {
                    cout << "Composite " << it->second->m_id << " "
                         << "contains nothing." << endl;
                }
            }

            pRoot->LinkEndChild(verTag);
        }

        void OutputXmlBase::WriteXmlDomain(TiXmlElement * pRoot)
        {
            // Write the <DOMAIN> subsection.
            TiXmlElement * domain = new TiXmlElement ("DOMAIN" );
            std::string list;
            CompositeMap::iterator it;
            
            for (it = m_mesh->m_composite.begin(); it != m_mesh->m_composite.end(); ++it)
            {
                if (it->second->m_items[0]->GetDim() == m_mesh->m_expDim)
                {
                    if (list.length() > 0)
                    {
                        list += ",";
                    }
                    list += boost::lexical_cast<std::string>(it->second->m_id);
                }
            }
            domain->LinkEndChild( new TiXmlText(" C[" + list + "] "));
            pRoot->LinkEndChild( domain );
        }

        void OutputXmlBase::WriteXmlExpansions(TiXmlElement * pRoot)
        {
            // Write a default <EXPANSIONS> section.
            TiXmlElement * expansions = new TiXmlElement ("EXPANSIONS");
            CompositeMap::iterator it;
            
            for (it = m_mesh->m_composite.begin(); it != m_mesh->m_composite.end(); ++it)
            {
                if (it->second->m_items[0]->GetDim() == m_mesh->m_expDim)
                {
                    TiXmlElement * exp = new TiXmlElement ( "E");
                    exp->SetAttribute("COMPOSITE", "C["
                        + boost::lexical_cast<std::string>(it->second->m_id)
                        + "]");
                    exp->SetAttribute("NUMMODES",4);
                    exp->SetAttribute("TYPE","MODIFIED");
                    
                    if (m_mesh->m_fields.size() == 0)
                    {
                        exp->SetAttribute("FIELDS","u");
                    }
                    else
                    {
                        string fstr;
                        for (int i = 0; i < m_mesh->m_fields.size(); ++i)
                        {
                            fstr += m_mesh->m_fields[i]+",";
                        }
                        fstr = fstr.substr(0,fstr.length()-1);
                        exp->SetAttribute("FIELDS", fstr);
                    }
                    
                    expansions->LinkEndChild(exp);
                }
            }
            pRoot->LinkEndChild(expansions);
        }
        
        void OutputXmlBase::WriteXmlConditions(TiXmlElement * pRoot)
        {
            TiXmlElement *conditions = 
                new TiXmlElement("CONDITIONS");
            TiXmlElement *boundaryregions = 
                new TiXmlElement("BOUNDARYREGIONS");
            TiXmlElement *boundaryconditions = 
                new TiXmlElement("BOUNDARYCONDITIONS");
            TiXmlElement *variables = 
                new TiXmlElement("VARIABLES");
            ConditionMap::iterator it;
            
            for (it = m_mesh->m_condition.begin(); it != m_mesh->m_condition.end(); ++it)
            {
                ConditionSharedPtr c = it->second;
                string tmp;
                
                // First set up boundary regions.
                TiXmlElement *b = new TiXmlElement("B");
                b->SetAttribute("ID", boost::lexical_cast<string>(it->first));
                
                for (int i = 0; i < c->m_composite.size(); ++i)
                {
                    tmp += boost::lexical_cast<string>(c->m_composite[i]) + ",";
                }
                
                tmp = tmp.substr(0, tmp.length()-1);

                TiXmlText *t0 = new TiXmlText("C["+tmp+"]");
                b->LinkEndChild(t0);
                boundaryregions->LinkEndChild(b);
                
                TiXmlElement *region = new TiXmlElement("REGION");
                region->SetAttribute(
                    "REF", boost::lexical_cast<string>(it->first));
                
                for (int i = 0; i < c->type.size(); ++i)
                {
                    string tagId;
                    
                    switch(c->type[i])
                    {
                        case eDirichlet:    tagId = "D"; break;
                        case eNeumann:      tagId = "N"; break;
                        case ePeriodic:     tagId = "P"; break;
                        case eHOPCondition: tagId = "N"; break;
                        default:                         break;
                    }
                    
                    TiXmlElement *tag = new TiXmlElement(tagId);
                    tag->SetAttribute("VAR", c->field[i]);
                    tag->SetAttribute("VALUE", c->value[i]);
                    
                    if (c->type[i] == eHOPCondition)
                    {
                        tag->SetAttribute("USERDEFINEDTYPE", "H");
                    }
                    
                    region->LinkEndChild(tag);
                }
                
                boundaryconditions->LinkEndChild(region);
            }

            for (int i = 0; i < m_mesh->m_fields.size(); ++i)
            {
                TiXmlElement *v = new TiXmlElement("V");
                v->SetAttribute("ID", boost::lexical_cast<std::string>(i));
                TiXmlText *t0 = new TiXmlText(m_mesh->m_fields[i]);
                v->LinkEndChild(t0);
                variables->LinkEndChild(v);
            }
            
            if (m_mesh->m_fields.size() > 0)
            {
                conditions->LinkEndChild(variables);
            }
            
            if (m_mesh->m_condition.size() > 0)
            {
                conditions->LinkEndChild(boundaryregions);
                conditions->LinkEndChild(boundaryconditions);
            }
            
            pRoot->LinkEndChild(conditions);
        }   
    }
}

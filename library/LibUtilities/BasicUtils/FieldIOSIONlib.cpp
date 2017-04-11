////////////////////////////////////////////////////////////////////////////////
//
//  File: FieldIOSIONlib.cpp
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
//  Description: I/O routines relating to Fields into SIONlib files.
//
////////////////////////////////////////////////////////////////////////////////

#include "sion.h"

#include <LibUtilities/BasicUtils/SIONlib.h>
#include <LibUtilities/BasicUtils/FieldIOSIONlib.h>
#include <boost/unordered_set.hpp>
#include <LibUtilities/Communication/CommMpi.h>


namespace berrc = boost::system::errc;

namespace Nektar
{
namespace LibUtilities
{ 

std::string FieldIOSIONlib::className =
    GetFieldIOFactory().RegisterCreatorFunction(
        "SIONlib", FieldIOSIONlib::create, "SIONlib-based output of field data.");

  
const unsigned int FieldIOSIONlib::FORMAT_VERSION = 1;
  
const std::string FieldIOSIONlib::ATTRNAME_FIELDS = "FIELDS";
const std::string FieldIOSIONlib::ATTRNAME_BASIS = "BASIS";
const std::string FieldIOSIONlib::ATTRNAME_SHAPE = "SHAPE";
const std::string FieldIOSIONlib::ATTRNAME_HOMOLENS = "HOMOGENEOUSLENGTHS";
const std::string FieldIOSIONlib::ATTRNAME_NUMMODES = "NUMMODESPERDIR";
const std::string FieldIOSIONlib::ATTRNAME_POINTSTYPE = "POINTSTYPE";
const std::string FieldIOSIONlib::ATTRNAME_NUMPOINTS = "NUMPOINTSPERDIR";
  
const std::string FieldIOSIONlib::ATTRVALUE_MIXORDER = "MIXORDER";
const std::string FieldIOSIONlib::ATTRVALUE_UNIORDER = "UNIORDER";

sion_int32 FieldIOSIONlib::block_size = -1;
sion_int64 FieldIOSIONlib::chunk_size = -1;
std::string FieldIOSIONlib::read_mode = "";
std::string FieldIOSIONlib::write_mode = "";

  
FieldIOSIONlib::FieldIOSIONlib(LibUtilities::CommSharedPtr pComm,
                         bool sharedFilesystem)
    : FieldIO(pComm, sharedFilesystem)
{
}


void FieldIOSIONlib::v_Init(const LibUtilities::SessionReaderSharedPtr session)
{
    block_size = -1;
    if (session->DefinesSolverInfo("IOBlockSize"))
    {
        std::string int_str = session->GetSolverInfo("IOBlockSize");
	std::istringstream(int_str) >> block_size;
    }

    sion_int64 blocks_per_chunk = 1;
    if (session->DefinesSolverInfo("IOBlocksPerChunk"))
    {
        std::string int_str = session->GetSolverInfo("IOBlocksPerChunk");
        std::istringstream(int_str) >> blocks_per_chunk;
    }

    chunk_size = ((sion_int64) block_size)*blocks_per_chunk;

    read_mode = "br";
    if (session->DefinesSolverInfo("IOReadMode"))
    {
        read_mode = session->GetSolverInfo("IOReadMode");
    }

    write_mode = "bw";
    if (session->DefinesSolverInfo("IOWriteMode"))
    {
        write_mode = session->GetSolverInfo("IOWriteMode");
    }
}


void FieldIOSIONlib::v_InitFromBenchmarker(const LibUtilities::IOSettingsSharedPtr iosettings)
{
    LibUtilities::IOSettings::iterator it;
    
    block_size = -1;
    it = iosettings->find("IOBlockSize");
    if (iosettings->end() != it)
    {
        std::istringstream(it->second) >> block_size;
    }

    sion_int64 blocks_per_chunk = 1;
    it = iosettings->find("IOBlocksPerChunk");
    if (iosettings->end() != it)
    {
        std::istringstream(it->second) >> blocks_per_chunk;
    }
    chunk_size = ((sion_int64) block_size)*blocks_per_chunk;

    read_mode = "br";
    it = iosettings->find("IOReadMode");
    if (iosettings->end() != it)
    {
        read_mode = it->second;
    }

    write_mode = "bw";
    it = iosettings->find("IOWriteMode");
    if (iosettings->end() != it)
    {
        write_mode = it->second;
    }
}

  
/**
 * Open a SIONlib file for writing by constructing and returning
 * a pointer to a SIONlib object, and also outputting the format version.
 *
 * @param outFilename  name of output file.
 * @return Pointer to SIONlib object.
 */
SIONlib::SIONFile *FieldIOSIONlib::OpenFileForWriting(const std::string &outFilename)
{
    CommMpi *commMpiPtr =  (CommMpi*) m_comm.get();
    SIONlib::SIONFile *fp = new SIONlib::SIONFile(outFilename, write_mode,
        1, chunk_size, block_size,
        m_comm->GetRank(), commMpiPtr->GetComm(), commMpiPtr->GetComm());

    ASSERTL0(NULL != fp,
             "Unable to construct SIONFile object for " + outFilename);

    SIONlib::SIONFile &f = *fp;
    f.open();
    
    ASSERTL0(-1 != f.getReturnCode(),
             "Unable to open SIONlib file " + outFilename);

    f << FORMAT_VERSION;
    
    return fp;
}

  
unsigned long FieldIOSIONlib::v_Write(const std::string &outFilePrefix,
    std::vector<FieldDefinitionsSharedPtr> &fielddefs,
    std::vector<std::vector<NekDouble> > &fielddata,
    const FieldMetaDataMap &fieldmetadatamap,
    const bool backup)
{
    // We make a number of assumptions in this code:
    //   1. All element ids have the same type: unsigned int
    //   2. All elements within a given field have the same number of values
    //   3. All element values have the same type, NekDouble

    ASSERTL1(fielddefs.size() == fielddata.size(),
             "fielddefs and fielddata have incompatible lengths.");

    size_t nFields = fielddefs.size();

    int homDim = -1;
    int varOrder = 0;

    for (size_t i = 0; i < nFields; ++i)
    {
        if (!fielddefs[i]->m_uniOrder)
        {
            varOrder = 1;
            break;
        }
    }

    m_comm->AllReduce(varOrder, LibUtilities::ReduceMax);

    SetUpOutput(outFilePrefix, false, backup);

    // Each MPI process iterates through its fields and outputs field
    // definitions and field data to the SIONlib file.
    SIONlib::SIONFile* fp = OpenFileForWriting(outFilePrefix);
    ASSERTL0(NULL != fp, "Cannot open SIONlib file.");
    SIONlib::SIONFile &f = *fp; 

    f << nFields;
    
    for (size_t i = 0; i < nFields; ++i)
    {
        ASSERTL1(fielddata[i].size() > 0,
            "fielddata vector must contain at least one value.");
        ASSERTL1(fielddata[i].size() ==
            fielddefs[i]->m_fields.size()*CheckFieldDefinition(fielddefs[i]),
            "fielddata vector has invalid size.");

        FieldDefinitionsSharedPtr def = fielddefs[i];

        size_t nAttrs = (def->m_numHomogeneousDir ? 5 : 4);
        if (def->m_pointsDef) nAttrs++;
        if (def->m_numPointsDef) nAttrs++;
        f << nAttrs;
          
        // FIELDS
        ////////////////////////////////////////////////////////////////
        f << ATTRNAME_FIELDS;
        
        size_t nAttrFields = def->m_fields.size();
        f << nAttrFields;
        
        for (size_t j = 0; j < nAttrFields; ++j)
        {
            f << def->m_fields[j];
        }
        ////////////////////////////////////////////////////////////////
                
        // BASIS
        ////////////////////////////////////////////////////////////////
        f << ATTRNAME_BASIS;
        
        nAttrFields = def->m_basis.size();
        f << nAttrFields;
        
        for (size_t j = 0; j < nAttrFields; ++j)
        {
            f << def->m_basis[j];
        }
        ////////////////////////////////////////////////////////////////
        
        // SHAPE
        ////////////////////////////////////////////////////////////////
        f << ATTRNAME_SHAPE;

        std::stringstream shapeStringStream;
        shapeStringStream << ShapeTypeMap[def->m_shapeType];
        
        if (def->m_numHomogeneousDir > 0)
        {
            if (homDim == -1)
            {
                homDim = def->m_numHomogeneousDir;
            }

            ASSERTL1(homDim == def->m_numHomogeneousDir,
                "HDF5 does not support variable homogeneous directions in "
                "the same file.");

            shapeStringStream << "-HomogenousExp"
                << def->m_numHomogeneousDir << "D";
        }
        
        if (def->m_homoStrips)
        {
            shapeStringStream << "-Strips";
        }
        
        f << shapeStringStream.str();
        ////////////////////////////////////////////////////////////////
        
        // Determine HOMOGENEOUS attributes
        if (def->m_numHomogeneousDir)
        {
            // HOMOGENEOUSLENGTHS
            ////////////////////////////////////////////////////////////////
            f << ATTRNAME_HOMOLENS;
            
            std::vector<NekDouble>& homoLens = def->m_homogeneousLengths;
            size_t nLens = homoLens.size();
            f << nLens;
            for (size_t j = 0; j < nLens; ++j)
            {
                f << homoLens[j];
            }
            ////////////////////////////////////////////////////////////////

            // homo IDs are streamed to separate files
            ////////////////////////////////////////////////////////////////
            std::vector<unsigned int>& homoYIDs = def->m_homogeneousYIDs;
            size_t nIDs = homoYIDs.size();
            f << nIDs;
            for (size_t j = 0; j < nIDs; ++j)
            {
                f << homoYIDs[j];
            }
            
            std::vector<unsigned int>& homoZIDs = def->m_homogeneousZIDs;
            nIDs = homoZIDs.size();
            f << nIDs;
            for (size_t j = 0; j < nIDs; ++j)
            {
                f << homoZIDs[j];
            }
            
            std::vector<unsigned int>& homoSIDs = def->m_homogeneousSIDs;
            nIDs = homoSIDs.size();
            f << nIDs;
            for (size_t j = 0; j < nIDs; ++j)
            {
                f << homoSIDs[j];
            }
            ////////////////////////////////////////////////////////////////
        }
        
        // NUMMODESPERDIR
        ////////////////////////////////////////////////////////////////
        f << ATTRNAME_NUMMODES;

        std::vector<unsigned int>& numModes = def->m_numModes;
        size_t nNumModes = def->m_basis.size();
        
        if (def->m_uniOrder && !varOrder)
        {
            f << ATTRVALUE_UNIORDER;
            stringstream numModesStringStream;
                
            for (size_t j = 0; j < nNumModes; ++j)
            {
                if (j > 0)
                {
                    numModesStringStream << ",";
                }
                numModesStringStream << numModes[j];
            }

            f << numModesStringStream.str();
        }
        else
        {
            f << ATTRVALUE_MIXORDER;
            f << nNumModes;
            
            for (size_t j = 0; j < nNumModes; ++j)
            {
                f << numModes[j];
            }
        }
        ////////////////////////////////////////////////////////////////

        // POINTSTYPE
        ////////////////////////////////////////////////////////////////
        if (def->m_pointsDef)
        {
            f << ATTRNAME_POINTSTYPE;

            stringstream pointsTypeStringStream;
            std::vector<LibUtilities::PointsType>& points = def->m_points;
            size_t nPoints = points.size();
            
            for (size_t j = 0; j < nPoints; ++j)
            {
                if (j > 0)
                {
                    pointsTypeStringStream << ",";
                }
                pointsTypeStringStream << kPointsTypeStr[points[j]];
            }

            f << pointsTypeStringStream.str();
        }
        ////////////////////////////////////////////////////////////////

        // NUMPOINTSPERDIR
        ////////////////////////////////////////////////////////////////
        if (def->m_numPointsDef)
        {
            f << ATTRNAME_NUMPOINTS;

            stringstream numPointsStringStream;
            std::vector<unsigned int>& numPoints = def->m_numPoints;
            size_t nNumPoints = numPoints.size();
            
            for (size_t j = 0; j < nNumPoints; ++j)
            {
                if (j > 0)
                {
                    numPointsStringStream << ",";
                }
                numPointsStringStream << numPoints[j];
            }

            f << numPointsStringStream.str();
        }
        ////////////////////////////////////////////////////////////////
                
        // elementIDs
        ////////////////////////////////////////////////////////////////
        std::vector<unsigned int>& elemIDs = def->m_elementIDs;
        size_t nElems = elemIDs.size();
        f << nElems;
        for (size_t j = 0; j < nElems; ++j)
        {
            f << elemIDs[j];
        }
        ////////////////////////////////////////////////////////////////

        // data
        ////////////////////////////////////////////////////////////////
        std::vector<NekDouble>& data = fielddata[i];
        size_t nVals = data.size();
        f << nVals;
        for (size_t j = 0; j < nVals; ++j)
        {
            f << data[j];
        }
        ////////////////////////////////////////////////////////////////
    }

    unsigned long nWritten = f.getBytesWritten();
    f.close();
    
    m_comm->Block();
    return nWritten;
}
  

/**
 * Open a SIONlib file for reading by constructing and returning
 * a pointer to a SIONlib object, and also read in the format version,
 * throwing an error if it is greater than FieldIOSIONlib::FORMAT_VERSION.
 *
 * @param inFilename  name of input file.
 * @return Pointer to SIONlib object.
 */
SIONlib::SIONFile *FieldIOSIONlib::OpenFileForReading(const std::string &inFilename)
{
    CommMpi *commMpiPtr =  (CommMpi*) m_comm.get();
    SIONlib::SIONFile *fp = new SIONlib::SIONFile(inFilename, read_mode,
        1, chunk_size, block_size,
        m_comm->GetRank(), commMpiPtr->GetComm(), commMpiPtr->GetComm());

    ASSERTL0(NULL != fp,
             "Unable to construct SIONFile object for " + inFilename);

    SIONlib::SIONFile &f = *fp;
    f.open();
    
    ASSERTL0(-1 != f.getReturnCode(),
             "Unable to open SIONlib file " + inFilename);

    unsigned int formatVersion = 0;
    f >> formatVersion;

    std::stringstream version;
    version << formatVersion;

    ASSERTL0(formatVersion <= FORMAT_VERSION,
             "File " + inFilename + " is at version " + version.str() +
             ", which is higher than supported in this version of Nektar++.");
    
    return fp;
}

unsigned long FieldIOSIONlib::v_Import(const std::string &infilename,
    std::vector<FieldDefinitionsSharedPtr> &fielddefs,
    std::vector<std::vector<NekDouble> > &fielddata,
    FieldMetaDataMap &fieldinfomap,
    const Array<OneD, int> &ElementIDs)  
{
    SIONlib::SIONFile* fp = OpenFileForReading(infilename);
    
    ASSERTL0(NULL != fp, "Cannot open SIONlib file.");

    SIONlib::SIONFile &f = *fp;    
    size_t nFields = 0;
        
    f >> nFields;
    
    for (size_t i = 0; i < nFields; ++i)
    {
        FieldDefinitionsSharedPtr def =
            MemoryManager<FieldDefinitions>::AllocateSharedPtr();
        
        std::string attrName, attrValue;
        size_t nAttrs, nAttrFields;

        f >> nAttrs;
        
        for (size_t j = 0; j < nAttrs; ++j)
        {         
            f >> attrName;

            if (attrName == ATTRNAME_FIELDS)
            {
                f >> nAttrFields;
                def->m_fields.resize(nAttrFields);

                for (size_t k = 0; k < nAttrFields; ++k)
                {
                    f >> def->m_fields[k];
                }
            }
            else if (attrName == ATTRNAME_BASIS)
            {
                f >> nAttrFields;
                def->m_basis.resize(nAttrFields);

                for (size_t k = 0; k < nAttrFields; ++k)
                {
                    f >> def->m_basis[k];
                }
            
                // check the basis is in range
                std::vector<BasisType>::const_iterator bIt  = def->m_basis.begin();
                std::vector<BasisType>::const_iterator bEnd = def->m_basis.end();
                for (; bIt != bEnd; ++bIt)
                {
                    BasisType bt = *bIt;
                    ASSERTL0(bt >= 0 && bt < SIZE_BasisType,
                        "Unable to correctly parse the basis types.");
                }
            }
            else if (attrName == ATTRNAME_SHAPE)
            {
                f >> attrValue;        
            
                // check to see if homogeneous expansion and if so
                // strip down the shapeString definition
                size_t loc;
                //---> this finds the first location of 'n'!
                if (attrValue.find("Strips") != string::npos)
                {
                    def->m_homoStrips = true;
                }

                if ((loc = attrValue.find_first_of("-")) != string::npos)
                {
                    if (attrValue.find("Exp1D") != string::npos)
                    {
                        def->m_numHomogeneousDir = 1;
                    }
                    else // HomogeneousExp1D
                    {
                        def->m_numHomogeneousDir = 2;
                    }

                    attrValue.erase(loc, attrValue.length());
                }

                // get the geometrical shape
                bool valid = false;
                for (unsigned int k = 0; k < SIZE_ShapeType; k++)
                {
                    if (ShapeTypeMap[k] == attrValue)
                    {
                        def->m_shapeType = (ShapeType) k;
                        valid = true;
                        break;
                    }
                }

                ASSERTL0(valid,
                    std::string("Unable to correctly parse the shape type: ")
                        .append(attrValue)
                        .c_str());
            }
            else if (attrName == ATTRNAME_HOMOLENS)
            {
                def->m_numHomogeneousDir = true;
          
                size_t nLens;
                f >> nLens;
                def->m_homogeneousLengths.resize(nLens);

                for (size_t k = 0; k < nLens; ++k)
                {
                    f >> def->m_homogeneousLengths[k];
                }

                size_t nIDs;
                f >> nIDs;
                def->m_homogeneousYIDs.resize(nIDs);

                for (size_t k = 0; k < nIDs; ++k)
                {
                    f >> def->m_homogeneousYIDs[k];
                }

                f >> nIDs;
                def->m_homogeneousZIDs.resize(nIDs);

                for (size_t k = 0; k < nIDs; ++k)
                {
                    f >> def->m_homogeneousZIDs[k];
                }

                f >> nIDs;
                def->m_homogeneousSIDs.resize(nIDs);

                for (size_t k = 0; k < nIDs; ++k)
                {
                    f >> def->m_homogeneousSIDs[k];
                }
            }
            else if (attrName == ATTRNAME_NUMMODES)
            {
                f >> attrValue;

                if (attrValue == ATTRVALUE_UNIORDER)
                {
                    def->m_uniOrder = true;

                    f >> attrValue;
                    bool valid = ParseUtils::GenerateOrderedVector(
                        attrValue.c_str(), def->m_numModes);

                    ASSERTL0(valid,
                        "Unable to correctly parse the number of modes.");
                }
                else if (attrValue == ATTRVALUE_MIXORDER)
                {
                    size_t nNumModes;
                
                    f >> nNumModes;
                    def->m_numModes.resize(nNumModes);
                
                    for (size_t k = 0; k < nNumModes; ++k)
                    {
                        f >> def->m_numModes[k];
                    }
                }
                else
                {
                    std::string errstr("Unknown " + attrName + " value: ");
                    errstr += attrValue;
                    ASSERTL1(false, errstr.c_str());
                }
            }
            else if (attrName == ATTRNAME_POINTSTYPE)
            {
                def->m_pointsDef = true;
            
                f >> attrValue;
            
                std::vector<std::string> pointsStrings;
                bool valid = ParseUtils::GenerateOrderedStringVector(
                    attrValue.c_str(), pointsStrings);

                ASSERTL0(valid,
                    "Unable to correctly parse the points types.");
                
                for (std::vector<std::string>::size_type k = 0;
                     k < pointsStrings.size();
                     ++k)
                {
                    valid = false;
                    for (unsigned int l = 0; l < SIZE_PointsType; ++l)
                    {
                        if (kPointsTypeStr[l] == pointsStrings[k])
                        {
                            def->m_points.push_back((PointsType)l);
                            valid = true;
                            break;
                        }
                    }

                    ASSERTL0(valid,
                        std::string(
                            "Unable to correctly parse the points type: ")
                            .append(pointsStrings[k])
                            .c_str());
                }
            }
            else if (attrName == ATTRNAME_NUMPOINTS)
            {
                def->m_numPointsDef = true;
            
                f >> attrValue;

                bool valid = ParseUtils::GenerateOrderedVector(
                    attrValue.c_str(), def->m_numPoints);
                ASSERTL0(valid,
                    "Unable to correctly parse the number of points.");
            }
            else
            {
                std::string errstr("Unknown field attribute: ");
                errstr += attrName;
                ASSERTL1(false, errstr.c_str());
            }
        }

        // element IDs
        ////////////////////////////////////////////////
        size_t nElems;
        f >> nElems;
        
        def->m_elementIDs.resize(nElems);        
        for (size_t j = 0; j < nElems; ++j)
        {
            f >> def->m_elementIDs[j];
        }
        ////////////////////////////////////////////////
        
        // element data
        ////////////////////////////////////////////////
        ASSERTL0(fielddata != NullVectorNekDoubleVector,
            "Null fielddata.");
        
        size_t nVals;
        f >> nVals;

        std::vector<NekDouble> data(nVals);
        for (size_t j = 0; j < nVals; ++j)
        {
            f >> data[j];
        }

        int datasize = CheckFieldDefinition(def);
        ASSERTL0(data.size() == datasize*def->m_fields.size(),
            "Input data is not the same length as header information.");
        ////////////////////////////////////////////////

        fielddefs.push_back(def);
        fielddata.push_back(data);
    }

    unsigned long nRead = f.getBytesRead();
    f.close();

    return nRead;
}

  
DataSourceSharedPtr FieldIOSIONlib::v_ImportFieldMetaData(
    const std::string &filename,
    FieldMetaDataMap &fieldmetadatamap)
{
    DataSourceSharedPtr ans = SIONlibDataSource::create(filename);
    
    return ans;
}


}
}

////////////////////////////////////////////////////////////////////////////////
//
//  File: SIONlib.cpp
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
//  Description: Minimal SIONlib wrapper
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SIONlib.h>


namespace Nektar
{
namespace LibUtilities
{
namespace SIONlib
{

const unsigned int SION_Base::LUSTRE_STRIPE_BASE_SIZE = 65536;  

// Common functions
char *SION_Base::getSionFileName() const
{
    return _sion_file_name;
}

void SION_Base::setMode(std::string mode)
{
    _mode = mode;
}

std::string SION_Base::getMode() const
{
    return _mode;
}


void SION_Base::setNumberOfFiles(int num_files)
{
    _num_files = num_files;
}

int SION_Base::getNumberOfFiles() const
{
    return _num_files;
}

void SION_Base::setNumberOfTasks(int num_tasks)
{
    _num_tasks = num_tasks;
}

int SION_Base::getNumberOfTasks() const
{
    return _num_tasks;
}

void SION_Base::setRank(int rank)
{
    _rank = rank;
}

int SION_Base::getRank() const
{
    return _rank;
}

void SION_Base::setChunkSize(sion_int64 chunk_size)
{
    _chunk_size = chunk_size;
}

sion_int64 SION_Base::getChunkSize() const
{
    return _chunk_size;
}

void SION_Base::setChunkSizes(sion_int64 *chunk_sizes)
{
    _chunk_sizes = chunk_sizes;
}

sion_int64 *SION_Base::getChunkSizes() const
{
    return _chunk_sizes;
}

void SION_Base::setGlobalRanks(int *global_ranks)
{
    _global_ranks = global_ranks;
}

int *SION_Base::getGlobalRanks() const
{
    return _global_ranks;
}

void SION_Base::setFileSystemBlockSize(sion_int32 fs_blk_size)
{
    _fs_blk_size = fs_blk_size;
}

sion_int32 SION_Base::getFileSystemBlockSize() const
{
    return _fs_blk_size;
}

int SION_Base::getNumberOfSuccessfulReadElements() const
{
    return _number_of_elements_sucessfully_read;
}

int SION_Base::getReturnCode() const
{
    return _return_code;
}

int SION_Base::getSid() const
{
    return _sid;
}


/* get information (with sion datatypes) */
int SION_Base::getFileEndianness() const
{
    return sion_get_file_endianness(_sid);
}

sion_int64 SION_Base::getBytesWritten() const
{
    return sion_get_bytes_written(_sid);
}

sion_int64 SION_Base::getBytesRead() const
{
    return sion_get_bytes_read(_sid);
}

sion_int64 SION_Base::getBytesAvailInBlock() const
{
    return sion_bytes_avail_in_block(_sid);
}

sion_int64 SION_Base::getBytesAvailInChunk() const
{
    return sion_bytes_avail_in_chunk(_sid);
}

sion_int64 SION_Base::getPosition() const
{
    return sion_get_position(_sid);
}

void SION_Base::seek()
{
    _return_code = sion_seek(_sid,
        SION_CURRENT_RANK, SION_CURRENT_BLK, SION_CURRENT_POS);
}

  

SIONFile::SIONFile(std::string sion_file_name, std::string mode,
    int num_files, sion_int64 chunk_size, sion_int32 block_size,
    int global_rank, MPI_Comm gComm, MPI_Comm lComm)
{
    size_t ncharacter = sion_file_name.length()+1;
    const char *tmp_sion_file_name = sion_file_name.c_str();
    _sion_file_name = new char[ncharacter];

    strncpy(_sion_file_name, tmp_sion_file_name, ncharacter);

    _mode = mode;
    _num_files = num_files;

    _global_rank = global_rank;
    
    _g_comm = gComm;
    _l_comm = lComm;

    const sion_int32 base_size = LUSTRE_STRIPE_BASE_SIZE;
    _fs_blk_size = -1;
    
    if (block_size > 0)
    {
        sion_int32 base_cnt = block_size / base_size;

	if (0 == base_cnt)
	{
	    // block_size can never be smaller than base_size
	    block_size = base_size;
	}
	else
	{
            if (0 != block_size % base_size)
            {
                // block_size must be a multiple of base_size
	        base_cnt += 1;
                block_size = base_size*base_cnt;  
            }

            if (base_cnt > 1 && 0 != base_cnt % 2)
            {
                // if block_size > base_size then
	        // block_size must be an even multiple of base_size
	        base_cnt += 1;
                block_size = base_size*base_cnt;
		// this ensures that Lustre file stripe size can be
	        // set equal to the SIONlib file system block size
	    }
        }

        _fs_blk_size = block_size;
    }

    if (chunk_size < base_size)
    {
        chunk_size = base_size;    
    }
    
    block_size = (_fs_blk_size > 0) ? _fs_blk_size : base_size;
    if (0 != chunk_size % block_size)
    {
        // chunk_size must be a multiple of block_size
        sion_int32 block_cnt = chunk_size / block_size + 1;
        chunk_size = block_size*block_cnt;
    }

    _chunk_size = chunk_size;
	
    _file_ptr = NULL;
    _new_sion_file_name = new char[255];
}

SIONFile::~SIONFile()
{
    delete [] _sion_file_name;
    _sion_file_name = NULL;
    delete [] _new_sion_file_name;
    _new_sion_file_name = NULL;
}


void SIONFile::setLocalCommunicator(MPI_Comm lComm)
{
    _l_comm = lComm;
}

MPI_Comm SIONFile::getLocalCommunicator() const
{
    return _l_comm;
}

void SIONFile::setGlobalCommunicator(MPI_Comm gComm)
{
    _g_comm = gComm;
}

MPI_Comm SIONFile::getGlobalCommunicator() const
{
    return _g_comm;
}

void SIONFile::setGlobalRank(int global_rank)
{
    _global_rank = global_rank;
}

int SIONFile::getGlobalRank() const
{
    return _global_rank;
}

char *SIONFile::getNewSionFileName() const
{
    return _new_sion_file_name;
}

void SIONFile::open()
{
    _sid = sion_paropen_mpi(_sion_file_name, _mode.c_str(), &_num_files,
        _g_comm, &_l_comm, &_chunk_size, &_fs_blk_size,
	&_global_rank, NULL, &_new_sion_file_name);
    
    _return_code = _sid;
}

void SIONFile::close()
{
    _return_code = sion_parclose_mpi (_sid);
}

    
void SIONFile::ensureFreeSpace(long numbytes)
{
    _return_code = sion_ensure_free_space(_sid, numbytes);
}

void SIONFile::endOfFile()
{
    _return_code = sion_feof(_sid);
}


template<class T>
void SIONFile::write(T data)
{
    _return_code = sion_fwrite(reinterpret_cast<char*>(&data), sizeof(data), 1, _sid);
}

SIONFile &operator<<(SIONFile &sf, const LibUtilities::BasisType &rhs)
{
    sf.write(rhs);
    return sf;
}
  
SIONFile &operator<<(SIONFile &sf, const unsigned int &rhs)
{
    sf.write(rhs);
    return sf;
}

SIONFile &operator<<(SIONFile &sf, const size_t &rhs)
{
    sf.write(rhs);
    return sf;
}

SIONFile &operator<<(SIONFile &sf, const double &rhs)
{
    sf.write(rhs);
    return sf;
}

/**
 * Write out a std:string object to file, but first write
 * out the string's length.
 *
 * @param data the string to be written to file.
 */
void SIONFile::write_string(const std::string &data)
{
    size_t tmp_size = data.size();

    _return_code = sion_fwrite(&tmp_size, sizeof(tmp_size), 1, _sid);
    if (1 == _return_code && tmp_size > 0)
    {
        _return_code = sion_fwrite(data.c_str(), tmp_size, 1, _sid);
    }
}
  
SIONFile &operator<<(SIONFile &sf, const std::string &rhs)
{
    sf.write_string(rhs);
    return sf;
}


template<class T>
void SIONFile::read(T *data)
{
    T tmp = *data;
    _return_code = sion_fread(data, sizeof(tmp), 1, _sid);
}

SIONFile &operator>>(SIONFile &sf, LibUtilities::BasisType &rhs)
{
    sf.read(&rhs);
    return sf;
}
  
SIONFile &operator>>(SIONFile &sf, unsigned int &rhs)
{
    sf.read(&rhs);
    return sf;
}

SIONFile &operator>>(SIONFile &sf, size_t &rhs)
{
    sf.read(&rhs);
    return sf;
}

SIONFile &operator>>(SIONFile &sf, double &rhs)
{
    sf.read(&rhs);
    return sf;
}

/**
 * Read in a string's length followed by the string's body.
 *
 * @param data the string object that will store the data
 * read from file.
 */
void SIONFile::read_string(std::string *data)
{
    size_t tmp_size = 0;
    
    _return_code = sion_fread(&tmp_size, sizeof(tmp_size), 1, _sid);

    if (1 == _return_code && tmp_size > 0)
    {
        char *tmp_buf = new char[tmp_size];
	if (NULL != tmp_buf)
	{
            _return_code = sion_fread(tmp_buf, tmp_size, 1, _sid);
	    if (1 == _return_code && NULL != data)
	    {
	        data->assign(tmp_buf, tmp_size);
	    }
	    delete [] tmp_buf;
	}
    }
}

SIONFile &operator>>(SIONFile &sf, std::string &rhs)
{
    sf.read_string(&rhs);
    return sf;
}
  

}
}
}

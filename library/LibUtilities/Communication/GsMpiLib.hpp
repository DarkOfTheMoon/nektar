///////////////////////////////////////////////////////////////////////////////
//
// File GsLib.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: MPI specific GS declarations
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_COMMUNICATION_GSMPILIB_HPP
#define NEKTAR_LIB_UTILITIES_COMMUNICATION_GSMPILIB_HPP

#include <iostream>
using namespace std;

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <mpi.h>
using namespace Nektar;

#include <LibUtilities/Communication/GsLib.hpp>
namespace Gs
{
    typedef enum { gs_double, gs_float, gs_int, gs_long, gs_dom_n } gs_dom;
    typedef enum { mode_plain, mode_vec, mode_many, mode_dry_run } gs_mode;

    typedef struct { void *ptr; size_t n,max; } array;
    typedef array buffer;
    typedef MPI_Comm comm_ext;
    typedef MPI_Request comm_req;

    struct comm {
      unsigned int id;
      unsigned int np;
      comm_ext c;
    };

    typedef struct {
      unsigned int n;      /* number of messages */
      unsigned int *p;     /* message source/dest proc */
      unsigned int *size;  /* size of message */
      unsigned int total;  /* sum of message sizes */
    } pw_comm_data;

    typedef struct {
      pw_comm_data comm[2];
      const unsigned int *map[2];
      comm_req *req;
      unsigned int buffer_size;
    } pw_data;

    typedef struct {
      const unsigned int *scatter_map, *gather_map;
      unsigned int size_r, size_r1, size_r2;
      unsigned int size_sk, size_s, size_total;
      unsigned int p1, p2;
      unsigned int nrecvn;
    } cr_stage;

    typedef struct {
      cr_stage *stage[2];
      unsigned int nstages;
      unsigned int buffer_size, stage_buffer_size;
    } cr_data;

    typedef struct {
      const unsigned int *map_to_buf[2], *map_from_buf[2];
      unsigned int buffer_size;
    } allreduce_data;

    typedef void exec_fun(
      void *data, gs_mode mode, unsigned vn, gs_dom dom, gs_op op,
      unsigned transpose, const void *execdata, const struct comm *comm, char *buf);
    typedef void fin_fun(void *data);

    typedef struct {
        unsigned int buffer_size, mem_size;
        void *data;
        exec_fun *exec;
        fin_fun *fin;
    } gs_remote;

    struct gs_data {
      struct comm comm;
      const unsigned int *map_local[2]; /* 0=unflagged, 1=all */
      const unsigned int *flagged_primaries;
      gs_remote r;
      unsigned int handle_size;
    };

    typedef enum {gs_auto, gs_pairwise, gs_crystal_router, gs_all_reduce} gs_method;

    extern "C"
    {
        void nektar_gs(void *u, gs_dom dom, gs_op op, unsigned transpose,
                gs_data *gsh, buffer *buf);
        gs_data *nektar_gs_setup(const long *id, unsigned int n, const struct comm *comm,
                                int unique, gs_method method, int verbose);
        void nektar_gs_free(gs_data *gsh);
        void nektar_gs_unique(long *id, unsigned int n, const struct comm *comm);
    }

    // The functions Init, Unique, Finalise and Gather are now part of
    // the Comm objects (in this case, in the CommMpi class).

}

#endif

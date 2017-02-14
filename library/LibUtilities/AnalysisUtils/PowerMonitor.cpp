/*

 Copyright (c) 2014 Harvey Richardson, Michael Bareford
 All rights reserved.

 See the LICENSE file elsewhere in this distribution for the
 terms under which this software is made available.

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <mpi.h>
#include <LibUtilities/AnalysisUtils/PowerMonitor.h>
  

namespace Nektar
{	
namespace LibUtilities
{

// initialise static constant attributes of PowerMonitor class
const char PowerMonitor::ver[] = "2.0.0";

const unsigned int PowerMonitor::MAX_FPATH_LEN = 128;
const unsigned int PowerMonitor::MAX_FLINE_LEN = 128;

const char PowerMonitor::sys_pm_cnt_dir[] = "/sys/cray/pm_counters/";
const char* PowerMonitor::cnt_fname[] = {
    "freshness",
    "power",
    "energy",
    "accel_power",
    "accel_energy",
    "startup",
    "power_cap",
    "accel_power_cap"
};
		
// initialise static attributes of PowerMonitor class
int PowerMonitor::rank(-1);
int PowerMonitor::min_node_rank(0);
int PowerMonitor::monitor_cnt(0);
int PowerMonitor::non_monitor_cnt(0);
bool PowerMonitor::first_record(false);
int PowerMonitor::mpi_comm_monitor(0);

FILE* PowerMonitor::cnt_fp[];
FILE* PowerMonitor::log_fp(NULL);
double PowerMonitor::tm0(0.0);
double PowerMonitor::entot0(0.0);
int PowerMonitor::last_nstep(0);
long int PowerMonitor::init_startup(0);

bool PowerMonitor::all_initialised(false);
        
        
// private methods
//////////////////////////////////////////////////////////////////////////////////
            
void PowerMonitor::OpenCounterFiles(void)
{
    char cnt_fpath[MAX_FPATH_LEN];
    unsigned int max_fname_len = MAX_FPATH_LEN - strlen(sys_pm_cnt_dir);

    strcpy(cnt_fpath, sys_pm_cnt_dir);
    strncat(cnt_fpath, cnt_fname[PM_COUNTER_FRESHNESS], max_fname_len);
    		
    // always open the freshness counter file first
    cnt_fp[PM_COUNTER_FRESHNESS] = fopen(cnt_fpath, "r");
    if (NULL == cnt_fp[PM_COUNTER_FRESHNESS])
    {
        fprintf(stderr, "PowerMonitor: failed to open %s!\n", cnt_fpath);
    }
    else
    {  	
  	// if the freshness counter has been opened successfully
  	// attempt to open the other counter files		
  	for (int i = 0; i < PM_NCOUNTERS; i++)
	{		
  	    if (PM_COUNTER_FRESHNESS == i)
	    {
  	        continue;
  	    }
  					
  	    strcpy(cnt_fpath, sys_pm_cnt_dir);
    	    strncat(cnt_fpath, cnt_fname[i], max_fname_len);

	    cnt_fp[i] = fopen(cnt_fpath, "r");
    	    if (NULL == cnt_fp[i] && !IsAcceleratorCounter(i))
	    {
      	        fprintf(stderr, "PowerMonitor: failed to open %s!\n", cnt_fpath);
      	    }			
  	}
    }
}

   	
void PowerMonitor::CloseCounterFiles(void)
{
    for (int i = 0; i < PM_NCOUNTERS; i++)
    {
        if (NULL != cnt_fp[i])
	{
      	    fclose(cnt_fp[i]);
      	    cnt_fp[i] = NULL;
    	}
    }
}
		
bool PowerMonitor::IsAcceleratorCounter(const unsigned int i)
{
    return (i == PM_COUNTER_ACCEL_POWER ||
      	    i == PM_COUNTER_ACCEL_ENERGY ||
  	    i == PM_COUNTER_ACCEL_POWER_CAP);
}
		
// return the first line from the counter file identified by i
void PowerMonitor::GetFirstLine(const unsigned int i, char* line, const unsigned int len)
{
    if (NULL != line)
    {
        memset(line, 0, len);
	if (i < PM_NCOUNTERS && NULL != cnt_fp[i])
	{
  	    rewind(cnt_fp[i]);
  	    while (NULL != fgets(line, len, cnt_fp[i]) || EAGAIN == errno);	
  	}
  	else
	{
  	    strcpy(line, "0");
  	}
    }
}
		
long int PowerMonitor::GetCounterValue(const unsigned int i)
{
    char line[MAX_FLINE_LEN];
    PowerMonitor::GetFirstLine(i, line, MAX_FLINE_LEN);
    return strtol(line, NULL, 10);
}

// determine the number of the node on which the process is running
int PowerMonitor::GetNodeNumber(void)
{
    int node_name_len, nn_i, nn_m, node_num(0);
    char node_name[MPI_MAX_PROCESSOR_NAME];
    
    MPI_Get_processor_name(node_name, &node_name_len);
    if (node_name_len > 0)
    {
        nn_i = node_name_len-1;
    	nn_m = 1;
    	node_num = 0;
    	while (nn_i > 0 && 0 != isdigit(node_name[nn_i]))
	{
      	    node_num = node_num + (node_name[nn_i]-'0')*nn_m;
      	    nn_m = 10*nn_m;
      	    nn_i = nn_i - 1;
    	}
    }

    return node_num;
}
		

// return true if PowerMonitor::Initialise has been called successfully
bool PowerMonitor::IsInitialised(void)
{
    bool ok = false;
  
    if (-1 != rank)
    {
        if (min_node_rank == rank)
	{
      	    ok = (monitor_cnt > 0);
      	    ok = (ok && (NULL != cnt_fp[PM_COUNTER_FRESHNESS]));
      	    ok = (ok && (NULL != cnt_fp[PM_COUNTER_POWER]));
      	    ok = (ok && (NULL != cnt_fp[PM_COUNTER_ENERGY]));
      				
      	    if (0 == rank)
	    {
                ok = (ok && (NULL != log_fp));
      	    }
    	}
    	else
	{
      	    ok = (non_monitor_cnt > 0);
    	}
    }

    return ok;
}


// public methods
//////////////////////////////////////////////////////////////////////////////////

// rank zero opens the output file
// allow the calling rank to self-identify as a monitoring process
// each monitoring process obtains the number of monitors
// call PowerMonitor::Record(-1,1)
void PowerMonitor::Initialise(const char* log_fpath)
{  
    int node_num = 0;
    MPI_Comm mpi_comm_node;
  
  
    if (PowerMonitor::IsInitialised())
    {
        return;
    }
  
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (0 == rank)
    {
        if (NULL != log_fp)
	{
      	    fclose(log_fp);
      	    log_fp = NULL;
    	}
    	// open power management counter log file
    	if (NULL == log_fp)
	{
      	    log_fp = fopen(log_fpath, "w"); 
    	}
    	if (NULL == log_fp)
	{
      	    log_fp = fopen("./pm_log.out", "w");
    	}
    }

    // determine the number of the node running this process
    node_num = GetNodeNumber();
      
    // determine if this rank is the minimum rank for the identified node number
    MPI_Comm_split(MPI_COMM_WORLD, node_num, rank, &mpi_comm_node);
    MPI_Allreduce(&rank, &min_node_rank, 1, MPI_INTEGER, MPI_MIN, mpi_comm_node);
 
    if (rank == min_node_rank)
    {
        // the minimum rank on a node is responsible for monitoring
        // the performance counters on that node  
    	MPI_Comm_split(MPI_COMM_WORLD, 1, rank, &mpi_comm_monitor);
    }
    else
    {
        MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &mpi_comm_monitor);
    }
  
    // ensure all monitor ranks have self-identified...
    MPI_Barrier(MPI_COMM_WORLD);
    // ...before each monitor rank obtains the total number of monitors
    if (min_node_rank == rank)
    {
        MPI_Comm_size(mpi_comm_monitor, &monitor_cnt);
    }
    else
    {
        // and other processes obtain the number of non-monitors
    	MPI_Comm_size(mpi_comm_monitor, &non_monitor_cnt);
    }
  
    if (min_node_rank == rank)
    {
        PowerMonitor::OpenCounterFiles();
	init_startup = PowerMonitor::GetCounterValue(PM_COUNTER_STARTUP);
    }
  
    bool ok = PowerMonitor::IsInitialised(), all_ok = false;
    MPI_Allreduce(&ok, &all_ok, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD);
    all_initialised = all_ok;
    if (all_initialised)
    {
        // do initial record, which ends with MPI_Barrier
        first_record = true;
        PowerMonitor::Record(-1, 1);
    }
    else
    {
        PowerMonitor::Finalise();
    }
  			
} // end of PowerMonitor::Initialise method


// read counter values if first rank on node,
// and output those values if rank zero
void PowerMonitor::Record(const int nstep, const int sstep) {
   
    if (!all_initialised)
    {
        return;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    if (min_node_rank == rank)
    {
        long int start_freshness, end_freshness;
    	double pmc_energy, tot_pmc_energy;
    	long int pmc_power, tot_pmc_power;
    
    	// get time
    	double tm = MPI_Wtime();
    	if (first_record)
	{
	    tm0 = tm;
      	    first_record = false;
    	}
    
    	// read the point-in-time power and accumulated energy counters
    	bool fresh(false);
    	while (!fresh)
	{
      	    start_freshness = PowerMonitor::GetCounterValue(PM_COUNTER_FRESHNESS);
      	    pmc_power = PowerMonitor::GetCounterValue(PM_COUNTER_POWER);
      	    pmc_energy = PowerMonitor::GetCounterValue(PM_COUNTER_ENERGY);
      	    end_freshness = PowerMonitor::GetCounterValue(PM_COUNTER_FRESHNESS);
      	    fresh = (end_freshness == start_freshness);
    	}
  
    	MPI_Reduce(&pmc_power, &tot_pmc_power, 1, MPI_LONG, MPI_SUM, 0, mpi_comm_monitor);
    	MPI_Reduce(&pmc_energy, &tot_pmc_energy, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_monitor);
         		
    	// output data
    	if (0 == rank)
	{
      	    if (tm0 == tm)
	    {
                // this function is being called by PowerMonitor::Initialise
        	entot0 = tot_pmc_energy;
        
        	if (NULL != log_fp)
		{
          	    fprintf(log_fp, "PowerMonitor v%s: time (s), step, substep, average power (W), energy used (J)\n", ver);
        	}
      	    }

	    tot_pmc_energy = round(tot_pmc_energy);
      	    double avg_pmc_power = (monitor_cnt > 0) ? ((double) tot_pmc_power)/((double) monitor_cnt) : 0.0;
      	    double dif_pmc_energy = tot_pmc_energy - entot0;
        
      	    if (NULL != log_fp)
	    {   
                // update counter data file   
        	fprintf(log_fp, "%f %d %d %f %f\n", tm-tm0, nstep, sstep, avg_pmc_power, dif_pmc_energy); 
      	    }
    	}
    } // end of <if (min_node_rank == rank)> clause
  
    last_nstep = nstep;
  
    MPI_Barrier(MPI_COMM_WORLD);
  
} // end of PowerMonitor::Record method


// close the files used to read and record counter data
void PowerMonitor::Finalise(void)
{
    if (all_initialised)
    {
    	// do the last record
    	PowerMonitor::Record(last_nstep+1, 1);
    }
  
    // if monitoring process (i.e., first process on node)    
    if (min_node_rank == rank)
    {
        long int final_startup = PowerMonitor::GetCounterValue(PM_COUNTER_STARTUP);
	if (final_startup != init_startup) {
	    fprintf(stderr, "PowerMonitor meaurements invalid! Blade-controller was restarted for node %d.\n", GetNodeNumber());
	}
	
        PowerMonitor::CloseCounterFiles();
        
        if (0 == rank && NULL != log_fp)
	{
    	    // close performance counter data file
            fclose(log_fp);
            log_fp = NULL;
    	}
    }	
  
    all_initialised = false;
    MPI_Barrier(MPI_COMM_WORLD);
  
} // end of PowerMonitor::Finalise method

}
}

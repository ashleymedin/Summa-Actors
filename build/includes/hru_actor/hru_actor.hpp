#pragma once

#include "caf/all.hpp"
#include "fortran_data_types.hpp"
#include "auxilary.hpp"
#include "timing_info.hpp"
#include "settings_functions.hpp"
#include <string>
#include "message_atoms.hpp"
#include "global.hpp"


struct Date {
  int y;
  int m;
  int d;
  int h;
};


/*********************************************
 * HRU Actor Fortran Functions
 *********************************************/
extern "C" {
  // Initialize HRU data_structures
  void initHRU(int* indxGRU, int* num_steps, void* hru_data, int* err);
  
  void setupHRUParam(int* indxGRU, int* indxHRU, void* hru_data, 
                     double* upArea, int* err);
  
  // Setup summa_readRestart File if this option has been chosen 
  void summa_readRestart(int* indxGRU, int* indxHRU, void* hru_data, 
                         double* dtInit, int* err);
  
  // Run the model for one timestep
  void RunPhysics(int* id, int* stepIndex, void* hru_data, double* dt, 
                  int* dt_int_factor, double* walltime_timestep, int* err);
  
  void hru_writeOutput(int* index_hru, int* index_gru, int* timestep, 
                       int* output_step, void* hru_data, int* y, int* m, 
                       int* d, int* h, int* err);
    
  int  hru_writeRestart(int* index_hru, int* index_gru, int* timestep, 
                        int* output_step, void* hru_data, int* err);

  void setTimeZoneOffset(int* iFile, void* hru_data, int* err);

  void HRU_readForcing(int* index_gru, int* iStep, int* iRead, int* iFile, 
                       void* hru_data,  int* err);
  
  // hru_writeOutput.f90
  void setFinalizeStatsFalse(int* indx_gru);

  void get_sundials_tolerances(void* hru_data, double* relTol, double* absTol);
  void set_sundials_tolerances(void* hru_data, double* relTol, double* absTol);

  void setIDATolerances(void* hru_data, double* relTolTempCas, 
                        double* absTolTempCas, double* relTolTempVeg, 
                        double* absTolTempVeg, double* relTolWatVeg, 
                        double* absTolWatVeg, double* relTolTempSoilSnow,
                        double* absTolTempSoilSnow, double* relTolWatSnow, 
                        double* absTolWatSnow, double* relTolMatric,
                        double* absTolMatric, double* relTolAquifr, 
                        double* absTolAquifr);
}

/*********************************************
 * HRU Actor state variables
 *********************************************/
struct hru_state {
  // Actor References
	caf::actor file_access_actor;
	caf::actor parent;

  int ref_gru;			  // gru_index from attributes/forcing files
  int indxGRU; 		  // index for gru part of derived types in FORTRAN
  int indxHRU; 	    // index for hru part of derived types in FORTRAN

    // Variables for forcing structures
	int stepsInCurrentFFile;        // number of time steps in current forcing file
  int num_steps_until_write;      // number of time steps until we pause for FA_Actor to write

  // HRU data structures (formerly summa_type)
  void *hru_data = new_handle_hru_type();

  // Serializable HRU data structure
  hru hru_data_serialized;

    // Misc Variables
  int     timestep = 1;	    // Current Timestep of HRU simulation
  int     forcingStep = 1;    // index of current time step in current forcing file
  int     num_steps = 0;      // number of time steps
  int     iFile = 1;              // index of current forcing file from forcing file list
  int     dt_init_factor = 1; // factor of dt_init (coupled_em)
  int     output_structure_step_index = 1; // index of current time step in output structure
  double  dt_init;            // used to initialize the length of the sub-step for each HRU
  double	upArea;             // area upslope of each HRU

  // Sundials variables
  double rtol = -9999; // -9999 uses default
  double atol = -9999; // -9999 uses default

  double walltime_timestep = 0.0; // walltime for the current timestep		

  // Checkpointing variables
  Date startDate =   {0,0,0,0}; // will be initalized when hru finishes first timestep
  Date currentDate = {0,0,0,0}; // will be initalized when hru finishes first timestep
  int  checkpoint= 0;
  int daysInMonth[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
  int restartFrequency; // 0=never(default) 1=hour 2=day 3=week? 4=month 5=year 
        

  // Settings
  HRU_Actor_Settings hru_actor_settings;
    
  ~hru_state() {
      delete_handle_hru_type(hru_data);
  }
};

caf::behavior hru_actor(caf::stateful_actor<hru_state>* self, int ref_gru, 
                        int indx_gru, HRU_Actor_Settings hru_actor_settings, 
                        caf::actor file_access_actor, caf::actor parent);

/*********************************************
 * Functions for the HRU Actor
 *********************************************/

/** Function to initalize the HRU for running */
void Initialize_HRU(caf::stateful_actor<hru_state>* self);

/** Function runs all of the hru time_steps */
int Run_HRU(caf::stateful_actor<hru_state>* self);

/** Function checks if the HRU is at a restart checkpoint */
bool isCheckpoint(caf::stateful_actor<hru_state>* self);

void serializeHru(caf::stateful_actor<hru_state>* self, hru& serialized_state);
void deserializeHru(caf::stateful_actor<hru_state>* self, hru& new_state);

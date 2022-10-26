#pragma once

#include "caf/all.hpp"
#include "output_manager.hpp"
#include "forcing_file_info.hpp"
#include "timing_info.hpp"

namespace caf {
struct file_access_state {
    // Variables set on Spwan
    caf::actor parent; 
    int startGRU;
    int numGRU;


    void *handle_forcing_file_info; // Handle for the forcing file information
    void *handle_ncid;              // output file ids
    OutputManager *output_manager;
    int num_vectors_in_output_manager;
    int num_steps;
    int outputStrucSize;
    int stepsInCurrentFile;
    int numFiles;
    int filesLoaded;
    int err;

    std::vector<Forcing_File_Info> forcing_file_list; // list of steps in file
    std::vector<bool> outputFileInitHRU;


    // vector of handles for attributes
    std::vector<void *> type_struct_handles;
    std::vector<void *> id_struct_handles;
    std::vector<void *> attr_struct_handles;

    // vector of handles for parameters
    std::vector<void *> mpar_struct_handles;
    std::vector<void *> bpar_struct_handles;


     // Timing Variables
    TimingInfo file_access_timing;
};

behavior file_access_actor(stateful_actor<file_access_state>* self, int startGRU, int numGRU, 
    int outputStrucSize, std::string configPath, actor parent);
void initalizeFileAccessActor(stateful_actor<file_access_state>* self);
int writeOutput(stateful_actor<file_access_state>* self, int indxGRU, int indxHRU, int numStepsToWrite, int returnMessage, caf::actor actorRef);
int readForcing(stateful_actor<file_access_state>* self, int currentFile);
int write(stateful_actor<file_access_state>* self, int listIndex);

// Read in the attributes for all HRUs that are in the run-domain
void readAttributes(stateful_actor<file_access_state>* self); 


} // end namespace
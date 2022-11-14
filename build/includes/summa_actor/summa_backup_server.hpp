#pragma once

#include "caf/all.hpp"
#include "caf/io/all.hpp"
#include "summa_server.hpp"
#include <string>
namespace caf {


// Inital behaviour that waits to connect to the lead server
behavior summa_backup_server_init(stateful_actor<summa_server_state>* self, Distributed_Settings distributed_settings, 
    Summa_Actor_Settings summa_actor_settings, File_Access_Actor_Settings file_access_actor_settings,
    Job_Actor_Settings job_actor_settings, HRU_Actor_Settings hru_actor_settings);

// Function that is called ot connect to the lead server
void connecting(stateful_actor<summa_server_state>* self, const std::string& host, uint16_t port);

// The behaviour of the backup server
behavior summa_backup_server(stateful_actor<summa_server_state>* self, const actor& summa_server);
}
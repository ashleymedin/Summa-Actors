#pragma once

#include "caf/all.hpp"
#include "caf/io/all.hpp"
#include "batch_manager.hpp"
#include <string>


namespace caf {

struct summa_server_state {
    int total_hru_count;
    int num_clients;
    int num_hru_per_batch;   
    std::vector<Batch*> batch_list;
    std::vector<Batch> solved_batches;
    std::vector<Batch> failed_batches;
};

behavior summa_server(stateful_actor<summa_server_state>* self, std::string config_path);
int parseSettings(stateful_actor<summa_server_state>* self, std::string config_path);
int assembleBatches(stateful_actor<summa_server_state>* self);
}
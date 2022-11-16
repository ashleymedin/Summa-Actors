#include "caf/all.hpp"
#include "caf/io/all.hpp"
#include "message_atoms.hpp"
#include "summa_backup_server.hpp"
#include "global.hpp"
#include <chrono>
#include <thread>
namespace caf {

behavior summa_backup_server_init(stateful_actor<summa_server_state>* self, Distributed_Settings distributed_settings, 
    Summa_Actor_Settings summa_actor_settings, File_Access_Actor_Settings file_access_actor_settings,
    Job_Actor_Settings job_actor_settings, HRU_Actor_Settings hru_actor_settings) {
    aout(self) << "Backup Server Started\n";
    char host[HOST_NAME_MAX];
    gethostname(host, HOST_NAME_MAX);
    self->state.hostname = host;
    
    self->state.distributed_settings = distributed_settings;
    self->state.summa_actor_settings = summa_actor_settings; 
    self->state.file_access_actor_settings = file_access_actor_settings;
    self->state.job_actor_settings = job_actor_settings;
    self->state.hru_actor_settings = hru_actor_settings;

    self->set_down_handler([=](const down_msg& dm){
        if(dm.source == self->state.current_server) {
            // aout(self) << "*** Lost Connection to Server" << std::endl;
            // uint16_t port = 4444;
            // std::string host = "a0449745d77d";
            // connecting(self, host, port);
        }
    });
    return {
        // Called by main to init the process
        [=](connect_atom, const std::string& host, uint16_t port) {
            connecting_backup(self, host, port);
        },
    };
}


void connecting_backup(stateful_actor<summa_server_state>* self, const std::string& host, uint16_t port) {
    self->state.current_server = nullptr;

    auto mm = self->system().middleman().actor_handle();
    self->request(mm, infinite, connect_atom_v, host, port)
        .await(
            [=](const node_id&, strong_actor_ptr serv,
                const std::set<std::string>& ifs) {
                if (!serv) {
                    aout(self) << R"(*** no server found at ")" << host << R"(":)" << port
                     << std::endl;
                    return;
                }
                if (!ifs.empty()) {
                    aout(self) << R"(*** typed actor found at ")" << host << R"(":)"
                        << port << ", but expected an untyped actor " << std::endl;
                    return;
                }
                aout(self) << "*** successfully connected to server" << std::endl;
                self->state.current_server = serv;
                auto hdl = actor_cast<actor>(serv);
                self->state.current_server_actor = hdl;
                self->monitor(hdl);
                self->become(summa_backup_server(self, hdl));

                },
            [=](const error& err) {
                aout(self) << R"(*** cannot connect to ")" << host << R"(":)" << port
                   << " => " << to_string(err) << std::endl;
        });
}

behavior summa_backup_server(stateful_actor<summa_server_state>* self, const actor& server_actor) {
    aout(self) << "We are the test behaviour\n";
    self->send(server_actor, connect_as_backup_v, self, self->state.hostname);

    return {

        [=](connect_as_backup) {
            aout(self) << "We are now connected to the lead server\n";
        },

        [=](update_with_current_state, Batch_Container& batch_container, Client_Container& client_container) {
            aout(self) << "Received the containers from the lead server\n";
            self->state.batch_container = &batch_container;
            self->state.client_container = &client_container;
        },

        // Client finished a batch and the lead server has sent an update
        [=](done_batch, actor client_actor, Batch& batch) {
            aout(self) << "Batch: " << batch.getBatchID() << " is done\n";
            self->state.batch_container->updateBatch_success(batch);
        },

        // Client has been assigned new batch by the lead server
        [=](new_assigned_batch, actor client_actor, Batch& batch) {
            aout(self) << "New Batch: " << batch.getBatchID() << " has been assigned\n";
            self->state.batch_container->setBatchAssigned(batch);
            self->state.client_container->setBatchForClient(client_actor, batch);
        },

        // Lead server has no more batches to distribute
        [=](no_more_batches, actor client_actor) {
            aout(self) << "No more batches to distribute\n";
            self->state.client_container->setBatchForClient(client_actor, {});

        },

        // Simulation has finished
        [=](time_to_exit) {
            aout(self) << "Received time to exit\n";
            self->quit();
        },
    };
}
}
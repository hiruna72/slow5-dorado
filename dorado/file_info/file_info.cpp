#include "file_info/file_info.h"

#include "utils/PostCondition.h"
#include "utils/fs_utils.h"
#include "utils/time_utils.h"

#include <pod5_format/c_api.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <iterator>
#include <set>
#include <stdexcept>
#include <string_view>
#include "slow5/slow5.h"

namespace {
using namespace std::literals;
constexpr std::string_view PROTOCOL_GROUP_ID_KEY = "protocol_group_id"sv;
}  // namespace

namespace dorado::file_info {

std::unordered_map<std::string, ReadGroup> load_read_groups(
        const std::vector<std::filesystem::directory_entry>& dir_files,
        int model_stride,
        const std::string& model_name,
        const std::string& modbase_model_names,
        TrimFlags trim_flags) {
    if (pod5_init() != POD5_OK) {
        throw std::runtime_error(
                fmt::format("Failed to initialise POD5: {}", pod5_get_error_string()));
    }

    std::unordered_map<std::string, ReadGroup> read_groups;
    for (const auto& entry : dir_files) {
        if (utils::has_pod5_extension(entry)) {
            // Open the file
            const auto file_path = entry.path().string();
            Pod5FileReader_t* file = pod5_open_file(file_path.c_str());
            if (!file) {
                spdlog::error("Failed to open file {}: {}", file_path, pod5_get_error_string());
                continue;
            }
            auto cleanup_file = utils::PostCondition([&file] {
                if (pod5_close_and_free_reader(file) != POD5_OK) {
                    spdlog::error("Failed to close and free POD5 reader: {}", pod5_get_error_string());
                }
            });

            // First get the run info count
            run_info_index_t run_info_count = 0;
            if (pod5_get_file_run_info_count(file, &run_info_count) != POD5_OK) {
                spdlog::error("Failed to query run info count: {}", pod5_get_error_string());
                continue;
            }

            for (run_info_index_t idx = 0; idx < run_info_count; idx++) {
                RunInfoDictData_t* run_info_data = nullptr;
                if (pod5_get_file_run_info(file, idx, &run_info_data) != POD5_OK) {
                    continue;
                }
                auto cleanup_run_info = utils::PostCondition([&run_info_data] {
                    if (pod5_free_run_info(run_info_data) != POD5_OK) {
                        spdlog::error("Failed to free run info: {}", pod5_get_error_string());
                    }
                });

                auto exp_start_time_ms = run_info_data->protocol_start_time_ms;
                auto acq_start_time_ms = run_info_data->acquisition_start_time_ms;
                std::string flowcell_id = run_info_data->flow_cell_id;
                std::string device_id = run_info_data->system_name;
                std::string run_id = run_info_data->protocol_run_id;
                std::string sample_id = run_info_data->sample_id;
                std::string position_id = run_info_data->sequencer_position;
                std::string experiment_id = run_info_data->experiment_name;

                if (experiment_id.empty()) {
                    for (size_t i = 0; i < run_info_data->tracking_id.size; ++i) {
                        if (run_info_data->tracking_id.keys[i] == PROTOCOL_GROUP_ID_KEY) {
                            experiment_id = run_info_data->tracking_id.values[i];
                            break;
                        }
                    }
                }

                std::string id = std::string(run_id).append("_").append(model_name);
                read_groups[id] = ReadGroup{
                        .run_id = std::move(run_id),
                        .basecalling_model = model_name,
                        .modbase_models = modbase_model_names,
                        .flowcell_id = std::move(flowcell_id),
                        .device_id = std::move(device_id),
                        .exp_start_time =
                                utils::get_string_timestamp_from_unix_time_ms(exp_start_time_ms),
                        .sample_id = std::move(sample_id),
                        .position_id = std::move(position_id),
                        .experiment_id = std::move(experiment_id),
                        .acq_start_time =
                                utils::get_string_timestamp_from_unix_time_ms(acq_start_time_ms),
                        .barcode_id = {},
                        .barcode_alias = {},
                        .model_stride = model_stride,
                        .trim_flags = trim_flags,
                };
            }
        } else if (utils::has_blow5_extension(entry)) {
            slow5_file_t *sp = slow5_open(entry.path().string().c_str(),"r");
            if(sp==NULL){
                spdlog::error("Error opening SLOW5/BLOW5 file {}", entry.path().string());
                continue;
            }
            int64_t read_group_count = sp->header->num_read_groups;
            for(int64_t j=0; j<read_group_count; j++){
                char* protocol_run_id_c = slow5_hdr_get("protocol_run_id", j, sp->header);
                std::string protocol_run_id = "";
                if(protocol_run_id_c){
                    protocol_run_id = std::string(protocol_run_id_c);
                } else{
                    spdlog::warn("protocol_run_id not found in read group {}", j);
                }

                char* flow_cell_id_c = slow5_hdr_get("flow_cell_id", j, sp->header);
                std::string flowcell_id = "";
                if(flow_cell_id_c){
                    flowcell_id = std::string(flow_cell_id_c);
                } else{
                    spdlog::warn("flow_cell_id not found in read group {}", j);
                }
                
                char* device_id_c = slow5_hdr_get("system_name", j, sp->header); //pod5 has started to use confusing variable name. device_id and sequencer_position refer to the same in the blow5 and pod5 file specifications.
                std::string device_id = "";
                if(!device_id_c){
                    device_id_c = slow5_hdr_get("host_product_serial_number", j, sp->header);
                }
                if(device_id_c){
                    device_id = std::string(device_id_c);
                } else{
                    spdlog::warn("device_id (system_name or host_product_serial_number) not found in read group {}", j);
                }

                char* experiment_id_c = slow5_hdr_get("experiment_name", j, sp->header);
                std::string experiment_id = "";
                if(experiment_id_c){
                    experiment_id = std::string(experiment_id_c);
                } else{
                    spdlog::warn("experiment_id not found in read group {}", j);
                }

                char* exp_start_time_ms_c = slow5_hdr_get("acquisition_start_time", j, sp->header);
                std::string exp_start_time_ms_str = "";
                if(!exp_start_time_ms_c){
                    exp_start_time_ms_c = slow5_hdr_get("exp_start_time", j, sp->header);
                }
                if(exp_start_time_ms_c) {
                    exp_start_time_ms_str = std::string(exp_start_time_ms_c);
                } else{
                    spdlog::warn("Neither acquisition_start_time nor exp_start_time found in read group {}", j);
                }
                // std::cerr << "exp_start_time_ms: " << exp_start_time_ms << std::endl;
                std::string normalized_acq_start_time;
                if (!exp_start_time_ms_str.empty()) {
                    try {
                        auto ms = utils::get_unix_time_ms_from_string_timestamp(exp_start_time_ms_str);
                        normalized_acq_start_time = utils::get_string_timestamp_from_unix_time_ms(ms);
                    } catch (const std::exception& e) {
                        spdlog::warn("Failed to parse start time '{}' in {}: {}", exp_start_time_ms_str, entry.path().string(), e.what());
                        normalized_acq_start_time.clear(); // leave empty if unparsable
                    }
                }


                char* protocol_start_time_c = slow5_hdr_get("protocol_start_time", j, sp->header);
                std::string protocol_start_time = "";
                if(!protocol_start_time_c){
                    spdlog::error("Failed to get protocol_start_time in {}", sp->meta.pathname);
                } else{
                    protocol_start_time = std::string(protocol_start_time_c);
                }
                auto protocol_start_time_ms = utils::get_unix_time_ms_from_string_timestamp(protocol_start_time);
                std::string normalized_protocol_start_time;
                if (!protocol_start_time.empty()) {
                    try {
                        normalized_protocol_start_time = utils::get_string_timestamp_from_unix_time_ms(protocol_start_time_ms);
                    } catch (const std::exception& e) {
                        spdlog::warn("Failed to parse protocol start time '{}' in {}: {}", protocol_start_time, entry.path().string(), e.what());
                        normalized_protocol_start_time.clear(); // leave empty if unparsable
                    }
                }

                char* sample_id_c = slow5_hdr_get("sample_id", j, sp->header);
                std::string sample_id = "";
                if(sample_id_c){
                    sample_id = std::string(sample_id_c);
                } else{
                    spdlog::warn("sample_id not found in read group {}", j);
                }

                char* position_id_c = slow5_hdr_get("sequencer_position", j, sp->header);
                std::string position_id = "";
                if(position_id_c){
                    position_id = std::string(position_id_c);
                } else{
                    spdlog::warn("position_id not found in read group {}", j);
                }

                std::string run_id = protocol_run_id;
                std::string id = std::string(run_id).append("_").append(model_name);
                read_groups[id] = ReadGroup{
                        .run_id = std::move(run_id),
                        .basecalling_model = model_name,
                        .modbase_models = modbase_model_names,
                        .flowcell_id = std::move(flowcell_id),
                        .device_id = std::move(device_id),
                        .exp_start_time = std::move(normalized_protocol_start_time),
                        .sample_id = std::move(sample_id),
                        .position_id = std::move(position_id),
                        .experiment_id = std::move(experiment_id),
                        .acq_start_time = std::move(normalized_acq_start_time),
                        .barcode_id = {},
                        .barcode_alias = {},
                        .model_stride = model_stride,
                        .trim_flags = trim_flags,
                };
            }
            slow5_close(sp);
        }
        
    }

    return read_groups;
}

size_t get_num_reads(const std::vector<std::filesystem::directory_entry>& dir_files,
                     std::optional<std::unordered_set<std::string>> read_list,
                     const std::unordered_set<std::string>& ignore_read_list) {
    return 1;
    if (pod5_init() != POD5_OK) {
        throw std::runtime_error(
                fmt::format("Failed to initialise POD5: {}", pod5_get_error_string()));
    }

    size_t num_reads = 0;
    for (const auto& entry : dir_files) {
        const auto ext = utils::get_extension(entry);

        if (ext == ".fast5") {
            throw std::runtime_error("FAST5 is not supported");
        } else if (ext != ".pod5") {
            continue;
        }

        // Open the file
        const auto file_path = entry.path().string();
        Pod5FileReader_t* file = pod5_open_file(file_path.c_str());
        if (!file) {
            spdlog::error("Failed to open file {}: {}", file_path, pod5_get_error_string());
            continue;
        }
        auto cleanup_file = utils::PostCondition([&file] {
            if (pod5_close_and_free_reader(file) != POD5_OK) {
                spdlog::error("Failed to close and free POD5 reader: {}", pod5_get_error_string());
            }
        });

        size_t read_count = 0;
        if (pod5_get_read_count(file, &read_count) != POD5_OK) {
            spdlog::error("Failed to query read count of {}: {}", file_path,
                          pod5_get_error_string());
            continue;
        }

        num_reads += read_count;
    }

    // Remove the reads in the ignore list from the total dataset read count.
    num_reads -= ignore_read_list.size();

    if (read_list) {
        // Get the unique read ids in the read list, since everything in the ignore
        // list will be skipped over.
        std::vector<std::string> final_read_list;
        std::set_difference(read_list->begin(), read_list->end(), ignore_read_list.begin(),
                            ignore_read_list.end(),
                            std::inserter(final_read_list, final_read_list.begin()));
        num_reads = std::min(num_reads, final_read_list.size());
    }

    return num_reads;
}

bool is_pod5_data_present(const std::vector<std::filesystem::directory_entry>& dir_files) {
    for (const auto& entry : dir_files) {
        if (utils::has_pod5_extension(entry) || utils::has_blow5_extension(entry)) {
            return true;
        }
    }
    return false;
}

uint16_t get_sample_rate(const std::vector<std::filesystem::directory_entry>& dir_files) {

    if (pod5_init() != POD5_OK) {
        throw std::runtime_error(
                fmt::format("Failed to initialise POD5: {}", pod5_get_error_string()));
    }


    for (const auto& entry : dir_files) {
        auto file_path = entry.path().string();
        if (utils::has_pod5_extension(entry)) {
            // Open the file
            Pod5FileReader_t* file = pod5_open_file(file_path.c_str());
            if (!file) {
                spdlog::error("Failed to open file {}: {}", file_path, pod5_get_error_string());
                continue;
            }
            auto cleanup_file = utils::PostCondition([&file, &file_path]() {
                if (pod5_close_and_free_reader(file) != POD5_OK) {
                    spdlog::error("Failed to close and free POD5 reader for file {}: {}", file_path,
                                pod5_get_error_string());
                }
            });

            // First get the run info count
            run_info_index_t run_info_count = 0;
            if (pod5_get_file_run_info_count(file, &run_info_count) != POD5_OK) {
                spdlog::error("Failed to fetch POD5 run info count for file {} : {}", file_path,
                            pod5_get_error_string());
                continue;
            }
            if (run_info_count > static_cast<run_info_index_t>(0)) {
                RunInfoDictData_t* run_info_data = nullptr;
                if (pod5_get_file_run_info(file, 0, &run_info_data) != POD5_OK) {
                    spdlog::error(
                            "Failed to fetch POD5 run info dict for file {} and run info index 0: {}",
                            file_path, pod5_get_error_string());
                    continue;
                }
                auto cleanup_run_info = utils::PostCondition([&run_info_data, &file_path] {
                    if (pod5_free_run_info(run_info_data) != POD5_OK) {
                        spdlog::error(
                                "Failed to free POD5 run info for file {} and run info index 0: {}",
                                file_path, pod5_get_error_string());
                    }
                });

                // Break out of loop if sample rate is found.
                return run_info_data->sample_rate;
            }
        } else if (utils::has_blow5_extension(entry)) {
            slow5_file_t *sp = slow5_open(file_path.c_str(),"r");
            if(sp==NULL){
                spdlog::error("Error opening SLOW5/BLOW5 file {}", file_path);
                continue;
            }
            // Assume sample rate is same across all read groups. So just check the first one. Similar to what is done in pod5.
            int j = 0;
            if (char* sr_c = slow5_hdr_get("sample_frequency", j, sp->header)) {
                try {
                    u_int16_t sample_rate = static_cast<u_int16_t>(std::lround(std::stod(sr_c)));
                    slow5_close(sp);
                    return sample_rate;
                } catch (const std::exception& e) {
                    spdlog::warn("Malformed sample_frequency ('{}') in {} RG {}: {}",
                                sr_c, file_path, j, e.what());
                }
            } else {
                spdlog::error("No sample_frequency found in {} (RG {})", file_path, j);
            }
        }
    }

    throw std::runtime_error("Unable to determine sample rate for data.");

}

static std::set<models::ChemistryKey> get_sequencing_chemistries(
        const std::vector<std::filesystem::directory_entry>& dir_files) {
    if (pod5_init() != POD5_OK) {
        throw std::runtime_error(
                fmt::format("Failed to initialise POD5: {}", pod5_get_error_string()));
    }

    std::set<models::ChemistryKey> chemistries;
    for (const auto& entry : dir_files) {
        const auto file_path = std::filesystem::path(entry).string();
        if (utils::has_pod5_extension(entry)) {

            // Open the file
            Pod5FileReader_t* file = pod5_open_file(file_path.c_str());
            if (!file) {
                spdlog::error("Failed to open file {}: {}", file_path, pod5_get_error_string());
                continue;
            }
            auto cleanup_file = utils::PostCondition([&file, &file_path] {
                if (pod5_close_and_free_reader(file) != POD5_OK) {
                    spdlog::error("Failed to close and free POD5 reader for file {}: {}", file_path,
                                pod5_get_error_string());
                }
            });

            // First get the run info count
            run_info_index_t run_info_count = 0;
            if (pod5_get_file_run_info_count(file, &run_info_count) != POD5_OK) {
                spdlog::error("Failed to fetch POD5 run info count for file {} : {}", file_path,
                            pod5_get_error_string());
                continue;
            }

            for (run_info_index_t ri_idx = 0; ri_idx < run_info_count; ri_idx++) {
                RunInfoDictData_t* run_info_data = nullptr;
                if (pod5_get_file_run_info(file, ri_idx, &run_info_data) != POD5_OK) {
                    spdlog::error(
                            "Failed to fetch POD5 run info dict for file {} and run info index {}: {}",
                            file_path, ri_idx, pod5_get_error_string());
                    continue;
                }
                auto cleanup_run_info = utils::PostCondition([&run_info_data, &file_path] {
                    if (pod5_free_run_info(run_info_data) != POD5_OK) {
                        spdlog::error("Failed to free POD5 run info for file {}: {}", file_path,
                                    pod5_get_error_string());
                    }
                });

                const auto chemistry_key = models::get_chemistry_key(
                        run_info_data->flow_cell_product_code, run_info_data->sequencing_kit,
                        run_info_data->sample_rate);
                spdlog::trace("POD5: {} {}", file_path, to_string(chemistry_key));
                chemistries.insert(chemistry_key);
            }
        } else if (utils::has_blow5_extension(entry)){
            slow5_file_t *sp = slow5_open(file_path.c_str(),"r");
            if(sp==NULL){
                spdlog::error("Error opening SLOW5/BLOW5 file {}", file_path);
                continue;
            }
            int64_t read_group_count = sp->header->num_read_groups;
            for(int64_t j=0; j<read_group_count; j++){
                char* flow_cell_product_code_c = slow5_hdr_get("flow_cell_product_code", j, sp->header);
                std::string flow_cell_product_code = "";
                if(!flow_cell_product_code_c){
                    spdlog::warn("No flow_cell_product_code found in {}. ({})", file_path.c_str(), "DataLoader::get_sequencing_chemistries");
                } else{
                    flow_cell_product_code = std::string(flow_cell_product_code_c);
                }
                char* sequencing_kit_c = slow5_hdr_get("sequencing_kit", j, sp->header);
                std::string sequencing_kit = "";
                if(!sequencing_kit_c){
                    spdlog::warn("No sequencing_kit found in {}. ({})", file_path.c_str(), "DataLoader::get_sequencing_chemistries");
                } else{
                    sequencing_kit = std::string(sequencing_kit_c);
                }

                int sample_rate = 0;
                if (char* sr_c = slow5_hdr_get("sample_frequency", j, sp->header)) {
                    try {
                        sample_rate = static_cast<int>(std::lround(std::stod(sr_c)));
                    } catch (const std::exception& e) {
                        spdlog::warn("Malformed sample_frequency ('{}') in {} RG {}: {}", sr_c, file_path, j, e.what());
                    }
                } else {
                    spdlog::error("No sample_frequency found in {} (RG {})", file_path, j);
                }

                const auto chemistry_key = models::get_chemistry_key(flow_cell_product_code, sequencing_kit, sample_rate);
                spdlog::trace("BLOW5: {} {}", file_path.c_str(), to_string(chemistry_key));
                chemistries.insert(chemistry_key);
            }
            slow5_close(sp);
            continue;            
        }
    }

    return chemistries;
}

models::Chemistry get_unique_sequencing_chemistry(
        const std::vector<std::filesystem::directory_entry>& dir_files) {
    std::set<models::ChemistryKey> data_chemistries = get_sequencing_chemistries(dir_files);

    if (data_chemistries.empty()) {
        throw std::runtime_error(
                "Failed to determine sequencing chemistry from data. Please select a model by "
                "path");
    }

    std::set<models::Chemistry> found;
    for (const auto& dc : data_chemistries) {
        const auto chemistry = models::get_chemistry(dc);
        if (chemistry == models::Chemistry::UNKNOWN) {
            spdlog::error("No supported chemistry found for {}", to_string(dc));
            spdlog::error(
                    "This is typically seen when using prototype kits. Please download an "
                    "appropriate model for your data and select it by model path");

            throw std::runtime_error("Could not resolve chemistry from data: Unknown chemistry");
        }
        found.insert(chemistry);
    }
    if (found.empty()) {
        throw std::runtime_error("Could not resolve chemistry from data: No data");
    }
    if (found.size() > 1) {
        spdlog::error("Multiple sequencing chemistries found in data");
        for (auto f : found) {
            spdlog::error("Found: {}", to_string(f));
        }

        throw std::runtime_error("Could not uniquely resolve chemistry from inhomogeneous data");
    }
    return *std::begin(found);
}

}  // namespace dorado::file_info

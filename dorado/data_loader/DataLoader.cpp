#include "data_loader/DataLoader.h"

#include "models/kits.h"
#include "read_pipeline/base/ReadPipeline.h"
#include "read_pipeline/base/messages/SimplexRead.h"
#include "utils/PostCondition.h"
#include "utils/fs_utils.h"
#include "utils/thread_utils.h"
#include "utils/time_utils.h"
#include "utils/types.h"

#include <ATen/ops/empty.h>
#include <pod5_format/c_api.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <ctime>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <vector>
#include <iostream>

#include "slow5/slow5.h"
#include "slow5_extra.h"
#include "slow5_thread.h"
#include <slow5/slow5_mt.h>

namespace dorado::data_loader {

namespace {
using namespace std::literals;

// 37 = number of bytes in UUID (32 hex digits + 4 dashes + null terminator)
const uint32_t POD5_READ_ID_LEN = 37;

constexpr std::string_view PROTOCOL_GROUP_ID_KEY = "protocol_group_id"sv;

std::vector<std::filesystem::directory_entry> collect_pod5_dataset(
        const std::vector<std::filesystem::directory_entry>& files) {
    std::vector<std::filesystem::directory_entry> pod5_entries;
    bool fast5_found = false;

    for (const auto& entry : files) {
        const auto ext = utils::get_extension(entry);
        if (ext == ".fast5") {
            fast5_found = true;
        } else if (ext == ".pod5") {
            pod5_entries.push_back(entry);
        } else if (ext == ".slow5" || ext == ".blow5") {
            pod5_entries.push_back(entry);
        }
    }

    if (fast5_found && pod5_entries.empty()) {
        spdlog::error(
                "FAST5 support in Dorado was removed in version 1.0.0. "
                "Please convert your dataset to POD5: "
                "https://pod5-file-format.readthedocs.io/en/latest/docs/"
                "tools.html#pod5-convert-fast5");
        throw std::runtime_error("FAST5 files are not supported.");
    }

    if (fast5_found && !pod5_entries.empty()) {
        spdlog::warn(
                "Skipping FAST5 files as support was dropped in Dorado version 1.0.0. "
                "Please convert your FAST5 dataset to POD5: "
                "https://pod5-file-format.readthedocs.io/en/latest/docs/"
                "tools.html#pod5-convert-fast5");
    }
    return pod5_entries;
}

// Parses pod5 run_info data into a ChemistryKey which is used to lookup the sequencing chemistry
models::ChemistryKey get_chemistry_key(const RunInfoDictData_t* const run_info_data) {
    return models::get_chemistry_key(run_info_data->flow_cell_product_code,
                                     run_info_data->sequencing_kit, run_info_data->sample_rate);
}

void issue_pod5_error(std::string_view err, std::string_view filename) {
    // "POD5 Failed to foo the bar - '<pod5 error>' @ '<filename>'."
    spdlog::error("POD5 {} - '{}' @ '{}'.", err, pod5_get_error_string(), filename);
}

void issue_pod5_error(std::string_view err,
                      std::string_view filename,
                      size_t batch_index,
                      size_t row) {
    // "POD5 Failed to foo the bar - '<pod5 error>' @ '<filename>' (<batch_index>,<row>)."
    spdlog::error("POD5 {} - '{}' @ '{}' ({},{}).", err, pod5_get_error_string(), filename,
                  batch_index, row);
}

void issue_pod5_error(std::string_view err, std::string_view filename, std::string_view read_id) {
    // "POD5 Failed to foo the bar - '<pod5 error>' @ '<filename>' [<read_id>]."
    spdlog::error("POD5 {} - '{}' @ '{}' [{}].", err, pod5_get_error_string(), filename, read_id);
}

bool should_process_pod5_row(const ReadBatchRowInfo_t& read_data,
                             const std::string& filename,
                             size_t batch_index,
                             size_t row,
                             const std::optional<std::unordered_set<std::string>>& allowed_read_ids,
                             const std::unordered_set<std::string>& ignored_read_ids) {
    char read_id_tmp[POD5_READ_ID_LEN]{};
    if (pod5_format_read_id(read_data.read_id, read_id_tmp) != POD5_OK) {
        issue_pod5_error("Failed to format read id", filename, batch_index, row);
        return false;
    }

    std::string read_id_str(read_id_tmp);
    bool read_in_ignore_list = ignored_read_ids.find(read_id_str) != ignored_read_ids.end();
    bool read_in_read_list =
            !allowed_read_ids || (allowed_read_ids->find(read_id_str) != allowed_read_ids->end());
    return read_in_read_list && !read_in_ignore_list;
}

SimplexReadPtr process_pod5_thread_fn(
        size_t row,
        size_t batch_index,
        const Pod5ReadRecordBatch* batch,
        const Pod5FileReader* file,
        const std::string& path,
        const std::unordered_map<int, std::vector<DataLoader::ReadSortInfo>>& reads_by_channel,
        const std::unordered_map<std::string, size_t>& read_id_to_index,
        const std::optional<std::unordered_set<std::string>>& allowed_read_ids,
        const std::unordered_set<std::string>& ignored_read_ids) {
    uint16_t read_table_version = 0;

    const std::string filename = std::filesystem::path(path).filename().string();

    ReadBatchRowInfo_t read_data{};
    if (pod5_get_read_batch_row_info_data(batch, row, READ_BATCH_ROW_INFO_VERSION, &read_data,
                                          &read_table_version) != POD5_OK) {
        issue_pod5_error("Failed to get read", filename, batch_index, row);
        return nullptr;
    }

    // Reading ReadBatchRowInfo_t is expensive, so we do the filtering here in the worker thread.
    if (!should_process_pod5_row(read_data, filename, batch_index, row, allowed_read_ids,
                                 ignored_read_ids)) {
        return nullptr;
    }

    //Retrieve global information for the run
    RunInfoDictData_t* run_info_data = nullptr;
    if (pod5_get_run_info(batch, read_data.run_info, &run_info_data) != POD5_OK) {
        issue_pod5_error("Failed to get Run Info", filename, batch_index, row);
        return nullptr;
    }
    auto cleanup = utils::PostCondition([&run_info_data, &filename, batch_index, row] {
        if (pod5_free_run_info(run_info_data) != POD5_OK) {
            issue_pod5_error("Failed to free Run Info", filename, batch_index, row);
        }
    });

    auto run_acquisition_start_time_ms = run_info_data->acquisition_start_time_ms;
    auto run_sample_rate = run_info_data->sample_rate;

    char read_id_tmp[POD5_READ_ID_LEN]{};
    if (pod5_format_read_id(read_data.read_id, read_id_tmp) != POD5_OK) {
        issue_pod5_error("Failed to format read id", filename, batch_index, row);
        return nullptr;
    }
    std::string read_id_str(read_id_tmp);

    auto options = at::TensorOptions().dtype(at::kShort);
    auto samples = at::empty(read_data.num_samples, options);

    if (pod5_get_read_complete_signal(file, batch, row, read_data.num_samples,
                                      samples.data_ptr<int16_t>()) != POD5_OK) {
        issue_pod5_error("Failed to get read signal", filename, read_id_str);
        return nullptr;
    }

    auto new_read = std::make_unique<SimplexRead>();
    new_read->read_common.raw_data = samples;
    new_read->read_common.attributes.sample_rate = run_sample_rate;

    auto start_time_ms = ((read_data.start_sample * 1000) / run_sample_rate);
    auto start_time = utils::get_string_timestamp_from_unix_time_ms(run_acquisition_start_time_ms +
                                                                    start_time_ms);
    new_read->run_acquisition_start_time_ms = run_acquisition_start_time_ms;
    new_read->read_common.start_time_ms = start_time_ms;
    new_read->scaling = read_data.calibration_scale;
    new_read->offset = read_data.calibration_offset;
    new_read->open_pore_level = read_data.open_pore_level;
    new_read->read_common.read_id = std::move(read_id_str);
    new_read->read_common.num_trimmed_samples = 0;
    new_read->read_common.attributes.read_number = read_data.read_number;
    new_read->read_common.attributes.filename = filename;
    new_read->read_common.attributes.mux = read_data.well;
    new_read->read_common.attributes.num_samples = read_data.num_samples;
    new_read->read_common.attributes.channel_number = read_data.channel;
    new_read->read_common.attributes.start_time = start_time;
    new_read->read_common.run_id = run_info_data->protocol_run_id;
    new_read->read_common.acquisition_id = run_info_data->acquisition_id;
    new_read->start_sample = read_data.start_sample;
    new_read->end_sample = read_data.start_sample + read_data.num_samples;
    new_read->read_common.flowcell_id = run_info_data->flow_cell_id;
    new_read->read_common.sequencing_kit = run_info_data->sequencing_kit;
    new_read->read_common.flow_cell_product_code = run_info_data->flow_cell_product_code;
    new_read->read_common.position_id = run_info_data->sequencer_position;
    new_read->read_common.sample_id = run_info_data->sample_id;
    new_read->read_common.protocol_start_time_ms = run_info_data->protocol_start_time_ms;
    new_read->read_common.is_duplex = false;

    // pod5s converted from fast5 don't necessarily have the experiment name
    // if it's missing, check the tracking info instead
    std::string experiment_name = run_info_data->experiment_name;
    if (experiment_name.empty()) {
        for (size_t i = 0; i < run_info_data->tracking_id.size; ++i) {
            if (run_info_data->tracking_id.keys[i] == PROTOCOL_GROUP_ID_KEY) {
                experiment_name = run_info_data->tracking_id.values[i];
                break;
            }
        }
    }

    new_read->read_common.experiment_id = std::move(experiment_name);
    new_read->read_common.num_minknow_events = read_data.num_minknow_events;

    // Get the condition_info from the run_info_data to determine if the sequencing kit
    // used has a rapid adapter and which one.
    const auto condition_info = models::ConditionInfo(get_chemistry_key(run_info_data));
    new_read->read_common.rapid_chemistry = condition_info.rapid_chemistry();
    new_read->read_common.chemistry = condition_info.chemistry();

    pod5_end_reason_t end_reason_value{POD5_END_REASON_UNKNOWN};
    char end_reason_string_value[200]{};
    size_t end_reason_string_value_size = sizeof(end_reason_string_value);
    {
        pod5_error_t pod5_ret =
                pod5_get_end_reason(batch, read_data.end_reason, &end_reason_value,
                                    end_reason_string_value, &end_reason_string_value_size);
        if (pod5_ret != POD5_OK) {
            issue_pod5_error("Failed to get end_reason", filename, read_id_str);
            return nullptr;
        } else if (end_reason_value == POD5_END_REASON_UNBLOCK_MUX_CHANGE ||
                   end_reason_value == POD5_END_REASON_MUX_CHANGE) {
            new_read->read_common.attributes.is_end_reason_mux_change = true;
        }
        new_read->read_common.attributes.end_reason = end_reason_string_value;
    }

    char pore_type_string_value[200]{};
    size_t pore_type_string_value_size = sizeof(pore_type_string_value);
    {
        pod5_error_t pod5_ret = pod5_get_pore_type(
                batch, read_data.pore_type, pore_type_string_value, &pore_type_string_value_size);
        if (pod5_ret != POD5_OK) {
            issue_pod5_error("Failed to get pore_type", filename, read_id_str);
            return nullptr;
        }
        new_read->read_common.attributes.pore_type = pore_type_string_value;
    }

    // Determine the time sorted predecessor of the read
    // if that information is available (primarily used for offline
    // duplex runs).
    if (reads_by_channel.find(read_data.channel) != reads_by_channel.end()) {
        auto& read_id = new_read->read_common.read_id;
        const auto& v = reads_by_channel.at(read_data.channel);
        auto read_id_iter = v.begin() + read_id_to_index.at(read_id);

        if (read_id_iter != v.begin()) {
            new_read->prev_read = std::prev(read_id_iter)->read_id;
        }
        if (std::next(read_id_iter) != v.end()) {
            new_read->next_read = std::next(read_id_iter)->read_id;
        }
    }
    return new_read;
}

SimplexReadPtr process_slow5_thread_fn(slow5_file_t *sp, slow5_rec_t * rec, const std::unordered_map<int, std::vector<DataLoader::ReadSortInfo>>& reads_by_channel,
        const std::unordered_map<std::string, size_t>& read_id_to_index, const std::unordered_set<std::string>* aux_fields){
    std::vector<int16_t> tmp(rec->raw_signal, rec->raw_signal + rec->len_raw_signal);
    char* run_id_c = slow5_hdr_get("run_id", rec->read_group, sp->header);
    std::string run_id = "";
    if(!run_id_c){
        spdlog::error("Failed to get run_id in {}", sp->meta.pathname);
    } else{
        run_id = std::string(run_id_c);
    }
    char* protocol_run_id_c = slow5_hdr_get("protocol_run_id", rec->read_group, sp->header);
    std::string protocol_run_id = "";
    if(!protocol_run_id_c){
        spdlog::error("Failed to get protocol_run_id in {}", sp->meta.pathname);
    } else{
        protocol_run_id = std::string(protocol_run_id_c);
    }
    char* flow_cell_id_c = slow5_hdr_get("flow_cell_id", rec->read_group, sp->header);
    std::string flow_cell_id = "";
    if(flow_cell_id_c){
        flow_cell_id = std::string(flow_cell_id_c);
    }
    char* flow_cell_product_code_c = slow5_hdr_get("flow_cell_product_code", rec->read_group, sp->header);
    std::string flow_cell_product_code = "";
    if(flow_cell_product_code_c){
        flow_cell_product_code = std::string(flow_cell_product_code_c);
    }
    char* position_id_c = slow5_hdr_get("sequencer_position", rec->read_group, sp->header);
    std::string position_id = "";
    if(position_id_c){
        position_id = std::string(position_id_c);
    }

    char* experiment_id_c = slow5_hdr_get("experiment_name", rec->read_group, sp->header);
    std::string experiment_id = "";
    if(experiment_id_c){
        experiment_id = std::string(experiment_id_c);
    }

    char* sequencing_kit_c = slow5_hdr_get("sequencing_kit", rec->read_group, sp->header);
    std::string sequencing_kit = "";
    if(sequencing_kit_c){
        sequencing_kit = std::string(sequencing_kit_c);
    }
    int ret = 0;
    uint64_t start_time = slow5_aux_get_uint64(rec, "start_time", &ret);
    if (ret != 0) {
        spdlog::error("Error in getting auxiliary attribute 'start_time' from the file.");
    }
    ret = 0;
    uint32_t mux = slow5_aux_get_uint8(rec, "start_mux", &ret);
    if (ret != 0) {
        spdlog::error("Error in getting auxiliary attribute 'start_mux' from the file.");
    }
    ret = 0;
    int32_t read_number = slow5_aux_get_int32(rec, "read_number", &ret);
    if (ret != 0) {
        spdlog::error("Error in getting auxiliary attribute 'read_number' from the file.");
    }
    ret = 0;
    uint64_t len;
    std::string channel_number_str = slow5_aux_get_string(rec, "channel_number", &len, &ret);
    if (ret != 0) {
        spdlog::error("Error in getting auxiliary attribute 'channel_number' from the file.");
    }
    int32_t channel_number = static_cast<int32_t>(std::stol(channel_number_str));

    char* exp_start_time_ms_c = slow5_hdr_get("acquisition_start_time", rec->read_group, sp->header);
    std::string exp_start_time_ms = "";
    if(!exp_start_time_ms_c){
        exp_start_time_ms_c = slow5_hdr_get("exp_start_time", rec->read_group, sp->header);
        if(!exp_start_time_ms_c) {
            spdlog::error("Neither acquisition_start_time nor exp_start_time found");
            return nullptr;
        }
    }
    exp_start_time_ms = std::string(exp_start_time_ms_c);
    auto run_acquisition_start_time_ms = utils::get_unix_time_ms_from_string_timestamp(exp_start_time_ms);
    // std::cerr << "run_acquisition_start_time_ms: " << run_acquisition_start_time_ms << std::endl;
    auto start_time_ms = ((start_time * 1000) /(uint64_t)rec->sampling_rate);
    auto start_time_str = utils::get_string_timestamp_from_unix_time_ms(run_acquisition_start_time_ms + start_time_ms);


    char* sample_id_c = slow5_hdr_get("sample_id", rec->read_group, sp->header);
    std::string sample_id = "";
    if(!sample_id_c){
        spdlog::error("Failed to get sample_id in {}", sp->meta.pathname);
        return nullptr;
    } else{
        sample_id = std::string(sample_id_c);
    }

    char* protocol_start_time_c = slow5_hdr_get("protocol_start_time", rec->read_group, sp->header);
    std::string protocol_start_time = "";
    if(!protocol_start_time_c){
        spdlog::error("Failed to get protocol_start_time in {}", sp->meta.pathname);
        return nullptr;
    } else{
        protocol_start_time = std::string(protocol_start_time_c);
    }
    auto protocol_start_time_ms = utils::get_unix_time_ms_from_string_timestamp(protocol_start_time);

    uint64_t num_minknow_events = 0;
    if (aux_fields && aux_fields->count("num_minknow_events")) {
    ret = 0;
        num_minknow_events = slow5_aux_get_uint64(rec, "num_minknow_events", &ret);
        if (ret != 0) {
            spdlog::error("Error in getting auxiliary attribute 'num_minknow_events' from the file.");
        }
    }

    float open_pore_level = nanf("");
    if (aux_fields && aux_fields->count("open_pore_level")) {
        ret = 0;
        open_pore_level = slow5_aux_get_float(rec, "open_pore_level", &ret);
        if (ret != 0 && ret != SLOW5_ERR_NOFLD) {
            spdlog::debug("Could not get 'open_pore_level' for read {}, err={}", rec->read_id, ret);
        }
    }

    auto new_read = std::make_unique<SimplexRead>();

    auto options = at::TensorOptions().dtype(at::kShort);
    // Previously: .clone().to(m_device_)
    new_read->read_common.raw_data = at::from_blob(tmp.data(), tmp.size(), options).clone();
    new_read->read_common.attributes.sample_rate = (uint64_t)rec->sampling_rate;
    new_read->run_acquisition_start_time_ms = run_acquisition_start_time_ms;
    new_read->read_common.start_time_ms = start_time_ms;
    new_read->scaling = rec->range / rec->digitisation;
    new_read->offset = rec->offset;
    new_read->open_pore_level = open_pore_level;
    new_read->read_common.read_id = std::string(rec->read_id);
    new_read->read_common.num_trimmed_samples = 0;
    new_read->read_common.attributes.read_number = read_number;
    new_read->read_common.attributes.filename = std::filesystem::path(sp->meta.pathname).filename().string();
    new_read->read_common.attributes.mux = mux;
    new_read->read_common.attributes.num_samples = rec->len_raw_signal;
    new_read->read_common.attributes.channel_number = channel_number;
    new_read->read_common.attributes.start_time = start_time_str;
    new_read->read_common.run_id = protocol_run_id; // not anymore run_id
    new_read->read_common.acquisition_id = run_id;
    new_read->start_sample = start_time;
    new_read->end_sample = start_time + rec->len_raw_signal;
    new_read->read_common.flowcell_id = flow_cell_id;
    new_read->read_common.sequencing_kit = sequencing_kit;
    new_read->read_common.flow_cell_product_code = flow_cell_product_code;
    new_read->read_common.position_id = position_id;
    new_read->read_common.sample_id = sample_id;
    new_read->read_common.protocol_start_time_ms = protocol_start_time_ms;
    new_read->read_common.is_duplex = false;
    new_read->read_common.experiment_id = experiment_id;
    new_read->read_common.num_minknow_events = num_minknow_events;
    // Get the condition_info from the run_info_data to determine if the sequencing kit
    // used has a rapid adapter and which one.

    RunInfoDictData_t run_info_data;
    run_info_data.flow_cell_product_code = flow_cell_product_code.c_str();
    run_info_data.sequencing_kit = sequencing_kit.c_str();
    run_info_data.sample_rate = rec->sampling_rate;

    const auto condition_info = models::ConditionInfo(get_chemistry_key(&run_info_data));
    new_read->read_common.rapid_chemistry = condition_info.rapid_chemistry();
    new_read->read_common.chemistry = condition_info.chemistry();

    /*ret = 0;
    uint8_t end_reason = slow5_aux_get_enum(rec,"end_reason",&ret);
    if(ret!=0){
        spdlog::error("Error in getting auxiliary attribute 'end_reason' from the file.");
    }
    uint8_t num_label = 0;
    char **labels = slow5_get_aux_enum_labels(sp->header, "end_reason", &num_label);
    if(labels==NULL){
        spdlog::error("Error in getting auxiliary attribute 'end_reason' labels from the file.");
    }
    if(end_reason != SLOW5_ENUM_NULL){
        // pod5_end_reason_t end_reason_value = static_cast<pod5_end_reason_t>(end_reason);
        // if (end_reason_value == POD5_END_REASON_UNBLOCK_MUX_CHANGE ||
        //        end_reason_value == POD5_END_REASON_MUX_CHANGE) {
        //     new_read->read_common.attributes.is_end_reason_mux_change = true;
        // }
        if(strcmp(labels[end_reason],"unblock_mux_change")==0 || strcmp(labels[end_reason],"mux_change")==0){
            new_read->read_common.attributes.is_end_reason_mux_change = true;
        }
    } else{
        spdlog::error("Failed to get read end_reason for read {}", new_read->read_common.read_id);
    }*/
    ret = 0;
    uint8_t end_reason = slow5_aux_get_enum(rec, "end_reason", &ret);
    if (ret != 0) {
        spdlog::error("Error getting aux 'end_reason' for read {}", new_read->read_common.read_id);
    } else if (end_reason != SLOW5_ENUM_NULL) {
        uint8_t num_label = 0;
        char **labels = slow5_get_aux_enum_labels(sp->header, "end_reason", &num_label);
        if (!labels) {
            spdlog::error("Failed to get enum labels for 'end_reason' for read {}", new_read->read_common.read_id);
        } else if (end_reason < num_label && labels[end_reason] != nullptr) {
            new_read->read_common.attributes.end_reason = labels[end_reason];
            if (std::strcmp(labels[end_reason], "unblock_mux_change") == 0 ||
                std::strcmp(labels[end_reason], "mux_change") == 0) {
                new_read->read_common.attributes.is_end_reason_mux_change = true;
            }
        } else {
            spdlog::error("end_reason index {} out of range (num_label={}) for read {}", static_cast<int>(end_reason), static_cast<int>(num_label), new_read->read_common.read_id);
        }
    } else {
        spdlog::debug("end_reason is SLOW5_ENUM_NULL for read {}", new_read->read_common.read_id);
    }

    // new_read->read_common.attributes.pore_type = not set


    // Determine the time sorted predecessor of the read
    // if that information is available (primarily used for offline
    // duplex runs).
    if (reads_by_channel.find(channel_number) != reads_by_channel.end()) {
        auto& read_id = new_read->read_common.read_id;
        const auto& v = reads_by_channel.at(channel_number);
        auto read_id_iter = v.begin() + read_id_to_index.at(read_id);

        if (read_id_iter != v.begin()) {
            new_read->prev_read = std::prev(read_id_iter)->read_id;
        }
        if (std::next(read_id_iter) != v.end()) {
            new_read->next_read = std::next(read_id_iter)->read_id;
        }
    }

    return new_read;

}

void process_read_data_slow5(core_t *core, db_t *db, int32_t i) {
    //
    struct slow5_rec *rec = NULL;
    if (slow5_rec_depress_parse(&db->mem_records[i], &db->mem_bytes[i], NULL, &rec, core->fp) != 0) {
        spdlog::error("Error in parsing read {}", i);
    } else {
        free(db->mem_records[i]);
    }
    auto new_read = process_slow5_thread_fn(core->fp, rec, core->reads_by_channel, core->read_id_to_index, core->slow5_aux_fields);
    //
    // db->read_data_ptrs[i] = new_read;
    db->read_data_ptrs[i] = std::move(new_read);
    slow5_rec_free(rec);
}

}  // namespace

void DataLoader::load_reads_by_channel(const std::vector<std::filesystem::directory_entry>& files) {
    std::string slow5_file_path = "";
    int slow5_file_count = 0;
    slow5_file_t *sp = NULL;
    for (const auto& entry : files) {
        if(utils::has_blow5_extension(entry)){
            slow5_file_path = entry.path();
            slow5_file_count++;
        }
    }
    if(slow5_file_count > 1){
        throw std::runtime_error("Please provide a single BLOW5 file path for duplex calling. Multiple files are not supported");
    } 
    if(slow5_file_count == 1){
        sp = slow5_open(slow5_file_path.c_str(), "r");
        if (sp == NULL) {
            throw std::runtime_error("Error in opening SLOW5/BLOW5 file");
        }
        int ret = slow5_idx_load(sp);
        if(ret<0){
            throw std::runtime_error("Error in loading index for SLOW5/BLOW5 file");
        }
        m_slow5_aux_fields.clear();
        uint64_t num_aux = 0;
        char **aux_names = slow5_get_aux_names(sp->header, &num_aux);
        if (aux_names != NULL) {
            for (uint64_t i = 0; i < num_aux; i++) {
                m_slow5_aux_fields.insert(aux_names[i]);
            }
        }
    }
    // If traversal in channel order is required, the following algorithm
    // is used -
    // 1. iterate through all the read metadata to collect channel information
    // across all pod5 files
    // 2. store the read list sorted by channel number
    spdlog::info("> Reading read channel info");
    if(slow5_file_count == 1){
        load_read_channels_slow5(sp, slow5_file_path);
    } else {
        load_read_channels(files);
    }
    spdlog::info("> Processed read channel info");
    // 3. for each channel, iterate through all files and in each iteration
    // only load the reads that correspond to that channel.
    for (int channel = 0; channel <= m_max_channel; channel++) {
        if (m_reads_by_channel.find(channel) != m_reads_by_channel.end()) {
            // Sort the read ids within a channel by its mux
            // and start time.
            spdlog::trace("Sort channel {}", channel);
            auto& reads = m_reads_by_channel.at(channel);
            std::sort(reads.begin(), reads.end(), [](ReadSortInfo& a, ReadSortInfo& b) {
                if (a.mux != b.mux) {
                    return a.mux < b.mux;
                } else {
                    return a.read_number < b.read_number;
                }
            });
            // Once sorted, create a hash table from read id
            // to index in the sorted list to quickly fetch the
            // read location and its neighbors.
            for (size_t i = 0; i < reads.size(); i++) {
                m_read_id_to_index[reads[i].read_id] = i;
            }
            spdlog::trace("Sorted channel {}", channel);
        }
        for (const auto& entry : files) {
            if (m_loaded_read_count == m_max_reads ||
                m_stop_loading.load(std::memory_order_relaxed)) {
                break;
            }

            if (utils::has_pod5_extension(entry)) {
                const auto path = std::filesystem::path(entry);
                auto& channel_to_read_ids = m_file_channel_read_order_map.at(path.string());
                auto& read_ids = channel_to_read_ids[channel];
                if (!read_ids.empty()) {
                    load_pod5_reads_from_file_by_read_ids(path.string(), read_ids);
                }
            } else if (utils::has_blow5_extension(entry) && sp != NULL) {
                auto& channel_to_read_ids = m_file_channel_read_order_map.at(slow5_file_path);
                auto& read_ids = channel_to_read_ids[channel];
                if (!read_ids.empty()) {
                    load_slow5_reads_from_file_by_read_ids(sp, read_ids);
                }
            } else {
                throw std::logic_error("Expected pod5/blow5 file");
            }

        }
        // Erase sorted list as it's not needed anymore.
        m_reads_by_channel.erase(channel);
    }
    if(slow5_file_count == 1){
        slow5_idx_unload(sp);
        slow5_close(sp);
    }
}

void DataLoader::load_reads_unrestricted(
        const std::vector<std::filesystem::directory_entry>& files) {
    for (const auto& entry : files) {
        if (m_loaded_read_count == m_max_reads || m_stop_loading.load(std::memory_order_relaxed)) {
            break;
        }
        spdlog::debug("Load reads from file {}", entry.path().string());
        std::string ext = std::filesystem::path(entry).extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (ext == ".pod5") {
            load_pod5_reads_from_file(entry.path().string());
        } else if (ext == ".slow5" || ext == ".blow5") {
            load_slow5_reads_from_file(entry.path().string());
        }
    }
}

void DataLoader::load_reads(const InputFiles& input_files, ReadOrder traversal_order) {
    if (pod5_init() != POD5_OK) {
        throw std::runtime_error(
                fmt::format("Failed to initialise POD5: {}", pod5_get_error_string()));
    }
    auto pod5_cleanup = utils::PostCondition([] { pod5_terminate(); });

    utils::start_busy_work();

    switch (traversal_order) {
    case ReadOrder::BY_CHANNEL:
        load_reads_by_channel(input_files.get());
        break;
    case ReadOrder::UNRESTRICTED:
        load_reads_unrestricted(input_files.get());
        break;
    default:
        throw std::runtime_error("Unsupported traversal order detected: " +
                                 dorado::to_string(traversal_order));
    }
}

void DataLoader::load_read_channels_slow5(::slow5_file* sp, const std::string& file_path_str) {
    // Use a std::map to store by sorted channel order.
	m_file_channel_read_order_map.emplace(file_path_str, channel_to_read_id_t());
	auto& channel_to_read_id = m_file_channel_read_order_map[file_path_str];

	slow5_rec_t **rec = NULL;
	int ret_batch=0;
	int ret=0;

	slow5_mt_t *mt = slow5_init_mt(slow5_threads,sp);
	slow5_batch_t *read_batch = slow5_init_batch(slow5_batchsize);
    while((ret_batch = slow5_get_next_batch(mt,read_batch,slow5_batchsize)) > 0){
	    rec = read_batch->slow5_rec;
	    for(int i=0;i<ret_batch;i++){
	        uint64_t len; //length of the array
	        char* channel_number = slow5_aux_get_string(rec[i], "channel_number", &len, &ret);
	        if(ret!=0){
                spdlog::error("Error in getting auxiliary attribute 'channel_number' from the file. Error code {}", ret);
	        }
	        if (channel_number == NULL){ //check if the field value exists and print the value
	            spdlog::error("channel_number is missing for the record {}", rec[i]->read_id);
	        } else{
	            int channel = atoi(channel_number);
	            // Update maximum number of channels encountered.
	            m_max_channel = std::max(m_max_channel, channel);
	            
	            // Store the read_id in the channel's list.
	            uint8_t  arr[16] = {0};
	            ret = sscanf(rec[i]->read_id, "%2hhx%2hhx%2hhx%2hhx-%2hhx%hhx-%2hhx%2hhx-%2hhx%2hhx-%2hhx%2hhx%2hhx%2hhx%2hhx%2hhx",
	                            &arr[0], &arr[1], &arr[2], &arr[3], &arr[4], &arr[5], &arr[6], &arr[7],
	                            &arr[8], &arr[9], &arr[10], &arr[11], &arr[12], &arr[13], &arr[14], &arr[15]);
	            if(ret !=16){
                    spdlog::error("Error in parsing uuid for the record {}. Return value: {}", rec[i]->read_id, ret);
	            }

	            ReadID read_id;
	            std::memcpy(read_id.data(), arr, POD5_READ_ID_SIZE);
	            channel_to_read_id[channel].push_back(std::move(read_id));

	            std::string rid(rec[i]->read_id);
	            ret = 0;
	            uint32_t mux = slow5_aux_get_uint8(rec[i], "start_mux", &ret);
	            if (ret != 0) {
	                throw std::runtime_error("Error in getting auxiliary attribute 'start_mux' from the file.");
	            }
	            ret = 0;
	            int32_t read_number = slow5_aux_get_int32(rec[i], "read_number", &ret);
	            if (ret != 0) {
	                throw std::runtime_error("Error in getting auxiliary attribute 'read_number' from the file.");
	            }
	            m_reads_by_channel[channel].push_back({rid, (int32_t)mux, (uint32_t)read_number});

	        }
	    }
	    if(ret_batch<slow5_batchsize){ //this indicates nothing left to read //need to handle errors
	        break;
	    }
	}
    slow5_free_batch(read_batch);
	slow5_free_mt(mt);
}

void DataLoader::load_read_channels(const std::vector<std::filesystem::directory_entry>& files) {
    for (const auto& entry : files) {
        auto file_path = std::filesystem::path(entry);
        std::string ext = file_path.extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (ext != ".pod5") {
            continue;
        }

        // Use a std::map to store by sorted channel order.
        const auto file_string = file_path.string();
        auto& channel_to_read_id = m_file_channel_read_order_map[file_string];

        // Open the file ready for walking:
        Pod5FileReader_t* file = pod5_open_file(file_string.c_str());
        if (!file) {
            issue_pod5_error("Failed to open file", file_string);
            continue;
        }
        auto cleanup_file = utils::PostCondition([&file, &file_string] {
            if (pod5_close_and_free_reader(file) != POD5_OK) {
                issue_pod5_error("Failed to close and free reader", file_string);
            }
        });

        std::size_t batch_count = 0;
        if (pod5_get_read_batch_count(&batch_count, file) != POD5_OK) {
            issue_pod5_error("Failed to query batch count", file_string);
            continue;
        }

        for (std::size_t batch_index = 0; batch_index < batch_count; ++batch_index) {
            Pod5ReadRecordBatch_t* batch = nullptr;
            if (pod5_get_read_batch(&batch, file, batch_index) != POD5_OK) {
                issue_pod5_error("Failed to get batch", file_string, batch_index, 0);
                continue;
            }
            auto cleanup_batch = utils::PostCondition([&batch, &file_string, batch_index] {
                if (pod5_free_read_batch(batch) != POD5_OK) {
                    issue_pod5_error("Failed to get batch", file_string, batch_index, 0);
                }
            });

            std::size_t batch_row_count = 0;
            if (pod5_get_read_batch_row_count(&batch_row_count, batch) != POD5_OK) {
                issue_pod5_error("Failed to get row count", file_string, batch_index, 0);
                continue;
            }

            for (std::size_t row = 0; row < batch_row_count; ++row) {
                uint16_t read_table_version = 0;
                ReadBatchRowInfo_t read_data{};
                if (pod5_get_read_batch_row_info_data(batch, row, READ_BATCH_ROW_INFO_VERSION,
                                                      &read_data, &read_table_version) != POD5_OK) {
                    issue_pod5_error("Failed to get read", file_string, batch_index, row);
                    continue;
                }

                int channel = read_data.channel;

                // Update maximum number of channels encountered.
                m_max_channel = std::max(m_max_channel, channel);

                // Store the read_id in the channel's list.
                ReadID read_id;
                std::memcpy(read_id.data(), read_data.read_id, POD5_READ_ID_SIZE);
                channel_to_read_id[channel].push_back(read_id);

                char read_id_tmp[POD5_READ_ID_LEN];
                if (pod5_format_read_id(read_data.read_id, read_id_tmp) != POD5_OK) {
                    issue_pod5_error("Failed to format read id", file_string, batch_index, row);
                }
                std::string rid(read_id_tmp);
                m_reads_by_channel[channel].push_back({rid, read_data.well, read_data.read_number});
            }
        }
    }
}

void DataLoader::load_pod5_reads_from_file_by_read_ids(const std::string& path,
                                                       const std::vector<ReadID>& read_ids) {
    // Open the file ready for walking:
    // TODO: The earlier implementation was caching the pod5 readers into a
    // map and re-using it during each iteration. However, we found a leak
    // in the pod5 traversal API which persists unless the reader is opened
    // and closed everytime. So the caching logic was reverted until the
    // leak is fixed in pod5 API.
    Pod5FileReader_t* file = pod5_open_file(path.c_str());
    if (!file) {
        issue_pod5_error("Failed to open file", path);
        return;
    }
    auto cleanup_file = utils::PostCondition([&file, &path] {
        if (pod5_close_and_free_reader(file) != POD5_OK) {
            issue_pod5_error("Failed to close and free reader", path);
        }
    });

    std::vector<uint8_t> read_id_array(POD5_READ_ID_SIZE * read_ids.size());
    for (size_t i = 0; i < read_ids.size(); i++) {
        std::memcpy(read_id_array.data() + POD5_READ_ID_SIZE * i, read_ids[i].data(),
                    POD5_READ_ID_SIZE);
    }

    std::size_t batch_count = 0;
    if (pod5_get_read_batch_count(&batch_count, file) != POD5_OK) {
        issue_pod5_error("Failed to query batch count", path);
        return;
    }

    std::vector<std::uint32_t> traversal_batch_counts(batch_count);
    std::vector<std::uint32_t> traversal_batch_rows(read_ids.size());
    size_t find_success_count;
    pod5_error_t err = pod5_plan_traversal(file, read_id_array.data(), read_ids.size(),
                                           traversal_batch_counts.data(),
                                           traversal_batch_rows.data(), &find_success_count);
    if (err != POD5_OK) {
        issue_pod5_error("Failed to plan traversal", path);
        return;
    }

    if (find_success_count != read_ids.size()) {
        spdlog::error(
                "POD5 failed to plan traversal of '{}' - Reads found by plan {}, reads in input {}",
                path, find_success_count, read_ids.size());
        throw std::runtime_error("Plan traversal didn't yield correct number of reads");
    }

    uint32_t row_offset = 0;
    for (std::size_t batch_index = 0; batch_index < batch_count; ++batch_index) {
        if (m_loaded_read_count == m_max_reads || m_stop_loading.load(std::memory_order_relaxed)) {
            break;
        }
        Pod5ReadRecordBatch_t* batch = nullptr;
        if (pod5_get_read_batch(&batch, file, batch_index) != POD5_OK) {
            issue_pod5_error("Failed to get batch", path, batch_index, 0);
            continue;
        }
        auto cleanup_batch = utils::PostCondition([&batch, &path, batch_index] {
            if (pod5_free_read_batch(batch) != POD5_OK) {
                issue_pod5_error("Failed to release batch", path, batch_index, 0);
            }
        });

        const std::size_t num_rows = traversal_batch_counts[batch_index];
        std::vector<std::future<SimplexReadPtr>> futures;
        futures.reserve(num_rows);
        for (std::size_t row_idx = 0; row_idx < num_rows; row_idx++) {
            uint32_t row = traversal_batch_rows[row_idx + row_offset];
            futures.push_back(m_thread_pool.push([row, batch_index, batch, file, &path, this] {
                return process_pod5_thread_fn(row, batch_index, batch, file, path,
                                              m_reads_by_channel, m_read_id_to_index,
                                              m_allowed_read_ids, m_ignored_read_ids);
            }));
        }
        row_offset += num_rows;

        wait_and_process_futures(std::move(futures));
    }
}

void DataLoader::load_pod5_reads_from_file(const std::string& path) {
    // Open the file ready for walking:
    Pod5FileReader_t* file = pod5_open_file(path.c_str());
    if (!file) {
        issue_pod5_error("Failed to open file", path);
        return;
    }
    auto file_cleanup = utils::PostCondition([&file, &path] {
        if (pod5_close_and_free_reader(file) != POD5_OK) {
            issue_pod5_error("Failed to close and free POD5 reader", path);
        }
    });

    std::size_t batch_count = 0;
    if (pod5_get_read_batch_count(&batch_count, file) != POD5_OK) {
        issue_pod5_error("Failed to query batch count", path);
        return;
    }

    for (std::size_t batch_index = 0; batch_index < batch_count; ++batch_index) {
        if (m_loaded_read_count == m_max_reads || m_stop_loading.load(std::memory_order_relaxed)) {
            break;
        }
        Pod5ReadRecordBatch_t* batch = nullptr;
        if (pod5_get_read_batch(&batch, file, batch_index) != POD5_OK) {
            issue_pod5_error("Failed to get batch", path, batch_index, 0);
            continue;
        }
        auto cleanup_batch = utils::PostCondition([&batch, &path, batch_index] {
            if (pod5_free_read_batch(batch) != POD5_OK) {
                issue_pod5_error("Failed to release batch", path, batch_index, 0);
            }
        });

        std::size_t batch_row_count = 0;
        if (pod5_get_read_batch_row_count(&batch_row_count, batch) != POD5_OK) {
            issue_pod5_error("Failed to get batch row count", path, batch_index, 0);
            continue;
        }
        batch_row_count = std::min(batch_row_count, m_max_reads - m_loaded_read_count);

        std::vector<std::future<SimplexReadPtr>> futures;

        futures.reserve(batch_row_count);
        for (std::size_t row = 0; row < batch_row_count; ++row) {
            futures.push_back(m_thread_pool.push([row, batch_index, batch, file, &path, this] {
                return process_pod5_thread_fn(row, batch_index, batch, file, path,
                                              m_reads_by_channel, m_read_id_to_index,
                                              m_allowed_read_ids, m_ignored_read_ids);
            }));
        }

        wait_and_process_futures(std::move(futures));
    }
}

void DataLoader::wait_and_process_futures(std::vector<std::future<SimplexReadPtr>> futures) {
    for (auto& v : futures) {
        auto read = v.get();
        if (!read) {
            // This was either a POD5 error, in which case the worker logged
            // an error, or a filtered read.
            continue;
        }
        if (!m_pipeline.is_running()) {
            // If the pipeline has finished early (--run-for) then stop processing
            // reads, but don't bail since we need to wait for all futures to finish.
            m_stop_loading.store(true, std::memory_order_relaxed);
            continue;
        }
        initialise_read(read->read_common);
        check_read(read);
        m_pipeline.push_message(std::move(read));
        m_loaded_read_count++;
    }
}

void DataLoader::initialise_read(ReadCommon& read_common) const {
    for (const auto& initialiser : m_read_initialisers) {
        initialiser(read_common);
    }
}

void DataLoader::check_read(const SimplexReadPtr& read) {
    if (read->read_common.chemistry == models::Chemistry::UNKNOWN &&
        m_log_unknown_chemistry.exchange(false)) {
        spdlog::warn(
                "Could not determine sequencing Chemistry from read data - "
                "some features might be disabled");
    }
}

static void on_worker_start() {
    // Setting the thread name is expensive, so do it once per thread at startup.
    utils::set_thread_name("dataloader_pod5");
}

void DataLoader::load_slow5_reads_from_file(const std::string &path) {    
    slow5_file_t *sp = slow5_open(path.c_str(), "r");
    if (sp == NULL) {
        spdlog::error("Error in opening file {}", path);
        return;
    }

    m_slow5_aux_fields.clear();
    uint64_t num_aux = 0;
    char **aux_names = slow5_get_aux_names(sp->header, &num_aux);
    if (aux_names != NULL) {
        for (uint64_t i = 0; i < num_aux; i++) {
            m_slow5_aux_fields.insert(aux_names[i]);
        }
    }

    int64_t batch_size = slow5_batchsize;
    int32_t num_threads = slow5_threads;

    while (1) {
        int flag_EOF = 0;
        db_t db = {0};
        db.mem_records = (char **) malloc(batch_size * sizeof *db.mem_records);
        db.mem_bytes = (size_t *) malloc(batch_size * sizeof *db.mem_bytes);

        int64_t record_count = 0;
        size_t bytes;
        char *mem;
        while (record_count < batch_size) {
            if (!(mem = (char *) slow5_get_next_mem(&bytes, sp))) {
                if (slow5_errno != SLOW5_ERR_EOF) {
                    throw std::runtime_error("Error in slow5_get_next_mem.");
                } else { //EOF file reached
                    flag_EOF = 1;
                    break;
                }
            } else {
                db.mem_records[record_count] = mem;
                db.mem_bytes[record_count] = bytes;
                record_count++;
            }
        }

        // Setup multithreading structures
        core_t core = {0};
        core.num_thread = (num_threads > record_count) ? record_count : num_threads;
        if (record_count == 0) {
            core.num_thread = 1;
        }
        core.fp = sp;
        // core.m_device_ = m_device;
        core.reads_by_channel = std::cref(m_reads_by_channel);
        core.read_id_to_index = std::cref(m_read_id_to_index);
        core.slow5_aux_fields  = &m_slow5_aux_fields;

        db.n_batch = record_count;
        db.read_data_ptrs = std::vector<dorado::SimplexReadPtr>(record_count);
        work_db(&core, &db, process_read_data_slow5);

        for (int64_t i = 0; i < record_count; i++) {
            if (!m_allowed_read_ids ||
                (m_allowed_read_ids->find(db.read_data_ptrs[i]->read_common.read_id) != m_allowed_read_ids->end())) {
                initialise_read(db.read_data_ptrs[i]->read_common);
                check_read(db.read_data_ptrs[i]);
                m_pipeline.push_message(std::move(db.read_data_ptrs[i]));
                m_loaded_read_count++;
            }
        }
        // Free everything
        free(db.mem_bytes);
        free(db.mem_records);

        if (flag_EOF == 1) {
            break;
        }
    }
    slow5_close(sp);
}

void DataLoader::load_slow5_reads_from_file_by_read_ids(::slow5_file* sp, const std::vector<ReadID>& read_ids) {
    int ret = 0;
    slow5_rec_t **rec = NULL;

    size_t num_rid = read_ids.size();
    size_t read_count = 0;
    char **rid = (char**)malloc(sizeof(char*)*slow5_batchsize);

    while(1){
        int local_batch_size = ((size_t)slow5_batchsize > (num_rid - read_count)) ? num_rid - read_count : (size_t)slow5_batchsize;
        for(size_t i=0; i<(size_t)local_batch_size; i++){
            char read_id_tmp[POD5_READ_ID_LEN];
            if (pod5_format_read_id(read_ids[read_count].data(), read_id_tmp) != POD5_OK) {
                spdlog::error("Failed to format read id");
            }
            std::string read_s(read_id_tmp);
            
            rid[i] = strdup(read_s.c_str());
            read_count++;
        }
        ret = slow5_get_batch_lazy(&rec, sp, rid, local_batch_size, slow5_threads);
        assert(ret==local_batch_size);
        for(int i=0;i<ret;i++){
            if (!m_allowed_read_ids ||
                (m_allowed_read_ids->find(std::string(rec[i]->read_id)) != m_allowed_read_ids->end())) {
                auto new_read = process_slow5_thread_fn(sp, rec[i], m_reads_by_channel, m_read_id_to_index, &m_slow5_aux_fields);
                spdlog::debug("read_id queued: {}",rec[i]->read_id);
                // m_pipeline.push_message(new_read);
                initialise_read(new_read->read_common);
                check_read(new_read);
                m_pipeline.push_message(std::move(new_read));
                m_loaded_read_count++;
            }
        }
        slow5_free_batch_lazy(&rec,ret);
        for(int i=0; i<local_batch_size; i++){
            free(rid[i]);
        }
        if(local_batch_size < slow5_batchsize){
            break;
        }
    }
    free(rid);
    
}

DataLoader::DataLoader(Pipeline& pipeline,
                       const std::string& device,
                       size_t num_worker_threads,
                       size_t max_reads,
                       std::optional<std::unordered_set<std::string>> read_list,
                       std::unordered_set<std::string> read_ignore_list,
                    int32_t slow5_threads_,
                    int64_t slow5_batchsize_)
        : m_pipeline(pipeline),
          m_device(device),
          m_thread_pool(num_worker_threads, on_worker_start),
          m_allowed_read_ids(std::move(read_list)),
          m_ignored_read_ids(std::move(read_ignore_list)) {
    // ReadID should be a drop-in replacement for read_id_t
    static_assert(sizeof(ReadID) == sizeof(read_id_t));

    slow5_threads = slow5_threads_;
    slow5_batchsize = slow5_batchsize_;
    m_max_reads = max_reads == 0 ? std::numeric_limits<decltype(m_max_reads)>::max() : max_reads;
    assert(m_thread_pool.n_threads() > 0);
}

InputFiles InputFiles::search_pod5s(const std::filesystem::path& path, bool recursive) {
    auto entries = collect_pod5_dataset(utils::fetch_directory_entries(path, recursive));

    // Intentionally returning a valid object even if there are 0 entries since duplex uses that
    // to differentiate between different modes of operation.
    InputFiles files;
    files.m_entries = std::move(entries);
    return files;
}

const std::vector<std::filesystem::directory_entry>& InputFiles::get() const { return m_entries; }

}  // namespace dorado::data_loader

#include "HtsWriter.h"

#include "htslib/bgzf.h"
#include "htslib/kroundup.h"
#include "htslib/sam.h"
#include "read_pipeline/ReadPipeline.h"
#include "utils/sequence_utils.h"

#include <indicators/progress_bar.hpp>
#include <spdlog/spdlog.h>

#include <stdexcept>
#include <string>
#include <unordered_set>

namespace dorado {

HtsWriter::HtsWriter(const std::string& filename, OutputMode mode, size_t threads, size_t num_reads)
        : MessageSink(10000), m_num_reads_expected(num_reads) {
    switch (mode) {
    case FASTQ:
        m_file = hts_open(filename.c_str(), "wf");
        break;
    case BAM:
        m_file = hts_open(filename.c_str(), "wb");
        break;
    case SAM:
        m_file = hts_open(filename.c_str(), "w");
        break;
    case UBAM:
        m_file = hts_open(filename.c_str(), "wb0");
        break;
    default:
        throw std::runtime_error("Unknown output mode selected: " + std::to_string(mode));
    }
    if (!m_file) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    if (m_file->format.compression == bgzf) {
        auto res = bgzf_mt(m_file->fp.bgzf, threads, 128);
        if (res < 0) {
            throw std::runtime_error("Could not enable multi threading for BAM generation.");
        }
    }

    if (m_num_reads_expected == 0) {
        m_progress_bar_interval = 100;
    } else {
        m_progress_bar_interval = m_num_reads_expected < 100 ? 1 : 100;
    }

    m_worker = std::make_unique<std::thread>(std::thread(&HtsWriter::worker_thread, this));
}

HtsWriter::~HtsWriter() {
    // Adding for thread safety in case worker thread throws exception.
    terminate();
    if (m_worker->joinable()) {
        join();
    }
    sam_hdr_destroy(header);
    hts_close(m_file);
}

HtsWriter::OutputMode HtsWriter::get_output_mode(std::string mode) {
    if (mode == "sam") {
        return SAM;
    } else if (mode == "bam") {
        return BAM;
    } else if (mode == "fastq") {
        return FASTQ;
    }
    throw std::runtime_error("Unknown output mode: " + mode);
}

void HtsWriter::join() { m_worker->join(); }

void HtsWriter::worker_thread() {
    std::unordered_set<std::string> processed_read_ids;
    size_t write_count = 0;

    // Initialize progress logging.
    if (m_num_reads_expected == 0) {
        std::cerr << "\r> Output records written: " << write_count;
    }

    Message message;
    while (m_work_queue.try_pop(message)) {
        auto aln = std::get<BamPtr>(std::move(message));
        write(aln.get());
        std::string read_id = bam_get_qname(aln.get());
        aln.reset();  // Free the bam alignment that's already written

        // For the purpose of estimating write count, we ignore duplex reads
        // these can be identified by a semicolon in their ID.
        // TODO: This is a hack, we should have a better way of identifying duplex reads.
        bool ignore_read_id = read_id.find(';') != std::string::npos;

        if (!ignore_read_id) {
            processed_read_ids.emplace(std::move(read_id));
        }

        if (m_num_reads_expected != 0) {
            write_count = processed_read_ids.size();
        } else {
            if (!ignore_read_id) {
                write_count++;
            }
        }

        if ((write_count % m_progress_bar_interval) == 0) {
            if ((write_count == 0) && !m_prog_bar_initialized) {
                m_progress_bar.set_progress(0.0f);
                m_prog_bar_initialized = true;
            }
            if (m_num_reads_expected != 0) {
                float progress = 100.f * static_cast<float>(write_count) / m_num_reads_expected;
                m_progress_bar.set_progress(progress);
#ifndef WIN32
                std::cerr << "\033[K";
#endif  // WIN32
            } else {
                std::cerr << "\r> Output records written: " << write_count;
            }
        }
    }
    // Clear progress information.
    if (m_num_reads_expected != 0 || write_count >= m_progress_bar_interval) {
        std::cerr << "\r";
    }
    spdlog::debug("Written {} records.", write_count);
}

int HtsWriter::write(bam1_t* record) {
    // track stats
    total++;
    if (record->core.flag & BAM_FUNMAP) {
        unmapped++;
    }
    if (record->core.flag & BAM_FSECONDARY) {
        secondary++;
    }
    if (record->core.flag & BAM_FSUPPLEMENTARY) {
        supplementary++;
    }
    primary = total - secondary - supplementary - unmapped;

    auto res = sam_write1(m_file, header, record);
    if (res < 0) {
        throw std::runtime_error("Failed to write SAM record, error code " + std::to_string(res));
    }
    return res;
}

int HtsWriter::write_header(const sam_hdr_t* hdr) {
    if (hdr) {
        header = sam_hdr_dup(hdr);
        return sam_hdr_write(m_file, header);
    }
    return 0;
}

stats::NamedStats HtsWriter::sample_stats() const { return stats::from_obj(m_work_queue); }

}  // namespace dorado
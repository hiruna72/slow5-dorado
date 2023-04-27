#include "TestUtils.h"
#include "htslib/sam.h"
#include "utils/bam_utils.h"

#include <catch2/catch.hpp>

#include <filesystem>

#define TEST_GROUP "[bam_utils][bam_reader]"

namespace fs = std::filesystem;

class BamWriterTestsFixture {
public:
    BamWriterTestsFixture() {
        fs::path aligner_test_dir = fs::path(get_data_dir("bam_reader"));
        m_in_sam = aligner_test_dir / "small.sam";
        m_out_bam = fs::temp_directory_path() / "out.bam";
    }

    ~BamWriterTestsFixture() { fs::remove(m_out_bam); }

protected:
    void generate_bam(int num_threads) {
        dorado::utils::BamReader reader(m_in_sam.string());
        dorado::utils::BamWriter writer(m_out_bam.string(), num_threads);

        dorado::utils::sq_t sequences;
        CHECK(sequences.size() == 0);  // No sequence information for this test.
        writer.write_header(reader.header, sequences);
        reader.read(writer, 1000);

        writer.join();
    }

private:
    fs::path m_in_sam;
    fs::path m_out_bam;
};

TEST_CASE_METHOD(BamWriterTestsFixture, "BamWriterTest: Write BAM", TEST_GROUP) {
    int num_threads = GENERATE(1, 10);
    REQUIRE_NOTHROW(generate_bam(num_threads));
}
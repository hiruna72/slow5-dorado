#pragma once
#include "read_pipeline/ReadPipeline.h"
#include "utils/stats.h"
#include "utils/types.h"

#include <atomic>
#include <string>
#include <string_view>
#include <vector>

namespace dorado {

namespace demux {

struct KitInfo {
    bool double_ends;
    bool ends_different;
    std::string top_front_flank;
    std::string top_rear_flank;
    std::string bottom_front_flank;
    std::string bottom_rear_flank;
    std::vector<std::string> barcodes;
};

struct AdapterSequence {
    std::vector<std::string> adapter;
    std::vector<std::string> adapter_rev;
    std::string top_primer;
    std::string top_primer_rev;
    std::string bottom_primer;
    std::string bottom_primer_rev;
    int top_primer_front_flank_len;
    int top_primer_rear_flank_len;
    int bottom_primer_front_flank_len;
    int bottom_primer_rear_flank_len;
    std::vector<std::string> adapter_name;
    std::string kit;
};

struct ScoreResults {
    float score = -1.f;
    float top_score = -1.f;
    float bottom_score = -1.f;
    float flank_score = -1.f;
    float top_flank_score = -1.f;
    float bottom_flank_score = -1.f;
    bool use_top = false;
    std::string adapter_name;
    std::string kit;
    int barcode_start = -1;
};

const ScoreResults UNCLASSIFIED = {-1.f, -1.f, -1.f,           -1.f,
                                   -1.f, -1.f, "unclassified", "unclassified"};

static const std::unordered_map<std::string, KitInfo> kit_info_map = {
        {"SQK-RBK004",
         {false,
          false,
          "GCTTGGGTGTTTAACC",
          "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
          "",
          "",
          {"BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11",
           "BC12"}}},
        {"SQK-RBK110-96",
         {false,
          false,
          "GCTTGGGTGTTTAACC",
          "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
          "",
          "",
          {"BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11",
           "BC12", "BC13", "BC14", "BC15", "BC16", "BC17", "BC18", "BC19", "BC20", "BC21", "BC22",
           "BC23", "BC24", "BC25", "BC26", "BC27", "BC28", "BC29", "BC30", "BC31", "BC32", "BC33",
           "BC34", "BC35", "BC36", "BC37", "BC38", "BC39", "BC40", "BC41", "BC42", "BC43", "BC44",
           "BC45", "BC46", "BC47", "BC48", "BC49", "BC50", "BC51", "BC52", "BC53", "BC54", "BC55",
           "BC56", "BC57", "BC58", "BC59", "BC60", "BC61", "BC62", "BC63", "BC64", "BC65", "BC66",
           "BC67", "BC68", "BC69", "BC70", "BC71", "BC72", "BC73", "BC74", "BC75", "BC76", "BC77",
           "BC78", "BC79", "BC80", "BC81", "BC82", "BC83", "BC84", "BC85", "BC86", "BC87", "BC88",
           "BC89", "BC90", "BC91", "BC92", "BC93", "BC94", "BC95", "BC96"}}},
        {"SQK-RBK114-24",
         {false,
          false,
          "C",
          "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
          "",
          "",
          {"BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08",
           "BC09", "BC10", "BC11", "BC12", "BC13", "BC14", "BC15", "BC16",
           "BC17", "BC18", "BC19", "BC20", "BC21", "BC22", "BC23", "BC24"}}},
        {"SQK-RBK114-96",
         {false,
          false,
          "C",
          "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
          "",
          "",
          {"BC01", "BC02", "BC03", "BC04",  "BC05", "BC06",  "BC07", "BC08",  "BC09",  "BC10",
           "BC11", "BC12", "BC13", "BC14",  "BC15", "BC16",  "BC17", "BC18",  "BC19",  "BC20",
           "BC21", "BC22", "BC23", "BC24",  "BC25", "RBK26", "BC27", "BC28",  "BC29",  "BC30",
           "BC31", "BC32", "BC33", "BC34",  "BC35", "BC36",  "BC37", "BC38",  "RBK39", "RBK40",
           "BC41", "BC42", "BC43", "BC44",  "BC45", "BC46",  "BC47", "RBK48", "BC49",  "BC50",
           "BC51", "BC52", "BC53", "RBK54", "BC55", "BC56",  "BC57", "BC58",  "BC59",  "RBK60",
           "BC61", "BC62", "BC63", "BC64",  "BC65", "BC66",  "BC67", "BC68",  "BC69",  "BC70",
           "BC71", "BC72", "BC73", "BC74",  "BC75", "BC76",  "BC77", "BC78",  "BC79",  "BC80",
           "BC81", "BC82", "BC83", "BC84",  "BC85", "BC86",  "BC87", "BC88",  "BC89",  "BC90",
           "BC91", "BC92", "BC93", "BC94",  "BC95", "BC96"}}},
        {"SQK-RPB004",
         {true,
          false,
          "CCGTGAC",
          "CGTTTTTCGTGCGCCGCTTC",
          "",
          "",
          {"BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11",
           "RLB12A"}}},
        {"SQK-PBK004",
         {true,
          true,
          "ATCGCCTACCGTGAC",
          "ACTTGCCTGTCGCTCTATCTTC",
          "ATCGCCTACCGTGAC",
          "TTTCTGTTGGTGCTGATATTGC",
          {"BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11",
           "BC12"}}},
        {"SQK-RAB204",
         {true,
          true,
          "ATCGCCTACCGTGAC",
          "AGAGTTTGATCMTGGCTCAG",
          "ATCGCCTACCGTGAC",
          "CGGTTACCTTGTTACGACTT",
          {"BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11",
           "BC12"}}},
        {"SQK-16S024",
         {true,
          true,
          "ATCGCCTACCGTGAC",
          "AGAGTTTGATCMTGGCTCAG",
          "ATCGCCTACCGTGAC",
          "CGGTTACCTTGTTACGACTT",
          {"BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08",
           "BC09", "BC10", "BC11", "BC12", "BC13", "BC14", "BC15", "BC16",
           "BC17", "BC18", "BC19", "BC20", "BC21", "BC22", "BC23", "BC24"}}},
        {"SQK-PCB109",
         {true,
          true,
          "ATCGCCTACCGTGAC",
          "ACTTGCCTGTCGCTCTATCTTC",
          "ATCGCCTACCGTGAC",
          "TTTCTGTTGGTGCTGATATTGC",
          {"BP01", "BP02", "BP03", "BP04", "BP05", "BP06", "BP07", "BP08",
           "BP09", "BP10", "BP11", "BP12", "BP13", "BP14", "BP15", "BP16",
           "BP17", "BP18", "BP19", "BP20", "BP21", "BP22", "BP23", "BP24"}}},
        {"SQK-PCB111-24",
         {true,
          true,
          "ATCGCCTACCGTGA",
          "TTGCCTGTCGCTCTATCTTC",
          "ATCGCCTACCGTGA",
          "TCTGTTGGTGCTGATATTGC",
          {"BP01", "BP02", "BP03", "BP04", "BP05", "BP06", "BP07", "BP08",
           "BP09", "BP10", "BP11", "BP12", "BP13", "BP14", "BP15", "BP16",
           "BP17", "BP18", "BP19", "BP20", "BP21", "BP22", "BP23", "BP24"}}},
        {"EXP-PBC096",
         {true,
          true,
          "GGTGCTG",
          "TTAACCTTTCTGTTGGTGCTGATATTGC",
          "GGTGCTG",
          "TTAACCTACTTGCCTGTCGCTCTATCTTC",
          {"BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11",
           "BC12", "BC13", "BC14", "BC15", "BC16", "BC17", "BC18", "BC19", "BC20", "BC21", "BC22",
           "BC23", "BC24", "BC25", "BC26", "BC27", "BC28", "BC29", "BC30", "BC31", "BC32", "BC33",
           "BC34", "BC35", "BC36", "BC37", "BC38", "BC39", "BC40", "BC41", "BC42", "BC43", "BC44",
           "BC45", "BC46", "BC47", "BC48", "BC49", "BC50", "BC51", "BC52", "BC53", "BC54", "BC55",
           "BC56", "BC57", "BC58", "BC59", "BC60", "BC61", "BC62", "BC63", "BC64", "BC65", "BC66",
           "BC67", "BC68", "BC69", "BC70", "BC71", "BC72", "BC73", "BC74", "BC75", "BC76", "BC77",
           "BC78", "BC79", "BC80", "BC81", "BC82", "BC83", "BC84", "BC85", "BC86", "BC87", "BC88",
           "BC89", "BC90", "BC91", "BC92", "BC93", "BC94", "BC95", "BC96"}}},
        {"SQK-NBD114-24",
         {true,
          true,
          "ATCGCCTACCGTGA",
          "TTGCCTGTCGCTCTATCTTC",
          "ATCGCCTACCGTGA",
          "TCTGTTGGTGCTGATATTGC",
          {"NB01", "NB02", "NB03", "NB04", "NB05", "NB06", "NB07", "NB08",
           "NB09", "NB10", "NB11", "NB12", "NB13", "NB14", "NB15", "NB16",
           "NB17", "NB18", "NB19", "NB20", "NB21", "NB22", "NB23", "NB24"}}},
        {"EXP-NBD104",
         {true,
          true,
          "AAGGTTAA",
          "CAGCACCT",
          "ATTGCTAAGGTTAA",
          "CAGCACC",
          {"NB01", "NB02", "NB03", "NB04", "NB05", "NB06", "NB07", "NB08", "NB09", "NB10", "NB11",
           "NB12"}}},
        {"EXP-NBD114",
         {true,
          true,
          "AAGGTTAA",
          "CAGCACCT",
          "ATTGCTAAGGTTAA",
          "CAGCACC",
          {"NB13", "NB14", "NB15", "NB16", "NB17", "NB18", "NB19", "NB20", "NB21", "NB22", "NB23",
           "NB24"}}},
        {"SQK-NBD114-96",
         {true,
          true,
          "AAGGTTAA",
          "CAGCACCT",
          "ATTGCTAAGGTTAA",
          "CAGCACC",
          {"NB01", "NB02", "NB03", "NB04", "NB05", "NB06", "NB07", "NB08", "NB09", "NB10", "NB11",
           "NB12", "NB13", "NB14", "NB15", "NB16", "NB17", "NB18", "NB19", "NB20", "NB21", "NB22",
           "NB23", "NB24", "NB25", "NB26", "NB27", "NB28", "NB29", "NB30", "NB31", "NB32", "NB33",
           "NB34", "NB35", "NB36", "NB37", "NB38", "NB39", "NB40", "NB41", "NB42", "NB43", "NB44",
           "NB45", "NB46", "NB47", "NB48", "NB49", "NB50", "NB51", "NB52", "NB53", "NB54", "NB55",
           "NB56", "NB57", "NB58", "NB59", "NB60", "NB61", "NB62", "NB63", "NB64", "NB65", "NB66",
           "NB67", "NB68", "NB69", "NB70", "NB71", "NB72", "NB73", "NB74", "NB75", "NB76", "NB77",
           "NB78", "NB79", "NB80", "NB81", "NB82", "NB83", "NB84", "NB85", "NB86", "NB87", "NB88",
           "NB89", "NB90", "NB91", "NB92", "NB93", "NB94", "NB95", "NB96"}}},
        {"EXP-NBD196",
         {true,
          true,
          "AAGGTTAA",
          "CAGCACCT",
          "ATTGCTAAGGTTAA",
          "CAGCACC",
          {"NB01", "NB02", "NB03", "NB04", "NB05", "NB06", "NB07", "NB08", "NB09", "NB10", "NB11",
           "NB12", "NB13", "NB14", "NB15", "NB16", "NB17", "NB18", "NB19", "NB20", "NB21", "NB22",
           "NB23", "NB24", "NB25", "NB26", "NB27", "NB28", "NB29", "NB30", "NB31", "NB32", "NB33",
           "NB34", "NB35", "NB36", "NB37", "NB38", "NB39", "NB40", "NB41", "NB42", "NB43", "NB44",
           "NB45", "NB46", "NB47", "NB48", "NB49", "NB50", "NB51", "NB52", "NB53", "NB54", "NB55",
           "NB56", "NB57", "NB58", "NB59", "NB60", "NB61", "NB62", "NB63", "NB64", "NB65", "NB66",
           "NB67", "NB68", "NB69", "NB70", "NB71", "NB72", "NB73", "NB74", "NB75", "NB76", "NB77",
           "NB78", "NB79", "NB80", "NB81", "NB82", "NB83", "NB84", "NB85", "NB86", "NB87", "NB88",
           "NB89", "NB90", "NB91", "NB92", "NB93", "NB94", "NB95", "NB96"}}},
};

static const std::unordered_map<std::string, std::string> barcodes = {
        // BC** barcodes.
        {"BC01", "AAGAAAGTTGTCGGTGTCTTTGTG"},
        {"BC02", "TCGATTCCGTTTGTAGTCGTCTGT"},
        {"BC03", "GAGTCTTGTGTCCCAGTTACCAGG"},
        {"BC04", "TTCGGATTCTATCGTGTTTCCCTA"},
        {"BC05", "CTTGTCCAGGGTTTGTGTAACCTT"},
        {"BC06", "TTCTCGCAAAGGCAGAAAGTAGTC"},
        {"BC07", "GTGTTACCGTGGGAATGAATCCTT"},
        {"BC08", "TTCAGGGAACAAACCAAGTTACGT"},
        {"BC09", "AACTAGGCACAGCGAGTCTTGGTT"},
        {"BC10", "AAGCGTTGAAACCTTTGTCCTCTC"},
        {"BC11", "GTTTCATCTATCGGAGGGAATGGA"},
        {"BC12", "CAGGTAGAAAGAAGCAGAATCGGA"},
        {"RLB12A", "GTTGAGTTACAAAGCACCGATCAG"},
        {"BC13", "AGAACGACTTCCATACTCGTGTGA"},
        {"BC14", "AACGAGTCTCTTGGGACCCATAGA"},
        {"BC15", "AGGTCTACCTCGCTAACACCACTG"},
        {"BC16", "CGTCAACTGACAGTGGTTCGTACT"},
        {"BC17", "ACCCTCCAGGAAAGTACCTCTGAT"},
        {"BC18", "CCAAACCCAACAACCTAGATAGGC"},
        {"BC19", "GTTCCTCGTGCAGTGTCAAGAGAT"},
        {"BC20", "TTGCGTCCTGTTACGAGAACTCAT"},
        {"BC21", "GAGCCTCTCATTGTCCGTTCTCTA"},
        {"BC22", "ACCACTGCCATGTATCAAAGTACG"},
        {"BC23", "CTTACTACCCAGTGAACCTCCTCG"},
        {"BC24", "GCATAGTTCTGCATGATGGGTTAG"},
        {"BC25", "GTAAGTTGGGTATGCAACGCAATG"},
        {"BC26", "CATACAGCGACTACGCATTCTCAT"},
        {"RBK26", "ACTATGCCTTTCCGTGAAACAGTT"},
        {"BC27", "CGACGGTTAGATTCACCTCTTACA"},
        {"BC28", "TGAAACCTAAGAAGGCACCGTATC"},
        {"BC29", "CTAGACACCTTGGGTTGACAGACC"},
        {"BC30", "TCAGTGAGGATCTACTTCGACCCA"},
        {"BC31", "TGCGTACAGCAATCAGTTACATTG"},
        {"BC32", "CCAGTAGAAGTCCGACAACGTCAT"},
        {"BC33", "CAGACTTGGTACGGTTGGGTAACT"},
        {"BC34", "GGACGAAGAACTCAAGTCAAAGGC"},
        {"BC35", "CTACTTACGAAGCTGAGGGACTGC"},
        {"BC36", "ATGTCCCAGTTAGAGGAGGAAACA"},
        {"BC37", "GCTTGCGATTGATGCTTAGTATCA"},
        {"BC38", "ACCACAGGAGGACGATACAGAGAA"},
        {"BC39", "CCACAGTGTCAACTAGAGCCTCTC"},
        {"RBK39", "TCTGCCACACACTCGTAAGTCCTT"},
        {"BC40", "TAGTTTGGATGACCAAGGATAGCC"},
        {"RBK40", "GTCGATACTGGACCTATCCCTTGG"},
        {"BC41", "GGAGTTCGTCCAGAGAAGTACACG"},
        {"BC42", "CTACGTGTAAGGCATACCTGCCAG"},
        {"BC43", "CTTTCGTTGTTGACTCGACGGTAG"},
        {"BC44", "AGTAGAAAGGGTTCCTTCCCACTC"},
        {"BC45", "GATCCAACAGAGATGCCTTCAGTG"},
        {"BC46", "GCTGTGTTCCACTTCATTCTCCTG"},
        {"BC47", "GTGCAACTTTCCCACAGGTAGTTC"},
        {"BC48", "CATCTGGAACGTGGTACACCTGTA"},
        {"RBK48", "GAGTCCGTGACAACTTCTGAAAGC"},
        {"BC49", "ACTGGTGCAGCTTTGAACATCTAG"},
        {"BC50", "ATGGACTTTGGTAACTTCCTGCGT"},
        {"BC51", "GTTGAATGAGCCTACTGGGTCCTC"},
        {"BC52", "TGAGAGACAAGATTGTTCGTGGAC"},
        {"BC53", "AGATTCAGACCGTCTCATGCAAAG"},
        {"BC54", "CAAGAGCTTTGACTAAGGAGCATG"},
        {"RBK54", "GGGTGCCAACTACATACCAAACCT"},
        {"BC55", "TGGAAGATGAGACCCTGATCTACG"},
        {"BC56", "TCACTACTCAACAGGTGGCATGAA"},
        {"BC57", "GCTAGGTCAATCTCCTTCGGAAGT"},
        {"BC58", "CAGGTTACTCCTCCGTGAGTCTGA"},
        {"BC59", "TCAATCAAGAAGGGAAAGCAAGGT"},
        {"BC60", "CATGTTCAACCAAGGCTTCTATGG"},
        {"RBK60", "GAACCCTACTTTGGACAGACACCT"},
        {"BC61", "AGAGGGTACTATGTGCCTCAGCAC"},
        {"BC62", "CACCCACACTTACTTCAGGACGTA"},
        {"BC63", "TTCTGAAGTTCCTGGGTCTTGAAC"},
        {"BC64", "GACAGACACCGTTCATCGACTTTC"},
        {"BC65", "TTCTCAGTCTTCCTCCAGACAAGG"},
        {"BC66", "CCGATCCTTGTGGCTTCTAACTTC"},
        {"BC67", "GTTTGTCATACTCGTGTGCTCACC"},
        {"BC68", "GAATCTAAGCAAACACGAAGGTGG"},
        {"BC69", "TACAGTCCGAGCCTCATGTGATCT"},
        {"BC70", "ACCGAGATCCTACGAATGGAGTGT"},
        {"BC71", "CCTGGGAGCATCAGGTAGTAACAG"},
        {"BC72", "TAGCTGACTGTCTTCCATACCGAC"},
        {"BC73", "AAGAAACAGGATGACAGAACCCTC"},
        {"BC74", "TACAAGCATCCCAACACTTCCACT"},
        {"BC75", "GACCATTGTGATGAACCCTGTTGT"},
        {"BC76", "ATGCTTGTTACATCAACCCTGGAC"},
        {"BC77", "CGACCTGTTTCTCAGGGATACAAC"},
        {"BC78", "AACAACCGAACCTTTGAATCAGAA"},
        {"BC79", "TCTCGGAGATAGTTCTCACTGCTG"},
        {"BC80", "CGGATGAACATAGGATAGCGATTC"},
        {"BC81", "CCTCATCTTGTGAAGTTGTTTCGG"},
        {"BC82", "ACGGTATGTCGAGTTCCAGGACTA"},
        {"BC83", "TGGCTTGATCTAGGTAAGGTCGAA"},
        {"BC84", "GTAGTGGACCTAGAACCTGTGCCA"},
        {"BC85", "AACGGAGGAGTTAGTTGGATGATC"},
        {"BC86", "AGGTGATCCCAACAAGCGTAAGTA"},
        {"BC87", "TACATGCTCCTGTTGTTAGGGAGG"},
        {"BC88", "TCTTCTACTACCGATCCGAAGCAG"},
        {"BC89", "ACAGCATCAATGTTTGGCTAGTTG"},
        {"BC90", "GATGTAGAGGGTACGGTTTGAGGC"},
        {"BC91", "GGCTCCATAGGAACTCACGCTACT"},
        {"BC92", "TTGTGAGTGGAAAGATACAGGACC"},
        {"BC93", "AGTTTCCATCACTTCAGACTTGGG"},
        {"BC94", "GATTGTCCTCAAACTGCCACCTAC"},
        {"BC95", "CCTGTCTGGAAGAAGAATGGACTT"},
        {"BC96", "CTGAACGGTCATAGAGTCCACCAT"},
        // BP** barcodes.
        {"BP01", "CAAGAAAGTTGTCGGTGTCTTTGTGAC"},
        {"BP02", "CTCGATTCCGTTTGTAGTCGTCTGTAC"},
        {"BP03", "CGAGTCTTGTGTCCCAGTTACCAGGAC"},
        {"BP04", "CTTCGGATTCTATCGTGTTTCCCTAAC"},
        {"BP05", "CCTTGTCCAGGGTTTGTGTAACCTTAC"},
        {"BP06", "CTTCTCGCAAAGGCAGAAAGTAGTCAC"},
        {"BP07", "CGTGTTACCGTGGGAATGAATCCTTAC"},
        {"BP08", "CTTCAGGGAACAAACCAAGTTACGTAC"},
        {"BP09", "CAACTAGGCACAGCGAGTCTTGGTTAC"},
        {"BP10", "CAAGCGTTGAAACCTTTGTCCTCTCAC"},
        {"BP11", "CGTTTCATCTATCGGAGGGAATGGAAC"},
        {"BP12", "CCAGGTAGAAAGAAGCAGAATCGGAAC"},
        {"BP13", "CAGAACGACTTCCATACTCGTGTGAAC"},
        {"BP14", "CAACGAGTCTCTTGGGACCCATAGAAC"},
        {"BP15", "CAGGTCTACCTCGCTAACACCACTGAC"},
        {"BP16", "CCGTCAACTGACAGTGGTTCGTACTAC"},
        {"BP17", "CACCCTCCAGGAAAGTACCTCTGATAC"},
        {"BP18", "CCCAAACCCAACAACCTAGATAGGCAC"},
        {"BP19", "CGTTCCTCGTGCAGTGTCAAGAGATAC"},
        {"BP20", "CTTGCGTCCTGTTACGAGAACTCATAC"},
        {"BP21", "CGAGCCTCTCATTGTCCGTTCTCTAAC"},
        {"BP22", "CACCACTGCCATGTATCAAAGTACGAC"},
        {"BP23", "CCTTACTACCCAGTGAACCTCCTCGAC"},
        {"BP24", "CGCATAGTTCTGCATGATGGGTTAGAC"},
        // NB** barcodes.
        {"NB01", "CACAAAGACACCGACAACTTTCTT"},
        {"NB02", "ACAGACGACTACAAACGGAATCGA"},
        {"NB03", "CCTGGTAACTGGGACACAAGACTC"},
        {"NB04", "TAGGGAAACACGATAGAATCCGAA"},
        {"NB05", "AAGGTTACACAAACCCTGGACAAG"},
        {"NB06", "GACTACTTTCTGCCTTTGCGAGAA"},
        {"NB07", "AAGGATTCATTCCCACGGTAACAC"},
        {"NB08", "ACGTAACTTGGTTTGTTCCCTGAA"},
        {"NB09", "AACCAAGACTCGCTGTGCCTAGTT"},
        {"NB10", "GAGAGGACAAAGGTTTCAACGCTT"},
        {"NB11", "TCCATTCCCTCCGATAGATGAAAC"},
        {"NB12", "TCCGATTCTGCTTCTTTCTACCTG"},
        {"NB13", "AGAACGACTTCCATACTCGTGTGA"},
        {"NB14", "AACGAGTCTCTTGGGACCCATAGA"},
        {"NB15", "AGGTCTACCTCGCTAACACCACTG"},
        {"NB16", "CGTCAACTGACAGTGGTTCGTACT"},
        {"NB17", "ACCCTCCAGGAAAGTACCTCTGAT"},
        {"NB18", "CCAAACCCAACAACCTAGATAGGC"},
        {"NB19", "GTTCCTCGTGCAGTGTCAAGAGAT"},
        {"NB20", "TTGCGTCCTGTTACGAGAACTCAT"},
        {"NB21", "GAGCCTCTCATTGTCCGTTCTCTA"},
        {"NB22", "ACCACTGCCATGTATCAAAGTACG"},
        {"NB23", "CTTACTACCCAGTGAACCTCCTCG"},
        {"NB24", "GCATAGTTCTGCATGATGGGTTAG"},
        {"NB25", "GTAAGTTGGGTATGCAACGCAATG"},
        {"NB26", "CATACAGCGACTACGCATTCTCAT"},
        {"NB27", "CGACGGTTAGATTCACCTCTTACA"},
        {"NB28", "TGAAACCTAAGAAGGCACCGTATC"},
        {"NB29", "CTAGACACCTTGGGTTGACAGACC"},
        {"NB30", "TCAGTGAGGATCTACTTCGACCCA"},
        {"NB31", "TGCGTACAGCAATCAGTTACATTG"},
        {"NB32", "CCAGTAGAAGTCCGACAACGTCAT"},
        {"NB33", "CAGACTTGGTACGGTTGGGTAACT"},
        {"NB34", "GGACGAAGAACTCAAGTCAAAGGC"},
        {"NB35", "CTACTTACGAAGCTGAGGGACTGC"},
        {"NB36", "ATGTCCCAGTTAGAGGAGGAAACA"},
        {"NB37", "GCTTGCGATTGATGCTTAGTATCA"},
        {"NB38", "ACCACAGGAGGACGATACAGAGAA"},
        {"NB39", "CCACAGTGTCAACTAGAGCCTCTC"},
        {"NB40", "TAGTTTGGATGACCAAGGATAGCC"},
        {"NB41", "GGAGTTCGTCCAGAGAAGTACACG"},
        {"NB42", "CTACGTGTAAGGCATACCTGCCAG"},
        {"NB43", "CTTTCGTTGTTGACTCGACGGTAG"},
        {"NB44", "AGTAGAAAGGGTTCCTTCCCACTC"},
        {"NB45", "GATCCAACAGAGATGCCTTCAGTG"},
        {"NB46", "GCTGTGTTCCACTTCATTCTCCTG"},
        {"NB47", "GTGCAACTTTCCCACAGGTAGTTC"},
        {"NB48", "CATCTGGAACGTGGTACACCTGTA"},
        {"NB49", "ACTGGTGCAGCTTTGAACATCTAG"},
        {"NB50", "ATGGACTTTGGTAACTTCCTGCGT"},
        {"NB51", "GTTGAATGAGCCTACTGGGTCCTC"},
        {"NB52", "TGAGAGACAAGATTGTTCGTGGAC"},
        {"NB53", "AGATTCAGACCGTCTCATGCAAAG"},
        {"NB54", "CAAGAGCTTTGACTAAGGAGCATG"},
        {"NB55", "TGGAAGATGAGACCCTGATCTACG"},
        {"NB56", "TCACTACTCAACAGGTGGCATGAA"},
        {"NB57", "GCTAGGTCAATCTCCTTCGGAAGT"},
        {"NB58", "CAGGTTACTCCTCCGTGAGTCTGA"},
        {"NB59", "TCAATCAAGAAGGGAAAGCAAGGT"},
        {"NB60", "CATGTTCAACCAAGGCTTCTATGG"},
        {"NB61", "AGAGGGTACTATGTGCCTCAGCAC"},
        {"NB62", "CACCCACACTTACTTCAGGACGTA"},
        {"NB63", "TTCTGAAGTTCCTGGGTCTTGAAC"},
        {"NB64", "GACAGACACCGTTCATCGACTTTC"},
        {"NB65", "TTCTCAGTCTTCCTCCAGACAAGG"},
        {"NB66", "CCGATCCTTGTGGCTTCTAACTTC"},
        {"NB67", "GTTTGTCATACTCGTGTGCTCACC"},
        {"NB68", "GAATCTAAGCAAACACGAAGGTGG"},
        {"NB69", "TACAGTCCGAGCCTCATGTGATCT"},
        {"NB70", "ACCGAGATCCTACGAATGGAGTGT"},
        {"NB71", "CCTGGGAGCATCAGGTAGTAACAG"},
        {"NB72", "TAGCTGACTGTCTTCCATACCGAC"},
        {"NB73", "AAGAAACAGGATGACAGAACCCTC"},
        {"NB74", "TACAAGCATCCCAACACTTCCACT"},
        {"NB75", "GACCATTGTGATGAACCCTGTTGT"},
        {"NB76", "ATGCTTGTTACATCAACCCTGGAC"},
        {"NB77", "CGACCTGTTTCTCAGGGATACAAC"},
        {"NB78", "AACAACCGAACCTTTGAATCAGAA"},
        {"NB79", "TCTCGGAGATAGTTCTCACTGCTG"},
        {"NB80", "CGGATGAACATAGGATAGCGATTC"},
        {"NB81", "CCTCATCTTGTGAAGTTGTTTCGG"},
        {"NB82", "ACGGTATGTCGAGTTCCAGGACTA"},
        {"NB83", "TGGCTTGATCTAGGTAAGGTCGAA"},
        {"NB84", "GTAGTGGACCTAGAACCTGTGCCA"},
        {"NB85", "AACGGAGGAGTTAGTTGGATGATC"},
        {"NB86", "AGGTGATCCCAACAAGCGTAAGTA"},
        {"NB87", "TACATGCTCCTGTTGTTAGGGAGG"},
        {"NB88", "TCTTCTACTACCGATCCGAAGCAG"},
        {"NB89", "ACAGCATCAATGTTTGGCTAGTTG"},
        {"NB90", "GATGTAGAGGGTACGGTTTGAGGC"},
        {"NB91", "GGCTCCATAGGAACTCACGCTACT"},
        {"NB92", "TTGTGAGTGGAAAGATACAGGACC"},
        {"NB93", "AGTTTCCATCACTTCAGACTTGGG"},
        {"NB94", "GATTGTCCTCAAACTGCCACCTAC"},
        {"NB95", "CCTGTCTGGAAGAAGAATGGACTT"},
        {"NB96", "CTGAACGGTCATAGAGTCCACCAT"},
};

static std::string barcode_kits_list_str() {
    return std::accumulate(kit_info_map.begin(), kit_info_map.end(), std::string(),
                           [](std::string& a, auto& b) -> std::string {
                               return a + (a.empty() ? "" : " ") + b.first;
                           });
}

class BarcodeClassifier {
public:
    BarcodeClassifier(const std::vector<std::string>& kit_names);
    ~BarcodeClassifier() = default;

    ScoreResults barcode(const std::string& seq);

private:
    std::vector<AdapterSequence> m_adapter_sequences;

    std::vector<AdapterSequence> generate_adapter_sequence(
            const std::vector<std::string>& kit_names);
    void calculate_adapter_score_different_double_ends(std::string_view read_seq,
                                                       const AdapterSequence& as,
                                                       std::vector<ScoreResults>& res);
    void calculate_adapter_score_double_ends(std::string_view read_seq,
                                             const AdapterSequence& as,
                                             std::vector<ScoreResults>& res);
    void calculate_adapter_score(std::string_view read_seq,
                                 const AdapterSequence& as,
                                 std::vector<ScoreResults>& res);
    ScoreResults find_best_adapter(const std::string& read_seq,
                                   std::vector<AdapterSequence>& adapter);
};

}  // namespace demux

}  // namespace dorado

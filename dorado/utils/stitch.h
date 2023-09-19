#pragma once
#include <memory>

namespace dorado {
class Read;
}  // namespace dorado

namespace dorado::utils {

// Given a read with unstitched chunks, stitch the chunks (accounting for overlap) and assign basecalled read and
// qstring to Read
void stitch_chunks(std::shared_ptr<Read> read);

}  // namespace dorado::utils
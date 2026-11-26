#include <fastbdf/fastbdf.hpp>

std::int64_t nanoseconds_in_one_second()
{
    return fastbdf::RecordingInfo{}.record_duration.count();
}

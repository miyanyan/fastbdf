#include <fastbdf/fastbdf.hpp>

std::int64_t ticks_in_one_second() {
  return fastbdf::Tick(std::chrono::seconds(1)).count();
}

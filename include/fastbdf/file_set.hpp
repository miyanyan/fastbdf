#pragma once

#include "reader.hpp"
#include "writer.hpp"
#include <algorithm>
#include <filesystem>
#include <list>
#include <variant>

namespace fastbdf {

// Optional application-owned handle table for clients migrating EDFlib file
// enumeration. References remain valid until their entry is closed.
class FileSet {
public:
  using Handle = std::uint64_t;

  Handle open(const std::filesystem::path &path,
              AnnotationMode mode = AnnotationMode::lazy) {
    auto key = available(path);
    files_.push_back({next_, key, Reader(path, mode)});
    return next_++;
  }

  Handle create(const std::filesystem::path &path, RecordingInfo info,
                WriterOptions options = {}) {
    auto key = available(path);
    files_.push_back(
        {next_, key, Writer::create(path, std::move(info), options)});
    return next_++;
  }

  Reader &reader(Handle handle) {
    auto &entry = *find(handle);
    if (auto *reader = std::get_if<Reader>(&entry.file))
      return *reader;
    throw Error("File handle is open for writing");
  }
  Writer &writer(Handle handle) {
    auto &entry = *find(handle);
    if (auto *writer = std::get_if<Writer>(&entry.file))
      return *writer;
    throw Error("File handle is open for reading");
  }
  void close(Handle handle) {
    auto entry = find(handle);
    if (auto *writer = std::get_if<Writer>(&entry->file))
      writer->finish();
    files_.erase(entry);
  }

  std::size_t size() const noexcept { return files_.size(); }
  bool contains(const std::filesystem::path &path) const {
    auto key =
        std::filesystem::weakly_canonical(std::filesystem::absolute(path));
    return std::any_of(files_.begin(), files_.end(), [&](const Entry &e) {
      std::error_code ignored;
      return e.path == key || std::filesystem::equivalent(e.path, key, ignored);
    });
  }
  std::vector<Handle> handles() const {
    std::vector<Handle> ids;
    ids.reserve(files_.size());
    for (const auto &entry : files_)
      ids.push_back(entry.id);
    return ids;
  }

private:
  struct Entry {
    Handle id;
    std::filesystem::path path;
    std::variant<Reader, Writer> file;
  };
  std::list<Entry> files_;
  Handle next_{1};

  std::filesystem::path available(const std::filesystem::path &path) const {
    if (next_ == 0)
      throw Error("File handle space exhausted");
    if (contains(path))
      throw Error("Path is already open in this FileSet");
    return std::filesystem::weakly_canonical(std::filesystem::absolute(path));
  }
  std::list<Entry>::iterator find(Handle handle) {
    auto it = std::find_if(files_.begin(), files_.end(),
                           [=](const Entry &e) { return e.id == handle; });
    if (it == files_.end())
      throw Error("Unknown file handle");
    return it;
  }
};
} // namespace fastbdf

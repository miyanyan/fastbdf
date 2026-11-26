#pragma once

#include <algorithm>
#include <cerrno>
#include <filesystem>
#include <limits>
#include <span>
#include <system_error>
#include <utility>
#include "../types.hpp"

#ifdef _WIN32
#    ifndef NOMINMAX
#        define NOMINMAX
#        define FASTBDF_UNDEF_NOMINMAX
#    endif
#    include <windows.h>
#    ifdef FASTBDF_UNDEF_NOMINMAX
#        undef NOMINMAX
#        undef FASTBDF_UNDEF_NOMINMAX
#    endif
#else
#    include <fcntl.h>
#    include <sys/stat.h>
#    include <unistd.h>
#endif

namespace fastbdf::detail
{

// Each request supplies an absolute offset. Independent reads on the same file
// never observe a shared file position, including on Windows.
class FileSource
{
public:
    explicit FileSource(const std::filesystem::path& path)
    {
#ifdef _WIN32
        handle_ = CreateFileW(path.c_str(), GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
                              nullptr, OPEN_EXISTING, FILE_FLAG_OVERLAPPED, nullptr);
        if (handle_ == INVALID_HANDLE_VALUE) fail("open", GetLastError());
        LARGE_INTEGER length;
        if (!GetFileSizeEx(handle_, &length)) {
            auto error = GetLastError();
            close();
            fail("get size", error);
        }
        size_ = static_cast<std::uint64_t>(length.QuadPart);
#else
        handle_ = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
        if (handle_ == -1) fail("open", errno);
        struct stat status{};
        if (::fstat(handle_, &status) < 0) {
            auto error = errno;
            close();
            fail("get size", error);
        }
        size_ = static_cast<std::uint64_t>(status.st_size);
#endif
    }

    FileSource(FileSource&& other) noexcept
        : handle_(std::exchange(other.handle_, invalid()))
        , size_(other.size_)
    {
    }
    FileSource& operator=(FileSource&& other) noexcept
    {
        if (this != &other) {
            close();
            handle_ = std::exchange(other.handle_, invalid());
            size_ = other.size_;
        }
        return *this;
    }
    FileSource(const FileSource&) = delete;
    FileSource& operator=(const FileSource&) = delete;
    ~FileSource() { close(); }

    std::uint64_t size() const noexcept { return size_; }

    void read_at(std::uint64_t offset, std::span<char> output) const
    {
        if (offset > size_ || output.size() > size_ - offset)
            throw Error("Read exceeds file length at byte " + std::to_string(offset), Errc::invalid_file);
        while (!output.empty()) {
            const auto count = (std::min)(output.size(), std::size_t{1} << 30);
#ifdef _WIN32
      // A separate event is required for each simultaneous OVERLAPPED
      // request. Reuse it across reads on this thread.
            struct Event
            {
                HANDLE value = CreateEventW(nullptr, TRUE, FALSE, nullptr);
                ~Event()
                {
                    if (value) CloseHandle(value);
                }
            };
            thread_local Event event;
            if (!event.value) fail("create read event", GetLastError());
            OVERLAPPED request{};
            request.Offset = static_cast<DWORD>(offset);
            request.OffsetHigh = static_cast<DWORD>(offset >> 32);
            request.hEvent = event.value;
            DWORD transferred{};
            if (!ReadFile(handle_, output.data(), static_cast<DWORD>(count), &transferred, &request)) {
                auto error = GetLastError();
                if (error != ERROR_IO_PENDING) fail("read", error);
                if (!GetOverlappedResult(handle_, &request, &transferred, TRUE)) fail("complete read", GetLastError());
            }
#else
            auto transferred = ::pread(handle_, output.data(), count, static_cast<off_t>(offset));
            if (transferred < 0) {
                if (errno == EINTR) continue;
                fail("read", errno);
            }
#endif
            if (!transferred)
                throw Error("Unexpected end of file at byte " + std::to_string(offset), Errc::invalid_file);
            output = output.subspan(static_cast<std::size_t>(transferred));
            offset += static_cast<std::uint64_t>(transferred);
        }
    }

private:
#ifdef _WIN32
    using Native = HANDLE;
    static Native invalid() noexcept { return INVALID_HANDLE_VALUE; }
    void close() noexcept
    {
        if (handle_ != invalid()) CloseHandle(handle_);
    }
#else
    using Native = int;
    static Native invalid() noexcept { return -1; }
    void close() noexcept
    {
        if (handle_ != invalid()) ::close(handle_);
    }
#endif
    Native handle_ = invalid();
    std::uint64_t size_{};
    [[noreturn]] static void fail(const char* operation, unsigned code)
    {
        throw Error(
            std::string("Cannot ") + operation + " file: " + std::system_category().message(static_cast<int>(code)),
            Errc::io);
    }
};
} // namespace fastbdf::detail

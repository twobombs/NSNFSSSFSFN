/* compile with:
 *
 * g++ -std=c++20 search_extqueries.cpp -o search_extqueries -W -Wall -O2
 *
 * Use the MAX_THREADS environment variable to cap the number of threads.
 */
#include <cctype>

#include <array>
#include <atomic>
#include <chrono>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

std::mutex results_mutex;
std::vector<std::pair<std::string, std::string>> extracted_results;
std::atomic<std::uintmax_t> total_bytes_processed {0};
std::atomic<std::size_t> total_extracted {0};

// Maps '0'-'9' to 0-9, and '-' to 10. Returns -1 for invalid characters.
constexpr int char_to_index(char c)
{
    if (c >= '0' && c <= '9')
        return c - '0';
    if (c == '-')
        return 10;
    return -1;
}

struct TrieNode {
    std::array<std::unique_ptr<TrieNode>, 11> children;
    bool is_end_of_word = false;
};

class NumericTrie
{
    std::unique_ptr<TrieNode> root;

  public:
    std::size_t size = 0;
    NumericTrie()
        : root(std::make_unique<TrieNode>())
    {
    }

    void insert(std::string const & key)
    {
        size++;
        TrieNode * current = root.get();
        for (char c: key) {
            int idx = char_to_index(c);
            if (idx == -1)
                continue; // Skip malformed characters

            if (!current->children[idx]) {
                current->children[idx] = std::make_unique<TrieNode>();
            }
            current = current->children[idx].get();
        }
        current->is_end_of_word = true;
    }

    TrieNode const * get_root() const { return root.get(); }
};

void process_chunk_with_trie(std::string const & filename,
                             std::uintmax_t start_offset,
                             std::uintmax_t end_offset,
                             NumericTrie const & trie)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file)
        return;

    file.seekg(start_offset);

    std::uintmax_t last_sync_pos = file.tellg();
    std::uintmax_t const update_threshold = 1024 * 1024;

    // Synchronize to the start of the next key/value pair
    if (start_offset > 0) {
        char c;
        while (file.get(c)) {
            if (c == ',' || c == '{')
                break;
        }
    }

    auto read_until = [&](char target) -> bool {
        char c;
        while (file.get(c)) {
            if (c == target)
                return true;
        }
        return false;
    };

    std::string key_buffer;
    std::string val_buffer;
    key_buffer.reserve(256);
    val_buffer.reserve(1024);

    while (file.tellg() != -1 &&
           static_cast<std::uintmax_t>(file.tellg()) < end_offset) {
        std::uintmax_t current_pos = file.tellg();
        if (current_pos - last_sync_pos >= update_threshold) {
            total_bytes_processed.fetch_add(current_pos - last_sync_pos,
                                            std::memory_order_relaxed);
            last_sync_pos = current_pos;
        }

        if (!read_until('"'))
            break;

        key_buffer.clear();
        TrieNode const * current_node = trie.get_root();
        bool is_match = true;
        char c;

        // Stream through the key character by character
        while (file.get(c) && c != '"') {
            if (is_match) {
                key_buffer += c; // Only append while we are matching
                int idx = char_to_index(c);

                if (idx != -1 && current_node->children[idx]) {
                    current_node = current_node->children[idx].get();
                } else {
                    is_match = false;
                }
            }
        }

        if (is_match && current_node && current_node->is_end_of_word) {
            if (!read_until(':'))
                break;

            val_buffer.clear();
            bool in_quotes = false;
            bool value_started = false;

            while (file.get(c)) {
                if (!value_started) {
                    if (std::isspace(c))
                        continue;
                    value_started = true;
                    if (c == '"') {
                        in_quotes = true;
                        continue;
                    }
                } else if (in_quotes) {
                    if (c == '"')
                        break;
                } else {
                    if (c == ',' || c == '}' || std::isspace(c)) {
                        file.unget();
                        break;
                    }
                }
                if (value_started) {
                    val_buffer += c;
                }
            }

            std::lock_guard<std::mutex> lock(results_mutex);
            extracted_results.emplace_back(key_buffer, val_buffer);
            total_extracted.fetch_add(1, std::memory_order_relaxed);
        }
    }

    std::uintmax_t final_pos;
    if (file.tellg() == -1) {
        final_pos = end_offset;
    } else {
        final_pos = file.tellg();
    }
    if (final_pos > last_sync_pos) {
        total_bytes_processed.fetch_add(final_pos - last_sync_pos,
                                        std::memory_order_relaxed);
    }
}

struct ProgressSnapshot {
    std::uintmax_t bytes;
    std::chrono::steady_clock::time_point time;
};

std::deque<ProgressSnapshot> history;

int main(int argc, char * argv[])
{
    if (argc != 4)
        throw std::runtime_error(
            "Usage: ./search_extqueries <json database> <file with "
            "list of queries> <output file>");

    std::string const input_filename = argv[1];
    std::string const queries_filename = argv[2];
    std::string const output_filename = argv[3];

    NumericTrie T;

    {
        std::ifstream queries(queries_filename);
        for (std::string s; std::getline(queries, s);)
            T.insert(s);
    }

    std::uintmax_t file_size = std::filesystem::file_size(input_filename);
    unsigned int num_threads = std::thread::hardware_concurrency();

    if (auto c = getenv("MAX_THREADS"); c != NULL) {
        unsigned int n = std::atoi(c);
        if (n && n < num_threads)
            num_threads = n;
    }

    if (num_threads == 0)
        num_threads = 4;

    std::vector<std::jthread> threads;
    std::uintmax_t chunk_size = file_size / num_threads;

    std::cout << "Processing " << file_size << " bytes across " << num_threads
              << " threads...\n";

    for (unsigned int i = 0; i < num_threads; ++i) {
        std::uintmax_t start = i * chunk_size;
        std::uintmax_t end =
            (i == num_threads - 1) ? file_size : (i + 1) * chunk_size;

        threads.emplace_back(process_chunk_with_trie, input_filename, start,
                             end, std::cref(T));
    }

    auto start_time = std::chrono::steady_clock::now();
    std::uintmax_t processed = 0;
    history.push_back({0, start_time});

    while (processed < file_size) {
        processed = total_bytes_processed.load(std::memory_order_relaxed);
        auto current_time = std::chrono::steady_clock::now();

        history.push_back({processed, current_time});

        if (history.size() > 11)
            history.pop_front();

        std::chrono::duration<double> total_elapsed = current_time - start_time;

        // Prevent division by zero on the very first quick iterations
        if (total_elapsed.count() > 0.1 && processed > 0) {
            int percentage = (processed * 100) / file_size;

            auto const & oldest = history.front();
            std::chrono::duration<double> window_elapsed =
                current_time - oldest.time;
            double window_elapsed_sec = window_elapsed.count();

            double bytes_per_sec = 0.0;
            if (window_elapsed_sec > 0) {
                bytes_per_sec = (processed - oldest.bytes) / window_elapsed_sec;
            }

            double speed_mbps = bytes_per_sec / (1024.0 * 1024.0);

            std::uintmax_t remaining_bytes = file_size - processed;
            double eta_sec = remaining_bytes / bytes_per_sec;

            int eta_mins = static_cast<int>(eta_sec) / 60;
            int eta_secs_remainder = static_cast<int>(eta_sec) % 60;

            std::cout << "\rProgress: " << percentage << "% "
                      << "| Speed: " << std::fixed << std::setprecision(1)
                      << speed_mbps << " MB/s "
                      << "| ETA: " << eta_mins << "m " << eta_secs_remainder
                      << "s    "
                      << "| Extracted: " << total_extracted.load() << "/"
                      << T.size << "  " << std::flush;
        }

        if (processed >= file_size)
            break;
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }

    std::cout << "\n";

    threads.clear();

    std::ofstream out_file(output_filename);
    if (out_file) {
        out_file << "{\n";
        bool first = true;
        for (auto const & [k, v]: extracted_results) {
            if (!first)
                out_file << ",\n";
            out_file << "  \"" << k << "\":"
                     << " \"" << v << "\"";
            first = false;
        }
        out_file << "\n}\n";
        std::cout << "Successfully wrote " << extracted_results.size()
                  << " pairs to " << output_filename << "\n";
    }

    {
        auto current_time = std::chrono::steady_clock::now();
        std::chrono::duration<double> total_elapsed = current_time - start_time;
        std::chrono::duration<double> elapsed = current_time - start_time;
        double elapsed_sec = elapsed.count();

        double bytes_per_sec = 0.0;
        if (elapsed_sec > 0) {
            bytes_per_sec = total_bytes_processed.load() / elapsed_sec;
        }
        double speed_mbps = bytes_per_sec / (1024.0 * 1024.0);
        std::cout << "Extracted " << total_extracted.load() << " pairs"
                  << std::fixed << std::setprecision(1)
                  << " in " << elapsed_sec << " s"
                  << " (" << speed_mbps << " MB/s)\n";
    }

    return 0;
}

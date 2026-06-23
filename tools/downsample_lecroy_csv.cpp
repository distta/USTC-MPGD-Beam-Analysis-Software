#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <regex>
#include <string>
#include <thread>
#include <vector>

using namespace std;
namespace fs = std::filesystem;

namespace {

bool ParseTime(const string& line, double& time) {
    const char* begin = line.c_str();
    char* end = nullptr;
    time = strtod(begin, &end);
    return end != begin && *end == ',';
}

bool Convert(const fs::path& path, double targetPitchSeconds, string& error) {
    ifstream input(path);
    if (!input) {
        error = "cannot open input";
        return false;
    }

    string firstHeader, segmentHeader;
    if (!getline(input, firstHeader) || !getline(input, segmentHeader)) {
        error = "incomplete header";
        return false;
    }

    smatch match;
    const regex headerPattern(R"(^Segments,([0-9]+),SegmentSize,([0-9]+)\r?$)");
    if (!regex_match(segmentHeader, match, headerPattern)) {
        error = "invalid segment header";
        return false;
    }
    const long long segments = stoll(match[1].str());
    const long long oldSegmentSize = stoll(match[2].str());
    if (segments <= 0 || oldSegmentSize < 2) {
        error = "invalid segment dimensions";
        return false;
    }

    vector<string> metadata;
    string line;
    bool foundDataHeader = false;
    while (getline(input, line)) {
        metadata.push_back(line);
        if (line.rfind("Time,Ampl", 0) == 0) {
            foundDataHeader = true;
            break;
        }
    }
    if (!foundDataHeader) {
        error = "missing Time,Ampl header";
        return false;
    }

    string sample0, sample1;
    if (!getline(input, sample0) || !getline(input, sample1)) {
        error = "missing initial samples";
        return false;
    }
    double time0 = 0.0, time1 = 0.0;
    if (!ParseTime(sample0, time0) || !ParseTime(sample1, time1) || !(time1 > time0)) {
        error = "invalid initial sample times";
        return false;
    }
    const double inputPitch = time1 - time0;
    const long long stride = max(1LL, llround(targetPitchSeconds / inputPitch));
    const long long newSegmentSize = (oldSegmentSize - 1) / stride + 1;
    if (stride == 1) return true;

    const fs::path temporary = path.string() + ".downsample.tmp";
    ofstream output(temporary, ios::trunc);
    if (!output) {
        error = "cannot create temporary output";
        return false;
    }
    output << firstHeader << '\n';
    output << "Segments," << segments << ",SegmentSize," << newSegmentSize << '\n';
    for (const string& metadataLine : metadata) output << metadataLine << '\n';

    const long long totalSamples = segments * oldSegmentSize;
    for (long long sample = 0; sample < totalSamples; ++sample) {
        if (sample == 0) {
            line = sample0;
        } else if (sample == 1) {
            line = sample1;
        } else if (!getline(input, line)) {
            output.close();
            fs::remove(temporary);
            error = "file ended before declared sample count";
            return false;
        }
        const long long localSample = sample % oldSegmentSize;
        if (localSample % stride == 0) output << line << '\n';
    }
    if (!output) {
        output.close();
        fs::remove(temporary);
        error = "failed while writing temporary output";
        return false;
    }
    output.close();
    fs::rename(temporary, path);
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2 || argc > 4) {
        cerr << "Usage: " << argv[0] << " CSV_DIRECTORY [TARGET_PITCH_NS=1] [THREADS=3]\n";
        return 2;
    }
    const fs::path directory = argv[1];
    const double targetPitchNs = argc >= 3 ? stod(argv[2]) : 1.0;
    const int threadCount = argc >= 4 ? stoi(argv[3]) : 3;
    if (!fs::is_directory(directory) || !(targetPitchNs > 0.0) || threadCount <= 0) return 2;

    const regex filePattern(R"(^C[234]Trace[0-9]+\.csv$)");
    vector<fs::path> files;
    for (const auto& entry : fs::directory_iterator(directory))
        if (entry.is_regular_file() && regex_match(entry.path().filename().string(), filePattern))
            files.push_back(entry.path());
    sort(files.begin(), files.end());

    atomic<size_t> next{0}, completed{0}, failed{0};
    mutex outputMutex;
    vector<thread> workers;
    for (int worker = 0; worker < threadCount; ++worker) {
        workers.emplace_back([&] {
            while (true) {
                const size_t index = next.fetch_add(1);
                if (index >= files.size()) break;
                string error;
                if (!Convert(files[index], targetPitchNs * 1.0e-9, error)) {
                    ++failed;
                    lock_guard<mutex> lock(outputMutex);
                    cerr << "\nFailed: " << files[index] << ": " << error << '\n';
                }
                const size_t done = ++completed;
                if (done % 25 == 0 || done == files.size()) {
                    lock_guard<mutex> lock(outputMutex);
                    cout << '\r' << done << '/' << files.size() << flush;
                }
            }
        });
    }
    for (auto& worker : workers) worker.join();
    cout << "\nConverted " << completed - failed << " files; failures=" << failed << '\n';
    return failed == 0 ? 0 : 1;
}

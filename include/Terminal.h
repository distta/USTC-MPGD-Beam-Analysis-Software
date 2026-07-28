#pragma once

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include <unistd.h>

namespace Terminal {

inline bool& VerboseFlag() {
    static bool enabled = false;
    return enabled;
}

inline bool& ColorOverrideFlag() {
    static bool enabled = true;
    return enabled;
}

inline bool Verbose() { return VerboseFlag(); }

inline bool Interactive() {
    static const bool interactive = ::isatty(STDOUT_FILENO) != 0;
    return interactive;
}

inline bool ColorEnabled() {
    const char* noColor = std::getenv("NO_COLOR");
    const char* term = std::getenv("TERM");
    return ColorOverrideFlag() && Interactive() && noColor == nullptr &&
           (!term || std::string(term) != "dumb");
}

inline std::string Paint(const std::string& text, const char* code) {
    if (!ColorEnabled()) return text;
    return "\033[" + std::string(code) + "m" + text + "\033[0m";
}

inline std::string Bold(const std::string& text) {
    return Paint(text, "1");
}

inline std::string Muted(const std::string& text) {
    return Paint(text, "2");
}

inline std::string Accent(const std::string& text) {
    return Paint(text, "1;36");
}

inline std::string Success(const std::string& text) {
    return Paint(text, "1;32");
}

inline std::string Warning(const std::string& text) {
    return Paint(text, "1;33");
}

inline std::string Failure(const std::string& text) {
    return Paint(text, "1;31");
}

template <typename Integer>
inline std::string Count(Integer value) {
    std::string text = std::to_string(value);
    for (int position = static_cast<int>(text.size()) - 3; position > 0;
         position -= 3) {
        text.insert(static_cast<size_t>(position), ",");
    }
    return text;
}

inline void StageStart(size_t index, size_t total, const std::string& name) {
    std::ostringstream number;
    number << '[' << index << '/' << total << ']';
    std::cout << Accent(number.str()) << ' ' << Bold(name) << '\n';
}

inline void StageDone(double seconds) {
    std::ostringstream elapsed;
    elapsed.setf(std::ios::fixed);
    elapsed.precision(1);
    elapsed << seconds << "s";
    std::cout << "      " << Success("✓") << ' '
              << Muted("completed · " + elapsed.str()) << "\n\n";
}

inline void Detail(const std::string& text) {
    std::cout << "      " << text << '\n';
}

inline void Note(const std::string& text) {
    std::cout << "      " << Warning("•") << ' ' << Muted(text) << '\n';
}

inline void ClearProgress() {
    if (Interactive()) std::cout << "\r\033[2K" << std::flush;
}

}  // namespace Terminal

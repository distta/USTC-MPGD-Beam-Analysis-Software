#include "Detector/DetectorFactory.h"
#include "Event/EventDisplayManager.h"
#include "Input/ConverterFactory.h"
#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptFactory.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include <nlohmann/json.hpp>

using json = nlohmann::json;
namespace fs = std::filesystem;

namespace {

constexpr int kUsageError = 2;
constexpr int kConfigurationError = 3;
constexpr int kInputError = 4;
constexpr int kInitializationError = 5;
constexpr int kExecutionError = 6;

void PrintUsage() {
   std::cout << "Usage:\n"
             << "  ./BeamAnalysis <base_dir> <run_id>\n\n"
             << "Example:\n"
             << "  ./BeamAnalysis /data/beam 1591\n";
}

json LoadConfig(const fs::path& path) {
   std::ifstream input(path);
   if (!input) {
      throw std::runtime_error("cannot open configuration file: " + path.string());
   }

   try {
      return json::parse(input);
   } catch (const json::parse_error& error) {
      throw std::runtime_error("invalid JSON in " + path.string() + ": " + error.what());
   }
}

void AddError(std::vector<std::string>& errors, const std::string& field,
              const std::string& message) {
   errors.push_back(field + ": " + message);
}

bool RequiredString(const json& object, const char* key,
                    const std::string& location,
                    std::vector<std::string>& errors) {
   if (!object.contains(key) || !object[key].is_string() ||
       object[key].get<std::string>().empty()) {
      AddError(errors, location + "." + key, "required string is missing");
      return false;
   }
   return true;
}

std::vector<std::string> ValidateConfig(const json& config) {
   std::vector<std::string> errors;
   const std::string mode =
       config.contains("mode") && config["mode"].is_string()
           ? config["mode"].get<std::string>()
           : "";

   if (!config.contains("mode") || !config["mode"].is_string()) {
      AddError(errors, "mode", "must be a string");
   } else {
      if (mode != "analysis" && mode != "eventDisplay") {
         AddError(errors, "mode",
                  "unsupported value '" + mode +
                      "' (expected 'analysis' or 'eventDisplay')");
      }
   }

   if (!config.contains("detectors") || !config["detectors"].is_array()) {
      AddError(errors, "detectors", "required array is missing");
   } else {
      for (size_t index = 0; index < config["detectors"].size(); ++index) {
         const auto& detector = config["detectors"][index];
         if (!detector.is_object()) {
            AddError(errors,
                     "detectors[" + std::to_string(index) + "]",
                     "must be an object");
         }
      }
   }

   if (config.contains("conversion")) {
      const auto& conversion = config["conversion"];
      if (!conversion.is_object()) {
         AddError(errors, "conversion", "must be an object");
      } else {
         if (conversion.contains("overwrite") &&
             !conversion["overwrite"].is_boolean()) {
            AddError(errors, "conversion.overwrite", "must be a boolean");
         }
         if (!conversion.contains("type") || !conversion["type"].is_string() ||
             conversion["type"].get<std::string>().empty()) {
            AddError(errors, "conversion.type",
                     "required string is missing");
         }
      }
   }

   size_t enabledCount = 0;
   if (!config.contains("scripts") || !config["scripts"].is_array()) {
      AddError(errors, "scripts", "required array is missing");
   } else if (mode == "analysis") {
      for (size_t index = 0; index < config["scripts"].size(); ++index) {
         const auto& script = config["scripts"][index];
         const std::string location =
             "scripts[" + std::to_string(index) + "]";
         if (!script.is_object()) {
            AddError(errors, location, "must be an object");
            continue;
         }
         if (script.contains("enabled") && !script["enabled"].is_boolean()) {
            AddError(errors, location + ".enabled", "must be a boolean");
            continue;
         }
         const bool enabled = script.value("enabled", true);
         if (!enabled || mode != "analysis") continue;

         RequiredString(script, "name", location, errors);
         const bool hasType = RequiredString(script, "type", location, errors);
         if (script.contains("config") && !script["config"].is_object()) {
            AddError(errors, location + ".config", "must be an object");
         }
         if (hasType) {
            const std::string type = script["type"].get<std::string>();
            if (!ScriptFactory::Instance().IsRegistered(type)) {
               AddError(errors, location + ".type",
                        "unknown script type '" + type + "'");
            }
         }
         if (enabled) ++enabledCount;
      }
   }

   if (mode == "analysis" && enabledCount == 0) {
      AddError(errors, "scripts",
               "analysis mode requires at least one enabled script");
   }
   return errors;
}

std::vector<const json*> EnabledScripts(const json& config) {
   std::vector<const json*> enabled;
   for (const auto& script : config["scripts"]) {
      if (!script.is_object() || !script.contains("name") ||
          !script["name"].is_string()) {
         continue;
      }
      if (!script.contains("enabled") ||
          (script["enabled"].is_boolean() && script["enabled"].get<bool>())) {
         enabled.push_back(&script);
      }
   }
   return enabled;
}

std::string DisplayPath(const fs::path& path) {
   std::error_code code;
   const fs::path absolute = fs::absolute(path, code).lexically_normal();
   if (code) return path.generic_string();

   const fs::path current = fs::current_path(code);
   if (code) return absolute.generic_string();

   const fs::path relative = fs::relative(absolute, current, code);
   if (!code && !relative.empty()) return relative.generic_string();
   return absolute.generic_string();
}

std::string FormatCount(Long64_t value) {
   std::string text = std::to_string(value);
   for (int position = static_cast<int>(text.size()) - 3; position > 0;
        position -= 3) {
      text.insert(static_cast<size_t>(position), ",");
   }
   return text;
}

std::string DetectorSummary(size_t trackers, size_t dut, size_t ignored) {
   return std::to_string(trackers) + " Tracker, " + std::to_string(dut) +
          " DUT, " + std::to_string(ignored) + " Ignored";
}

void PrintRunPlan(const std::string& runID, const json& config,
                  const fs::path& configFile, const fs::path& baseDir,
                  const fs::path& outputDirectory,
                  const std::string& inputSummary) {
   size_t trackers = 0;
   size_t dut = 0;
   size_t ignored = 0;
   for (const auto& detector : config["detectors"]) {
      const std::string role = detector.value("role", "Tracker");
      if (role == "Tracker")
         ++trackers;
      else if (role == "DUT")
         ++dut;
      else
         ++ignored;
   }
   std::cout << "BeamAnalysis run " << runID << '\n'
             << "Mode     : " << config.value("mode", "analysis") << '\n'
             << "Base     : " << DisplayPath(baseDir) << '\n'
             << "Config   : " << DisplayPath(configFile) << '\n'
             << "Input    : " << inputSummary << '\n'
             << "Output   : " << DisplayPath(outputDirectory) << '\n'
             << "Geometry : " << DetectorSummary(trackers, dut, ignored)
             << '\n'
             << "Scripts  : ";
   const auto scripts = EnabledScripts(config);
   if (scripts.empty()) {
      std::cout << "(none)\n";
   } else {
      for (size_t index = 0; index < scripts.size(); ++index) {
         if (index > 0) std::cout << ", ";
         std::cout << (*scripts[index])["name"].get<std::string>();
      }
      std::cout << '\n';
   }
   std::cout << '\n'
             << std::flush;
}

bool EnsureDirectory(const fs::path& directory, std::string& error) {
   std::error_code code;
   fs::create_directories(directory, code);
   if (code || !fs::is_directory(directory)) {
      error = "cannot create directory " + directory.string();
      if (code) error += ": " + code.message();
      return false;
   }
   return true;
}

bool RunAnalysis(const json& config,
                 const std::shared_ptr<RawDataParser>& parser,
                 const fs::path& outputDirectory,
                 std::string& failure) {
   const auto enabled = EnabledScripts(config);
   for (size_t index = 0; index < enabled.size(); ++index) {
      const auto& scriptConfig = *enabled[index];
      const std::string name = scriptConfig["name"].get<std::string>();
      const std::string type = scriptConfig["type"].get<std::string>();
      const auto started = std::chrono::steady_clock::now();
      std::cout << "------------------------------------------------------------\n"
                << "Script " << index + 1 << '/' << enabled.size() << ": "
                << name << '\n'
                << "------------------------------------------------------------\n";

      std::shared_ptr<IScript> script;
      bool finalizeAttempted = false;
      bool success = false;
      std::string scriptError;
      const auto finalize = [&] {
         if (!script || finalizeAttempted) return;
         finalizeAttempted = true;
         script->Finalize();
      };
      try {
         script = ScriptFactory::Instance().CreateScript(
             type, scriptConfig.value("config", json::object()));
         if (!script) {
            throw std::runtime_error("factory returned a null script");
         }
         script->SetParser(parser);
         script->SetConfig(config);
         script->SetOutputDir(outputDirectory.string() + "/");

         if (!script->Validate()) {
            throw std::runtime_error("configuration validation failed");
         }
         if (!script->Initialize()) {
            throw std::runtime_error("initialization failed");
         }
         success = script->Execute();
         finalize();
         if (!success) scriptError = "Execute() returned false";
      } catch (const std::exception& error) {
         success = false;
         scriptError = error.what();
         try {
            finalize();
         } catch (const std::exception& finalizeError) {
            scriptError += "; Finalize() failed: ";
            scriptError += finalizeError.what();
         } catch (...) {
            scriptError += "; Finalize() failed with an unknown exception";
         }
      } catch (...) {
         success = false;
         scriptError = "unknown exception";
         try {
            finalize();
         } catch (const std::exception& finalizeError) {
            scriptError += "; Finalize() failed: ";
            scriptError += finalizeError.what();
         } catch (...) {
            scriptError += "; Finalize() failed with an unknown exception";
         }
      }

      const double elapsedSeconds =
          std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                        started)
              .count();
      if (!success) {
         failure = name + ": " + scriptError;
         std::cerr << "\nScript " << index + 1 << '/' << enabled.size()
                   << " failed: " << name << " after " << std::fixed
                   << std::setprecision(1) << elapsedSeconds << " s: "
                   << scriptError << "\n"
                   << "------------------------------------------------------------\n";
         return false;
      }
      std::cout << "\nScript " << index + 1 << '/' << enabled.size()
                << " complete: " << name << " in " << std::fixed
                << std::setprecision(1) << elapsedSeconds << " s\n"
                << "------------------------------------------------------------\n\n";
   }
   return true;
}

bool RunEventDisplay(const std::shared_ptr<RawDataParser>& parser,
                     const fs::path& outputDirectory,
                     std::string& error) {
   EventDisplayManager manager(parser, outputDirectory.string());
   if (!manager.Initialize()) {
      error = "EventDisplayManager initialization failed";
      return false;
   }
   manager.RunInteractive();
   return true;
}

}  // namespace

int main(int argc, char* argv[]) {
   if (argc == 2 && std::string(argv[1]) == "--help") {
      PrintUsage();
      return 0;
   }
   if (argc != 3) {
      PrintUsage();
      return kUsageError;
   }

   const std::string runID = argv[2];
   if (runID.empty() ||
       !std::all_of(runID.begin(), runID.end(),
                    [](unsigned char c) { return std::isdigit(c); })) {
      std::cerr << "Configuration error:\n"
                << "  run_id: must be a non-empty decimal number\n";
      return kConfigurationError;
   }

   const auto runStarted = std::chrono::steady_clock::now();

   const fs::path baseDir = fs::absolute(argv[1]).lexically_normal();
   const fs::path rawDirectory = baseDir / "raw";
   const fs::path rootConfigFile = baseDir / "config.json";
   const fs::path rawConfigFile = rawDirectory / "config.json";
   const fs::path configFile =
       fs::is_regular_file(rootConfigFile) ? rootConfigFile : rawConfigFile;
   const fs::path processedDirectory = baseDir / "processed";
   const fs::path outputDirectory = baseDir / "result" / runID;

   if (!fs::is_directory(baseDir)) {
      std::cerr << "Input/output error: base directory does not exist: "
                << baseDir << '\n';
      return kInputError;
   }

   std::string error;
   if (!EnsureDirectory(rawDirectory, error) ||
       !EnsureDirectory(processedDirectory, error) ||
       !EnsureDirectory(outputDirectory, error)) {
      std::cerr << "Input/output error: " << error << '\n';
      return kInputError;
   }

   json config;
   try {
      config = LoadConfig(configFile);
   } catch (const std::exception& exception) {
      std::cerr << "Configuration error: " << exception.what() << '\n';
      return kConfigurationError;
   }
   if (!config.is_object()) {
      std::cerr << "Configuration error: root must be an object\n";
      return kConfigurationError;
   }
   if (!config.contains("mode")) config["mode"] = "analysis";

   const auto configurationErrors = ValidateConfig(config);
   if (!configurationErrors.empty()) {
      std::cerr << "Configuration error:\n";
      for (const auto& error : configurationErrors) {
         std::cerr << "  " << error << '\n';
      }
      return kConfigurationError;
   }

   const fs::path processedInput = processedDirectory / ("run" + runID + ".root");
   std::string inputSummary;
   std::unique_ptr<IRawDataConverter> converter;
   bool shouldConvert = false;
   const bool conversionEnabled =
       config.contains("conversion") &&
       config["conversion"].value("enabled", true);
   if (conversionEnabled) {
      const std::string converterType =
          config["conversion"].value("type", "APV25SRS");
      converter = ConverterFactory::Create(converterType);
      if (!converter) {
         std::cerr << "Input/conversion error: unknown converter type '"
                   << converterType << "'\n";
         return kInputError;
      }
      converter->Configure(config["conversion"]);
      const bool overwrite = config["conversion"].value("overwrite", false);
      shouldConvert = overwrite || !fs::is_regular_file(processedInput);
      inputSummary = shouldConvert
                         ? DisplayPath(rawDirectory) + " -> " +
                               DisplayPath(processedInput) + " (" +
                               converterType + " conversion)"
                         : DisplayPath(processedInput) + " (reused)";
   } else {
      if (!fs::is_regular_file(processedInput)) {
         std::cerr << "Input/conversion error: processed input does not exist: "
                   << processedInput << '\n';
         return kInputError;
      }
      inputSummary = DisplayPath(processedInput) + " (conversion disabled)";
   }
   PrintRunPlan(runID, config, configFile, baseDir, outputDirectory,
                inputSummary);

   if (shouldConvert) {
      const auto started = std::chrono::steady_clock::now();
      const std::string converterType =
          config["conversion"].value("type", "APV25SRS");
      if (!converter->AcquireRawData(rawDirectory, runID, error)) {
         std::cerr << "Input/conversion error: " << error << '\n';
         return kInputError;
      }
      std::cout << "\n============================================================\n"
                << "  Conversion started: " << converterType << '\n'
                << "============================================================\n";
      if (!converter->Convert(processedInput.string())) {
         std::cerr << "Input/conversion error: conversion failed\n";
         return kInputError;
      }
      if (!fs::is_regular_file(processedInput)) {
         std::cerr << "Input/conversion error: converter reported success but "
                      "output is missing: "
                   << processedInput << '\n';
         return kInputError;
      }
      std::cout << "============================================================\n"
                << "  Conversion completed in " << std::fixed
                << std::setprecision(1)
                << std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - started)
                       .count()
                << " s\n"
                << "============================================================\n";
   }

   auto& detectorFactory = DetectorFactory::GetInstance();
   detectorFactory.Clear();
   if (!detectorFactory.Initialize(config)) {
      return kInitializationError;
   }

   auto parser = std::make_shared<RawDataParser>(processedInput.string());
   if (!parser->Initialize()) {
      return kInitializationError;
   }

   bool success = false;
   try {
      if (config["mode"] == "analysis") {
         success = RunAnalysis(config, parser, outputDirectory, error);
      } else {
         success = RunEventDisplay(parser, outputDirectory, error);
      }
   } catch (const std::exception& exception) {
      error = exception.what();
      std::cerr << "Execution error: " << error << '\n';
      success = false;
   } catch (...) {
      error = "unknown exception";
      std::cerr << "Execution error: " << error << '\n';
      success = false;
   }

   const double elapsed =
       std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                     runStarted)
           .count();
   std::cout << "\nSummary\n"
             << "Status  : " << (success ? "OK" : "FAILED") << '\n'
             << "Run ID  : " << runID << '\n'
             << "Output  : " << DisplayPath(outputDirectory) << '\n'
             << "Elapsed : " << std::fixed << std::setprecision(1)
             << elapsed << " s\n";
   return success ? 0 : kExecutionError;
}

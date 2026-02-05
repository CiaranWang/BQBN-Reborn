#include <iostream>
#include <ctime>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <optional>
#include <cstdlib>

#ifdef _WIN32
#include <windows.h>
#endif

#include "BQBN-Reborn.h"
#include "helptext.h"

namespace fs = std::filesystem;
using namespace std;

// ======================================================
// Version
// ======================================================
static const std::string PROGRAM_VERSION = "1.0.0";

static void print_version()
{
    std::cout << "BQBN-Reborn Version: " << PROGRAM_VERSION << std::endl;
}

static void print_help()
{
    std::cout << HELP_TEXT << std::endl;
}

// ======================================================
// Utility: executable location
// ======================================================
fs::path get_exe_path(char* argv0)
{
    return fs::absolute(argv0);
}

// ======================================================
// Windows binary updater
// ======================================================
static void run_update(char* argv0)
{
#ifdef _WIN32
    fs::path exe_path = get_exe_path(argv0);
    fs::path exe_dir = exe_path.parent_path();

    fs::path new_exe = exe_dir / "BQBN_new.exe";
    fs::path updater = exe_dir / "BQBN_updater.exe";

    std::cout << "[UPDATE] Downloading new version...\n";

        std::string cmd =
            "curl -L -o \"" + new_exe.string() +
            "\" https://github.com/CiaranWang/BQBN-Reborn/releases/latest/download/BQBN-Reborn.exe";

    if (system(cmd.c_str()) != 0) {
        std::cout << "[UPDATE] Download failed.\n";
        return;
    }

    std::cout << "[UPDATE] Launching updater...\n";

    std::string args =
        "\"" + exe_path.string() + "\" \"" + new_exe.string() + "\"";

    ShellExecuteA(
        nullptr,
        "open",
        updater.string().c_str(),
        args.c_str(),
        exe_dir.string().c_str(),
        SW_SHOWNORMAL
    );

    exit(0); // Must exit so updater can replace the exe
#endif
}

int main(int argc, char* argv[])
{
    // ---------- Early check for updates ----------
    if (argc > 1) {
        std::string arg1 = argv[1];
        if (arg1 == "--version" || arg1 == "-v" || arg1 == "-V") {
            print_version();
            return 0;
        }
        if (arg1 == "--update" || arg1 == "-u" || arg1 == "-U") {
            run_update(argv[0]);
            return 0;
        }
        if (arg1 == "--help" || arg1 == "-h" || arg1 == "-H") {
            print_help();
            return 0;
        }
    }

    optional<int> seed; // <<< optional seed

    // ----- Parse command line arguments -----
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) {
            input_file = argv[++i];
        }
        else if (arg == "-o" && i + 1 < argc) {
            output_folder = argv[++i];
        }
        else if (arg == "-n" && i + 1 < argc) {
            project = argv[++i];
        }
        else if (arg == "--seed" && i + 1 < argc) {
            seed = std::stoi(argv[++i]);  // <<< store in optional
        }   
        else if (arg == "--step" && i + 1 < argc) {
            cycles = std::stoi(argv[++i]);
        }
        else if (arg == "--checkpoint" && i + 1 < argc) {
            checkpoint = std::stoi(argv[++i]);
        }
        else if (arg == "--tunnel_radias" && i + 1 < argc) {
            tunnel_r = std::stoi(argv[++i]);
        }
        else if (arg == "--mirror") {
            if (i + 1 < argc) {
                string next = argv[i + 1];
                if (next == "1" || next == "true" || next == "TRUE") {
                    mirror = true;
                    ++i; // skip the value
                }
                else if (next == "0" || next == "false" || next == "FALSE") {
                    mirror = false;
                    ++i; // skip the value
                }
                else {
                    // next argument is not a valid value → treat --mirror as flag
                    mirror = true;
                }
            }
            else {
                // no value → treat --mirror as flag
                mirror = true;
            }
        }
        else {
            cerr << "Unknown or incomplete argument: " << arg << "\n";
            return 1;
        }
    }

    // Ensure the input file exists
    if (input_file.empty()) {
        cerr << "Error: -i input_file is required\n";
        return 1;
    }

    // ----- Set default output path if not provided -----
    if (output_folder.empty()) {
        output_folder = input_file.parent_path();
    }
    // Ensure the output folder exists
    if (!fs::exists(output_folder)) {
        try {
            fs::create_directories(output_folder); // creates all missing intermediate directories
            std::cout << "[INFO] Created output folder: " << output_folder << std::endl;
        } catch (const fs::filesystem_error& e) {
            std::cerr << "Error: Failed to create output folder '"
                << output_folder << "': " << e.what() << std::endl;
            return 1;
        }
    }

    //print running info on screen
    cout << "Input file: \t" << input_file.string() << "\n";
    cout << "Output folder: \t" << output_folder.string() << "\n";
    cout << "Project name: \t" << project << "\n";
    cout << "Cycle number: \t" << cycles << "\n";
    cout << "Check every: \t" << checkpoint << " steps" << "\n";

    if (seed.has_value()) {
        cout << "Seed:        " << seed.value() << " (using deterministic per-thread RNG)\n";
    } else {
        cout << "Seed:        (not provided, using time-based RNG per thread)\n";
    }

    if (tunnel_r > 0)
        std::cout << "Tunnel enabled, radius = " << tunnel_r << std::endl;
    else
        std::cout << "No tunnel" << std::endl;

    std::cout << "Mirror: " << (mirror ? "ON" : "OFF") << std::endl;

    unsigned int rng_seed = seed.value_or(static_cast<unsigned int>(time(nullptr)));
    std::srand(rng_seed);  // for std::rand() usage

    if (!load_parameters(input_file.string()))
    {
        cerr << "Error: Failed to load parameters from" << input_file.string() << "\n";
    }
    print_parameters();

    GlobalConst();

    OutInfo(seed);

    LBM();

    return 0;
}

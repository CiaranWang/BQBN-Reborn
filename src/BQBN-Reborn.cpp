#include <iostream>
#include <ctime>
#include <windows.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <optional>
#include <cstdlib>

#include "BQBN-Reborn.h"
#include "helptext.h"

namespace fs = std::filesystem;
using namespace std;

// ======================================================
// Version
// ======================================================
static const std::string PROGRAM_VERSION = "0.1.1";

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
        "\" https://github.com/CiaranWang/BQBN-Reborn/releases/download/v0.1.1/BQBN-Reborn.exe";

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
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        if (arg == "--version" || arg == "-v" || arg == "-V")
        {
            print_version();
            return 0;
        }
        else if (arg == "--help" || arg == "-h" || arg == "-H")
        {
            print_help();
            return 0;
        }
        else if (arg == "--update" || arg == "-u" || arg == "-U")
        {
            run_update(argv[0]);
            return 0;
        }
    }

	cout << "Hello World." << endl;
	return 0;
}

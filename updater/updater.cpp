#include <windows.h>
#include <filesystem>

namespace fs = std::filesystem;

int main(int argc, char* argv[])
{
    if (argc < 3)
        return 1;

    fs::path oldExe = argv[1];
    fs::path newExe = argv[2];

    // Give main app time to exit
    Sleep(1000);

    // Replace executable
    if (!MoveFileExA(
        newExe.string().c_str(),
        oldExe.string().c_str(),
        MOVEFILE_REPLACE_EXISTING))
    {
        return 2;
    }

    // Restart application
    ShellExecuteA(
        nullptr,
        "open",
        oldExe.string().c_str(),
        nullptr,
        nullptr,
        SW_SHOWNORMAL
    );

    return 0;
}

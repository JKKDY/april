#include <april/april.h>
#include <iostream>
#include <ankerl/unordered_dense.h>

namespace april {
    void hello() {
        std::cout << "Hello from APRIL!" << std::endl;

        ankerl::unordered_dense::map<std::string, int> wordCount;

        wordCount["april"] = 42;
        wordCount["acropolis"] = 1337;

        for (const auto& [word, count] : wordCount) {
            std::cout << word << ": " << count << '\n';
        }

        #if defined(_WIN32)
                std::cout << "Running on: Windows\n";
        #elif defined(__linux__)
                std::cout << "Running on: Linux\n";
        #elif defined(__APPLE__)
                std::cout << "Running on: macOS\n";
        #else
                std::cout << "Running on: Unknown platform\n";
        #endif
    }
}
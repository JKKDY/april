#include "april/env/environment.h"
#include "april/env/particle.h"
#include "april/env/interaction.h"


#pragma once

// Define APRIL_API for DLL import/export
// #ifdef _WIN32
//   #ifdef APRIL_BUILD_DLL
//     #define APRIL_API __declspec(dllexport)
//   #else
//     #define APRIL_API __declspec(dllimport)
//   #endif
// #else
//   #define APRIL_API
// #endif

#define APRIL_API


namespace april {

    APRIL_API void hello();
}
#pragma once

#ifndef NDEBUG

#include <iostream>
#include <stdexcept>

#define AP_ASSERT(Expr, Msg) april::utils::debug::assert(#Expr, (Expr), __FILE__, __LINE__, (Msg))

namespace april::utils::debug {

        [[noreturn]] inline void assert(const char* expr_str, bool expr, const char* file, int line, const char* msg) {
            if (!expr) {
                std::cerr << "Assert failed:\t" << msg << "\n"
                    << "Expected:\t" << expr_str << "\n"
                    << "Source:\t\t" << file << ", line " << line << "\n";
                throw std::runtime_error(msg);
            }
        }

} // namespace april::utils::debug

#else

#define AP_ASSERT(Expr, Msg)

#endif

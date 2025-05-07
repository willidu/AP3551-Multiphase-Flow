#pragma once

#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_matrix.h>
#include <sunlinsol/sunlinsol_dense.h>

using real_t = sunrealtype;

#include <memory>
#include <limits>
#include <iomanip>
#include <sstream>
#include <functional>

#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>


#ifdef __GNUC__  // Suppress warnings from matplotlibcpp
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#ifndef WITHOUT_NUMPY
    #define WITHOUT_NUMPY
#endif

#include <matplotlibcpp.h>
namespace plt = matplotlibcpp;

#ifdef __GNUC__ // Restore warnings
    #pragma GCC diagnostic pop
#endif


namespace CMF
{

class Log
{
public:
    static void Init();
    inline static void SetLevel(spdlog::level::level_enum level) { s_Logger->set_level(level); }
    inline static std::shared_ptr<spdlog::logger>& getLogger() { return s_Logger; }

private:
    static std::shared_ptr<spdlog::logger> s_Logger;
};


// Solver for the linear system Au = b
// Currently the solver does not support restarts
class LinearSolver
{
public:
    LinearSolver(size_t N);
    ~LinearSolver();

    [[nodiscard]] inline SUNMatrix Matrix() noexcept { return m_Matrix;   }
    [[nodiscard]] inline N_Vector  RHS()    noexcept { return m_RHS;      }
    [[nodiscard]] inline N_Vector  SolVec() noexcept { return m_Solution; }
    [[nodiscard]] std::vector<real_t> getSolutionVector() const;

    void solve(real_t tol = 1e-10);

private:
    const size_t m_Vars;
    SUNContext m_Context;
    SUNMatrix m_Matrix;
    SUNLinearSolver m_Solver;
    N_Vector m_RHS;
    N_Vector m_Solution;
};

std::vector<real_t> channelMesh1D(size_t N, real_t H, real_t s = 1.0);

std::vector<real_t> logspace(real_t start, real_t end, size_t num_points);

/**
 * @brief Van Driest damping function.
 * 
 * @param yplus   [--] y+ at the wall
 * @param aplus   [--] A+ parameter. Default is 26.0
 * @return real_t [--] Van Driest damping function value
 */
real_t vanDriest(real_t yplus, real_t aplus = 26.0);


/**
 * @brief Calculate U+ given some y+ and a roughness height ks+.
 * 
 * @param yplus   [--] y+ at the wall
 * @param ksplus  [--] ks+ roughness height. Default is 0.0
 * @return real_t [--] U+ at the wall
 */
real_t uplus(real_t yplus, real_t ksplus = 0.0);

} // namespace CMF

#define LOG_TRACE(...)    ::CMF::Log::getLogger()->trace(__VA_ARGS__)
#define LOG_INFO(...)     ::CMF::Log::getLogger()->info(__VA_ARGS__)
#define LOG_WARN(...)     ::CMF::Log::getLogger()->warn(__VA_ARGS__)
#define LOG_ERROR(...)    ::CMF::Log::getLogger()->error(__VA_ARGS__)
#define LOG_CRITICAL(...) ::CMF::Log::getLogger()->critical(__VA_ARGS__)

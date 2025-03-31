#include "core.hpp"

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>


namespace CMF
{

std::shared_ptr<spdlog::logger> Log::s_Logger;

void Log::Init()
{
    spdlog::set_pattern("%^[%T] [%l] %n: %v%$");
    s_Logger = spdlog::stdout_color_mt("APP");
    s_Logger->set_level(spdlog::level::trace);
}


std::vector<real_t> channelMesh1D(size_t N, real_t H, real_t s)
{
    std::vector<real_t> mesh(N);

    for (size_t i = 0; i < N; ++i)
    {
        // Transform xi = [-1, 1] to y = (0, H) (not including endpoints)
        const real_t xi = -1.0 + 2.0 * (i + 1) / (N + 2);
        mesh[i] = (H / 2) * (1.0 + std::tanh(s * xi) / std::tanh(s));
    }

    return mesh;
}


LinearSolver::LinearSolver(size_t N)
  : m_Vars(N)
{
    if (SUNErrCode flag = SUNContext_Create(SUN_COMM_NULL, &m_Context); flag != 0)
    {
        LOG_ERROR("SUNContext_Create returned error code {0}", flag);
        throw std::runtime_error("Failed to create SUNDIALS context");
    }

    m_Matrix   = SUNDenseMatrix(N, N, m_Context);
    m_Solution = N_VNew_Serial(N, m_Context);
    m_RHS      = N_VNew_Serial(N, m_Context);
    m_Solver   = SUNLinSol_Dense(m_Solution, m_Matrix, m_Context);
}


LinearSolver::~LinearSolver()
{
    N_VDestroy(m_RHS);
    N_VDestroy(m_Solution);
    SUNLinSolFree_Dense(m_Solver);
    SUNMatDestroy(m_Matrix);
    SUNContext_Free(&m_Context);
}


std::vector<real_t> LinearSolver::getSolutionVector() const
{
    const real_t* u_data = N_VGetArrayPointer_Serial(m_Solution);
    return std::vector<real_t>(u_data, u_data + m_Vars);
}


void LinearSolver::solve(real_t tol)
{
    if (SUNErrCode flag = SUNLinSolSetup(m_Solver, m_Matrix); flag != 0)
    {
        LOG_ERROR("SUNLinSolSetup returned error code {0}", flag);
        throw std::runtime_error("Failed to setup linear solver");
    }

    if (SUNErrCode flag = SUNLinSolSolve(m_Solver, m_Matrix, m_Solution, m_RHS, tol);
        flag != 0)
    {
        LOG_ERROR("SUNLinSolSolve returned error code {0}", flag);
        throw std::runtime_error("Failed to solve linear system");
    }
}


real_t vanDriest(real_t yplus, real_t aplus)
{
    return 1.0 - std::exp(-yplus / aplus);
}


real_t uplus(real_t yplus, real_t ksplus)
{
    // TODO - Add roughness
    static constexpr real_t kappa = 0.41;    // von Kármán constant
    static constexpr real_t B = 5.0;         // Empirical constant

    assert(yplus > 0.0 && "Invalid y+");

    if (yplus < 5.0)
    {
        // Viscous sublayer: Linear law
        return yplus;
    }
    else if (yplus < 30.0)
    {
        // Blended transition
        real_t t = (yplus - 5.0) / (30.0 - 5.0);
        real_t U_viscous = 5.0;                            // U+ at y+ = 5
        real_t U_log = (1.0 / kappa) * std::log(30.0) + B; // Log-law at y+ = 30
        return (1 - t) * U_viscous + t * U_log;
    }

    // Log-law region
    return (1.0 / kappa) * std::log(yplus) + B;
}


real_t utau(real_t y, real_t U, real_t nu)
{
    static constexpr real_t kappa = 0.41;    // von Kármán constant
    static constexpr real_t B = 5.0;         // Empirical constant - TODO add roughness
    static constexpr real_t tol = 1e-3;      // Convergence tolerance
    static constexpr size_t maxIter = 1000;  // Max iterations for Newton's method

    // Blended law of the wall residual function
    auto f = [&](real_t u_tau)
    {
        const real_t y_plus = y * u_tau / nu;
        // TODO - add roughness
        return uplus(y_plus, 0.0) - (U / u_tau);
    };

    // Initial guess: Blended approach
    real_t u_tau = std::sqrt(nu * U / y);

    // Newton's method to solve for u_tau
    for (size_t iter = 0; iter < maxIter; ++iter)
    {
        const real_t y_plus = y * u_tau / nu;
        const real_t df_du_tau = (y / nu) / ((1.0 / kappa) * (1.0 / y_plus)); // Approximate derivative
        const real_t u_tau_new = u_tau - f(u_tau) / df_du_tau;

        if (std::abs((u_tau_new - u_tau)/u_tau) < tol)
            break;

        u_tau = u_tau_new;

        if (iter == maxIter - 1)
        {
            LOG_WARN("u_tau calculation did not converge after {0} iterations.", maxIter);
            LOG_TRACE("Final y+ {0:.2e}, u_tau {1:.2e}", y_plus, u_tau);
        }
    }

    if (std::isnan(u_tau))
        LOG_WARN("u_tau calculation resulted in u_tau = NaN. "
                 "This warning may be ignored if the final y+ is valid.");

    return u_tau;
}

} // namespace CMF

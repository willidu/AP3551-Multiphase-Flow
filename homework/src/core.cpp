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
        const real_t xi = -1.0 + 2.0 * i / (N - 1);
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

} // namespace CMF

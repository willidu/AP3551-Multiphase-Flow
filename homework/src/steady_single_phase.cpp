#include "steady_single_phase.hpp"

#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_matrix.h>
#include <sunlinsol/sunlinsol_dense.h>

namespace CMF
{

namespace T1
{

void fillSystem(
    SUNMatrix A,
    N_Vector b,
    const std::vector<real_t>& mesh,
    const std::function<real_t(real_t)>& viscosity,
    const BC& bc)
{
    // TODO: Check signs of the forcing terms
    // TODO: Check indexes (i, i+1/2, i-1/2 etc.)

    const size_t N = mesh.size();

    for (size_t i = 1; i < N - 1; ++i)
    {
        // Forward and backward spacing
        const real_t dy_f = mesh.at(i+1) - mesh.at(i);
        const real_t dy_b = mesh.at(i) - mesh.at(i-1);
        const real_t dy = (dy_f + dy_b) / 2.0;

        const real_t mu_ip = viscosity((mesh.at(i) + mesh.at(i+1)) / 2);
        const real_t mu_im = viscosity((mesh.at(i) + mesh.at(i-1)) / 2);

        SM_ELEMENT_D(A, i, i - 1) = -mu_im / (dy_b * dy);
        SM_ELEMENT_D(A, i, i)     =  mu_ip / (dy_f * dy)
                                  +  mu_im / (dy_b * dy);
        SM_ELEMENT_D(A, i, i + 1) = -mu_ip / (dy_f * dy);

        NV_Ith_S(b, i) = -bc.global.second;  // Forcing term
    }

    switch (bc.lowerWall.first)
    {
        case WallBC::Velocity:
        {
            const real_t dy = mesh.at(1) - mesh.at(0);
            const real_t mu0 = viscosity(mesh.at(0) + dy / 2.0);
            SM_ELEMENT_D(A, 0, 0) = 3 * mu0 / (dy * dy);
            SM_ELEMENT_D(A, 0, 1) = -mu0 / (dy * dy);
            NV_Ith_S(b, 0) = bc.global.second 
                           + (2 * mu0 * bc.lowerWall.second / (dy * dy));
            break;
        }

        case WallBC::VelocityGradient:
        {
            const real_t dy = mesh.at(1) - mesh.at(0);
            const real_t mu0 = viscosity(mesh.at(0) + dy / 2.0);
            SM_ELEMENT_D(A, 0, 0) = 2 * mu0 / (dy * dy);
            SM_ELEMENT_D(A, 0, 1) = -2 * mu0 / (dy * dy);
            NV_Ith_S(b, 0) = bc.global.second
                           + bc.lowerWall.second / dy;
            break;
        }

        default:
        {
            const char* msg = "Unsupported lower wall BC";
            LOG_ERROR(msg);
            throw std::runtime_error(msg);
        }
    }

    switch (bc.upperWall.first)
    {
        case WallBC::Velocity:
        {
            const real_t dy = mesh.at(N-1) - mesh.at(N-2);
            const real_t muN = viscosity(mesh.at(N-1) - dy / 2.0);
            SM_ELEMENT_D(A, N - 1, N - 2) = -muN / (dy * dy);
            SM_ELEMENT_D(A, N - 1, N - 1) = 3 * muN / (dy * dy);
            NV_Ith_S(b, N - 1) = bc.global.second
                               + (2 * muN * bc.upperWall.second / (dy * dy));
            break;
        }

        case WallBC::VelocityGradient:
        {
            const real_t dy = mesh.at(N-1) - mesh.at(N-2);
            const real_t muN = viscosity(mesh.at(N-1) - dy / 2.0);
            SM_ELEMENT_D(A, N - 1, N - 2) = -2 * muN / (dy * dy);
            SM_ELEMENT_D(A, N - 1, N - 1) = 2 * muN / (dy * dy);
            NV_Ith_S(b, N - 1) = bc.global.second
                               + bc.upperWall.second / dy;
            break;
        }

        default:
        {
            const char* msg = "Unsupported upper wall BC";
            LOG_ERROR(msg);
            throw std::runtime_error(msg);
        }
    }
}


std::vector<real_t> steadyChannelFlow(
    const std::vector<real_t>& mesh,
    std::function<real_t(real_t)> viscosity,
    BC bc
)
{
    if (bc.global.first != GlobalBC::PressureGradient)
    {
        const char* msg = "Only PressureGradient BC is supported";
        LOG_ERROR(msg);
        throw std::runtime_error(msg);
    }
    const size_t N = mesh.size();

    // Step 1: Create the SUNContext
    SUNContext sunctx;
    if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
    {
        const char* msg = "Failed to create SUNContext";
        LOG_ERROR(msg);
        throw std::runtime_error(msg);
    }

    // Step 2: Create dense matrix and RHS vector
    SUNMatrix A = SUNDenseMatrix(N, N, sunctx);
    N_Vector b = N_VNew_Serial(N, sunctx);
    N_Vector u = N_VNew_Serial(N, sunctx);

    // Step 3: Fill system matrix and RHS vector
    fillSystem(A, b, mesh, viscosity, bc);

    // Step 3: Solve the system using SUNLINSOL_DENSE
    SUNLinearSolver LS = SUNLinSol_Dense(u, A, sunctx);
    SUNLinSolSetup(LS, A);
    SUNLinSolSolve(LS, A, u, b, 1e-10);

    // Step 4: Copy solution to std::vector
    const real_t* u_data = N_VGetArrayPointer_Serial(u);
    const std::vector<real_t> u_result(u_data, u_data + N);

    // Step 5: Cleanup
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(b);
    N_VDestroy(u);
    SUNContext_Free(&sunctx);

    return u_result;
}

} // namespace T1
} // namespace CMF

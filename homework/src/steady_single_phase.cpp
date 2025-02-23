#include "steady_single_phase.hpp"

#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_matrix.h>
#include <sunlinsol/sunlinsol_dense.h>

namespace CMF
{

namespace T1
{

std::vector<real_t> steadyChannelFlow(
    const std::vector<real_t>& mesh,
    std::function<real_t(real_t)> viscosity,
    BC bc
)
{
    // Temporary, need to implement all types of BCs
    const real_t dPdx = bc.global.second;
    const real_t U0 = bc.wall.second;
    const real_t UH = bc.wall.second;

    const size_t N = mesh.size();

    // Step 1: Create the SUNContext
    SUNContext sunctx;
    if (SUNContext_Create(NULL, &sunctx))
    {
        throw std::runtime_error("Failed to create SUNContext");
    }

    // Step 2: Create dense matrix and RHS vector
    SUNMatrix A = SUNDenseMatrix(N, N, sunctx);
    N_Vector b = N_VNew_Serial(N, sunctx);
    N_Vector u = N_VNew_Serial(N, sunctx);

    // Step 3: Fill system matrix
    real_t dy = mesh[1] - mesh[0]; // TODO: Implement non-uniform mesh
    for (int i = 1; i < N - 1; ++i)
    {
        real_t mu_ip = viscosity(mesh[i] + dy / 2);
        real_t mu_im = viscosity(mesh[i] - dy / 2);

        SM_ELEMENT_D(A, i, i - 1) = -mu_im / (dy * dy);
        SM_ELEMENT_D(A, i, i) = (mu_ip + mu_im) / (dy * dy);
        SM_ELEMENT_D(A, i, i + 1) = -mu_ip / (dy * dy);

        NV_Ith_S(b, i) = -dPdx;  // Forcing term
    }

    // Step 4: Apply Ghost Node BC at y=0
    real_t mu0 = viscosity(mesh[0] + dy / 2);
    SM_ELEMENT_D(A, 0, 0) = 3 * mu0 / (dy * dy);
    SM_ELEMENT_D(A, 0, 1) = -mu0 / (dy * dy);
    NV_Ith_S(b, 0) = dPdx + (2 * mu0 * U0 / (dy * dy));

    // Step 5: Apply Ghost Node BC at y=H
    real_t muN = viscosity(mesh[N - 1] - dy / 2);
    SM_ELEMENT_D(A, N - 1, N - 2) = -muN / (dy * dy);
    SM_ELEMENT_D(A, N - 1, N - 1) = 3 * muN / (dy * dy);
    NV_Ith_S(b, N - 1) = dPdx + (2 * muN * UH / (dy * dy));

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

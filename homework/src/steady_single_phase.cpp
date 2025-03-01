#include "steady_single_phase.hpp"


namespace CMF
{

namespace T1
{

bool validEntry(
    const std::vector<real_t>& mesh,
    const std::function<real_t(real_t)>& viscosity,
    const BC& bc)
{
    if (mesh.size() < 3)
    {
        LOG_ERROR("Mesh must have at least 3 nodes");
        return false;
    }

    if (bc.global.first != GlobalBC::PressureGradient)
    {
        LOG_ERROR("Only PressureGradient BC is supported");
        return false;
    }

    if (bc.lowerWall.first == WallBC::VelocityGradient
        && bc.upperWall.first == WallBC::VelocityGradient)
    {
        LOG_ERROR("Both walls cannot have VelocityGradient BC");
        return false;
    }

    return true;
}


void fillSystem(
    SUNMatrix A,
    N_Vector b,
    const std::vector<real_t>& mesh,
    const std::function<real_t(real_t)>& viscosity,
    const BC& bc)
{
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

        SM_ELEMENT_D(A, i, i - 1) = + mu_im / (dy_b * dy);
        SM_ELEMENT_D(A, i, i)     = - mu_ip / (dy_f * dy)
                                    - mu_im / (dy_b * dy);
        SM_ELEMENT_D(A, i, i + 1) = + mu_ip / (dy_f * dy);

        NV_Ith_S(b, i) = bc.global.second;
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
                           - bc.lowerWall.second / dy;
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
                               - bc.upperWall.second / dy;
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
    if (!validEntry(mesh, viscosity, bc))
    {
        const char* msg = "Invalid input";
        LOG_ERROR(msg);
        throw std::runtime_error(msg);
    }

    const size_t N = mesh.size();
    LinearSolver solver(N);

    fillSystem(solver.Matrix(), solver.RHS(), mesh, viscosity, bc);
    
    solver.solve();
    return solver.getSolutionVector();
}

} // namespace T1
} // namespace CMF

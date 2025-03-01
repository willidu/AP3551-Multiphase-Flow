#include "steady_single_phase.hpp"


namespace CMF
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


size_t findIndex(
    const std::vector<real_t>& mesh,
    real_t y
)
{
    const auto it = std::lower_bound(mesh.begin(), mesh.end(), y);

    if (it == mesh.end())
        return mesh.size() - 1;

    if (it == mesh.begin())
        return 0;

    return std::distance(mesh.begin(), it);
}


real_t velocityGradient(
    const std::vector<real_t>& mesh,
    const std::vector<real_t>& u,
    size_t i
)
{
    if (i == 0)
        return (u.at(1) - u.at(0)) / (mesh.at(1) - mesh.at(0));
    
    if (i == mesh.size() - 1)
        return (u.at(i) - u.at(i-1)) / (mesh.at(i) - mesh.at(i-1));

    return (u.at(i+1) - u.at(i-1)) / (mesh.at(i+1) - mesh.at(i-1));
}


std::vector<real_t> steadyChannelFlow(
    const std::vector<real_t>& mesh,
    std::function<real_t(real_t)> viscosity_mol,
    std::function<real_t(real_t)> mixingLength,
    BC bc,
    size_t maxIter,
    real_t tol
)
{
    if (!validEntry(mesh, viscosity_mol, bc))
    {
        const char* msg = "Invalid input";
        LOG_ERROR(msg);
        throw std::runtime_error(msg);
    }

    const size_t N = mesh.size();

    // We will fist assume that the eddy viscosity is zero
    std::vector<real_t> u;
    {
        LinearSolver solver(N);
        fillSystem(solver.Matrix(), solver.RHS(), mesh, viscosity_mol, bc);
        solver.solve();
        u = solver.getSolutionVector();
    }

    // We will now iterate to find a solution with non-zero eddy viscosity
    for (size_t iter = 0; iter < maxIter; ++iter)
    {
        LinearSolver solver(N);
        std::function<real_t(real_t)> viscosity = [&](real_t y) -> real_t
        {
            const real_t mu_mol = viscosity_mol(y);
            const size_t i = findIndex(mesh, y);
            const real_t dudy = velocityGradient(mesh, u, i);
            const real_t mu_t = std::pow(mixingLength(y), 2.0) * std::abs(dudy);

            return mu_mol + mu_t;
        };

        fillSystem(solver.Matrix(), solver.RHS(), mesh, viscosity, bc);
        solver.solve();
        std::vector<real_t> u_new = solver.getSolutionVector();

        // Convergence with L2 norm
        const real_t error = std::sqrt(std::inner_product(
            u.begin(), u.end(), u_new.begin(), 0.0, std::plus<real_t>(),
            [](real_t a, real_t b) { return std::pow(a - b, 2.0); }
        ));

        u = u_new;

        if (error < tol)
        {
            LOG_INFO("Converged after {0} iterations", iter);
            break;
        }

        if (iter == maxIter - 1)
        {
            LOG_WARN("Did not converge after {0} iterations", maxIter);
        }
    }

    return u;
}

} // namespace CMF

#include "steady_single_phase.hpp"


namespace CMF
{

bool validEntry(
    const Mesh& mesh,
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
    const Mesh& mesh,
    const BC& bc)
{
    for (size_t i = 1; i < mesh.size() - 1; ++i)
    {
        const auto& P = mesh[i];    // Center
        const auto& N = mesh[i+1];  // North
        const auto& S = mesh[i-1];  // South

        // Linear interpolation of viscosity at faces
        const real_t mu_n = (1.0 - P.width / N.width) * P.viscosity
                                 + P.width / N.width  * N.viscosity;
        const real_t mu_s = (1.0 - P.width / S.width) * P.viscosity
                                 + P.width / S.width  * S.viscosity;

        SM_ELEMENT_D(A, i, i + 1) = + mu_n * 2.0 / (P.width + N.width);
        SM_ELEMENT_D(A, i, i)     = - mu_n * 2.0 / (P.width + N.width)
                                    - mu_s * 2.0 / (P.width + S.width);
        SM_ELEMENT_D(A, i, i - 1) = + mu_s * 2.0 / (P.width + S.width);

        NV_Ith_S(b, i) = bc.global.second * P.width;
    }

    switch (bc.lowerWall.first)
    {
        case WallBC::Velocity:
        {
            const auto& P = mesh[0];
            const auto& N = mesh[1];

            const real_t mu_n = (1.0 - P.width / N.width) * P.viscosity
                                     + P.width / N.width  * N.viscosity;
            const real_t mu_s = P.viscosity; // TODO - Double check

            SM_ELEMENT_D(A, 0, 0) = + 2.0 * mu_n / (P.width + N.width)
                                    + 2.0 * mu_s / (P.width);
            SM_ELEMENT_D(A, 0, 1) = - 2.0 * mu_n / (P.width + N.width);
            NV_Ith_S(b, 0) = bc.global.second * P.width
                           + 2.0 * mu_s * bc.lowerWall.second / P.width;
            break;
        }

        case WallBC::VelocityGradient:
        {
            const auto& P = mesh[0];
            const auto& N = mesh[1];

            const real_t mu_n = (1.0 - P.width / N.width) * P.viscosity
                                     + P.width / N.width  * N.viscosity;
            const real_t mu_s = P.viscosity; // TODO - Double check

            SM_ELEMENT_D(A, 0, 0) = - 2.0 * mu_n / (P.width + N.width);
            SM_ELEMENT_D(A, 0, 1) = + 2.0 * mu_n / (P.width + N.width);
            NV_Ith_S(b, 0) = bc.global.second * P.width
                           + 2.0 * mu_s * bc.lowerWall.second / P.width;
            break;
        }

        case WallBC::WallFunction:
        {
            SM_ELEMENT_D(A, 0, 0) = 1.0;
            NV_Ith_S(b, 0) = bc.lowerWall.second;
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
            const size_t N = mesh.size() - 1;
            const auto& P = mesh[N];
            const auto& S = mesh[N - 1];

            const real_t mu_n = P.viscosity; // TODO - Double check
            const real_t mu_s = (1.0 - P.width / S.width) * P.viscosity
                                     + P.width / S.width  * S.viscosity;

            SM_ELEMENT_D(A, N, N)   = + 2.0 * mu_n / (P.width)
                                      + 2.0 * mu_s / (P.width + S.width);
            SM_ELEMENT_D(A, N, N-1) = - 2.0 * mu_s / (P.width + S.width);
            NV_Ith_S(b, N) = bc.global.second * P.width
                           + 2.0 * mu_n * bc.upperWall.second / P.width;
            break;
        }

        case WallBC::VelocityGradient:
        {
            const size_t N = mesh.size() - 1;
            const auto& P = mesh[N];
            const auto& S = mesh[N - 1];

            const real_t mu_n = P.viscosity; // TODO - Double check
            const real_t mu_s = (1.0 - P.width / S.width) * P.viscosity
                                     + P.width / S.width  * S.viscosity;

            SM_ELEMENT_D(A, N, N)   = - 2.0 * mu_n / (P.width + S.width);
            SM_ELEMENT_D(A, N, N-1) = + 2.0 * mu_s / (P.width + S.width);
            NV_Ith_S(b, N) = bc.global.second * P.width
                           - 2.0 * mu_n * bc.upperWall.second / P.width;
            break;
        }

        case WallBC::WallFunction:
        {
            const size_t N = mesh.size() - 1;
            SM_ELEMENT_D(A, N, N) = 1.0;
            NV_Ith_S(b, N) = bc.upperWall.second;
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


void steadyChannelFlow(
    Mesh& mesh,
    BC bc
)
{
    if (!validEntry(mesh, bc))
    {
        const char* msg = "Invalid input";
        LOG_ERROR(msg);
        throw std::runtime_error(msg);
    }

    const size_t N = mesh.size();
    LinearSolver solver(N);

    fillSystem(solver.Matrix(), solver.RHS(), mesh, bc);
    
    solver.solve();
    mesh.setVelocityProfile(solver.getSolutionVector());
}


real_t velocityGradient(
    const Mesh& mesh,
    const BC& bc,
    size_t i
)
{
    // Backward difference in upper half
    if (i > mesh.size() / 2)
        return (mesh[i].u - mesh[i-1].u) / mesh[i].width;
    
    // Forward difference in lower half
    if (i < mesh.size() / 2)
        return (mesh[i+1].u - mesh[i].u) / mesh[i].width;

    // Central difference in middle
    return (mesh[i+1].u - mesh[i-1].u) / (mesh[i+1].y - mesh[i-1].y);
}


void steadyChannelFlow(
    Mesh& mesh,
    const BC& bc,
    real_t viscosity_mol,  // Also stored in mesh, but will be overwritten
    real_t density,
    std::function<real_t(real_t)> mixingLength,
    size_t maxIter,
    real_t relax,
    real_t tol
)
{
    assert(viscosity_mol > 0.0 && "Invalid molecular viscosity");
    assert(density > 0.0 && "Invalid density");
    assert(mixingLength && "Invalid mixing length function");
    assert(bc.lowerWall.first == WallBC::Velocity
        && bc.upperWall.first == WallBC::Velocity
        && "Only Velocity BC is supported"
    );
    assert(relax > 0.0 && relax < 1.0 && "Invalid relaxation factor");

    if (!validEntry(mesh, bc))
    {
        const char* msg = "Invalid input";
        LOG_ERROR(msg);
        throw std::runtime_error(msg);
    }
    
    // Initial guess with no eddy viscosity
    std::vector<real_t> mu_effective(mesh.size(), viscosity_mol);
    mesh.setViscosityProfile(mu_effective);
    steadyChannelFlow(mesh, bc);
    std::vector<real_t> u = mesh.getSolution().second;
    std::vector<real_t> u_new;

    // Initial guess
    real_t u_tau = std::sqrt(mesh.height() * std::abs(bc.global.second) / (4.0 * density));

    // We will now iterate to find a solution with non-zero eddy viscosity
    for (size_t iter = 0; iter < maxIter; ++iter)
    {
        for (size_t i = 0; i < mesh.size(); ++i)
        {
            const real_t wallDistance = std::min(mesh[i].y, mesh.height() - mesh[i].y);
            const real_t yplus = density * u_tau * wallDistance / viscosity_mol;
            const real_t mu_eff = viscosity_mol + density
                            * std::pow(mixingLength(wallDistance), 2) 
                            * std::abs(velocityGradient(mesh, bc, i))
                            * vanDriest(yplus);

            mu_effective.at(i) = (1.0 - relax) * mesh[i].viscosity + relax * mu_eff;
            mesh[i].yplus = yplus;
        }
        mesh.setViscosityProfile(mu_effective);

        steadyChannelFlow(mesh, bc);
        u_new = mesh.getSolution().second;

        // Convergence with L2 norm
        const real_t error = std::sqrt(std::inner_product(
            u.begin(), u.end(), u_new.begin(), 0.0, std::plus<real_t>(),
            [](real_t a, real_t b) { return std::pow(a - b, 2.0); }
        ));

        if (error < tol)
        {
            LOG_INFO("Converged after {0} iterations", iter);
            break;
        }
        if (iter == maxIter - 1)
        {
            LOG_WARN(
                "Did not converge after {0} iterations. Final residual: {1:.2e}",
                maxIter, error);
        }
        u = u_new;
        u_tau = (1.0 - relax) * u_tau + relax * (
            u_new.at(0) / uplus(mesh[0].yplus, 0.0)
        );
    }
    LOG_TRACE("Final y+ at first node: {0:.2f}", mesh[0].yplus);
}


void steadyChannelFlow(
    Mesh& mesh,
    BC bc,
    std::function<real_t(real_t)> mixingLength,
    real_t viscosity_mol,  // Also stored in mesh, but will be overwritten
    real_t density,
    real_t averageRoughness,
    size_t maxIter,
    real_t relax,
    real_t tol
)
{
    assert(viscosity_mol > 0.0 && "Invalid molecular viscosity");
    assert(density > 0.0 && "Invalid density");
    assert(averageRoughness >= 0.0 && "Invalid average roughness");
    assert(mixingLength && "Invalid mixing length function");
    assert(bc.lowerWall.first == WallBC::Velocity
        && bc.upperWall.first == WallBC::Velocity
        && "Only Velocity BC is supported"
    );
    assert(relax > 0.0 && relax <= 1.0 && "Invalid relaxation factor");

    if (!validEntry(mesh, bc))
    {
        const char* msg = "Invalid input";
        LOG_ERROR(msg);
        throw std::runtime_error(msg);
    }
    
    std::vector<real_t> mu_effective(mesh.size(), viscosity_mol);
    std::vector<real_t> u = mesh.getSolution().second;
    std::vector<real_t> u_new;

    // Initial guess
    real_t u_tau = std::sqrt(mesh.height() * std::abs(bc.global.second) / (4.0 * density));

    static constexpr real_t kappa = 0.41;  // TODO - remove?

    // We will now iterate to find a solution with non-zero eddy viscosity
    for (size_t iter = 0; iter < maxIter; ++iter)
    {
        const real_t delta_plus = density * u_tau * 0.5 * mesh.height() / viscosity_mol;

        for (size_t i = 0; i < mesh.size(); ++i)
        {
            const real_t wallDistance = std::min(mesh[i].y, mesh.height() - mesh[i].y);
            const real_t yplus = density * u_tau * wallDistance / viscosity_mol;
            real_t mu_eff = viscosity_mol;
            
            if (yplus < 30.0 && i == 0)
            {
                LOG_WARN(
                    "Invalid y+ at iteration {0}: {1:.2f}. "
                    "Ignore this if final solution is converged.",
                    iter, yplus);
            }
            else if (yplus < 0.2 * delta_plus)
            {
                mu_eff += density * kappa * wallDistance * u_tau;
            }
            else if (yplus < 0.3 * delta_plus)
            {
                // Apply blending between log-law region and wake
                const real_t a = 0.5 * (1.0 + std::cos(M_PI * (yplus - 0.2 * delta_plus) / (0.1 * delta_plus)));
                mu_eff += density * kappa * wallDistance * u_tau * a;
                mu_eff += density
                        * std::pow(mixingLength(wallDistance), 2) 
                        * std::abs(velocityGradient(mesh, bc, i))
                        * (1.0 - a);
            }
            else
            {
                mu_eff += density
                        * std::pow(mixingLength(wallDistance), 2) 
                        * std::abs(velocityGradient(mesh, bc, i));
            }

            mu_effective.at(i) = (1.0 - relax) * mesh[i].viscosity + relax * mu_eff;
            mesh[i].yplus = yplus;
        }
        mesh.setViscosityProfile(mu_effective);

        // TODO - Add roughness
        const real_t slipVelocity = u_tau * uplus(mesh[0].yplus, 0.0);
        bc.lowerWall = {WallBC::WallFunction, slipVelocity};
        bc.upperWall = {WallBC::WallFunction, slipVelocity};

        steadyChannelFlow(mesh, bc);
        u_new = mesh.getSolution().second;

        // Convergence with L2 norm
        const real_t error = std::sqrt(std::inner_product(
            u.begin(), u.end(), u_new.begin(), 0.0, std::plus<real_t>(),
            [](real_t a, real_t b) { return std::pow(a - b, 2.0); }
        ));

        if (error < tol)
        {
            LOG_INFO("Converged after {0} iterations", iter);
            LOG_TRACE("Final y+ at first node: {0:.2f}", mesh[0].yplus);
            LOG_TRACE("Channel BL delta+ {0:.2f}", delta_plus);
            break;
        }
        if (iter == maxIter - 1)
        {
            LOG_WARN(
                "Did not converge after {0} iterations. Final residual: {1:.2e}",
                maxIter, error);
            LOG_TRACE("Final y+ at first node: {0:.2f}", mesh[0].yplus);
            LOG_TRACE("Channel BL delta+ {0:.2f}", delta_plus);
        }

        u = u_new;
        u_tau = (1.0 - relax) * u_tau + relax * (
                u_new.at(0) / uplus(mesh[0].yplus, 0.0)
        );
    }
}

} // namespace CMF

#include "timeStep.hpp"

std::pair<double, double> adaptTimestep(std::vector<double> &distributions,
                                        std::vector<double> &density, std::vector<double> &mass,
                                        std::vector<flag_t> &flags, double oldTau,
                                        double oldTimeStep, std::array<double, 3> &gravitation) {
    // We resize the timestep if the maximum velocity is too large.
    double maximumVelocityNorm = 0.0;

#pragma omp parallel for reduction(max : maximumVelocityNorm)
    for (size_t i = 0; i < flags.size(); ++i) {
        if (flags[i] != flag_t::FLUID && flags[i] != flag_t::INTERFACE)
            continue;

        const auto fieldIndex = i * Q;
        const auto curDensity = density[i];
        std::array<double, 3> velocity;
        computeVelocity(&distributions[fieldIndex], curDensity, velocity.data());
        const double curNorm = std::sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] +
                                         velocity[2] * velocity[2]);
        maximumVelocityNorm = std::max(maximumVelocityNorm, curNorm);
    }
    // Critical is half of the velocity for which numerics become unstable.
    const double criticalVelocity = 0.5 * C_S * C_S;
    const double multiplier = 4.0 / 5.0; // Safety factor.
    const double upperLimit = criticalVelocity / multiplier;
    const double lowerLimit = criticalVelocity * multiplier;

    double newTimeStep = oldTimeStep;
    if (maximumVelocityNorm > upperLimit) {
        // Decrease time step
        newTimeStep *= multiplier;
    } else if (maximumVelocityNorm < lowerLimit) {
        // Increase time step
        newTimeStep /= multiplier;
    } else {
        // Current time step works
        return std::pair<double, double>(oldTau, oldTimeStep);
    }
    const double timeRatio = newTimeStep / oldTimeStep;
    const double newTau = timeRatio * (oldTau - 0.5) + 0.5;
    const double tauRatio = timeRatio * ((1.0 / oldTau) / newTau);
    if (newTau < (1.0 / 1.99)) {
        // Time step would be too small!
        return std::pair<double, double>(oldTau, oldTimeStep);
    }

    // With the new time step, we need to rescale nearly everything.
    for (auto &g : gravitation) {
        g *= timeRatio * timeRatio;
    }

    double totalFluidVolume = 0.0;
    double totalMass = 0.0;
#pragma omp parallel for reduction(+ : totalFluidVolume, totalMass)
    for (size_t i = 0; i < flags.size(); ++i) {
        // TODO: Do we need to calculate this for boundary cells as well?
        if (flags[i] == flag_t::FLUID) {
            ++totalFluidVolume;
            totalMass += density[i];
        } else if (flags[i] == flag_t::INTERFACE) {
            totalMass += mass[i];
            totalFluidVolume += mass[i] / density[i];
        }
    }
    const double medianDensity = totalFluidVolume / totalMass;

    // Now we can rescale the distributions and the densities
    std::array<double, 3> oldVelocity;
    std::array<double, 3> newVelocity;
    std::array<double, Q> oldFeq;
    std::array<double, Q> newFeq;
#pragma omp parallel for private(oldVelocity, newVelocity, oldFeq, newFeq)
    for (size_t i = 0; i < density.size(); ++i) {
        // No need to rescale distributions of boundary cells, they are generated in each time step
        // anyway.
        if (flags[i] != flag_t::FLUID && flags[i] != flag_t::INTERFACE)
            continue;

        const auto fieldIndex = i * Q;

        const double oldDensity = density[i];
        const double newDensity = timeRatio * (density[i] - medianDensity) + medianDensity;

        computeVelocity(&distributions[fieldIndex], density[i], oldVelocity.data());
        newVelocity = oldVelocity;
        for (auto &v : newVelocity) {
            v *= timeRatio;
        }

        computeFeq(oldDensity, oldVelocity.data(), oldFeq.data());
        computeFeq(newDensity, newVelocity.data(), newFeq.data());

        for (int j = 0; j < Q; ++j) {
            const double feqRatio = newFeq[j] / oldFeq[j];
            // Rescale off-eq. parts of distributions.
            distributions[fieldIndex + j] =
                feqRatio * (oldFeq[j] + tauRatio * (distributions[fieldIndex + j] - oldFeq[j]));
        }

        if (flags[i] == flag_t::INTERFACE) {
            // We also need to track the change in mass!
            mass[i] = mass[i] * (oldDensity / newDensity);
        }
        density[i] = newDensity;
    }

    return std::pair<double, double>(newTau, newTimeStep);
}

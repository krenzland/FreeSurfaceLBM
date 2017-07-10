#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

#include <array>
#include <math.h>
#include <unordered_map>

using coord_t = std::array<int, 3>;

// Hash function for our coord_type.
namespace std
{
    template<> struct hash<coord_t>
    {
        using argument_type = coord_t;
        using result_type = size_t;
        result_type operator()(argument_type const& c) const
        {
            // We need to combine the hashes, this uses the algorithm used by boost's hash_combine.
            // See e.g. http://www.boost.org/doc/libs/1_35_0/doc/html/boost/hash_combine_id241013.html
            result_type h = 0;
            for (const auto i : c){
                // Black magic.
                h ^= std::hash<int>{}(i) + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            return h;
        }
    };
}

enum class flag_t {
    FLUID = 0,
    NO_SLIP = 1,
    MOVING_WALL = 2,
    FREE_SLIP = 3,
    INFLOW = 4,
    OUTFLOW = 5,
    PRESSURE_IN = 6,
    PARALLEL_BOUNDARY = 7,
    EMPTY = 8,
    INTERFACE = 9,
};

// This is needed to parse the config file.
const auto stringToFlag =
    std::unordered_map<std::string, flag_t>({{"FLUID", flag_t::FLUID},
                                             {"NO_SLIP", flag_t::NO_SLIP},
                                             {"MOVING_WALL", flag_t::MOVING_WALL},
                                             {"FREE_SLIP", flag_t::FREE_SLIP},
                                             {"INFLOW", flag_t::INFLOW},
                                             {"OUTFLOW", flag_t::OUTFLOW},
                                             {"PRESSURE_IN", flag_t::PRESSURE_IN},
                                             {"PARALLEL_BOUNDARY", flag_t::PARALLEL_BOUNDARY}});

struct boundary_t {
    flag_t bcLeft;
    flag_t bcRight;
    flag_t bcTop;
    flag_t bcBottom;
    flag_t bcFront;
    flag_t bcBack;
    double velocityWall[3];
    double velocityIn[3];
    double pressureIn;
};

static const std::array<std::array<int, 3>, 19> LATTICEVELOCITIES = {{
    {0, -1, -1}, {-1, 0, -1}, {0, 0, -1}, {1, 0, -1}, {0, 1, -1}, {-1, -1, 0}, {0, -1, 0},
    {1, -1, 0},  {-1, 0, 0},  {0, 0, 0},  {1, 0, 0},  {-1, 1, 0}, {0, 1, 0},   {1, 1, 0},
    {0, -1, 1},  {-1, 0, 1},  {0, 0, 1},  {1, 0, 1},  {0, 1, 1}}};

static const double LATTICEWEIGHTS[19] = {(1. / 36), (1. / 36), (2. / 36), (1. / 36), (1. / 36),
                                          (1. / 36), (2. / 36), (1. / 36), (2. / 36), (12. / 36),
                                          (2. / 36), (1. / 36), (2. / 36), (1. / 36), (1. / 36),
                                          (1. / 36), (2. / 36), (1. / 36), (1. / 36)};

static const double C_S = 1.0 / std::sqrt(3);

static const int Q = 19;

// otherDims[i] gives all dims that are not i.
const static std::array<std::array<int, 2>, 3> otherDims = {{{1, 2}, {0, 2}, {0, 1}}};

#endif

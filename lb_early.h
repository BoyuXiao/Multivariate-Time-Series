#ifndef LB_EARLY_H
#define LB_EARLY_H

#include <cmath>
#include <limits>

inline bool lb_can_early_stop(double best_dist) {
    return std::isfinite(best_dist) && best_dist < 1e100;
}

inline double lb_threshold_sq(double best_dist) {
    return best_dist * best_dist;
}

#endif

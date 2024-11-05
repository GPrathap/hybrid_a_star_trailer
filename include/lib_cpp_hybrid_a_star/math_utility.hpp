#ifndef MATH_UTILITY
#define MATH_UTILITY

#include <cmath>
#include <limits>

namespace math_utility
{
    inline double pi_to_pi(double angle){
        while (angle > M_PI){
            angle -= 2.0*M_PI;
        }
        while (angle < -M_PI)
        {
            angle += 2.0*M_PI;
        }
        return angle;
    }

    constexpr inline double deg2rad(double degrees) {
        return degrees * (M_PI / 180.0);
    }

    inline  double mod2pi(double x){
        double v = std::fmod(x, (2.0*M_PI));
        if (v < -M_PI){
            v += 2.0*M_PI;
        }else{
            if(v > M_PI){
                v -= 2.0*M_PI;
            }
        }
        return v;
    }

    inline void polar(double x, double y, double& r, double& theta){
        r = sqrt(pow(x, 2) + pow(y, 2));
        theta = atan2(y, x);
    }


    
} // namespace math_utility


#endif
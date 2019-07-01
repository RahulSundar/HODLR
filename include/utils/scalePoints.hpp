#ifndef __scalePoints_hpp__
#define __scalePoints_hpp__

// Maps the points from a domain having center c1 and radius r1
// To a domain having center c2 and radius r2
// Basically mapping from [c1 - r1 / 2, c1 + r1 / 2] --> [c2 - r2 / 2, c2 + r2 / 2]
Vec scalePoints(const double c1, const double r1, const Vec &x1,
                 const double c2, const double r2)
{   
    return (c2 + (r2 * (x1.array() - c1)).array() / r1);
}

#endif

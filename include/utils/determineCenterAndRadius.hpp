#ifndef __determineCenterAndRadius_hpp__
#define __determineCenterAndRadius_hpp__

void determineCenterAndRadius(const Vec &input_array, double &center, double &radius)
{
    double max_val, min_val;
    max_val = input_array.maxCoeff();
    min_val = input_array.minCoeff();

    center = 0.5 * (max_val + min_val);
    radius = 0.5 * (max_val - min_val);
}

#endif

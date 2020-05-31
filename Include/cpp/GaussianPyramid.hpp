#ifndef GAUSSIAN_PYRAMID_H
#define GAUSSIAN_PYRAMID_H
 
 
#include "Image.hpp"
 
class gaussian_pyramid
{
private:
    DImage* image_pyramid;
    int pyramid_levels;
public:
    gaussian_pyramid(void);
    ~gaussian_pyramid(void);
    void form_image_pyramid(const DImage& image, double ratio=0.75, int min_width=30);
    
    inline int get_pyramid_levels() const {return pyramid_levels;};
    inline DImage& Image(int index) {return image_pyramid[index];};
};
 
#endif

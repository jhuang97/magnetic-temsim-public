#ifndef VZAVGLUT_HPP
#define VZAVGLUT_HPP

#include "slicelib.hpp"
#include "cfpix.hpp"
#include "newD.hpp"

class VzPxavgMaker {
    public:
        const double rmax = 3.0;
        double spacing;
        bool LUT_allocated;
        bool LUT_calculated;
        VzPxavgMaker();
        ~VzPxavgMaker();
        void allocate_LUT(double _spacing, vectori _Zs, double _px_x, double _px_y);
        void calculate_LUT();
        void fill_vz(const vectorf &x, const vectorf &y, const vectorf &occ,
            const vectori &Znum, const int natom, const int istart,
            const float ax, const float by,
            cfpix &trans, const int nx, const int ny);
        void check_sums();
        void check_map(int Z_idx, int map_idx);
    private:
        vectori Zs;
        double px_x, px_y;
        int map_radius_x_px, map_radius_y_px, map_size_x, map_size_y;
        int num_maps_x_radius, num_maps_y_radius, num_maps_x, num_maps_y;
        double spacing_x, spacing_y;
        int num_elems, num_maps, map_num_pixels;
        float*** vz_LUT = nullptr;
};


#endif // VZAVGLUT_HPP

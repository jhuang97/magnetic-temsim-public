#include "vzavgLUT.hpp"
#include "slicelib.hpp"

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;

VzPxavgMaker::VzPxavgMaker(  ) {
    LUT_allocated = false;
    LUT_calculated = false;
    vz_LUT = nullptr;
}

void VzPxavgMaker::allocate_LUT(double _spacing, vectori _Zs, double _px_x, double _px_y) {
    spacing = _spacing;
    Zs = _Zs;
    px_x = _px_x;
    px_y = _px_y;
    map_radius_x_px = static_cast<int>(floor(rmax/px_x + 0.5));
    map_radius_y_px = static_cast<int>(floor(rmax/px_y + 0.5));
    map_size_x = 2*map_radius_x_px + 1;
    map_size_y = 2*map_radius_y_px + 1;

    num_maps_x_radius = static_cast<int>(ceil(px_x/2.0/spacing));
    num_maps_y_radius = static_cast<int>(ceil(px_y/2.0/spacing));
    spacing_x = px_x/2.0/num_maps_x_radius;
    spacing_y = px_y/2.0/num_maps_y_radius;
    num_maps_x = 2*num_maps_x_radius + 1;
    num_maps_y = 2*num_maps_y_radius + 1;

    num_elems = Zs.size();
    map_num_pixels = map_size_x * map_size_y;
    num_maps = num_maps_x * num_maps_y;

    printf("Map spacing: %f %f, map size (px): %d %d, # maps: %d %d\n", spacing_x, spacing_y, map_size_x, map_size_y, num_maps_x, num_maps_y);
    printf("Allocating array with dimensions %d x %d x %d = %d\n", num_elems, num_maps, map_num_pixels, num_elems*num_maps*map_num_pixels);
    vz_LUT = new3D<float>( num_elems, num_maps, map_num_pixels, "vz_LUT" );
    LUT_allocated = true;
}

void VzPxavgMaker::calculate_LUT() {
    double rmax2 = rmax*rmax;
    for (size_t iZ = 0; iZ < Zs.size(); ++iZ) {
        int Z = Zs[iZ];
        // (x0, y0) is the position of the atom relative to the central pixel
        for (int ix0 = 0; ix0 < num_maps_x; ++ix0) {
            double x0 = (ix0 - num_maps_x_radius)*spacing_x;
            int ix0_offset = ix0*num_maps_y;
            for (int iy0 = 0; iy0 < num_maps_y; ++iy0) {
                double y0 = (iy0 - num_maps_y_radius)*spacing_y;
                int map_idx = ix0_offset + iy0;
                if (map_idx % 10 == 0) {
                    printf("%d/%d ", map_idx, num_maps);
                    fflush(stdout);
                }
#pragma omp parallel for
                for (int ixp = 0; ixp < map_size_x; ++ixp) {
                    double dx = (ixp - map_radius_x_px) * px_x - x0;
                    int ixp_offset = ixp*map_size_y;
                    for (int iyp = 0; iyp < map_size_y; ++iyp) {
                        double dy = (iyp - map_radius_y_px) * px_y - y0;
                        if (dx*dx + dy*dy > rmax2) {
                            vz_LUT[iZ][map_idx][ixp_offset + iyp] = 0.0;
                        } else {
                            vz_LUT[iZ][map_idx][ixp_offset + iyp] =
                                (float) vzatom_pxavg(Z, dx, dy, px_x, px_y);
                        }
                    }
                }
            }
        }
        cout << endl << "Generated images for Z = " << Z << endl;
    }
    LUT_calculated = true;
}

void VzPxavgMaker::fill_vz(const vectorf &x, const vectorf &y, const vectorf &occ,
            const vectori &Znum, const int natom, const int istart,
            const float ax, const float by,
            cfpix &trans, const int nx, const int ny) {
    for (int i = istart; i < (istart + natom); ++i) {
        float x_px = (float) (x[i] / px_x);
        float y_px = (float) (y[i] / px_y);
        float x_px_i = roundf(x_px);
        float y_px_i = roundf(y_px);
        // figure out which map(s) to use
        float ixo = x_px-x_px_i;
        float iyo = y_px-y_px_i;
//        printf("%f %f; %f %f; %f %f\n", x[i], y[i], x_px, y_px, ixo, iyo);
        ixo = (ixo + 0.5f) * 2.0f * num_maps_x_radius;
        iyo = (iyo + 0.5f) * 2.0f * num_maps_y_radius;
        float ixo_idx, iyo_idx;
        ixo = modf(ixo, &ixo_idx);
        iyo = modf(iyo, &iyo_idx);
        int ixo_int = static_cast<int>(ixo_idx);
        int iyo_int = static_cast<int>(iyo_idx);
//        printf("%f %f; %f %f; %f %f; %d %d\n", x[i], y[i], x_px, y_px, ixo, iyo, ixo_int, iyo_int);

        int num_maps = 0;
        int   m1, m2, m3, m4;  // map indices
        float w1, w2, w3, w4;  // weights
        if ((ixo == 0) && (iyo == 0)) {
            num_maps = 1;
            m1 = ixo_int * num_maps_y + iyo_int;
        } else if (ixo == 0) {
            num_maps = 2;
            m1 = ixo_int * num_maps_y + iyo_int;
            m2 = m1 + 1;
            w1 = 1.0f - iyo;
            w2 = iyo;
        } else if (iyo == 0) {
            num_maps = 2;
            m1 = ixo_int * num_maps_y + iyo_int;
            m2 = (ixo_int+1) * num_maps_y + iyo_int;
            w1 = 1.0f - ixo;
            w2 = ixo;
        } else {
            num_maps = 4;
            m1 = ixo_int * num_maps_y + iyo_int;
            m2 = (ixo_int+1) * num_maps_y + iyo_int;
            m3 = (ixo_int+1) * num_maps_y + iyo_int+1;
            m4 = ixo_int * num_maps_y + iyo_int+1;
            w1 = (1.0f - ixo) * (1.0f - iyo);
            w2 = ixo * (1.0f - iyo);
            w3 = ixo * iyo;
            w4 = (1.0f - ixo) * iyo;
//            printf("%f %f; %d %d; %f %f %f %f\n", ixo, iyo, ixo_int, iyo_int, w1, w2, w3, w4);
        }

        // where map starts
        int ix_start = static_cast<int>(x_px_i) - map_radius_x_px;
        int iy_start = static_cast<int>(y_px_i) - map_radius_y_px;
        if (ix_start < 0)   ix_start += nx;
        if (ix_start >= nx) ix_start -= nx;
        if (iy_start < 0)   iy_start += ny;
        if (iy_start >= ny) iy_start -= ny;

        // which element
        auto it = find(Zs.begin(), Zs.end(), Znum[i]);
        int Zidx;
        if (it != Zs.end()) {
            Zidx = it - Zs.begin();
        } else {
            cout << "oh no, element " << Znum[i] << " not in table" << endl;
        }

        // fill in
        int ixw, iyw;
        for (int ixp = 0; ixp < map_size_x; ++ixp) {
            ixw = ix_start + ixp;
            if (ixw >= nx) ixw -= nx;

            int ixp_offset = ixp*map_size_y;
            for (int iyp = 0; iyp < map_size_y; ++iyp) {
                iyw = iy_start + iyp;
                if (iyw >= ny) iyw -= ny;
                int idxp = ixp_offset + iyp;

                if (num_maps == 1) {
                    trans.re(ixw, iyw) += occ[i]*vz_LUT[Zidx][m1][idxp];
                } else if (num_maps == 2) {
                    trans.re(ixw, iyw) += occ[i]*
                    (w1*vz_LUT[Zidx][m1][idxp] + w2*vz_LUT[Zidx][m2][idxp]);
                } else if (num_maps == 4) {
                    trans.re(ixw, iyw) += occ[i]*
                    (w1*vz_LUT[Zidx][m1][idxp] + w2*vz_LUT[Zidx][m2][idxp]
                   + w3*vz_LUT[Zidx][m3][idxp] + w4*vz_LUT[Zidx][m4][idxp]);
                }
            }
        }
    }
}

void VzPxavgMaker::check_sums() {
    if (!LUT_allocated) {
        cout << "LUT not yet allocated" << endl;
        return;
    }
    for (int iZ = 0; iZ < num_elems; ++iZ) {
        cout << "Z = " << Zs[iZ] << ":" << endl;
        vectord vzsums;
        for (int im = 0; im < num_maps; ++im) {
            double vzsum = 0.0;
            for (int ip = 0; ip < map_num_pixels; ++ip) {
                vzsum += vz_LUT[iZ][im][ip];
            }
            vzsums.push_back(vzsum * px_x * px_y);
        }
        double sum = accumulate(vzsums.begin(), vzsums.end(), 0.0);
        double mean = sum / vzsums.size();

        vector<double> diff(vzsums.size());
        transform(vzsums.begin(), vzsums.end(), diff.begin(), [mean](double x) { return x - mean; });
        double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        double stdev = sqrt(sq_sum / vzsums.size());
        printf("mean: %f, std dev: %f, %e\n", mean, stdev, stdev/mean);
    }
}

void VzPxavgMaker::check_map(int Z_idx, int map_idx) {
    if (!LUT_allocated) {
        cout << "LUT not yet allocated" << endl;
        return;
    }
    for (int iy = 0; iy < map_size_y; ++iy) {
        for (int ix = 0; ix < map_size_x; ++ix) {
            printf("%8.4e ", vz_LUT[Z_idx][map_idx][map_size_y*ix + iy]);
        }
        printf("\n");
    }
}

VzPxavgMaker::~VzPxavgMaker() {
    if (LUT_allocated) {
        delete3D(vz_LUT, num_elems, num_maps);
    }
}

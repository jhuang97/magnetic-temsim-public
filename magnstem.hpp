#pragma once

#ifndef MAGNSTEM_HPP   // only include this file if its not already
#define MAGNSTEM_HPP   // remember that this has been included

#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <string>   // STD string class
#include <vector>

#include <iostream>

using namespace std;

#include "cfpix.hpp"       // complex image handler with FFT 
#include "slicelib.hpp"    // misc. routines for multislice 
#include "newD.hpp"      //  for 2D and 3D arrays
#include "paramgroup.hpp"
#include "vzavgLUT.hpp"

#include <algorithm>
#include "floatTIFF.hpp"   // file I/O routines in TIFF format
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

// modes of collector
enum { ADF = 0, CONFOCAL = 1, ADF_SEG = 2, TOTAL = 3 };

class magnstem {

public:

    magnstem();         // constructor functions
    ~magnstem();        //  destructor function
    int lwobble, lxzimage, lpacbed, lverbose;
    int lrepeatz;
    // removed l1d.  assume never do 1D mode.

    int calculate(vectorf& param, int multiMode, int natomin, 
        double xi, double xf, double yi, double yf, int nxout, int nyout,
        vectori& Znum, vectorf& xa, vectorf& ya, vectorf& za, vectorf& occ, vectorf& wobble,
        vector<cfpix>& mag_phase, unsigned long iseed,
        vectori& ThickSave, int nThick,
        vectord& almin, vectord& almax, vectori& collectorMode, int ndetect,
        vectord& phiMin, vectord& phiMax,
        float*** pixr, float** rmin, float** rmax,
        float** pacbedPix, opts_4dsave opts4d);
    
    // miscellaneous information that may be used in calling program
    long nbeamt;
    double totmin, totmax;

    void magnstem::trlayer(const vectorf& x, const vectorf& y, const vectorf& occ,
        const vectori& Znum, const int natom, const int istart,
        const float ax, const float by, const float kev,
        cfpix& trans, cfpix& mag_image, const long nx, const long ny,
        double* phirms, long* nbeams, const float k2max);

    void magnstem::trlayer_pxavg(const vectorf& x, const vectorf& y, const vectorf& occ,
        const vectori& Znum, const int natom, const int istart,
        const float ax, const float by, const float kev,
        cfpix& trans, cfpix& mag_image, const long nx, const long ny,
        double* phirms, long* nbeams, const float k2max,
        const int it, vectorf& param, opts_4dsave opts4d);

    void magnstem::trlayer_pxavg_LUT(const vectorf& x, const vectorf& y, const vectorf& occ,
        const vectori& Znum, const int natom, const int istart,
        const float ax, const float by, const float kev,
        cfpix& trans, cfpix& mag_image, const long nx, const long ny,
        double* phirms, long* nbeams, const float k2max, VzPxavgMaker& maker);
        //const int it, vectorf& param, opts_4dsave opts4d);

    void magnstem::CountBeams(vectorf& param, int& nbeamp, int& nbeampo, float& res, float& almax);

private:
    int NZMAX, NRMAX;
    float BW;
    double ABERR, RMIN, RMAX;
    enum { xFALSE = 0, xTRUE = 1 };
    double pi;
    double subpixel_spacing;
    int use_LUT_when_no_phonons;

    // structure
    vectorf xa, ya, za, occ, wobble;
    vectorf xa2, ya2, za2, occ2;
    vectori Znum, Znum2;
    int slice_per_cell;
    vectori slice_start_idx;
    vectori slice_end_idx;  // inclusive

    // from params
    int nx, ny, nxprobe, nyprobe, nslice;
    float ax, by, cz;                   //  specimen dimensions
    double keV, df, apert1, apert2, deltaz;
    int natom;

    double wavlen;

    /* extra for confocal mode */
    int doConfocal;
    float dfa2C, dfa2phiC, dfa3C, dfa3phiC;   /* astigmatism parameters */
    double* collectMin, * collectMax;
    double Cs3C, Cs5C, dfC, apert1C, apert2C; /* aberrations of collector lens */

    std::string sbuffer;
    void messageAST(std::string& smsg, int level = 0);  // common error message handler

    vectorf kx, ky, kx2, ky2, kxp, kyp, kxp2, kyp2;
    double k2maxp, k2maxt;
    vectorf xp, yp;
    vectord k2max, k2min;


    cfpix cprop;           // complex propagator in Fourier space

    // complex probe and transmission functions
    cfpix* probe;
    vector<cfpix> trans;
    VzPxavgMaker vz_maker;

    double periodic(double pos, double size);
    void magnstem::STEMsignals(vectord& x, vectord& y, int npos, vectorf& p,
        vector<cfpix>& trans, 
        int multiMode, double*** detect, int ndetect,
        vectori& ThickSave, int nThick, vectord& sum, vectori& collectorMode,
        vectord& phiMin, vectord& phiMax, opts_4dsave opts4d);

    void magnstem::STEMsignals_wobble(double xp, double yp, int pos_ix, int pos_iy, int nwobble,
        vectorf& p, int multiMode, 
        vectori& Znum, vectorf& xa, vectorf& ya, vectorf& za, vectorf& occ, vectorf& wobble,
        vector<unsigned long>& iseeds, vector<cfpix>& mag_phase,
        vectori& ThickSave, float** ca_cbed,
        double** detect, int ndetect, double& sum, vectori& collectorMode, VzPxavgMaker& maker,
        vectord& phiMin, vectord& phiMax, opts_4dsave opts4d);

    void magnstem::invert2D(float** pix, long nx, long ny);

    void magnstem::write_cbed_from_probe(cfpix& probe, const string path, 
        vectorf& param, opts_4dsave opts4d);
    void magnstem::write_cbed(float** pix, const string path,
        vectorf& param, opts_4dsave opts4d);
    void magnstem::write_wvfn(cfpix& probe, const string path, vectorf& param);
    void magnstem::write_real_space_real(cfpix& probe, const string path, vectorf& param);

    void magnstem::record_probe_intensity_distributions(cfpix& probe, floatTIFF& xz_img, floatTIFF& yz_img, int z_idx);

};

#endif
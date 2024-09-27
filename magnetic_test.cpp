#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

#include <string>
#include <iostream>  //  C++ stream IO
#include <fstream>
#include <iomanip>   //  to format the output
#include <vector>   // STD vector class

using namespace std;

#include "cfpix.hpp"       // complex image handler with FFT
#include "slicelib.hpp"    // misc. routines for multislice
#include "floatTIFF.hpp"   // file I/O routines in TIFF format
#include "probe.hpp"       // probe calculation

//  use only one;  the calculation part
//#include "autoslic_cuda.hpp"    //  header for cuda nvcc version of this program
#include "autoslic.hpp"    //  header for C++ version of this class
#include "magnstem.hpp"
#include "paramgroup.hpp"

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#ifdef USE_OPENMP
#include <omp.h>
/*  get wall time for benchmarking openMP */
#define walltim() ( omp_get_wtime() )
double walltimer;
#endif

// only works for windows?
#include <direct.h>
#define GetCurrentDir _getcwd

enum { TRUE = 1, FALSE = 0 };

std::string get_current_dir() {
    char buff[FILENAME_MAX]; //create string buffer to hold path
    GetCurrentDir(buff, FILENAME_MAX);
    string current_working_dir(buff);
    return current_working_dir;
}

// image size in pixels nx, ny
vectorf get_params(float apert_mrad, int nx, int ny, int n_phonons) {
    floatTIFF myFile;

    cout << get_current_dir() << endl;

    cout << "Earl Kirkland's multislice, modified by Jeff\n";
    double pi = 4.0 * atan(1.0);

    vectorf param;
    /*  many of the parameters will be passed in as a vector of floats */
    int i;
    int NPARAM = myFile.maxParam();
    param.resize(NPARAM);
    for (i = 0; i < NPARAM; i++) param[i] = 0.0F;

    //int mode = MPWAVE;
    //string probe_fileout = "magnetic_test_images\\probe.tif";

    param[pNXPRB] = (float)nx;
    param[pNYPRB] = (float)ny;

    // size of output image in angstroms
    param[pNX] = (float)nx;
    param[pNY] = (float)ny;

    //double ax = 50.7387;
    //double by = 46.8772;
    //float dx, dy; // size of pixel
    //param[pAX] = (float)ax;
    //param[pBY] = (float)by;
    //param[pDX] = dx = (float)(ax / nx);
    //param[pDY] = dy = (float)(by / ny);
    // inverse of output image size??
    //double rx = 1.0 / ax;
    //double rx2 = rx * rx;
    //double ry = 1.0 / by;
    //double ry2 = ry * ry;

    int ixmid = nx / 2;
    int iymid = ny / 2;

    // probe parameters
    double keV = 300; // V0 (kV)
    double v0 = keV;
    double Cs3 = 0.0; // mm
    double Cs5 = 0.0; // mm
    double df = 0.0; // angstroms
    //double apert = 3.4; // mrad
    double apert = (double)apert_mrad;
    param[pENERGY] = (float)keV;
    param[pCS] = (float)(Cs3 * 1.0e7);
    param[pCS5] = (float)(Cs5 * 1.0e7);
    param[pDEFOCUS] = (float)df;
    param[pOAPERT] = (float)(apert / 1000.0);
    //param[pYCTILT] = 0.0014; // tilt in radians

    // Spacing between slices in angstroms
    //     This depends on the structure. For crystals it is best to let the slices correspond to planes of atoms.
    param[pDELTAZ] = 2.013F;

    if (n_phonons > 1) {
        param[pTEMPER] = 300.0;
        param[pNWOBBLE] = (float)n_phonons;
    }

    return param;
}

int magnstem_test(string mdir, unsigned long iseed) {
    magnstem mstem;
    double walltimer = walltim();  /* wall time for openMP */

    const int NSMAX = 1000;   // max number of slices
    const int NCMAX = 1024;   // max characters in file names
    const int NZMAX = 103;    // max atomic number Z
    const float tiny = 1.0e-4F;

    string file_coords = "In2Se3_m100.xyz";
    int ncellx = 16;
    int ncelly = 4;
    int ncellz = 1; // should be 1.  If want thicker sample, handle it using pNCELLREPZ

    float apert_mrad = 1.267; // convergence semiangle
    int nx = 2048; int ny = 2048; // image size in pixels

    // read in specimen coords
    float ax_f, by_f, cz_f;
    vectori Znum;
    vectorf x, y, z, occ, wobble;
    string description;
    int natom = ReadXYZcoord(file_coords.c_str(), ncellx, ncelly, ncellz,
        &ax_f, &by_f, &cz_f, Znum, x, y, z, occ, wobble,
        description);

    if (natom < 1) {
        cout << "no atom coord. read in from file " << file_coords << endl;
        cout << "can not continue...." << endl;
        return EXIT_SUCCESS;
    }

    if ((ax_f < tiny) || (by_f < tiny) || (cz_f < tiny)) {
        cout << "Bad unit cell size ax,by,cz = " << ax_f << ", " << by_f << ", " << cz_f << endl;
        exit(0);
    }

    //for (int i = 0; i < natom; i++) {
    //    cout << Znum[i] << " " << x[i] << " " << y[i] << " " << z[i] << " " << occ[i] << " " << wobble[i] << "\n";
    //}

    cout << natom << " atomic coordinates read in" << endl;
    cout << description << endl;

    //cout << "Size in pixels Nx, Ny= " << nx << " x " << ny << " = " << nx * ny
    //    << " beams" << endl;
    cout << "Lattice constant a,b = " << ax_f << ", " << by_f << endl;

    int n_phonons = 10;  // calculation with phonons enabled if n_phonons > 1
    vectorf params = get_params(apert_mrad, nx, ny, n_phonons);
    
    // thickness
    int ncellrepz = 300;  // Assuming lrepeatz == 1 in magnstem, ncellrepz is the number of times the slices are repeated in the z direction
    vectori ThickSave;  // indices at which to save CBED patterns, ranging from 0 to ncellrepz-1 inclusive
    //ThickSave.push_back(1);
    for (int k = 0; k < 300; k += 3) {
        ThickSave.push_back(k);
    }
    /*for (int k = 101; k < 150; k += 3) {
        ThickSave.push_back(k);
    }*/
    //ThickSave.push_back(0);

    ThickSave.push_back(ncellrepz - 1);
    int nThick = ThickSave.size();

    params[pNCELLREPZ] = ncellrepz;
    
    //float mag_phase_multiplier = 0;

    // Probe positions
    int num_probex = -1; // was 15
    int num_probey = -1; // was 45

    double xi = 0; double xf = 0; double yi = 0; double yf = 0;
    //xi = ax_f * 0.35;
    //xf = ax_f * 0.65;
    //yi = by_f * 0.35;
    //yf = by_f * 0.65;
    //xi = 14.1188;
    //xf = 36.6199;
    //yi = 12.1881;
    //yf = 34.6891;

    int probe_pos_version = 5;
    if (probe_pos_version == 1 || probe_pos_version == 2) {
        if (probe_pos_version == 1) {
            num_probex = 25; num_probey = 25;
        }
        else if (probe_pos_version == 2) {
            num_probex = 15;
            num_probey = 45;
        }
        float cell_sizex = ax_f / ncellx;
        float cell_sizey = by_f / ncelly;
        float scan_step_sizex;
        float scan_step_sizey = cell_sizey / num_probey; // two magnetic unit cells along c-axis
        if (probe_pos_version == 2)
            scan_step_sizex = 7.248392 / num_probex; // two unit cells along a-axis
        else scan_step_sizex = scan_step_sizey;
        int xni = -num_probex / 2;
        int xnf = num_probex / 2 + (num_probex % 2 != 0) - 1;
        int yni = -num_probey / 2;
        int ynf = num_probey / 2 + (num_probey % 2 != 0) - 1;
        xi = double(xni) * scan_step_sizex + ax_f / 2;
        xf = double(xnf) * scan_step_sizex + ax_f / 2;
        yi = double(yni) * scan_step_sizey + by_f / 2;
        yf = double(ynf) * scan_step_sizey + by_f / 2;
    }
    else if (probe_pos_version == 5) {
        num_probex = 1;
        num_probey = 9;
        float x_cell_length = 4.026;
        float y_cell_length = 9.583333333;
        float scan_step_sizex = x_cell_length / 2 / num_probex; // by symmetry, only need to cover half the u.c. in x-direction
        float scan_step_sizey = y_cell_length / num_probey;
        int xni = 0;
        int xnf = num_probex - 1;
        int yni = 0;
        int ynf = num_probey - 1;
        xi = double(xni) * scan_step_sizex + ax_f / 2;
        xf = double(xnf) * scan_step_sizex + ax_f / 2;
        yi = double(yni) * scan_step_sizey + by_f / 2;
        yf = double(ynf) * scan_step_sizey + by_f / 2;
    }
    else if (probe_pos_version == 3) {
        num_probex = 1;
        num_probey = 1;
        xi = ax_f / 2;
        xf = ax_f / 2;
        yi = by_f / 2;
        yf = by_f / 2;
    }
    else if (probe_pos_version == 4 || probe_pos_version == 6) {
        num_probex = 15;
        num_probey = 45;
        float cell_sizex = ax_f / ncellx;
        float cell_sizey = by_f / ncelly;
        float scan_step_sizex = 7.248392 / num_probex; // two unit cells along a-axis
        float scan_step_sizey = cell_sizey / num_probey; // two magnetic unit cells along c-axis
        int xni = -num_probex / 2;
        int xnf = num_probex / 2 + (num_probex % 2 != 0) - 1;
        int yni = -num_probey / 2;
        int ynf = num_probey / 2 + (num_probey % 2 != 0) - 1;
        xi = double(xni) * scan_step_sizex + ax_f / 2;
        xf = double(xnf) * scan_step_sizex + ax_f / 2;
        yi = double(yni) * scan_step_sizey + by_f / 2;
        yf = double(ynf) * scan_step_sizey + by_f / 2;

        double dx = (xf - xi) / ((double)(num_probex - 1));
        double dy = (yf - yi) / ((double)(num_probey - 1));
        if (probe_pos_version == 4) {
            // choose only one column of probe positions (only one x value)
            num_probex = 1;
            int x_index = 3;
            xi = xi + dx * ((double)x_index);
            xf = xi;
        }
        else if (probe_pos_version == 6) {
            num_probey = 1;
            yi = by_f / 2;
            yf = yi;

            num_probex = 8;
            // xi unchanged
            xf = xi + dx * ((double)num_probex - 1);
        }
    }

    cout << "Probes: " << endl;
    cout << "  x: " << xi << " .. " << xf << " (" << num_probex << ")" << endl;
    cout << "  y: " << yi << " .. " << yf << " (" << num_probey << ")" << endl;
    
    //double ax = 50.7387;
    //double by = 46.8772;
    float dx, dy; // size of pixel
    params[pAX] = ax_f;
    params[pBY] = by_f;
    params[pDX] = dx = (float)(ax_f / nx);
    params[pDY] = dy = (float)(by_f / ny);
     // inverse of output image size??
    double rx = 1.0 / ax_f;
    double rx2 = rx * rx;
    double ry = 1.0 / by_f;
    double ry2 = ry * ry;

    params[pCZ] = cz_f;

    vector<cfpix> mag_phase;
    // Here, if you want, you could add code to read magnetic phase shift from file and load it into mag_phase

    // specify detectors
    vectori collectorMode;
    collectorMode.push_back(ADF);
    vectord almin;
    vectord almax;
    almin.push_back(60.0 * 0.001); // almin, almax need to be in radians for ADF detector?
    almax.push_back(200.0 * 0.001);

    int ndetect = 1;
    int lpacbed = TRUE;
    mstem.lpacbed = lpacbed;
    if (n_phonons > 1) {
        mstem.lwobble = 1;
    }

    /* init the min/max record of total integrated intensity */
    float **rmin = new2D<float>(nThick, ndetect, "rmin");
    float **rmax = new2D<float>(nThick, ndetect, "rmax");
    
    float*** pixr = new3D<float>(ndetect * nThick, num_probex, num_probey, "pixr");

    float** pacbedPix;
    if (lpacbed == TRUE) {
        pacbedPix = new2D<float>(nx, ny, "pacbedPix");
        for (int ix = 0; ix < nx; ix++) for (int iy = 0; iy < ny; iy++)
            pacbedPix[ix][iy] = 0;
    }
    else pacbedPix = NULL;

    vectord phiMin, phiMax;
    phiMin.resize(ndetect);
    phiMax.resize(ndetect);

    int multiMode = 0;

    int nbeamp, nbeampo;
    float res, thetamax;
    mstem.CountBeams(params, nbeamp, nbeampo, res, thetamax);

    params[pMODE] = mMAGNSTEM;  // save mode = autostem

    fs::path mdir_path(mdir);
    if (!fs::create_directory(mdir_path)) {
        cout << "directory " << mdir_path << " already exists?" << endl;
    }

    opts_4dsave opts4d;
    opts4d.save_cbed = true; // whether to save CBED patterns
    opts4d.save_wvfn = false; // whether to save wavefunctions
    opts4d.save_probe_int_distr = false; // whether to save probe intensity distributions

    // how many pixels in CBED pattern
    opts4d.cbed_nxout = 200;
    opts4d.cbed_nyout = 200;

    // save CBED and wavefunctions in subfolders
    opts4d.cbed_dir = mdir + "/cbed";
    opts4d.wvfn_dir = mdir + "/wvfn";
    if (opts4d.save_cbed) {
        fs::path dir(opts4d.cbed_dir);
        if (!fs::create_directory(dir)) {
            cout << "directory " << dir << " already exists?" << endl;
        }
    }
    if (opts4d.save_wvfn || opts4d.save_probe_int_distr) {
        fs::path dir(opts4d.wvfn_dir);
        if (!fs::create_directory(dir)) {
            cout << "directory " << dir << " already exists?" << endl;
        }
    }

    // run the multislice calculation
    int magn_status = mstem.calculate(params, multiMode, natom, 
        xi, xf, yi, yf, num_probex, num_probey,
        Znum, x, y, z, occ, wobble, mag_phase, iseed,
        ThickSave, nThick,
        almin, almax, collectorMode, ndetect,
        phiMin, phiMax,
        pixr, rmin, rmax,
        pacbedPix, opts4d);

    int nslice = (int)round(cz_f / params[pDELTAZ]) * params[pNCELLREPZ];   // may be off by 1 or 2 with wobble
    long nbeamt = mstem.nbeamt;   //  ??? get beam count - should do this better
    double totmin = mstem.totmin;
    double totmax = mstem.totmax;

    // output images
    double img_dx = (xf - xi) / ((double)(num_probex - 1));  // pixels size for image output
    double img_dy = (yf - yi) / ((double)(num_probey - 1));

    /*  directory file listing parameters for each image file */
    ofstream fp;
    string fileoutpre = "magnstem_test";
    string fileout = (mdir_path / (fileoutpre + ".txt")).make_preferred().string();
    string pacbedFile = (mdir_path / (fileoutpre + "_pacbed.tif")).make_preferred().string();
    fp.open(fileout.c_str());
    if (fp.bad()) {
        cout << "Cannot open output file " << fileout << endl;
        exit(0);
    }
    fp << "C" << endl;
    fp << "C   output of magnetic_test with magnstem " << endl;
    fp << "C" << endl;
    fp << "C   nslice= " << nslice << endl;
    fp << "C deltaz= " << params[pDELTAZ] << ", file in= " << file_coords << endl;
    fp << "C V0= " << params[pENERGY] << ", Cs3= " << params[pCS] << ", Cs5= " << params[pCS5] << ", df= " << params[pDEFOCUS] << endl;
    fp << "C Apert= " << params[pOAPMIN] * 1000.0 << " mrad to " << params[pOAPERT] * 1000.0 << " mrad" << endl;
    fp << "C Crystal tilt x,y= " << params[pXCTILT] << ", " << params[pXCTILT] << endl;

    double pi = 4.0 * atan(1.0);
    for (int idetect = 0; idetect < ndetect; idetect++) {
        if (ADF == collectorMode[idetect])
            fp << "C Detector " << idetect << ", Almin= " << almin[idetect] * 1000.0
            << " mrad, Almax= " << almax[idetect] * 1000.0 << " mrad" << endl;
        else if (CONFOCAL == collectorMode[idetect])
            fp << "C Detector " << idetect << ", cmin= " << almin[idetect]
            << " Angst, cmax= " << almax[idetect] << " Angst." << endl;
        else if (ADF_SEG == collectorMode[idetect])
            fp << "C Detector " << idetect << ", Almin= " << almin[idetect] * 1000.0
            << " mrad, Almax= " << almax[idetect] * 1000.0 << " mrad, "
            << "phimin= " << phiMin[idetect] * 180.0 / pi << ", phimax= "
            << phiMax[idetect] * 180.0 / pi << " deg." << endl;
    }

    fp << "C ncellx= " << ncellx << endl;
    fp << "C ncelly= " << ncelly << endl;
    fp << "C ax= " << ax_f << " A, by= " << by_f << " A, cz= " << cz_f << " A" << endl;
    fp << "C xi= " << xi << endl;
    fp << "C xf= " << xf << endl;
    fp << "C yi= " << yi << endl;
    fp << "C yf= " << yf << endl;
    fp << "C num_probex= " << num_probex << endl;
    fp << "C num_probey= " << num_probey << endl;
    fp << "C ncellrepz= " << ncellrepz << endl;
    fp << "C thickness= " << cz_f * float(ncellrepz) << endl;
    fp << "C nx= " << nx << endl;
    fp << "C ny= " << ny << endl;
    fp << "C Number of symmetrical anti-aliasing "
        << "beams in probe wave function= " << nbeamp << endl;
    fp << "C with a resolution (in Angstroms) = " << res << endl;
    if (mstem.lwobble == 1) {
        fp << "C Number of thermal configurations = " << n_phonons << endl;
        //fp << "C Source size = " << sourceFWHM << " Ang. (FWHM)" << endl;
    }
    //if (TRUE == labErr) {
    //    fp << "C  add pi/4 aberr. tuning errors" << endl;
    //}
    fp << endl;

    /*  store params plus min and max */
    params[pIMAX] = 0.0F;
    params[pIMIN] = 0.0F;
    params[pDX] = (float)dx;
    params[pDY] = (float)dy;
    params[pWAVEL] = (float)wavelength(params[pENERGY]);
    params[pNSLICES] = (float)nslice;
    params[pNXOUT] = (float) num_probex;  // size of output (in pixels)
    params[pNYOUT] = (float) num_probey;

    floatTIFF myFile;
    int NPARAM = myFile.maxParam();
    for (int i = 0; i < NPARAM; i++) myFile.setParam(i, params[i]);
    myFile.setnpix(1);
    myFile.resize(num_probex, num_probey);
    float aimin, aimax;
    aimin = aimax = 0.0F;

    for (int it = 0; it < nThick; it++)
        for (int i = 0; i < ndetect; i++) {
            //sprintf( fileout, "%s%03d%03d.tif", fileoutpre, i, it );   // plain C
            fileout = (mdir_path / (fileoutpre + toString(i) + "_" + toString(it) + ".tif")).make_preferred().string();
            cout << fileout << ": output pix range : " << rmin[it][i] << " to "
                << rmax[it][i] << endl;
            myFile.setParam(pRMAX, rmax[it][i]);
            myFile.setParam(pRMIN, rmin[it][i]);
            myFile.setParam(pMINDET, 0);
            myFile.setParam(pMAXDET, 0);
            myFile.setParam(pPMINDET, 0);
            myFile.setParam(pPMAXDET, 0);
            myFile.setParam(pCMINDET, 0);
            myFile.setParam(pCMAXDET, 0);
            if (ADF == collectorMode[i]) {
                myFile.setParam(pMINDET, (float)(almin[i]));
                myFile.setParam(pMAXDET, (float)(almax[i]));
            }
            else if (CONFOCAL == collectorMode[i]) {
                myFile.setParam(pCMINDET, (float)almin[i]);
                myFile.setParam(pCMAXDET, (float)almax[i]);
            }
            else if (ADF_SEG == collectorMode[i]) {
                myFile.setParam(pMINDET, (float)(almin[i]));
                myFile.setParam(pMAXDET, (float)(almax[i]));
                myFile.setParam(pPMINDET, (float)(phiMin[i]));
                myFile.setParam(pPMAXDET, (float)(phiMax[i]));
            }
            for (int ix = 0; ix < num_probex; ix++) for (int iy = 0; iy < num_probey; iy++)
                myFile(ix, iy) = pixr[i + it * ndetect][ix][iy];
            if (myFile.write(fileout.c_str(), rmin[it][i], rmax[it][i], aimin, aimax,
                (float)dx, (float)dy) != 1) {
                cout << "Cannot write output file " << fileout << endl;
            }

            if (ADF == collectorMode[i])
                fp << "file: " << fileout
                << ", detector= " << almin[i] * 1000.0 << " to " << almax[i] * 1000.0 << " mrad, "
                << "thickness #= " << ThickSave[it] << ", range= " << rmin[it][i]
                << " to " << rmax[it][i] << endl;
            else if (CONFOCAL == collectorMode[i])
                fp << "file: " << fileout
                << ", detector= " << almin[i] << " to " << almax[i] << " Angst., "
                << "thickness #= " << ThickSave[it] << ", range= " << rmin[it][i]
                << " to " << rmax[it][i] << endl;
            else if (ADF_SEG == collectorMode[i])
                fp << "file: " << fileout
                << ", detector= " << almin[i] * 1000.0 << " to " << almax[i] * 1000.0 << " mrad, "
                << "thickness #= " << ThickSave[it] << ", range= " << rmin[it][i]
                << " to " << rmax[it][i] << endl;
        }  /*  end for(i=... */

    fp.close();

    /*   save pos. aver. CBED if needed */
    if (lpacbed == TRUE) {
        myFile.setnpix(1);
        int nxprobe = (int)round(params[pNXPRB]);
        int nyprobe = (int)round(params[pNYPRB]);
        int nx1 = nxprobe / 6;
        int nx2 = (5 * nxprobe) / 6;  /*  cut out center portion without anti-aliasing zeros */
        int ny1 = nyprobe / 6;
        int ny2 = (5 * nyprobe) / 6;
        int nxout2 = nx2 - nx1 + 1;
        int nyout2 = ny2 - ny1 + 1;
        if ((nxout2 < 1) || (nyout2 < 1)) {
            nx1 = ny1 = 0;
            nx2 = nxout2 = nxprobe;
            ny2 = nyout2 = nyprobe;
        }
        myFile.resize(nxout2, nyout2);
        aimin = aimax = 0.0F;
        float scalef = (float)(1.0 / (nxout2 * nyout2));
        int ixo = 0;
        for (int ix2 = nx1; ix2 <= nx2; ix2++) {
            int iyo = 0;
            for (int iy2 = ny1; iy2 <= ny2; iy2++) {
                myFile(ixo, iyo++) = scalef * pacbedPix[ix2][iy2];
            }  ixo++;
        }
        float rmin0 = myFile.min(0);
        float rmax0 = myFile.max(0);
        double dxp = ax_f * ((double)nxprobe) / nx;
        double dyp = by_f * ((double)nyprobe) / ny;
        dxp = 1.0 / dxp;
        dyp = 1.0 / dyp;
        myFile.setParam(pRMAX, rmax0);
        myFile.setParam(pRMIN, rmin0);
        myFile.setParam(pDX, (float)dxp);
        myFile.setParam(pDY, (float)dyp);
        cout << "pos. averg. CBED (unaliased) size " << nxout2 << " x " << nyout2 << " pixels\n"
            << " and range (arb. units): " << rmin0 << " to " << rmax0 << endl;
        if (myFile.write(pacbedFile.c_str(), rmin0, rmax0, aimin, aimax,
            (float)dxp, (float)dyp) != 1) {
            cout << "Cannot write output file " << pacbedFile << endl;
        }
    }   /*  end if( lpacbed.... */

    // delete arrays that were allocated
    delete2D<float>(rmin, nThick);
    delete2D<float>(rmax, nThick);
    delete3D<float>(pixr, ndetect * nThick, num_probex);
    if (lpacbed == TRUE) {
         delete2D<float>(pacbedPix, nx);
    }

    cout << "wall time = " << walltim() - walltimer << " sec." << endl;

    return 0;
}

int main() {
    magnstem_test("./inse7/", 10000ul);

    return 0;
}
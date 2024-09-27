#include <ctime>  /* ANSI C libraries used */

#include <sstream>      // string streams
#include <cstdio>

#include "magnstem.hpp"

#define USE_OPENMP      // define to use openMP
#define MANY_ABERR      //  define to include many aberrations 

#ifdef USE_OPENMP
#include <omp.h>
#endif

magnstem::magnstem()
{

    BW = (2.0F / 3.0F);  // bandwidth limit
    ABERR = 1.0e-4;    // max error for a,b

    NZMAX = 103;    // max atomic number Z

    NRMAX = 100;    // number in look-up-table in vzatomLUT
    RMIN = 0.01;   // r (in Ang) range of LUT for vzatomLUT()
    RMAX = 5.0;
    subpixel_spacing = 0.005; // 0.005, maybe 0.02 for testing
    use_LUT_when_no_phonons = 1;

    pi = 4.0 * atan(1.0);

    // init control flags
    lwobble = lxzimage = lpacbed = 0;
    lrepeatz = 1;
    doConfocal = xFALSE;

    sbuffer = "new magnstem instance";
	messageAST( sbuffer, 0);
	return;
}

magnstem::~magnstem()
{}

/* -------------------  messageAST() -------------------
   message output
   direct all output message here to redirect to the command line
   or a GUI status line or message box when appropriate

   msg[] = character string with message to disply
   level = level of seriousness
        0 = simple status message
        1 = significant warning
        2 = possibly fatal error
*/
void magnstem::messageAST(std::string& smsg, int level)
{
    messageSL(smsg.c_str(), level);  //  just call slicelib version for now

}  // end autostem::messageAST()

/*------------------------ periodic() ---------------------*/
/*
     make probe positions periodic in the supercell
     in case some wobble off the edge with source size of user excess

    pos = input position (x or y);
    size = supercell size ( 0 to size)

    return positive value  0 <= x < size
*/
double magnstem::periodic(double pos, double size)
{
    double x = pos;
    while (x < 0) x += size;
    x = fmod(x, size);
    return(x);
}

/*
  mag_phase = phase shift to be multiplied to the real space wavefunction after every unit cell.
    Dimensions of mag_phase must match that of the transmission function (nx, ny) or mag_phase will be ignored.
  pixr = the images produced by all the detectors
*/
int magnstem::calculate(vectorf& param, int multiMode, int natomin,
    double xi, double xf, double yi, double yf, int nxout, int nyout,
    vectori& Znum, vectorf& xa, vectorf& ya, vectorf& za, vectorf& occ, vectorf& wobble,
    vector<cfpix>& mag_phase, unsigned long iseed,
    vectori& ThickSave, int nThick,
    vectord& almin, vectord& almax, vectori& collectorMode, int ndetect,
    vectord& phiMin, vectord& phiMax,
    float*** pixr, float** rmin, float** rmax,
    float** pacbedPix, opts_4dsave opts4d) {
    int ix, iy, i, idetect, iwobble, nwobble,
        nprobes, ip, it, nbeamp, nbeampo, ix2, iy2;

    float prr, pri, temp, temperature;

    double scale, sum, wx, w,
        tctx, tcty, dx, dy, ctiltx, ctilty, k2maxa, k2maxb, k2;

    vectord x, y, sums;

    // ---- get setup parameters from param[]
    ax = param[pAX];
    by = param[pBY];
    cz = param[pCZ];
    nx = ToInt(param[pNX]);         // size of transmission function (in pixels)
    ny = ToInt(param[pNY]);
    keV = param[pENERGY];             // electron beam energy in keV
    df = param[pDEFOCUS];               // defocus
    ctiltx = param[pXCTILT];          // crystal tilt
    ctilty = param[pYCTILT];
    apert1 = param[pOAPMIN];          // obj. apert size in radians (min, max for annular apert)
    apert2 = param[pOAPERT];
    temperature = param[pTEMPER];             // temperature
    nwobble = ToInt(param[pNWOBBLE]);       // number config. to average
    deltaz = param[pDELTAZ];                  // slice thickness

    nxprobe = ToInt(param[pNXPRB]);  // probe size in pixels
    nyprobe = ToInt(param[pNYPRB]);

    //  confocal collector lens parameters
    dfC = param[pCDF];                // collector defocus
    dfa2C = param[pCDFA2];            // collector astig 2nd order
    dfa2phiC = param[pCDFA2PHI];
    dfa3C = param[pCDFA3];            // collector astig 2nd order
    dfa3phiC = param[pCDFA3PHI];
    Cs3C = param[pCCS3];              // collector spherical aberr.
    Cs5C = param[pCCS5];
    apert1C = param[pCCAPMIN];        //  collector apert. in radians
    apert2C = param[pCCAPMAX];

    natom = natomin;
    wavlen = wavelength(keV);

    //  check specified nxout, nyout, deltaz
    // nxout, nyout are number of probe positions
    nprobes = (lwobble == 0) ? nyout : nwobble;
    if (nxout < 1) {
        sbuffer = "nxout must be > 1 in autostem but it is " + toString(nxout);
        messageAST(sbuffer, 2);
        return(-1);
    }
    if (nyout < 1) {
        sbuffer = "nyout must be > 1 in autostem but it is " + toString(nyout);
        messageAST(sbuffer, 2);
        return(-2);
    }

    if (deltaz < 0.1F) {
        sbuffer = "delta z is too small; it is " + toString(deltaz);
        messageAST(sbuffer, 2);
        return(-3);
    }

    // initialize iseeds (use a separate seed for each probe for reproducibility, since the
    // probes may be propagated in parallel)
    vector<unsigned long> iseeds(nprobes);
    for (int i = 0; i < nprobes; i++) {
        iseeds[i] = iseed + i;
    }

    //  see if confocal needed
    doConfocal = xFALSE;
    for (i = 0; i < ndetect; i++)
        if (CONFOCAL == collectorMode[i]) doConfocal = xTRUE;

    vectori Znum_unique;
    for (vectori::iterator it = begin(Znum); it != end(Znum); ++it) {
        if (find(Znum_unique.begin(), Znum_unique.end(), *it) == Znum_unique.end()) {
            Znum_unique.push_back(*it);
        }
    }
    cout << "Structure has Z values: ";
    for (auto& zi : Znum_unique) {
        cout << zi << " ";
    }
    cout << endl;

#ifdef USE_OPENMP
    /*  force LUT init. to avoid redundant init in parallel form */
    double rsq = 0.5;  /* arbitrary position */
    double vz;
    for (auto& zi : Znum_unique) {
        vz = vzatomLUT(zi, rsq);
        vz = vzatom_pxavg(zi, rsq, rsq, rsq, rsq);
    }
#endif

    // structure
    /* to add random offsets */
    xa2.resize(natom);
    ya2.resize(natom);
    za2.resize(natom);
    Znum2.resize(natom);
    occ2.resize(natom);

    if (lrepeatz == 1) {
        cout << cz / deltaz << endl;
        slice_per_cell = (int) round(cz / deltaz);
        deltaz = cz / (double)slice_per_cell;
        sbuffer = "Using " + toString(slice_per_cell) + " slice(s) per cell, deltaz = " + toString(deltaz) + ".";
        messageAST(sbuffer, 0);

        if (0 != lverbose) {
            sbuffer = "Sorting atoms by depth...";
            messageAST(sbuffer, 0);
        }
        sortByZ(xa, ya, za, occ, Znum, wobble, natom);

        if (za[0] < 0.0 || za[natom - 1] >= cz) {
            sbuffer = "Atom depth out of bounds. Need to satisfy 0 <= z < cz.";
            messageAST(sbuffer, 2);
            return -1;
        }

        // compute indices of each slice
        slice_start_idx.assign(slice_per_cell, -1);
        slice_end_idx.assign(slice_per_cell, -1);
        int prevIdx = 0;
        int currIdx;
        for (int aidx = 0; aidx < natom; aidx++) {
            //cout << za[aidx] / deltaz << " ";
            currIdx = (int)floor(za[aidx] / deltaz + 0.001); // this is a terrible hack

            if (currIdx < prevIdx) {
                sbuffer = "Atoms out of order";
                messageAST(sbuffer, 2);
                return -1;
            }
            else if (currIdx > prevIdx) {
                slice_start_idx[currIdx] = aidx;
                slice_end_idx[prevIdx] = aidx - 1;
                if (aidx == 0) {
                    slice_start_idx[0] = -1;
                }
                for (int midx = prevIdx + 1; midx < currIdx; midx++) {
                    slice_start_idx[midx] = -1;
                    slice_end_idx[midx] = -1;
                }
                prevIdx = currIdx;
            }
            else if (aidx == 0) {
                slice_start_idx[0] = 0;
            }
            if (aidx == natom - 1) {
                slice_end_idx[currIdx] = aidx;
            }
        }
        for (int sidx = currIdx + 1; sidx < slice_per_cell; sidx++) {
            slice_start_idx[sidx] = -1;
            slice_end_idx[sidx] = -1;
        }
        sbuffer = "start indices: ";
        for (vectori::const_iterator i = slice_start_idx.begin(); i != slice_start_idx.end(); ++i) {
            sbuffer += toString(*i) + ", ";
        }
        sbuffer += "; end indices: ";
        for (vectori::const_iterator i = slice_end_idx.begin(); i != slice_end_idx.end(); ++i) {
            sbuffer += toString(*i) + ", ";
        }
        messageAST(sbuffer, 0);

    }
    else { // not really implemented
        slice_per_cell = 1;
        if (lwobble == 0) { // if lwobble == 1, then first add wobble, then sort atoms
            if (0 != lverbose) {
                sbuffer = "Sorting atoms by depth...";
                messageAST(sbuffer, 0);
            }
            sortByZ(xa, ya, za, occ, Znum, natom);
        }
    }


    /*  check that requested probe size is not bigger
    than transmission function size (or too small)
    */
    if ((nxprobe > nx) || (nxprobe < 2)) {
        nxprobe = nx;
        sbuffer = "Probe size reset to nx= " + toString(nxprobe);
        messageAST(sbuffer, 0);
    }

    if ((nyprobe > ny) || (nyprobe < 2)) {
        nyprobe = ny;
        sbuffer = "probe size reset to ny= " + toString(nyprobe);
        messageAST(sbuffer, 0);
    }

    /*  calculate spatial frequencies for future use
        (one set for transmission function and one for probe
        wavefunction)
    NOTE: zero freg is in the bottom left corner and
        expands into all other corners - not in the center
        this is required for FFT - don't waste time rearranging

    remember : the x values are the same for both sets

    x2, y2 are used for confocal

    */

    kx.resize(nx);
    ky.resize(ny);
    kx2.resize(nx);
    ky2.resize(ny);
    xp.resize(nx);
    yp.resize(ny);

    freqn(kx, kx2, xp, nx, ax);
    freqn(ky, ky2, yp, ny, by);

    kxp.resize(nxprobe);
    kyp.resize(nyprobe);
    kxp2.resize(nxprobe);
    kyp2.resize(nyprobe);

    freqn(kxp, kxp2, xp, nxprobe, ax * ((double)nxprobe) / nx);
    freqn(kyp, kyp2, yp, nyprobe, by * ((double)nyprobe) / ny);

    //for (vectorf::const_iterator i = kxp.begin(); i != kxp.end(); ++i) {
    //    cout << *i << ' ';
    //}
    //cout << endl;

    // impose anti-aliasing bandwidth limit on transmission functions

    sum = ((double)nx) / (2.0 * ax);
    k2maxp = ((double)ny) / (2.0 * by);
    if (sum < k2maxp) k2maxp = sum;
    k2maxt = k2maxp;  //  for trans. function potential in cuAtompot()
    k2maxp = BW * k2maxp;
    //printf("Bandwidth limited to a real space resolution of %f Angstroms\n",   //???
    //                 1.0F/k2maxp);
    //printf("   (= %.2f mrad)  for symmetrical anti-aliasing.\n",
    //     wavlen*k2maxp*1000.0F);
    k2maxp = k2maxp * k2maxp;  // for probe 
    k2maxt = k2maxt * k2maxt;  // for trans function potential

    /*  initialize propagator */
    cprop.resize(nxprobe, nyprobe);
    tctx = 2.0 * tan(ctiltx);
    tcty = 2.0 * tan(ctilty);
    scale = pi * deltaz;
    for (ix = 0; ix < nxprobe; ix++) {
        wx = (kxp2[ix] * wavlen - kxp[ix] * tctx);
        for (iy = 0; iy < nyprobe; iy++) {
            if (((double) kxp2[ix] + (double) kyp2[iy]) < k2maxp) {
                w = scale * (wx + kyp2[iy] * wavlen - kyp[iy] * tcty);
                cprop.re(ix, iy) = (float)cos(w);
                cprop.im(ix, iy) = (float)-sin(w);
            }
            else {
                cprop.re(ix, iy) = 0.0F;
                cprop.im(ix, iy) = 0.0F;
            }  //  end if( kx2[ix]... 
        } // end for(iy..) 
    } // end for(ix..) 

    /*   calculate number of pixels in the probe and obj. apert. */
    k2maxa = apert1 / wavlen;
    k2maxa = k2maxa * k2maxa;
    k2maxb = apert2 / wavlen;
    k2maxb = k2maxb * k2maxb;

    nbeamp = nbeampo = 0;
    for (iy = 0; iy < nyprobe; iy++)
        for (ix = 0; ix < nxprobe; ix++) {
            k2 = (double) kyp2[iy] + (double) kxp2[ix];
            if (k2 < k2maxp) nbeamp++;
            if ((k2 >= k2maxa) && (k2 <= k2maxb)) nbeampo++;
        }

    //  output this in outer calling program if wanted - not every time here
    sbuffer = "Number of symmetrical anti-aliasing beams in probe= "
           + toString(nbeamp);
    messageAST( sbuffer, 0 );
    sbuffer = "Number of beams in probe aperture= " + toString(nbeampo);
    messageAST( sbuffer, 0 );

    if (nbeamp < 200) {
        sbuffer = "WARNING: the probe is under sampled, this is a bad idea...";
        messageAST(sbuffer, 1);
    }
    if (nbeampo < 100) {
        sbuffer = "WARNING: the probe aperture is under sampled, this is a bad idea...";
        messageAST(sbuffer, 2);
        //exit(EXIT_FAILURE);
    }

    /*  convert aperture dimensions */

    k2min.resize(ndetect);
    k2max.resize(ndetect);

    for (idetect = 0; idetect < ndetect; idetect++) {
        if ((ADF == collectorMode[idetect]) || (ADF_SEG == collectorMode[idetect])) {
            k2max[idetect] = almax[idetect] / wavlen;
            k2max[idetect] = k2max[idetect] * k2max[idetect];
            k2min[idetect] = almin[idetect] / wavlen;
            k2min[idetect] = k2min[idetect] * k2min[idetect];
        }
        else if (CONFOCAL == collectorMode[idetect]) {
            k2max[idetect] = almax[idetect] * almax[idetect];
            k2min[idetect] = almin[idetect] * almin[idetect];
        }
    }

    /*  init the min/max record of total integrated intensity */

    totmin = 10.0;
    totmax = -10.0;
    sums.resize(nprobes);

    // allocate probe wave function and transmission function arrays
    probe = new cfpix[nprobes];
    if (NULL == probe) {
        sbuffer = "Cannot allocate probe array";
        messageAST(sbuffer, 2);
        exit(EXIT_FAILURE);
    }
    ix = probe[0].resize(nxprobe, nyprobe);
    if (ix < 0) {
        sbuffer = "Cannot allocate probe array storage";
        messageAST(sbuffer, 2);
        exit(EXIT_FAILURE);
    }
    probe[0].init();
    if (nprobes > 1) for (ip = 1; ip < nprobes; ip++) {
        ix = probe[ip].resize(nxprobe, nyprobe);
        if (ix < 0) {
            sbuffer = "Cannot allocate probe array storage";
            messageAST(sbuffer, 2);
            exit(EXIT_FAILURE);    //  should do something better here ??
        }
        probe[ip].copyInit(probe[0]);
    }

    int num_trans = (lwobble == 0) ? slice_per_cell : nwobble;
    trans.reserve(num_trans);
    trans.push_back(cfpix());
    ix = trans[0].resize(nx, ny);
    if (ix < 0) {
        sbuffer = "Cannot allocate transmission function storage";
        messageAST(sbuffer, 2);
        exit(EXIT_FAILURE);
    }
    trans[0].init();
    if (num_trans > 1) for (ip = 1; ip < num_trans; ip++) {
        trans.push_back(cfpix());
        ix = trans[ip].resize(nx, ny);
        if (ix < 0) {
            sbuffer = "Cannot allocate transmission function array storage";
            messageAST(sbuffer, 2);
            exit(EXIT_FAILURE);    //  should do something better here ??
        }
        trans[ip].copyInit(trans[0]);
    }

    // allocate pacbed image
    if (lpacbed == xTRUE) {
        for (ix = 0; ix < nxprobe; ix++) for (iy = 0; iy < nyprobe; iy++)
            pacbedPix[ix][iy] = 0;
    }

    /* ------------- start here for a full image output -------------- */
    /*
      do one whole line at once NOT the whole image (which may be huge)
    */

    sbuffer = "output file size in pixels is " + toString(nxout) + " x "
        + toString(nyout);
    messageAST(sbuffer, 0);
    if (lwobble == 0 && nprobes != nyout) {
        sbuffer = "Error, nprobes= " + toString(nprobes) + " must be the same as nyout= "
            + toString(nyout) + ", in image mode.";
        messageAST(sbuffer, 2);
        exit(0);
    }
    /* double up first index to mimic a 4D array */
    for (i = 0; i < (nThick * ndetect); i++) {
        for (ix = 0; ix < nxout; ix++)
            for (iy = 0; iy < nyout; iy++)
                pixr[i][ix][iy] = 0.0F;
    }

    /*  iterate the multislice algorithm proper for each
        position of the focused probe */

    if (nxout > 1) dx = (xf - xi) / ((double)(nxout - 1));
    else dx = 1.0;
    if (nyout > 1) dy = (yf - yi) / ((double)(nyout - 1));
    else dy = 1.0;
    x.resize(nprobes);
    y.resize(nprobes);

    if (lrepeatz == 1) {
        if (lwobble == 0) {
            opts4d.pos_ix.resize(nprobes);
            opts4d.pos_iy.resize(nprobes);
            double*** detect = new3D<double>(nThick, ndetect, nprobes,
                "detect");

            nwobble = 1;

            // calculate transmission functions
            if (use_LUT_when_no_phonons == 1) {
                if (!vz_maker.LUT_allocated) {
                    vz_maker.allocate_LUT(subpixel_spacing, Znum_unique, ax / nx, by / ny);
                }
                if (!vz_maker.LUT_calculated) {
                    vz_maker.calculate_LUT();
                }
            }

            for (it = 0; it < slice_per_cell; it++) {
                if (slice_start_idx[it] >= 0) {
                    double phirms;
                    int natom = slice_end_idx[it] - slice_start_idx[it] + 1;

                    if (use_LUT_when_no_phonons == 1) {
                        trlayer_pxavg_LUT(xa, ya, occ,
                            Znum, natom, slice_start_idx[it], (float)ax, (float)by, (float)keV,
                            trans[it], mag_phase[it], (long)nx, (long)ny, &phirms, &nbeamt, (float)k2maxp, vz_maker);
                            //it, param, opts4d);
                    }
                    else {
                        trlayer_pxavg(xa, ya, occ,
                            Znum, natom, slice_start_idx[it], (float)ax, (float)by, (float)keV,
                            trans[it], mag_phase[it], (long)nx, (long)ny, &phirms, &nbeamt, (float)k2maxp,
                            it, param, opts4d);
                    }
                    sbuffer = "calculated transmission function " + toString(it) + ", " + toString(natom) + " atoms.";
                    messageAST(sbuffer, 0);
                }
                else {
                    trans[it] = 1.0F;
                    sbuffer = "warning: no atoms in slice " + toString(it);
                    messageAST(sbuffer, 1);
                }
            }

            //exit(0);

            for (ix = 0; ix < nxout; ix++) {

                for (iy = 0; iy < nyout; iy++) {
                    x[iy] = xi + dx * ((double)ix);
                    //  + sourcesize * rangauss(iseed);  - does not converge well
                    y[iy] = yi + dy * ((double)iy);
                    //  + sourcesize * rangauss(iseed);  - does not converge well
                    x[iy] = periodic(x[iy], ax);   /* put back in supercell */
                    y[iy] = periodic(y[iy], by);   /* if necessary */
                    opts4d.pos_ix[iy] = ix;
                    opts4d.pos_iy[iy] = iy;
                }

                sbuffer = "calculating line " + toString(ix);
                messageAST(sbuffer, 0);

                STEMsignals(x, y, nyout, param,
                    trans,
                    multiMode, detect, ndetect,
                    ThickSave, nThick, sums, collectorMode,
                    phiMin, phiMax, opts4d);

                for (iy = 0; iy < nyout; iy++) {
                    if (sums[iy] < totmin) totmin = sums[iy];
                    if (sums[iy] > totmax) totmax = sums[iy];
                    for (it = 0; it < nThick; it++) {
                        for (idetect = 0; idetect < ndetect; idetect++)
                            pixr[idetect + it * ndetect][ix][iy] += (float)
                            (detect[it][idetect][iy] / ((double)nwobble));
                    }
                    if (sums[iy] < 0.9) {
                        sbuffer = "Warning integrated intensity too small, = "
                            + toString(sums[iy]) + " at " + toString(x[iy]) + ", " + toString(y[iy]);
                        messageAST(sbuffer, 0);
                    }
                    if (sums[iy] > 1.1) {
                        sbuffer = "Warning integrated intensity too large, = "
                            + toString(sums[iy]) + " at " + toString(x[iy]) + ", " + toString(y[iy]);
                        messageAST(sbuffer, 0);
                    }
                }

                /*   sum position averaged CBED if requested
                     - assume probe still left from stemsignal()  */
                if (lpacbed == xTRUE) {
                    for (iy = 0; iy < nyout; iy++) {
                        for (ix2 = 0; ix2 < nxprobe; ix2++)
                            for (iy2 = 0; iy2 < nyprobe; iy2++) {
                                prr = probe[iy].re(ix2, iy2);
                                pri = probe[iy].im(ix2, iy2);
                                pacbedPix[ix2][iy2] += (prr * prr + pri * pri);
                            }
                    }
                }
            }

            if (lpacbed == xTRUE) {
                invert2D(pacbedPix, nxprobe, nyprobe);  /*  put zero in middle */
            }
            delete3D<double>(detect, nThick, ndetect);
        }
        else {  // Handle the case of phonons
            opts4d.pos_ix.resize(1);
            opts4d.pos_iy.resize(1);
            double** detect = new2D<double>(nThick, ndetect, "detect");
            double xp, yp;  // probe position
            double sum = 0;

            float** ca_cbed = new2D<float>(nxprobe, nyprobe, "ca_cbed"); // CBED averaged over multiple configurations at one scan position

            if (!vz_maker.LUT_allocated) {
                vz_maker.allocate_LUT(subpixel_spacing, Znum_unique, ax/nx, by/ny);
            }
            if (!vz_maker.LUT_calculated) {
                vz_maker.calculate_LUT();
            }
            
            for (ix = 0; ix < nxout; ix++) {
                opts4d.pos_ix[0] = ix;
                xp = periodic(xi + dx * ((double)ix), ax);

                for (iy = 0; iy < nyout; iy++) {
                    sbuffer = "calculating line " + toString(ix) + ", column " + toString(iy);
                    messageAST(sbuffer, 0);
                    opts4d.pos_iy[0] = iy;
                    yp = periodic(yi + dy * ((double)iy), by);

                    // calculate for a single probe position
                    STEMsignals_wobble(xp, yp, ix, iy, nwobble,
                        param, multiMode,
                        Znum, xa, ya, za, occ, wobble,
                        iseeds, mag_phase,
                        ThickSave, ca_cbed, 
                        detect, ndetect, sum, collectorMode, vz_maker,
                        phiMin, phiMax, opts4d);

                    if (sum < totmin) totmin = sum;
                    if (sum > totmax) totmax = sum;
                    for (it = 0; it < nThick; it++) {
                        for (idetect = 0; idetect < ndetect; idetect++)
                            pixr[idetect + it * ndetect][ix][iy] += (float)detect[it][idetect];
                    }

                    if (sum < 0.9) {
                        sbuffer = "Warning integrated intensity too small, = "
                            + toString(sum) + " at " + toString(xp) + ", " + toString(yp);
                        messageAST(sbuffer, 0);
                    }
                    if (sum > 1.1) {
                        sbuffer = "Warning integrated intensity too large, = "
                            + toString(sum) + " at " + toString(xp) + ", " + toString(yp);
                        messageAST(sbuffer, 0);
                    }

                    if (lpacbed == xTRUE) {
                        for (ix2 = 0; ix2 < nxprobe; ix2++)
                            for (iy2 = 0; iy2 < nyprobe; iy2++) {
                                pacbedPix[ix2][iy2] += ca_cbed[ix2][iy2];
                            }
                    }
                    // no need to shift pacbedPix; ca_cbed was already shifted
                }
            }

            delete2D<double>(detect, nThick);
            delete2D<float>(ca_cbed, nxprobe);
        }

        /*  find range to output data files  */
        for (it = 0; it < nThick; it++)
            for (i = 0; i < ndetect; i++) {
                rmin[it][i] = rmax[it][i] = pixr[i + it * ndetect][0][0];
                for (ix = 0; ix < nxout; ix++)
                    for (iy = 0; iy < nyout; iy++) {
                        temp = pixr[i + it * ndetect][ix][iy];
                        if (temp < rmin[it][i])rmin[it][i] = (float)temp;
                        if (temp > rmax[it][i])rmax[it][i] = (float)temp;
                    }
            }
    }
    else {
        sbuffer = "Non-perfect crystal mode is not implemented.";
        messageAST(sbuffer, 2);
    }

    if (nbeampo < 100) {
        sbuffer = "WARNING: the probe aperture was under sampled, this was a bad idea...";
        messageAST(sbuffer, 2);
        //exit(EXIT_FAILURE);
    }

    delete[] probe;

    return 0;
}

/*------------------------ STEMsignals() ---------------------*/
/*

  adapted from STEMsignals() in autostem.cpp

  subroutine to calculate the stem signal arising from a given
  probe position

  iterate the multislice algorithm proper for each position of
  the focused probe

  This version uses massive amounts of memory to avoid
  recalculating the transmission functions more than necessary

   note zero freg is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for fft - don't waste time rearranging

     change to propagate thru whole unit cell not just
     range of atoms on 14-sep-2006 ejk
     add multipole aberrations 9-may-2011 ejk
     change to cfpix for probe and trans 10-nov-2012 ejk

  x[],y[]     = real positions of the incident probe
  npos        = int number of positions
  param[]     = parameters of probe
  multiMode   = flag to add multipole aberrations
  detect[][][]= real array to get signal into each detector
            for each probe position and thickness
  ndetect     = number of detector geometries
  ThickSave[] = thicknesses at which to save data (other than the last)
  nThick      = number of thickness levels (including the last)
  sum         = real total integrated intensity

  modifications: (nov-2020)
  * use precomputed transmission functions
  * save wavefunctions or CBED patterns
  * apply an extra phase shift

  the assumed global variables are:

  nxprobe,nyprobe   = int size of probe wavefunction in pixels
  nx,ny         = int size of transmission function in pixels
  trans                  = float complex transmission function
  propxr[][], propxi[][] = float real,imag propagator vs x
  propyr[][], propyi[][] = float real,imag propagator vs y
  ax,by,cz      = float unit cell size in Angs
  kxp[], kyp[]      = float spatial frequencies vs x, y
  kxp2[], kyp2[]    = float square of kxp[], kyp[]
  xp[], yp[]        = float real space positions in probe (confocal)
  apert1, apert2    = double min,max objective aperture (in radians)
  k2maxp            = double max spatial freq of probe squared
  pi                = double constant PI
  wavlen            = double electron wavelength in Angs
  df                = double defocus (in Ang)
  Cs3,Cs5           = double spherical aberration (in Ang)

  xa[],ya[],za[]    = atom coordinates
  occ[]         = atomic occupancy
  Znum[]        = atomic numbers
  natom         = number of atoms
  deltaz        = slice thickness
  v0            = beam energy
  nbeamt        = number of beams in transmission function
  zmin, zmax    = range of z coord. of the atoms
  nslice        = number of slices
  doConfocal    = flag indicating confocal is needed
  k2min[], k2max[] = detector angles in 1/Ang

    NOTE:  too many thing come in as globals, but...

*/
void magnstem::STEMsignals(vectord& x, vectord& y, int npos, vectorf& p,
    vector<cfpix>& trans, 
    int multiMode, double*** detect, int ndetect,
    vectori& ThickSave, int nThick, vectord& sum, vectori& collectorMode,
    vectord& phiMin, vectord& phiMax, opts_4dsave opts4d)
{
    int ix, iy, ixt, iyt, idetect, ixmid, iymid;
    int istart, na, ip, i, it;

    long nxl, nyl;

    float scale, prr, pri, tr, ti;

    double  chi0, chi1, k2maxa, k2maxb,
        w, k2, phi, phirms, alx, aly;
    double sum0, sum1, delta, zslice, totalz;

    vectori ixoff, iyoff;
    vectord xoff, yoff;

    /* extra for confocal */
    float hr, hi;
    double chi2C, chi3C, k2maxaC, k2maxbC, r2, rx2;
    cfpix cpix;            /* complex confocal image */
    vector<unique_ptr<floatTIFF>> distr_xz;
    vector<unique_ptr<floatTIFF>> distr_yz;

    /* ------ make sure x,y are ok ------ */

    for (ip = 0; ip < npos; ip++) {
        if ((x[ip] < 0.0) || (x[ip] > ax) ||
            (y[ip] < 0.0) || (y[ip] > by)) {
            sum[ip] = -1.2345;
            sbuffer = "bad x,y in STEMsignals() = \n" + toString(x[ip]) + toString(y[ip]);
            messageAST(sbuffer, 0);
            return;
        }

        if (opts4d.save_probe_int_distr) {
            distr_xz.push_back(make_unique<floatTIFF>());
            distr_yz.push_back(make_unique<floatTIFF>());
            int NPARAM = distr_xz[ip]->maxParam();
            for (i = 0; i < NPARAM; i++) {
                distr_xz[ip]->setParam(i, p[i]);
                distr_yz[ip]->setParam(i, p[i]);
            }
            distr_xz[ip]->setnpix(1);
            distr_xz[ip]->resize(nxprobe, p[pNCELLREPZ] + 1);
            distr_yz[ip]->setnpix(1);
            distr_yz[ip]->resize(nyprobe, p[pNCELLREPZ] + 1);
        }
    }

    /*  generate a probe at position x,y which must be inside the
        0->ax and 0->by region
        normalize to have unit integrated intensity
        probe position (x,y) = (0,0) = lower left corner

    NOTE: The probe wavefunction is nxprobe x nyprobe which may
        be smaller than the transmission function.
        Position the probe near the center of the nxprobe x nyprobe
        region because it does not wrap around properly
        (i.e. manually wrap around in the nx x ny region if
         the probe is near a boundary of the nx x ny region))
    */

    ixmid = nxprobe / 2;
    iymid = nyprobe / 2;
    chi1 = pi * wavlen;
    k2maxa = apert1 / wavlen;
    k2maxa = k2maxa * k2maxa;
    k2maxb = apert2 / wavlen;
    k2maxb = k2maxb * k2maxb;

    /* extra for confocal */
    chi2C = 0.5 * Cs3C * wavlen * wavlen;
    chi3C = Cs5C * wavlen * wavlen * wavlen * wavlen / 3.0;
    k2maxaC = apert1C / wavlen;
    k2maxaC = k2maxaC * k2maxaC;
    k2maxbC = apert2C / wavlen;
    k2maxbC = k2maxbC * k2maxbC;

    ixoff.resize(npos);
    iyoff.resize(npos);
    xoff.resize(npos);
    yoff.resize(npos);

    /* ------- calculate all of the probe wave functions at once ------
        to reuse the transmission functions which takes a long
        time to calculate*/

        /*  paralleling this loop has little effect */
#pragma omp parallel for private(ix,iy,sum0,k2,w,chi0,scale,tr,ti,alx,aly) 
    for (ip = 0; ip < npos; ip++) {
        ixoff[ip] = (int)floor(x[ip] * ((double)nx) / ax) - ixmid;
        xoff[ip] = x[ip] - ax * ((double)ixoff[ip]) / ((double)nx);

        iyoff[ip] = (int)floor(y[ip] * ((double)ny) / by) - iymid;
        yoff[ip] = y[ip] - by * ((double)iyoff[ip]) / ((double)ny);

        sum0 = 0.0;
        for (ix = 0; ix < nxprobe; ix++) {
            alx = wavlen * kxp[ix];  /* x component of angle alpha */
            for (iy = 0; iy < nyprobe; iy++) {
                aly = wavlen * kyp[iy];  /* y component of angle alpha */
                k2 = kxp2[ix] + kyp2[iy];
                if ((k2 >= k2maxa) && (k2 <= k2maxb)) {
                    w = 2. * pi * (xoff[ip] * kxp[ix] + yoff[ip] * kyp[iy]);
                    chi0 = (2.0 * pi / wavlen) * chi(p, alx, aly, multiMode);
                    chi0 = -chi0 + w;
                    probe[ip].re(ix, iy) = tr = (float)cos(chi0);
                    probe[ip].im(ix, iy) = ti = (float)sin(chi0);
                    sum0 += (double)(tr * tr + ti * ti);
                }
                else {
                    probe[ip].re(ix, iy) = 0.0F;
                    probe[ip].im(ix, iy) = 0.0F;
                }
            }
        }  /* end for( ix... */

        scale = (float)(1.0 / sqrt(sum0));
        probe[ip] *= scale;

    }  /* end for( ip...) */

    /* -------- transmit thru nslice layers ------------------------
        beware ixoff,iyoff must be with one nx,ny
            of the array bounds
    */


    nxl = (long)nx;
    nyl = (long)ny;

    scale = 1.0F / (((float)nx) * ((float)ny));

    zslice = 0.75 * deltaz;  /*  start a little before top of unit cell */
    istart = 0;
    nslice = 0;
    it = 0;     /* thickness level index */

    bool save_slice_wvfn = false;

    for (nslice = 0; nslice < p[pNCELLREPZ]; nslice++) {
        save_slice_wvfn = opts4d.save_wvfn && find(ThickSave.begin(), ThickSave.end(), nslice - 1) != ThickSave.end();

        //sbuffer = "starting nslice " + toString(nslice); 
        //messageAST(sbuffer, 0);

        /*----- multislice trans/prop through one unit cell for all probes ---- */
#pragma omp parallel for private(ix,iy,ixt,iyt,prr,pri)
        for (ip = 0; ip < npos; ip++) {
            /* range of unit cell */
            for (int itr = 0; itr < slice_per_cell; itr++) {
                probe[ip].ifft();
                if (itr == 0) {
                    if (save_slice_wvfn) {
                        fs::path w_path = (fs::path(opts4d.wvfn_dir) / ("layer_" + to_string(nslice - 1) + "_wvfn_ix_" + to_string(opts4d.pos_ix[ip])
                            + "_iy_" + to_string(opts4d.pos_iy[ip]) + ".tif")).make_preferred();
                        //cout << "need to save " << w_path.string() << endl;
                        write_wvfn(probe[ip], w_path.string(), p);
                    }
                    if (opts4d.save_probe_int_distr) {
                        record_probe_intensity_distributions(probe[ip], *distr_xz[ip], *distr_yz[ip], nslice);
                    }
                }

                /* apply transmission function without checking if there are atoms in this slice */
                // the probe is not shifted; the transmission function is shifted
                for (ix = 0; ix < nxprobe; ix++) {
                    ixt = ix + ixoff[ip];
                    if (ixt >= nx) ixt = ixt - nx;
                    else if (ixt < 0) ixt = ixt + nx;
                    for (iy = 0; iy < nyprobe; iy++) {
                        iyt = iy + iyoff[ip];
                        if (iyt >= ny) iyt = iyt - ny;
                        else if (iyt < 0) iyt = iyt + ny;
                        prr = probe[ip].re(ix, iy);
                        pri = probe[ip].im(ix, iy);
                        probe[ip].re(ix, iy) = prr * trans[itr].re(ixt, iyt)
                            - pri * trans[itr].im(ixt, iyt);  // real
                        probe[ip].im(ix, iy) = prr * trans[itr].im(ixt, iyt)
                            + pri * trans[itr].re(ixt, iyt);  // imag
                    } /* end for(iy...) */
                }  /* end for(ix...) */

                probe[ip].fft();

                /*  multiplied by the propagator function */
                probe[ip] *= cprop;

            }  /* end  for( itr... */

            if (ip == 0) {
                float sum_probe0 = 0;
                for (ix = 0; ix < nxprobe; ix++) for (iy = 0; iy < nyprobe; iy++) {
                    prr = probe[0].re(ix, iy);
                    pri = probe[0].im(ix, iy);
                    sum_probe0 += prr * prr + pri * pri;
                }
                printf("%d %.4f, ", nslice, sum_probe0);
            }
        }  /* end for( ip...) */

           /*  if this is a good thickness level then save the ADF or confocal signals
               - remember that the last level may be off by one layer with
                  thermal displacements so special case it
            */

            /*  look at all values because they may not be in order */
        for (it = 0; it < nThick; it++)
            if (ThickSave[it] == nslice) {

                if (0 != lverbose) {
                    sbuffer = "save ADF/confocal signals, thickness level " + toString(ThickSave[it]);  // diagnostic
                    messageAST(sbuffer, 0);
                }

                /*  loop over all probes again */
#pragma omp parallel for private(ix,iy,idetect,prr,pri,delta,k2,cpix,phi,chi0,hr,hi,sum0,sum1,rx2,r2)
                for (ip = 0; ip < npos; ip++) {

                    /*  zero sum count */
                    sum[ip] = 0.0;
                    for (ix = 0; ix < ndetect; ix++) detect[it][ix][ip] = 0.0;

                    /*  sum intensity incident on the ADF detector
                            and calculate total integrated intensity
                        - changed detector limits to >= min and < max
                           so many concentric ADF detectors sum correctly
                           7-jul-2011
                    */
                    for (ix = 0; ix < nxprobe; ix++) {
                        for (iy = 0; iy < nyprobe; iy++) {

                            prr = probe[ip].re(ix, iy);
                            pri = probe[ip].im(ix, iy);
                            delta = prr * prr + pri * pri;
                            sum[ip] += delta;

                            k2 = kxp2[ix] + kyp2[iy];
                            phi = atan2(kyp[iy], kxp[ix]);  //  for ADF_SEG detector

                            for (idetect = 0; idetect < ndetect; idetect++) {
                                if (ADF == collectorMode[idetect]) {
                                    if ((k2 >= k2min[idetect]) &&
                                        (k2 < k2max[idetect]))
                                        detect[it][idetect][ip] += delta;
                                }
                                else if (ADF_SEG == collectorMode[idetect]) {
                                    if ((k2 >= k2min[idetect]) && (k2 < k2max[idetect])
                                        && (phi >= phiMin[idetect]) && (phi < phiMax[idetect]))
                                        detect[it][idetect][ip] += delta;
                                }
                            }
                        } /* end for(iy..) */
                    }  /* end for(ix...) */

                    if (opts4d.save_cbed) {
                        fs::path c_path = (fs::path(opts4d.cbed_dir) / ("layer_" + to_string(nslice) + "_cbed_ix_" + to_string(opts4d.pos_ix[ip]) 
                            + "_iy_" + to_string(opts4d.pos_iy[ip]) + ".tif")).make_preferred();
                        //cout << "need to save " << c_path.string() << endl;
                        write_cbed_from_probe(probe[ip], c_path.string(), p, opts4d);
                    }

                    /*  transform back if confocal needed
                        - use copy of probe so original can continue in use  */
                    if (doConfocal == xTRUE) {
                        /*  allocate/deallocate here so openMP will work
                            otherwise have to allocate nxprobe cpix arrays
                            - a littel slow but a lot less memory */
                        cpix.resize(nxprobe, nyprobe);
                        cpix.copyInit(probe[0]);
                        sum0 = 0;
                        for (ix = 0; ix < nxprobe; ix++) {
                            for (iy = 0; iy < nyprobe; iy++) {
                                k2 = kxp2[ix] + kyp2[iy];
                                if ((k2 >= k2maxaC) && (k2 <= k2maxbC)) {
                                    phi = atan2(ky[iy], kx[ix]);
                                    /*  offset defocus by zslice so both lens referenced to
                                       entrance surface of specimen  */
                                    chi0 = chi1 * k2 * ((chi2C + chi3C * k2) * k2 - dfC + zslice
                                        + dfa2C * sin(2.0 * (phi - dfa2phiC))
                                        + 2.0F * dfa3C * wavlen * sqrt(k2) *
                                        sin(3.0 * (phi - dfa3phiC)) / 3.0);
                                    chi0 = -chi0;
                                    hr = (float)cos(chi0);
                                    hi = (float)sin(chi0);
                                    prr = probe[ip].re(ix, iy);  // real
                                    pri = probe[ip].im(ix, iy);  // imag
                                    cpix.re(ix, iy) = prr * hr - pri * hi;
                                    cpix.im(ix, iy) = prr * hi + pri * hr;
                                    sum0 += prr * prr + pri * pri;
                                }
                                else {
                                    cpix.re(ix, iy) = 0.0F;
                                    cpix.im(ix, iy) = 0.0F;
                                }
                            }  /*  end for( iy... )  */
                        }  /*  end for( ix... )  */

                        cpix.ifft();

                        /* find normalization constant
                          i.e. correct for constants in the FFT */
                        sum1 = 0.0;
                        for (ix = 0; ix < nxprobe; ix++) {
                            for (iy = 0; iy < nyprobe; iy++) {
                                prr = cpix.re(ix, iy);
                                pri = cpix.im(ix, iy);
                                sum1 += prr * prr + pri * pri;
                            }
                        }

                        /* integrate over real space detector and normalize */
                        for (ix = 0; ix < nxprobe; ix++) {
                            rx2 = xoff[ip] - xp[ix];
                            rx2 = rx2 * rx2;
                            for (iy = 0; iy < nyprobe; iy++) {
                                r2 = yoff[ip] - yp[iy];
                                r2 = rx2 + r2 * r2;
                                prr = cpix.re(ix, iy);
                                pri = cpix.im(ix, iy);
                                delta = prr * prr + pri * pri;
                                for (idetect = 0; idetect < ndetect; idetect++) {
                                    if (CONFOCAL == collectorMode[idetect]) {
                                        if ((r2 >= k2min[idetect]) &&
                                            (r2 < k2max[idetect]))
                                            detect[it][idetect][ip] += delta * (sum0 / sum1);
                                    }
                                }
                            }  /* end for( iy... )*/
                        }  /*  end for( ix....) */

                    }  /* end if( doConfocal==TRUE) */

                    if ((opts4d.save_wvfn || opts4d.save_probe_int_distr)
                        && nslice == p[pNCELLREPZ] - 1) {
                        probe[ip].ifft();
                        if (opts4d.save_wvfn) {
                            fs::path w_path = (fs::path(opts4d.wvfn_dir) / ("layer_" + to_string(nslice) + "_wvfn_ix_" + to_string(opts4d.pos_ix[ip])
                                + "_iy_" + to_string(opts4d.pos_iy[ip]) + ".tif")).make_preferred();
                            //cout << "need to save " << w_path.string() << endl;
                            write_wvfn(probe[ip], w_path.string(), p);
                        }

                        if (opts4d.save_probe_int_distr) {
                            record_probe_intensity_distributions(probe[ip], *distr_xz[ip], *distr_yz[ip], p[pNCELLREPZ]);
                        }
                        probe[ip].fft();
                    }

                }  /* end for( ip.. */

            }  /* end if( ((it...*/
    }

    if (opts4d.save_probe_int_distr)
        for (ip = 0; ip < npos; ip++) {
            int status;
            fs::path xz_path = (fs::path(opts4d.wvfn_dir) / ("xz_probeintdistr_ix_" + to_string(opts4d.pos_ix[ip])
                + "_iy_" + to_string(opts4d.pos_iy[ip]) + ".tif")).make_preferred();
            status = distr_xz[ip]->write(xz_path.string().c_str(), distr_xz[ip]->min(0), distr_xz[ip]->max(0),
                0.0F, 0.0F, 1, 1);
            if (status != 1) cout << "cannot write TIF file " << xz_path << endl;

            fs::path yz_path = (fs::path(opts4d.wvfn_dir) / ("yz_probeintdistr_ix_" + to_string(opts4d.pos_ix[ip])
                + "_iy_" + to_string(opts4d.pos_iy[ip]) + ".tif")).make_preferred();
            status = distr_yz[ip]->write(yz_path.string().c_str(), distr_yz[ip]->min(0), distr_yz[ip]->max(0),
                0.0F, 0.0F, 1, 1);
            if (status != 1) cout << "cannot write TIF file " << yz_path << endl;
        }

    return;

}// end magnstem::STEMsignals() - openMP version

 /*
 * calculates for a single probe position, with multiple phonon configurations
 * 
 * globals that are used: slice_start_idx, slice_end_idx, ...
*/
void magnstem::STEMsignals_wobble(double xp, double yp, int pos_ix, int pos_iy, int nwobble,
    vectorf& p, int multiMode,
    vectori& Znum, vectorf& xa, vectorf& ya, vectorf& za, vectorf& occ, vectorf& wobble,
    vector<unsigned long>& iseeds, vector<cfpix>& mag_phase,
    vectori& ThickSave, float** ca_cbed, 
    double** detect, int ndetect, double& sum, vectori& collectorMode, VzPxavgMaker& maker,
    vectord& phiMin, vectord& phiMax, opts_4dsave opts4d) {
    // generate probes, all at same probe position
    if ((xp < 0.0) || (xp > ax) ||
        (yp < 0.0) || (yp > by)) {
        sbuffer = "bad x,y in STEMsignals() = \n" + toString(xp) + toString(yp);
        messageAST(sbuffer, 0);
        return;
    }
    int ixmid = nxprobe / 2;
    int iymid = nyprobe / 2;
    double chi1 = pi * wavlen;
    double k2maxa = apert1 / wavlen;
    k2maxa = k2maxa * k2maxa;
    double k2maxb = apert2 / wavlen;
    k2maxb = k2maxb * k2maxb;
    /* extra for confocal */
    double chi2C = 0.5 * Cs3C * wavlen * wavlen;
    double chi3C = Cs5C * wavlen * wavlen * wavlen * wavlen / 3.0;
    double k2maxaC = apert1C / wavlen;
    k2maxaC = k2maxaC * k2maxaC;
    double k2maxbC = apert2C / wavlen;
    k2maxbC = k2maxbC * k2maxbC;

    int ixoff = (int)floor(xp * ((double)nx) / ax) - ixmid;
    double xoff = xp - ax * ((double)ixoff) / ((double)nx);
    int iyoff = (int)floor(yp * ((double)ny) / by) - iymid;
    double yoff = yp - by * ((double)iyoff) / ((double)ny);
    double sum0 = 0.0;
    for (int ix = 0; ix < nxprobe; ix++) {
        double alx = wavlen * kxp[ix];  /* x component of angle alpha */
        for (int iy = 0; iy < nyprobe; iy++) {
            double aly = wavlen * kyp[iy];  /* y component of angle alpha */
            double k2 = kxp2[ix] + kyp2[iy];
            if ((k2 >= k2maxa) && (k2 <= k2maxb)) {
                double w = 2.0 * pi * (xoff * kxp[ix] + yoff * kyp[iy]);
                double chi0 = (2.0 * pi / wavlen) * chi(p, alx, aly, multiMode);
                chi0 = -chi0 + w;
                float tr = (float)cos(chi0);
                probe[0].re(ix, iy) = tr;
                float ti = (float)sin(chi0);
                probe[0].im(ix, iy) = ti;
                sum0 += (double)(tr * tr + ti * ti);
            }
            else {
                probe[0].re(ix, iy) = 0.0F;
                probe[0].im(ix, iy) = 0.0F;
            }
        }
    }  /* end for( ix... */

    float scale = (float)(1.0 / sqrt(sum0));
    probe[0] *= scale;
    for (int ip = 1; ip < nwobble; ip++) {
        probe[ip] = probe[0];
    }

    float temperature = p[pTEMPER];
    float temp_scale = (float) sqrt(temperature / 300.0);
    for (nslice = 0; nslice < p[pNCELLREPZ]; nslice++) {
#pragma omp parallel for
        for (int ip = 0; ip < nwobble; ip++) {
            // add random thermal displacements for atoms in unit cell
            vectorf xa2, ya2;
            xa2.resize(natom);
            ya2.resize(natom);
            //double rand_num = 0.0;
            for (int i = 0; i < natom; i++) {
                //rand_num = rangauss(&iseeds[ip]);
                xa2[i] = xa[i] +
                    (float)(wobble[i] * rangauss(&iseeds[ip]) * temp_scale);
                //rand_num = rangauss(&iseeds[ip]);
                ya2[i] = ya[i] +
                    (float)(wobble[i] * rangauss(&iseeds[ip]) * temp_scale);
                //cout << rand_num << " ";
                //if (i == 100) {
                //    cout << xa[i] << "_" << xa2[i] << " ";
                //}
            }

            for (int itr = 0; itr < slice_per_cell; itr++) {
                probe[ip].ifft();

                // calculate transmission function
                // use separate transmission function for each probe
                double phirms;
                int natom_slice = slice_end_idx[itr] - slice_start_idx[itr] + 1;
                trlayer_pxavg_LUT(xa2, ya2, occ,
                    Znum, natom_slice, slice_start_idx[itr], (float)ax, (float)by, (float)keV,
                    trans[ip], mag_phase[itr], (long)nx, (long)ny, &phirms, &nbeamt, (float)k2maxp, maker);
                    //pos_ix*100000 + pos_iy*10000 + nslice*1000 + ip*10 + itr, p, opts4d);

                // apply transmission function
                // (transmission function shifted for probe position)
                for (int ix = 0; ix < nxprobe; ix++) {
                    int ixt = ix + ixoff;
                    if (ixt >= nx) ixt = ixt - nx;
                    else if (ixt < 0) ixt = ixt + nx;
                    for (int iy = 0; iy < nyprobe; iy++) {
                        int iyt = iy + iyoff;
                        if (iyt >= ny) iyt = iyt - ny;
                        else if (iyt < 0) iyt = iyt + ny;
                        float prr = probe[ip].re(ix, iy);
                        float pri = probe[ip].im(ix, iy);
                        probe[ip].re(ix, iy) = prr * trans[ip].re(ixt, iyt)
                            - pri * trans[ip].im(ixt, iyt);  // real
                        probe[ip].im(ix, iy) = prr * trans[ip].im(ixt, iyt)
                            + pri * trans[ip].re(ixt, iyt);  // imag
                    } /* end for(iy...) */
                }  /* end for(ix...) */

                probe[ip].fft();

                /*  multiplied by the propagator function */
                probe[ip] *= cprop;
                
                //cout << " " << ip << " " ;
            }
        } // end for( ip...)

        vectori::iterator iter = find(ThickSave.begin(), ThickSave.end(), nslice);
        if (iter != ThickSave.end()) {
            sum = 0.0;
            int it = distance(ThickSave.begin(), iter);

            // average CBED pattern at this scan position and thickness
            for (int ix = 0; ix < nxprobe; ix++) {
                for (int iy = 0; iy < nyprobe; iy++) {
                    ca_cbed[ix][iy] = 0;
                }
            }
            for (int ip = 0; ip < nwobble; ip++) {
                for (int ix = 0; ix < nxprobe; ix++) {
                    for (int iy = 0; iy < nyprobe; iy++) {
                        float prr = probe[ip].re(ix, iy);
                        float pri = probe[ip].im(ix, iy);
                        ca_cbed[ix][iy] += (prr * prr + pri * pri) / (float)nwobble;
                    }
                }
            }

            // sum intensity incident on detectors
            for (int id = 0; id < ndetect; id++) detect[it][id] = 0.0;
            for (int ix = 0; ix < nxprobe; ix++) {
                for (int iy = 0; iy < nyprobe; iy++) {
                    double delta = (double)(ca_cbed[ix][iy]);
                    sum += delta;
                    double k2 = (double)kxp2[ix] + (double)kyp2[iy];
                    double phi = atan2(kyp[iy], kxp[ix]);
                    for (int idetect = 0; idetect < ndetect; idetect++) {
                        if (ADF == collectorMode[idetect]) {
                            if ((k2 >= k2min[idetect]) &&
                                (k2 < k2max[idetect]))
                                detect[it][idetect] += delta;
                        }
                        else if (ADF_SEG == collectorMode[idetect]) {
                            if ((k2 >= k2min[idetect]) && (k2 < k2max[idetect])
                                && (phi >= phiMin[idetect]) && (phi < phiMax[idetect]))
                                detect[it][idetect] += delta;
                        }
                    }
                }
            }

            // don't bother handling confocal (not sure what that is)

            // save cropped CBED
            invert2D(ca_cbed, nxprobe, nyprobe);
            if (opts4d.save_cbed) {
                fs::path c_path = (fs::path(opts4d.cbed_dir) / ("layer_" + to_string(nslice) + "_cbed_ix_" + to_string(opts4d.pos_ix[0])
                    + "_iy_" + to_string(opts4d.pos_iy[0]) + ".tif")).make_preferred();
                write_cbed(ca_cbed, c_path.string(), p, opts4d);
            }
        }
    }

} // end magnstem::STEMsignals_wobble()

/* -------------------  CountBeams() -------------------
   copied from autostem.cpp

   count the number of beams (Fourier coefficients) in the
   transmission function and the probe for informational purposes
   - to communicate to the main calling program
   -  must be recalculated here (should be identical to that
      used in calculate()
  input:
        param[] = array of parameters

  output:
        nbeamp  = number of symmetrical beams in the probe
        nbeampo = number of beams in the probe aperture
        res  = bandwidth limited resolution in Angstroms
        almax = alpha max in radians
*/
void magnstem::CountBeams(vectorf& param, int& nbeamp, int& nbeampo, float& res, float& almax)
{
    //  use all local array variables so I don't accidentally
    //     disturb the main calculation
    int ix, iy, nx1, ny1, nxprobe1, nyprobe1;
    float ax, by, wavl, apert1, apert2;
    double sum, k2maxp, k2maxa, k2maxb, k2;

    vectorf xp1, yp1, kxp1, kyp1, kxp21, kyp21;

    ax = param[pAX];                  // supercell size in Ang.
    by = param[pBY];

    nx1 = ToInt(param[pNX]);        // size of transmission function (in pixels)
    ny1 = ToInt(param[pNY]);

    nxprobe1 = ToInt(param[pNXPRB]);        // probe size in pixels
    nyprobe1 = ToInt(param[pNYPRB]);

    wavl = (float)wavelength(param[pENERGY]);      //  electron wavelength in Ang.

    apert1 = param[pOAPMIN];  // obj. apert size in radians (min, max for annular apert)
    apert2 = param[pOAPERT];

    //  transmission function sampling
    // impose anti-aliasing bandwidth limit on transmission functions
    sum = ((double)nx1) / (2.0 * ax);
    k2maxp = ((double)ny1) / (2.0 * by);
    if (sum < k2maxp) k2maxp = sum;
    k2maxp = BW * k2maxp;
    res = (float)(1.0 / k2maxp);
    almax = (float)(wavl * k2maxp);
    k2maxp = k2maxp * k2maxp;

    //  probe sampling
    //  the only way to do this is to do part of the calculation and
    //    throw away the results (kx,ky etc.) - but only a small bit wasted
    xp1.resize(nxprobe1);
    yp1.resize(nyprobe1);
    kxp1.resize(nxprobe1);
    kyp1.resize(nyprobe1);
    kxp21.resize(nxprobe1);
    kyp21.resize(nyprobe1);

    freqn(kxp1, kxp21, xp1, nxprobe1, ax * ((double)nxprobe1) / nx1);
    freqn(kyp1, kyp21, yp1, nyprobe1, by * ((double)nyprobe1) / ny1);

    /*   calculate number of pixels in the probe and obj. apert. */
    k2maxa = apert1 / wavl;
    k2maxa = k2maxa * k2maxa;
    k2maxb = apert2 / wavl;
    k2maxb = k2maxb * k2maxb;

    nbeamp = nbeampo = 0;
    for (iy = 0; iy < nyprobe1; iy++)
        for (ix = 0; ix < nxprobe1; ix++) {
            k2 = kyp21[iy] + kxp21[ix];
            if (k2 < k2maxp) nbeamp++;
            if ((k2 >= k2maxa) && (k2 <= k2maxb)) nbeampo++;
        }

};  // end magnstem::CountBeams()

/*--------------------- trlayer() -----------------------*/
/*  copied from autostem.cpp

  Calculate complex specimen transmission function
  for one layer using real space projected atomic potentials

  x[],y[] = real array of atomic coordinates
  occ[]   = real array of occupancies
  Znum[]  = array of atomic numbers
  istart  = starting index of atom coord.
  natom   = number of atoms
  ax, by  = size of transmission function in Angstroms
  kev     = beam energy in keV
  trans   = 2D array to get complex specimen
        transmission function
  nx, ny  = dimensions of transmission functions
  *phirms = average phase shift of projected atomic potential
  *nbeams = will get number of Fourier coefficients
  k2max   = square of max k = bandwidth limit

  convert to cfpix class for trans 10-nov-2012 ejk
  convert arrays to vector<> 25-dec-2017 ejk
  add option to return just potential if k2max<0 28-sep-2018 ejk

*/
void magnstem::trlayer(const vectorf& x, const vectorf& y, const vectorf& occ,
    const vectori& Znum, const int natom, const int istart,
    const float ax, const float by, const float kev,
    cfpix& trans, cfpix& mag_image, const long nx, const long ny,
    double* phirms, long* nbeams, const float k2max)
{
    int idx, idy, i, ixo, iyo, ix, iy, ixw, iyw, nx1, nx2, ny1, ny2;
    float k2;
    double r, rx2, rsq, vz, rmin, rmin2, sum, scale, scalex, scaley;

    const double rmax = 3.0, rmax2 = rmax * rmax; /* max atomic radius in Angstroms */

    scale = sigma(kev) / 1000.0;  /* in 1/(volt-Angstroms) */

    scalex = ax / nx;
    scaley = by / ny;

    /* min radius to avoid  singularity */
    rmin = ax / ((double)nx);
    r = by / ((double)ny);
    rmin = 0.25 * sqrt(0.5 * (rmin * rmin + r * r));
    rmin2 = rmin * rmin;

    idx = (int)(nx * rmax / ax) + 1;
    idy = (int)(ny * rmax / by) + 1;

    for (ix = 0; ix < nx; ix++) {
        for (iy = 0; iy < ny; iy++)
            trans.re(ix, iy) = 0.0F;
    }

    bool do_magnetism = (mag_image.nx() == nx) && (mag_image.ny() == ny);
    //if (do_magnetism) {
    //    cout << "adding magnetism phase to transmission function" << endl;
    //}

    /*  paralleling this loop has little effect   */
    /*#pragma omp parallel for private(ix,iy,ixo,iyo,nx1,nx2,ny1,ny2,rx2,r,ixw,iyw,vz,rsq)*/
//#pragma omp parallel for private(ix,iy,ixo,iyo,nx1,nx2,ny1,ny2,rx2,ixw,iyw,vz,rsq)
    for (i = istart; i < (istart + natom); i++) {
        ixo = (int)(x[i] / scalex);
        iyo = (int)(y[i] / scaley);
        nx1 = ixo - idx;
        nx2 = ixo + idx;
        ny1 = iyo - idy;
        ny2 = iyo + idy;

        /* add proj. atomic potential at a local region near its center
           taking advantage of small range of atomic potential */

        for (ix = nx1; ix <= nx2; ix++) {
            rx2 = x[i] - ((double)ix) * scalex;
            rx2 = rx2 * rx2;
            ixw = ix;
            while (ixw < 0) ixw = ixw + nx;
            ixw = ixw % nx;
            for (iy = ny1; iy <= ny2; iy++) {
                rsq = y[i] - ((double)iy) * scaley;
                rsq = rx2 + rsq * rsq;
                if (rsq <= rmax2) {
                    iyw = iy;
                    while (iyw < 0) iyw = iyw + ny;
                    iyw = iyw % ny;
                    if (rsq < rmin2) rsq = rmin2;
                    /*r = sqrt( rx2 + r*r );
                    vz = occ[i] * vzatom( Znum[i], r ); slow */
                    vz = occ[i] * vzatomLUT(Znum[i], rsq);
//#pragma omp critical
                    {
                        trans.re(ixw, iyw) += (float)vz;
                    }
                }
            } /* end for(iy... */
        }  /* end for(ix... */

    } /* end for(i=0... */

    if (k2max > 0.0) {
        /* convert phase to a complex transmission function */
        sum = 0;
        for (ix = 0; ix < nx; ix++) {
            for (iy = 0; iy < ny; iy++) {
                vz = scale * trans.re(ix, iy);
                if (do_magnetism) {
                    vz += mag_image.re(ix, iy);
                }
                sum += vz;
                trans.re(ix, iy) = (float)cos(vz);
                trans.im(ix, iy) = (float)sin(vz);
            }
        }

        /* bandwidth limit the transmission function */
        *nbeams = 0;
        trans.fft();
        for (ix = 0; ix < nx; ix++) {
            for (iy = 0; iy < ny; iy++) {
                k2 = ky2[iy] + kx2[ix];
                if (k2 < k2max) *nbeams += 1;
                else trans.re(ix, iy) = trans.im(ix, iy) = 0.0F;
            }
        }
        trans.ifft();

    }
    else {

        //  just scale to sigma * Vz and return real valued potential
        sum = 0;
        for (ix = 0; ix < nx; ix++)
            for (iy = 0; iy < ny; iy++) {
                trans.re(ix, iy) *= (float)scale;
                sum += trans.re(ix, iy);
            }
    }

    *phirms = sum / (((double)nx) * ((double)ny));

    return;

};  /* end magnstem::trlayer() */

/*--------------------- trlayer_pxavg() -----------------------*/
/*

  Calculate complex specimen transmission function
  for one layer using real space projected atomic potentials
  Integrates over the region of each pixel.

  x[],y[] = real array of atomic coordinates
  occ[]   = real array of occupancies
  Znum[]  = array of atomic numbers
  istart  = starting index of atom coord.
  natom   = number of atoms
  ax, by  = size of transmission function in Angstroms
  kev     = beam energy in keV
  trans   = 2D array to get complex specimen
        transmission function
  nx, ny  = dimensions of transmission functions
  *phirms = average phase shift of projected atomic potential
  *nbeams = will get number of Fourier coefficients
  k2max   = square of max k = bandwidth limit

  2021-12

*/
void magnstem::trlayer_pxavg(const vectorf& x, const vectorf& y, const vectorf& occ,
    const vectori& Znum, const int natom, const int istart,
    const float ax, const float by, const float kev,
    cfpix& trans, cfpix& mag_image, const long nx, const long ny,
    double* phirms, long* nbeams, const float k2max,
    const int it, vectorf& param, opts_4dsave opts4d)
{
    int idx, idy, i, ixo, iyo, ix, iy, ixw, iyw, nx1, nx2, ny1, ny2;
    float k2;
    double r, rx2, rsq, vz, rmin, rmin2, sum, scale, scalex, scaley;
    double dx, dy;

    const double rmax = 3.0, rmax2 = rmax * rmax; /* max atomic radius in Angstroms */

    scale = sigma(kev) / 1000.0;  /* in 1/(volt-Angstroms) */

    scalex = ax / nx;
    scaley = by / ny;

    /* min radius to avoid  singularity */
    rmin = ax / ((double)nx);
    r = by / ((double)ny);
    rmin = 0.25 * sqrt(0.5 * (rmin * rmin + r * r));
    rmin2 = rmin * rmin;

    idx = (int)(nx * rmax / ax) + 1;
    idy = (int)(ny * rmax / by) + 1;

    for (ix = 0; ix < nx; ix++) {
        for (iy = 0; iy < ny; iy++)
            trans.re(ix, iy) = 0.0F;
    }

    bool do_magnetism = (mag_image.nx() == nx) && (mag_image.ny() == ny);
    //if (do_magnetism) {
    //    cout << "adding magnetism phase to transmission function" << endl;
    //}

    /*  parallelize this loop because vzatom_pxavg is slow and no parallel loops at higher level */
    /*#pragma omp parallel for private(ix,iy,ixo,iyo,nx1,nx2,ny1,ny2,rx2,r,ixw,iyw,vz,rsq)*/
#pragma omp parallel for private(ix,iy,ixo,iyo,nx1,nx2,ny1,ny2,rx2,ixw,iyw,vz,rsq,dx,dy)
    for (i = istart; i < (istart + natom); i++) {
        //if (i % 100 == 0) {
        //    printf("%d/%d ", i, istart+natom);
        //    fflush(stdout);
        //}
        ixo = (int)(x[i] / scalex);
        iyo = (int)(y[i] / scaley);
        nx1 = ixo - idx;
        nx2 = ixo + idx;
        ny1 = iyo - idy;
        ny2 = iyo + idy;

        /* add proj. atomic potential at a local region near its center
           taking advantage of small range of atomic potential */

        for (ix = nx1; ix <= nx2; ix++) {
            dx = x[i] - ((double)ix) * scalex;
            rx2 = dx * dx;
            ixw = ix;
            while (ixw < 0) ixw = ixw + nx;
            ixw = ixw % nx;
            for (iy = ny1; iy <= ny2; iy++) {
                dy = y[i] - ((double)iy) * scaley;
                rsq = rx2 + dy * dy;
                if (rsq <= rmax2) {
                    iyw = iy;
                    while (iyw < 0) iyw = iyw + ny;
                    iyw = iyw % ny;
                    // if (rsq < rmin2) rsq = rmin2;
                    /*r = sqrt( rx2 + r*r );
                    vz = occ[i] * vzatom( Znum[i], r ); slow */
                    vz = occ[i] * vzatom_pxavg(Znum[i], dx, dy, scalex, scaley);
#pragma omp critical
                    {
                        trans.re(ixw, iyw) += (float)vz;
                    }
                }
            } /* end for(iy... */
        }  /* end for(ix... */

    } /* end for(i=0... */

    //fs::path w_path = (fs::path(opts4d.cbed_dir) / ("layer_" + to_string(it) + "_vz_pxavg_nx_" + to_string(nxprobe)
    //    + "_ny_" + to_string(nyprobe) + ".tif")).make_preferred();
    ////cout << "need to save " << w_path.string() << endl;
    //write_real_space_real(trans, w_path.string(), param);

    if (k2max > 0.0) {
        /* convert phase to a complex transmission function */
        sum = 0;
        for (ix = 0; ix < nx; ix++) {
            for (iy = 0; iy < ny; iy++) {
                vz = scale * trans.re(ix, iy);
                if (do_magnetism) {
                    vz += mag_image.re(ix, iy);
                }
                sum += vz;
                trans.re(ix, iy) = (float)cos(vz);
                trans.im(ix, iy) = (float)sin(vz);
            }
        }

        /* bandwidth limit the transmission function */
        *nbeams = 0;
        trans.fft();
        for (ix = 0; ix < nx; ix++) {
            for (iy = 0; iy < ny; iy++) {
                k2 = ky2[iy] + kx2[ix];
                if (k2 < k2max) *nbeams += 1;
                else trans.re(ix, iy) = trans.im(ix, iy) = 0.0F;
            }
        }
        trans.ifft();

    }
    else {

        //  just scale to sigma * Vz and return real valued potential
        sum = 0;
        for (ix = 0; ix < nx; ix++)
            for (iy = 0; iy < ny; iy++) {
                trans.re(ix, iy) *= (float)scale;
                sum += trans.re(ix, iy);
            }
    }

    *phirms = sum / (((double)nx) * ((double)ny));

    return;

};  /* end magnstem::trlayer_pxavg() */

/*--------------------- trlayer_pxavg_LUT() -----------------------*/
/*

  Calculate complex specimen transmission function
  for one layer using real space projected atomic potentials
  Integrates over the region of each pixel.
  Uses an existing lookup table in a VzPxavgMaker instance.

  x[],y[] = real array of atomic coordinates
  occ[]   = real array of occupancies
  Znum[]  = array of atomic numbers
  istart  = starting index of atom coord.
  natom   = number of atoms
  ax, by  = size of transmission function in Angstroms
  kev     = beam energy in keV
  trans   = 2D array to get complex specimen
        transmission function
  nx, ny  = dimensions of transmission functions
  *phirms = average phase shift of projected atomic potential
  *nbeams = will get number of Fourier coefficients
  k2max   = square of max k = bandwidth limit

  2021-12

*/
void magnstem::trlayer_pxavg_LUT(const vectorf& x, const vectorf& y, const vectorf& occ,
    const vectori& Znum, const int natom, const int istart,
    const float ax, const float by, const float kev,
    cfpix& trans, cfpix& mag_image, const long nx, const long ny,
    double* phirms, long* nbeams, const float k2max, VzPxavgMaker& maker)
    //const int it, vectorf& param, opts_4dsave opts4d)
{
    int idx, idy, ixo, iyo, ix, iy, ixw, iyw, nx1, nx2, ny1, ny2;
    float k2;
    double r, rx2, rsq, vz, rmin, rmin2, sum, scale, scalex, scaley;
    double dx, dy;

    const double rmax = 3.0, rmax2 = rmax * rmax; /* max atomic radius in Angstroms */

    scale = sigma(kev) / 1000.0;  /* in 1/(volt-Angstroms) */

    scalex = ax / nx;
    scaley = by / ny;

    /* min radius to avoid  singularity */
    rmin = ax / ((double)nx);
    r = by / ((double)ny);
    rmin = 0.25 * sqrt(0.5 * (rmin * rmin + r * r));
    rmin2 = rmin * rmin;

    idx = (int)(nx * rmax / ax) + 1;
    idy = (int)(ny * rmax / by) + 1;

    for (ix = 0; ix < nx; ix++) {
        for (iy = 0; iy < ny; iy++)
            trans.re(ix, iy) = 0.0F;
    }

    maker.fill_vz(x, y, occ, Znum, natom, istart, ax, by, trans, nx, ny);

    //if (it < 100) {
    //    fs::path w_path = (fs::path(opts4d.cbed_dir) / ("trans_fn_" + to_string(it) + "_vz_pxavg_nx_" + to_string(nxprobe)
    //        + "_ny_" + to_string(nyprobe) + ".tif")).make_preferred();
    //    //cout << "need to save " << w_path.string() << endl;
    //    write_real_space_real(trans, w_path.string(), param);

    //    fs::path w_path2 = (fs::path(opts4d.cbed_dir) / ("coords_" + to_string(it) + "_vz_pxavg_nx_" + to_string(nxprobe)
    //        + "_ny_" + to_string(nyprobe) + ".txt")).make_preferred();
    //    FILE *coord_file;
    //    coord_file = fopen(w_path2.string().c_str(), "w");
    //    fprintf(coord_file, "Z x y\n");
    //    for (int i = istart; i < (istart + natom); ++i) {
    //        fprintf(coord_file, "%d %.6f %.6f\n", Znum[i], x[i], y[i]);
    //    }
    //    fclose(coord_file);
    //}

    if (k2max > 0.0) {
        /* convert phase to a complex transmission function */
        sum = 0;
        for (ix = 0; ix < nx; ix++) {
            for (iy = 0; iy < ny; iy++) {
                vz = scale * trans.re(ix, iy);
                sum += vz;
                trans.re(ix, iy) = (float)cos(vz);
                trans.im(ix, iy) = (float)sin(vz);
            }
        }

        /* bandwidth limit the transmission function */
        *nbeams = 0;
        trans.fft();
        for (ix = 0; ix < nx; ix++) {
            for (iy = 0; iy < ny; iy++) {
                k2 = ky2[iy] + kx2[ix];
                if (k2 < k2max) *nbeams += 1;
                else trans.re(ix, iy) = trans.im(ix, iy) = 0.0F;
            }
        }
        trans.ifft();

    }
    else {

        //  just scale to sigma * Vz and return real valued potential
        sum = 0;
        for (ix = 0; ix < nx; ix++)
            for (iy = 0; iy < ny; iy++) {
                trans.re(ix, iy) *= (float)scale;
                sum += trans.re(ix, iy);
            }
    }

    *phirms = sum / (((double)nx) * ((double)ny));

    return;

};  /* end magnstem::trlayer_pxavg() */

/*------------------------- invert2D() ----------------------*/
/*
    (copied from autostem.cpp)
        rearrange pix with corners moved to center (a la FFT's)

         pix[ix][iy] = real array with image
         nx,ny = range of image 0<ix<(nx-1) and 0<iy<(ny-1)

*/
void magnstem::invert2D(float** pix, long nx, long ny)
{
#define SWAP(a,b)       {t=a; a=b; b=t;}

    long ix, iy, ixmid, iymid;
    float t;

    ixmid = nx / 2;
    iymid = ny / 2;

    for (ix = 0; ix < nx; ix++)
        for (iy = 0; iy < iymid; iy++)
            SWAP(pix[ix][iy], pix[ix][iy + iymid]);

    for (ix = 0; ix < ixmid; ix++)
        for (iy = 0; iy < ny; iy++)
            SWAP(pix[ix][iy], pix[ix + ixmid][iy]);

#undef SWAP
};  // end autostem::invert2D()


/*
    Takes in a probe wavefunction in k-space and outputs CBED to file
*/
void magnstem::write_cbed_from_probe(cfpix& probe, const string path, 
    vectorf& param, opts_4dsave opts4d) {
    float** cbed_pix = new2D<float>(nxprobe, nyprobe, "cbed_pix");
    float prr, pri;
    for (int ix2 = 0; ix2 < nxprobe; ix2++)
        for (int iy2 = 0; iy2 < nyprobe; iy2++) {
            prr = probe.re(ix2, iy2);
            pri = probe.im(ix2, iy2);
            cbed_pix[ix2][iy2] = (prr * prr + pri * pri);
        }

    invert2D(cbed_pix, nxprobe, nyprobe);
    write_cbed(cbed_pix, path, param, opts4d);
    delete2D<float>(cbed_pix, nxprobe);
}

void magnstem::write_cbed(float** pix, const string path,
    vectorf& param, opts_4dsave opts4d) {
    int nxout2 = 0;
    int nyout2 = 0;
    int nx1, nx2, ny1, ny2;
    if (opts4d.cbed_nxout > 1) {
        nx1 = nxprobe / 2 - opts4d.cbed_nxout / 2;
        nx2 = nxprobe / 2 + opts4d.cbed_nxout / 2;
        if (opts4d.cbed_nxout % 2 == 0) {
            nx2--;
        }
        if (nx1 < 0) nx1 = 0;
        if (nx2 >= nxprobe) nx2 = nxprobe - 1;
        nxout2 = nx2 - nx1 + 1;
    }
    if (opts4d.cbed_nyout > 1) {
        ny1 = nyprobe / 2 - opts4d.cbed_nyout / 2;
        ny2 = nyprobe / 2 + opts4d.cbed_nyout / 2;
        if (opts4d.cbed_nyout % 2 == 0) {
            ny2--;
        }
        if (ny1 < 0) ny1 = 0;
        if (ny2 >= nyprobe) ny2 = nyprobe - 1;
        nyout2 = ny2 - ny1 + 1;
    }
    if ((nxout2 < 1) || (nyout2 < 1)) {
        nx1 = ny1 = 0;
        nx2 = nxprobe - 1;
        nxout2 = nxprobe;
        ny2 = nyprobe - 1;
        nyout2 = nyprobe;
    }
    floatTIFF myFile;
    int NPARAM = myFile.maxParam();
    for (int i = 0; i < NPARAM; i++) myFile.setParam(i, param[i]);
    myFile.setnpix(1);
    myFile.resize(nxout2, nyout2);
    float aimin = 0.0F;
    float aimax = 0.0F;
    //float scalef = (float)(1.0 / (nxout2 * nyout2));
    int ixo = 0;
    for (int ix2 = nx1; ix2 <= nx2; ix2++) {
        int iyo = 0;
        for (int iy2 = ny1; iy2 <= ny2; iy2++) {
            myFile(ixo, iyo++) = pix[ix2][iy2];
        }  ixo++;
    }

    float rmin0 = myFile.min(0);
    float rmax0 = myFile.max(0);
    double dxp = ax * ((double)nxprobe) / nx;
    double dyp = by * ((double)nyprobe) / ny;
    dxp = 1.0 / dxp;
    dyp = 1.0 / dyp;
    myFile.setParam(pRMAX, rmax0);
    myFile.setParam(pRMIN, rmin0);
    myFile.setParam(pDX, (float)dxp);
    myFile.setParam(pDY, (float)dyp);
    //cout << "CBED (unaliased) size " << nxout2 << " x " << nyout2 << " pixels and range (arb. units): " << rmin0 << " to " << rmax0 << endl;
    if (myFile.write(path.c_str(), rmin0, rmax0, aimin, aimax,
        (float)dxp, (float)dyp) != 1) {
        cout << "Cannot write output file " << path << endl;
    }
}

/*
    Takes in a probe wavefunction in real space and outputs wavefunction to file
*/
void magnstem::write_wvfn(cfpix& probe, const string path, vectorf& param) {
    floatTIFF myFile;
    float rmin, rmax, aimin, aimax, dx, dy;

    int NPARAM = myFile.maxParam();
    for (int i = 0; i < NPARAM; i++) myFile.setParam(i, param[i]);
    param[pDX] = dx = (float)(ax / ((float)nx));
    param[pDY] = dy = (float)(by / ((float)ny));    
    
    myFile.setnpix(2);
    myFile.resize(2 * nxprobe, nyprobe);
   
    probe.findRange(rmin, rmax, aimin, aimax);
    for (int ix = 0; ix < nxprobe; ix++) for (int iy = 0; iy < nyprobe; iy++) {
        myFile(ix, iy) = probe.re(ix, iy);
        myFile(ix + nx, iy) = probe.im(ix, iy);
    }

    int i = myFile.write(path.c_str(), rmin, rmax, aimin, aimax, dx, dy);
    cout << "wvfn pix range " << rmin << " to " << rmax << " real, " << aimin << " to " << aimax << " imag" << endl;
    if (i != 1) cout << "cannot write TIF file " << path << endl;
}

/*
    Outputs a real-valued real space image to file
*/
void magnstem::write_real_space_real(cfpix& probe, const string path, vectorf& param) {
    floatTIFF myFile;
    float rmin, rmax, aimin, aimax, dx, dy;

    int NPARAM = myFile.maxParam();
    param[pDX] = dx = (float)(ax / ((float)nx));
    param[pDY] = dy = (float)(by / ((float)ny));
    for (int i = 0; i < NPARAM; i++) myFile.setParam(i, param[i]);

    myFile.setnpix(1);
    myFile.resize(nxprobe, nyprobe);

    //probe.findRange(rmin, rmax, aimin, aimax);
    for (int ix = 0; ix < nxprobe; ix++) for (int iy = 0; iy < nyprobe; iy++) {
        myFile(ix, iy) = probe.re(ix, iy);
    }

    float rmin0 = myFile.min(0);
    float rmax0 = myFile.max(0);

    int i = myFile.write(path.c_str(), rmin0, rmax0, 0.0F, 0.0F, dx, dy);
    cout << "real space pix range " << rmin0 << " to " << rmax0 << " real" << endl;
    if (i != 1) cout << "cannot write TIF file " << path << endl;
}

void magnstem::record_probe_intensity_distributions(cfpix& probe, floatTIFF& xz_img, floatTIFF& yz_img, int z_idx) {
    float sum2;
    int ix, iy;
    float prr, pri;
    for (ix = 0; ix < nxprobe; ix++) {
        sum2 = 0;
        for (iy = 0; iy < nyprobe; iy++) {
            prr = probe.re(ix, iy);
            pri = probe.im(ix, iy);
            sum2 += prr * prr + pri * pri;
        }
        xz_img(ix, nslice) = sum2;
    }
    for (iy = 0; iy < nyprobe; iy++) {
        sum2 = 0;
        for (ix = 0; ix < nxprobe; ix++) {
            prr = probe.re(ix, iy);
            pri = probe.im(ix, iy);
            sum2 += prr * prr + pri * pri;
        }
        yz_img(iy, nslice) = sum2;
    }
}
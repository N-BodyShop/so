#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include "kd2.h"

/* The following 2 functions come from Kitayama & Suto 1996, ApJ, 469, 480 */

#define sqr(x) ((x)*(x))
double Omegaf(double Omega0, double Lambda0, double z)
{
    double answer,zplus12,zplus13;

    zplus12 = sqr(1.0+z);
    zplus13 = zplus12*(1.0+z);
    answer = Omega0*zplus13/(Omega0*zplus13+(1.-Omega0-Lambda0)*zplus12+Lambda0);
    /*    fprintf(stderr,"Omegaf: %lg\n",answer); */
    return (answer);
}

double rhovir_over_rhobar(double Omega0, int Lambda_OPT, double z)
{
    double etaf,wf,Lambda0,answer;

    if (Omega0 == 1.0) {
	return(178.0);
    }
    if (Lambda_OPT) {
	Lambda0 = 1.0-Omega0;
	wf = 1./Omegaf(Omega0,Lambda0,z) - 1.;
	answer = 18.*sqr(M_PI)*(1.+0.4093*pow(wf,0.9052));
    } else {
	etaf = acosh(2.0/Omegaf(Omega0,0.0,z) - 1.0);
	answer =  4.0*sqr(M_PI)/(pow(sinh(etaf)-etaf,2));
	answer *= pow(cosh(etaf)-1,3);
    }
    /*    fprintf(stderr,"r/r: %lg\n",answer); */
    return(answer);
}

void usage(void)
{
    fprintf(stderr,"USAGE:\n");
    fprintf(stderr,"so -i <FOF .gtp file> [-o <outfilebase>] [([-dark] [-gas] [-star]) || [-all])]\n");
    fprintf(stderr,"          [-mark <markfile>]  [-std] [-list <File containing group indexes>]\n");
    fprintf(stderr,"          [-rho <fThreshold>] [-M <fMinGTPMass>] [-m <mMinSOMembers>]\n");
    fprintf(stderr,"          [-O <fOmega0>]  [-L]  [-z <fRedshift>]\n");
    fprintf(stderr,"          [-p <xyzPeriod>]  [-c <xyzCenter>]\n");
    fprintf(stderr,"          [-cx <xCenter>]  [-cy <yCenter>]  [-cz <zCenter>]\n");
    fprintf(stderr,"          [-u <fMassUnit> <fMpcUnit>]\n\n");

    fprintf(stderr,"  Finds smallest spherical radius R within which the density is less than\n");
    fprintf(stderr,"  the value <fThreshold> and calculates enclosed mass, half mass and quarter\n");
    fprintf(stderr,"  mass radii, maximum circular velocity and max Vcirc radius, as well as\n");
    fprintf(stderr,"  Vcirc in 1/4R increments from 1/4R to 2R.\n");
    fprintf(stderr,"   -o: Main output will be written to file <outfilebase>.sovcirc.\n");
    fprintf(stderr,"       default = so.sovcirc.\n");
    fprintf(stderr,"   -dark -gas -star: Outputs radial profile of dark, gas, and/or star particles\n");
    fprintf(stderr,"       in a separate file <outfilebase>.sodark / .sogas / .sostar\n");
    fprintf(stderr,"   -all: same as '-dark -gas -star'.\n");
    fprintf(stderr,"   -mark: Outputs radial profile of marked particles given in <markfile> to\n");
    fprintf(stderr,"       to a separate file <outfilebase>.somark\n.");
    fprintf(stderr,"   -std: Read and write standard TIPSY binaries.\n");
    fprintf(stderr,"   -rho: By default, <fThreshold> is calculated automatically from cosmological\n");
    fprintf(stderr,"       parameters to be the virial density, but may be overridden using the\n");
    fprintf(stderr,"       -rho option.\n");
    fprintf(stderr,"   -L: option will set Lambda0 = 1-Omega0.\n");
    fprintf(stderr,"   -p: THIS VERSION ASSUMES PERIODIC BOUNDRY CONDITIONS.  Default period = 1\n");
    fprintf(stderr,"   -z: <fRedshift> will be automatically calculated as 1/h.time-1 where h.time\n");
    fprintf(stderr,"       is read from the TIPSY particle file header.  -z will override this.\n");
    fprintf(stderr,"   -M: <fMinGTPMass> specifies the minimum mass a group given in the .gtp file\n");
    fprintf(stderr,"       must have in order to be considered.\n");
    fprintf(stderr,"   -list: names an ASCII file which lists indexes of groups in the .gtp file\n");
    fprintf(stderr,"       to be considered.  This may be used in conjunction with -M.\n");
    fprintf(stderr,"   -m: <fMinSOMembers> specifies the minimum number of particles a group must\n");
    fprintf(stderr,"       have for its parameters to be accurately calculated.\n");
    fprintf(stderr,"   -u: Output is in system units unless simulation units are specified with -u\n");
    fprintf(stderr,"       option, whereupon output will be in the form:\n");
    fprintf(stderr,"       Mass [Msol], Distance [kpc], Velocity [km/s].\n");
    fprintf(stderr,"  See output file header for output format.\n");
    fprintf(stderr,"  Groupwise error codes are returned in Mvir and Rvir columns as follows:\n");
    fprintf(stderr,"     -1.0:  Fewer than nMembers particles within 1.2 times group's .grp radius\n");
    fprintf(stderr,"     -2.0:  Density < <fThreshold> within radius enclosing no more than\n");
    fprintf(stderr,"            <nMembers> particles.\n");
    fprintf(stderr,"     -3.0:  Density > <fThreshold> beyond 3 times group's .grp radius\n");
    exit(1);
}


int main(int argc,char **argv)
{
	KD kd;

	int i,j;
	int bThreshold, bStandard, bLambda, bPeriodic, bRedshift;
	int bDark, bGas, bStar, bMark;
	int nBucket, nMembers, nSmooth, sec, usec;
	float fOmega, fLambda, fRedshift, fThreshold, fMinMass, fPeriod[3], fCenter[3];
	float fMassUnit, fMpcUnit;
	float G, H0;
	char *achGTPFile, *achListFile, *achOutFileBase, *achMarkFile;
	char achDefOutBase[3], achLongWord[128];
	time_t timeRun;
	FILE *fpOutFile;

	fprintf(stderr,"SO v1.1: Jeff Gardner, Jan 2000\n");
	strcpy(achDefOutBase,"so");
	/*
	 ** Bucket size set to 16, user cannot affect this!
	 */
   	nBucket = 16;
	/*
	 * Default tipsy native format.
	 */
	bStandard = 0;
	/* 
	 * Default to find virial density
	 */
	bThreshold = 0;
	/* 
	 * Default minimum FOF group mass to consider
	 */ 
	fMinMass = 0.0;
	/* 
	 * Default minimum particles within Rvir
	 */ 
	nMembers = 8;
	/*
	 ** Default Cosmological parameters.
	 */
	fRedshift=-9.9999;  /* This should be set later */
	bRedshift=0;
	fMassUnit=-9.9;  /* Less than 0 is the "not set" value */
	fMpcUnit=-9.9;
	fOmega = 1.0;
	fLambda = 0.0;
	bLambda = 0;
	bPeriodic = 1.0;
	for (j=0;j<3;++j) {
	    fPeriod[j] = 1.0;
	    fCenter[j] = 0.0;
	}
	/* G and H0 are fixed (they are unused) */
	G = 1.0;
	H0 = 2.8944;
	/*
	 ** Default group finding parameters (nSmooth is initial group size).
	 */
	nSmooth = 1028;
	/* 
	 * Default output params
	 */ 
	bDark=0;
	bStar=0;
	bGas=0;
	bMark=0;
	

	achGTPFile=NULL;
	achListFile=NULL;
	achOutFileBase=NULL;
	achMarkFile=NULL;
	/*
	 ** Now get the command line arguments!
	 */
	i = 1;
	while (i < argc) {
	    if (!strcmp(argv[i],"-i")) {
		++i;
		if (i >= argc) usage();
		achGTPFile = argv[i];
		++i;
	    }
	    else if (!strcmp(argv[i],"-o")) {
		++i;
		if (i >= argc) usage();
		achOutFileBase = argv[i];
		++i;
	    }
	    else if (!strcmp(argv[i],"-z")) {
		++i;
		if (i >= argc) usage();
		bRedshift = 1;
		fRedshift = atof(argv[i]);
		++i;
	    }
	    else if (!strcmp(argv[i],"-O")) {
		++i;
		if (i >= argc) usage();
		fOmega = atof(argv[i]);
		++i;
	    }
	    else if (!strcmp(argv[i],"-L")) {
		++i;
		bLambda=1;
	    }
	    else if (!strcmp(argv[i],"-s")) {
		++i;
		if (i >= argc) usage();
		nSmooth = atoi(argv[i]);
		++i;
	    }
	    else if (!strcmp(argv[i],"-rho")) {
		++i;
		if (i >= argc) usage();
		fThreshold = atof(argv[i]);
		bThreshold = 1;
		++i;
	    }
	    else if (!strcmp(argv[i],"-m")) {
		++i;
		if (i >= argc) usage();
		nMembers = atoi(argv[i]);
		++i;
	    }
	    else if (!strcmp(argv[i],"-p")) {
		++i;
		if (i >= argc) usage();
		fPeriod[0] = atof(argv[i]);
		fPeriod[1] = atof(argv[i]);
		fPeriod[2] = atof(argv[i]);
		bPeriodic = 1;
		++i;
	    }
	    else if (!strcmp(argv[i],"-c")) {
		++i;
		if (i >= argc) usage();
		fCenter[0] = atof(argv[i]);
		fCenter[1] = atof(argv[i]);
		fCenter[2] = atof(argv[i]);
		++i;
	    }
	    else if (!strcmp(argv[i],"-cx")) {
		++i;
		if (i >= argc) usage();
		fCenter[0] = atof(argv[i]);
		++i;
	    }
	    else if (!strcmp(argv[i],"-cy")) {
		++i;
		if (i >= argc) usage();
		fCenter[1] = atof(argv[i]);
		++i;
	    }
	    else if (!strcmp(argv[i],"-cz")) {
		++i;
		if (i >= argc) usage();
		fCenter[2] = atof(argv[i]);
		++i;
	    }
	    else if (!strcmp(argv[i],"-std")) {
		bStandard = 1;
		++i;
	    }
	    else if (!strcmp(argv[i],"-M")) {
		++i;
		if (i >= argc) usage();
		fMinMass = atof(argv[i]);
		++i;
	    }
	    else if (!strcmp(argv[i],"-u")) {
		++i;
		if (i >= argc) usage();
		fMassUnit = atof(argv[i]);
		++i;
		if (i >= argc) usage();
		fMpcUnit = atof(argv[i]);
		++i;
	    }
	    else if (!strcmp(argv[i],"-list")) {
		++i;
		if (i >= argc) usage();
		achListFile = argv[i];
		++i;
	    }
	    else if (!strcmp(argv[i],"-mark")) {
		++i;
		if (i >= argc) usage();
		achMarkFile = argv[i];
		bMark=1;
		++i;
	    }
	    else if (!strcmp(argv[i],"-dark")) {
		bDark = 1;
		++i;
	    }
	    else if (!strcmp(argv[i],"-gas")) {
		bGas = 1;
		++i;
	    }
	    else if (!strcmp(argv[i],"-star")) {
		bStar = 1;
		++i;
	    }
	    else if (!strcmp(argv[i],"-all")) {
		bDark = 1;
		bGas = 1;
		bStar = 1;
		++i;
	    }
	    else usage();
	}

	if (achGTPFile==NULL) usage();
	if (achOutFileBase==NULL) 
	    achOutFileBase=achDefOutBase;
      

	if (bLambda)
	    fLambda=1.0-fOmega;

	kdInit(&kd,nBucket,fPeriod,fCenter,0,nMembers,bPeriodic,bDark,bGas,bStar,bMark);
	/*
	 * Read TIPSY file to get redshift 
	 */
	i=kdReadTipsy(kd,stdin, bStandard);
	fprintf(stderr,"Read %d particles from TIPSY file.\n",i);

	/*
	 * Read in mark file if needed
	 */
	if (bMark) {
	    assert(achMarkFile != NULL);
	    i=kdReadMark(kd,achMarkFile);
	    fprintf(stderr,"%d mark particles read from %s\n",i,achMarkFile);
	}

	/* Set Redshift according to file header if not set by user */
	if (!bRedshift) {   
	    fRedshift = (1.0/kd->fTime)-1.0;
	}

	/*
	 * Find virial density in simulation units
	 */
	if (!bThreshold) {
	    fThreshold = rhovir_over_rhobar(fOmega,bLambda,fRedshift) * fOmega;
	}

	/*
	 * Output various parameters for log purposes
	 */
	sprintf(achLongWord,"%s.sovcirc",achOutFileBase);
	fpOutFile = fopen(achLongWord,"w");
	assert(fpOutFile != NULL);
	time(&timeRun);
	fprintf(fpOutFile,"#SO v1.0: Jeff Gardner, Jan 2000\n");
	fprintf(fpOutFile,"# Run on %s",ctime(&timeRun));
	fprintf(fpOutFile,"# Input .gtp file: %s\n",achGTPFile);
	if (achListFile != NULL)
	    fprintf(fpOutFile,"# Groups list from file: %s\n",achListFile);
	if (bThreshold) {
	    fprintf(fpOutFile,"# fThreshold = %g  (set by user)\n",fThreshold);
	} else {
	    fprintf(fpOutFile,"# fThreshold = %g  (VIRIAL DENSITY)\n",fThreshold);
	}
	fprintf(fpOutFile,"# fRedshift: %g   fOmega: %g   fLambda: %g\n",fRedshift,fOmega,fLambda);
	fprintf(fpOutFile,"# bPeriodic: %d  fPeriod[i]: %g %g %g   fCenter[i]: %g %g %g\n",
	      bPeriodic,fPeriod[0],fPeriod[1],fPeriod[2],fCenter[0],fCenter[1],fCenter[2]);
	fprintf(fpOutFile,"# fMinMass: %g  nMembers: %d\n",fMinMass,nMembers);
	if (fMassUnit < 0.0) {
	    fprintf(fpOutFile,"# fMassUnit: UNSPECIFIED  fMpcUnit: UNSPECIFIED\n#\n");
	} else {
	    fprintf(fpOutFile,"# fMassUnit: %g  fMpcUnit: %g\n#\n",fMassUnit,fMpcUnit);
	}
	
	
	kdSetUniverse(kd,G,fOmega,fLambda,H0,fRedshift,fMassUnit,fMpcUnit);
	kdBuildTree(kd);

	/*
	 * Get list of groups centers to gather around
	 */
	i=kdReadGTPList(kd,achGTPFile,achListFile,fMinMass,bStandard);
	fprintf(stderr,"Read %d groups to process.\n",i);

	/*
	 * SO all groups in list
	 */
	kdTime(kd,&sec,&usec);
	kdSO(kd,fThreshold,nSmooth);
	kdTime(kd,&sec,&usec);

	/*
	 * Output
	 */
	if (bDark)
	    kdWriteProfile(kd,achOutFileBase,timeRun,fpOutFile,DARK);
	if (bGas)
	    kdWriteProfile(kd,achOutFileBase,timeRun,fpOutFile,GAS);
	if (bStar)
	    kdWriteProfile(kd,achOutFileBase,timeRun,fpOutFile,STAR);
	if (bMark)
	    kdWriteProfile(kd,achOutFileBase,timeRun,fpOutFile,MARK);
	kdWriteOut(kd,fpOutFile);
	fclose(fpOutFile);

	fprintf(stderr,"SO CPU Time:");
	fprintf(stderr,"   %d.%06d\n\n",sec,usec);

}




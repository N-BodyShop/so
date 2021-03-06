/*
 * SO v1.7 (May 2003)
 *   -- Fixed bug of the handling of the following relatively rare case:
 *          Imagine SKID groups 1 and 2.  SKID 1 is less massive than SKID 2 and
 *          is therefore SO-ed first.  However, the SO virial radius of group 1 is
 *          much larger than group 2.  Therefore group 2 is actually subsumed by
 *          group 1.  Previously, this would have meant that all of the particles
 *          in group 2 would have been left at their previous group 1 tag, but
 *          group 2 would never have been zeroed in the .sovcirc file.  Now,
 *          this event is handled properly as a "slurp".
 */

/* SO v1.61 (April 2002)
 *   -- Fixed initialization bug in determination of group center of mass
 *      (only affects group velocities in .gtp output on alphas)
 */

/* SO v1.6: (May 2001)
 *   -- Added extra definition NMASSPROFILE which is the number of bins output for the 
 *      mass profiles invoked by -gas,-dark,-star,-all.  This is currently set to 16
 *      It can be different from NVCIRC (to which it was set previously) which is currently 8.
 */

/* SO v1.5:
 *   -- Improved file error handling.
 *   -- -rho option removed and replaced with -delta in light of Jenkins et al (2001) results.
 *   -- The biggest change is that the halos are now processed in increasing order of input mass.
 *      This allows smaller halos to be subsumed by larger ones.  Before now, there was no
 *      deterministic behavior the virial radii of two halos overlapped and they "shared"
 *      particles.  Now, if two halos are vying for the same particle, the following occurs:
 *          1) If the smaller-mass halo is centered *inside* the SO radius of the larger-mass halo,
 *             the larger-mass halos completely subsumes all particles of the the smaller one
 *  	       which are inside its (the larger group's) SO radius.  Particles outside
 *             the larger group's SO radius are set to group ID 0.  The SO radius of the smaller
 *	       group is set to -10.0*(larger-group).  The virial mass is multiplied by -1.
 *	    2) If the smaller-mass halo is centered *outside* th SO radius of the larger-mass halo,
 *             It retains ownership of all of its particles.  The larger halo will still use
 *             the smaller halo's particles when calculating it's own mass, radius, etc., but
 *             will not claim those particles as its own in .sogrp files.  This means that
 *             The total mass in SO groups can be greater than the total mass of particles
 *             with group ID > 0 in the .sogrp file.  This is usually not a larger effect.
 *      This change allows SKID .gtp files to be used as the input rather than FOF!
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <time.h>
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
    fprintf(stderr,"so -i <SKID .gtp file> [-o <outfilebase>] [([-dark] [-gas] [-star]) || [-all])]\n");
    fprintf(stderr,"      [-mark <markfile>]  [-std]  [-grp] [-gtp] [-subsumed] [-ignored]\n");
    fprintf(stderr,"      [-list <File containing group indexes>]\n");
    fprintf(stderr,"      [-pot || -stat <SKID .stat file containing most-bound-particle positions>]\n");
    fprintf(stderr,"      [-delta <fThreshold>] [-M <fMinGTPMass>] [-m <mMinSOMembers>]\n");
    fprintf(stderr,"      [-O <fOmega0>]  [-L]  [-z <fRedshift>]\n");
    fprintf(stderr,"      [-p <xyzPeriod>]  [-c <xyzCenter>]\n");
    fprintf(stderr,"      [-cx <xCenter>]  [-cy <yCenter>]  [-cz <zCenter>]\n");
    fprintf(stderr,"      [-u <fMassUnit> <fMpcUnit>]\n\n");

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
    fprintf(stderr,"       to a separate file <outfilebase>.somark.\n");
    fprintf(stderr,"   -std: Read and write standard TIPSY binaries.\n");
    fprintf(stderr,"   -grp -gtp: Causes .sogrp and/or .sogtp file to be written.  The group\n");
    fprintf(stderr,"       numbers given in these output files are the SAME as the group numbers\n");
    fprintf(stderr,"       in the input .gtp file.  Groups in the .sogtp file that are not examined\n");
    fprintf(stderr,"       at all will have 0 mass and 0 radius.  Groups for which SO tries and\n");
    fprintf(stderr,"       fails to find radius R will have mass 0 and '-errorcode' radius.\n");
    fprintf(stderr,"   -pot: Center on the most bound particle of the group.  Default is to center\n");
    fprintf(stderr,"       on the group position as given in the input .gtp file.  In this case,\n");
    fprintf(stderr,"       the most bound particle is the particle within the radius of the input\n");
    fprintf(stderr,"       groups that has the smallest Phi field given in the TIPSY file.\n");
    fprintf(stderr,"   -stat: Center on the most bound particle as reported by the SKID .stat file.\n");
    fprintf(stderr,"       This method is preferable to to -pot as SKID constructs potentials from\n");
    fprintf(stderr,"       only those particles in the group and is not subject to the global\n");
    fprintf(stderr,"       potential field.  NOTE: this requires SKID revisions from 2/25/00.\n");
    fprintf(stderr,"   -delta: By default, <fThreshold> is calculated automatically from cosmological\n");
    fprintf(stderr,"       parameters to be the virial density, but may be overridden using the\n");
    fprintf(stderr,"       -delta option by specifying the overdensity of the threshold.\n");
    fprintf(stderr,"   -L: option will set Lambda0 = 1-Omega0.\n");
    fprintf(stderr,"   -p: THIS VERSION ASSUMES PERIODIC BOUNDRY CONDITIONS.  Default period = 1\n");
    fprintf(stderr,"   -z: <fRedshift> will be automatically calculated as 1/h.time-1 where h.time\n");
    fprintf(stderr,"       is read from the TIPSY particle file header.  -z will override this.\n");
    fprintf(stderr,"   -M: <fMinGTPMass> specifies the minimum mass a group given in the input .gtp\n");
    fprintf(stderr,"        file must have in order to be considered.\n");
    fprintf(stderr,"   -list: names an ASCII file which lists indexes of groups in the input .gtp\n");
    fprintf(stderr,"        file to be considered.  This may be used in conjunction with -M.\n");
    fprintf(stderr,"   -m: <fMinSOMembers> specifies the minimum number of particles a group must\n");
    fprintf(stderr,"       have for its parameters to be accurately calculated.\n");
    fprintf(stderr,"   -u: Output is in system units unless simulation units are specified with -u\n");
    fprintf(stderr,"       option, whereupon output will be in the form:\n");
    fprintf(stderr,"       Mass [Msol], Distance [kpc], Velocity [km/s].\n");
    fprintf(stderr,"   -subsumed/-ignored: causes .sosub/.soign TIPSY array files to be written\n");
    fprintf(stderr,"       which give the number of times a particle was subsumed into a larger\n");
    fprintf(stderr,"       group or ignored by a larger group because it belonged to a smaller one.\n");
    fprintf(stderr,"\n  See output file header for output format.\n\n");
    fprintf(stderr,"  Groupwise error codes are returned in Mvir and Rvir columns as follows:\n");
    fprintf(stderr,"     -1.0:  Fewer than nMembers particles within 1.2 times group's .grp radius\n");
    fprintf(stderr,"     -2.0:  Density < <fThreshold> within radius enclosing no more than\n");
    fprintf(stderr,"            <nMembers> particles.\n");
    fprintf(stderr,"     -3.0:  Density > <fThreshold> beyond 3 times group's .grp radius\n");
    fprintf(stderr,"     -Mvir: Group was 'subsumed' or 'slurped' by the group 'groupID' that is given\n");
    fprintf(stderr,"            in the Rvir columns as -10*groupID.  If the Vc columns are set with\n");
    fprintf(stderr,"            meaningful numbers, then it was 'subsumed'.  If the numbers are 0,\n");
    fprintf(stderr,"            then it was 'slurped'.\n");
    fprintf(stderr,"\n'Subsuming,' 'Slurping,' and 'Retaining' behavior:\n");
    fprintf(stderr,"    The groups are SO'd in increasing order of SKID or FOF mass.  This is an attempt\n");
    fprintf(stderr,"    to produce the results whereby the most massive SO groups take precidence over\n");
    fprintf(stderr,"    the smaller SO groups during a group collision.  After SO finds the virial mass\n");
    fprintf(stderr,"    and radius of group 'A', it examines the member particles of A to see if any of them\n");
    fprintf(stderr,"    are members of other groups.  When it encounters a particle that is already a member\n");
    fprintf(stderr,"    of another group 'B' it checks as follows:\n");
    fprintf(stderr,"        1. If the center of B is within the virial radius of A:\n");
    fprintf(stderr,"                - Zero all particles that were members of B.\n");
    fprintf(stderr,"                - Mark Mvir and Rvir of B according to case 4 (-Mvir:) above.\n");
    fprintf(stderr,"                - Tag all particles in the virial radius of A as belonging to A.\n");
    fprintf(stderr,"                - Group B is now 'subsumed'.\n");
    fprintf(stderr,"        2. ELSE If the center of A is within the virial radius of B:\n");
    fprintf(stderr,"                - Zero all particles that were members of A.\n");
    fprintf(stderr,"                - Mark Mvir and Rvir of A according to case 4 (-Mvir:) above.\n");
    fprintf(stderr,"                - Tag all particles in the virial radius of B as belonging to B.\n");
    fprintf(stderr,"                - Group A is now 'slurped' by B.  Thus, slurping is an odd case where\n");
    fprintf(stderr,"                  the Rvir of the smaller SKID or FOF group is *larger* than the\n");
    fprintf(stderr,"                  Rvir of the larger SKID or FOF group.\n");
    fprintf(stderr,"                - Particles in Group A will be tabulated as 'subsumed' even though the\n");
    fprintf(stderr,"                  group itself will be tabulated as 'slurped.'\n");
    fprintf(stderr,"        2. ELSE Niether A nor B are within each other's radius:\n");
    fprintf(stderr,"                - Group B retains all of its particles.\n");
    fprintf(stderr,"                - However, all of the particles that fall within Rvir of A,\n");
    fprintf(stderr,"                  regardless of membership, are used in calculating group A's\n");
    fprintf(stderr,"                  quantities.  Thus, there will be a discrepancy in the total mass\n");
    fprintf(stderr,"                  of all groups and the total mass of all particles in groups.\n");
    fprintf(stderr,"                  This discrepancy is reported at the end of the run.\n");
    fprintf(stderr,"                - Group B particles that would have belonged to Group A are reported\n");
    fprintf(stderr,"                  as having been 'retained in the face of adversity.'\n\n");
    

    
    exit(1);
}


int main(int argc,char **argv)
{
	KD kd;

	int i,j;
	int bThreshold, bStandard, bLambda, bPeriodic, bRedshift;
	int bDark, bGas, bStar, bMark, bGrp, bGtp, bPot, bSubsumed, bIgnored;
	int nBucket, nMembers, nSmooth, sec, usec;
	float fOmega, fLambda, fRedshift, fThreshold, fMinMass, fPeriod[3], fCenter[3];
	float fMassUnit, fMpcUnit;
	float G, H0;
	char *achGTPFile, *achListFile, *achOutFileBase, *achMarkFile, *achStatFile;
	char achDefOutBase[3], achLongWord[128];
	time_t timeRun;
	FILE *fpOutFile;

	fprintf(stderr,"SO Release 1.7: Jeff Gardner, May 2003\n");
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
	bGrp=0;
	bGtp=0;
	bPot=0;
        bSubsumed=0;
        bIgnored=0;

	achGTPFile=NULL;
	achListFile=NULL;
	achOutFileBase=NULL;
	achMarkFile=NULL;
	achStatFile=NULL;
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
                fprintf(stderr,"-rho option is no longer availible.  Use -delta instead.\n");
                usage();
		++i;
	    }
	    else if (!strcmp(argv[i],"-delta")) {
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
	    else if (!strcmp(argv[i],"-grp")) {
		bGrp = 1;
		++i;
	    }
	    else if (!strcmp(argv[i],"-gtp")) {
		bGtp = 1;
		++i;
	    }
	    else if (!strcmp(argv[i],"-pot")) {
		bPot = 1;
		++i;
		if (achStatFile != NULL) usage();
	    }
	    else if (!strcmp(argv[i],"-subsumed")) {
		bSubsumed = 1;
		++i;
	    }
	    else if (!strcmp(argv[i],"-ignored")) {
		bIgnored = 1;
		++i;
	    }
	    else if (!strcmp(argv[i],"-stat")) {
		++i;
		if (i >= argc) usage();
		achStatFile = argv[i];
		++i;
		if (bPot) usage();
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

	kdInit(&kd,nBucket,fPeriod,fCenter,0,nMembers,bPeriodic,bDark,bGas,bStar,bMark,bPot);
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
	} else {  /* If set by user, convert from overdensity to density */
            fThreshold *= fOmega;
        }
        

	/*
	 * Output various parameters for log purposes
	 */
	sprintf(achLongWord,"%s.sovcirc",achOutFileBase);
	fpOutFile = fopen(achLongWord,"w");
	assert(fpOutFile != NULL);
	time(&timeRun);
	fprintf(fpOutFile,"#SO v1.61: Jeff Gardner, April 2002\n");
	fprintf(fpOutFile,"# Run on %s",ctime(&timeRun));
	fprintf(fpOutFile,"# Input .gtp file: %s\n",achGTPFile);
	if (achListFile != NULL)
	    fprintf(fpOutFile,"# Groups list from file: %s\n",achListFile);
	if (achStatFile != NULL)
	    fprintf(fpOutFile,"# Group potential centers from file: %s\n",achStatFile);
	if (bThreshold) {
	    fprintf(fpOutFile,"# fThreshold = %g  (set by user)\n",fThreshold);
	} else {
	    fprintf(fpOutFile,"# fThreshold = %g  (VIRIAL DENSITY)\n",fThreshold);
	}
	fprintf(fpOutFile,"# fRedshift: %g   fOmega: %g   fLambda: %g\n",fRedshift,fOmega,fLambda);
	fprintf(fpOutFile,"# bPeriodic: %d  fPeriod[i]: %g %g %g   fCenter[i]: %g %g %g\n",
	      bPeriodic,fPeriod[0],fPeriod[1],fPeriod[2],fCenter[0],fCenter[1],fCenter[2]);
	fprintf(fpOutFile,"# fMinMass: %g  nMembers: %d  bPot: %d\n",fMinMass,nMembers,bPot);
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
	 * Replace centers with position of most bound particles from .stat file
	 * if requested
	 */
	if (achStatFile != NULL) {
	    j=kdReadStat(kd,achStatFile);
	    fprintf(stderr,"Replaced %d group centers.\n",j);
	    if (i != j ) {
		fprintf(stderr,"ERROR in reading .stat file!\n");
		exit(1);
	    }
	}

	/*
	 * SO all groups in list
	 */
	kdTime(kd,&sec,&usec);
	kdSO(kd,fThreshold,nSmooth);
	kdTime(kd,&sec,&usec);

        /*
         * Stats
         */
        kdOutStats(kd,fpOutFile);
        

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
	if (bGrp)
	    kdWriteArray(kd,achOutFileBase);
	if (bGtp)
	    kdWriteGTP(kd,achOutFileBase,bStandard);
        if (bSubsumed)
            kdWriteConflict(kd,achOutFileBase,KD_SUBSUMED);
        if (bIgnored)
            kdWriteConflict(kd,achOutFileBase,KD_IGNORED);
        

	fprintf(stderr,"SO CPU Time:");
	fprintf(stderr,"   %d.%06d\n\n",sec,usec);

}




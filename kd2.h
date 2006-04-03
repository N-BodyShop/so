#ifndef KD_HINCLUDED
#define KD_HINCLUDED

#include "cosmo.h"
#include <sys/time.h>

#define sqr(x) ((x)*(x))

/* Number of Vcirc bins */
#define NVCIRC 8
/* Number of mass profile bins */
#define NMASSPROFILE 16

#define FLOAT float

#define ROOT		1
#define LOWER(i)	(i<<1)
#define UPPER(i)	((i<<1)+1)
#define PARENT(i)	(i>>1)
#define SIBLING(i) 	((i&1)?i-1:i+1)
#define SETNEXT(i)\
{\
	while (i&1) i=i>>1;\
	++i;\
	}

#define DARK	1
#define GAS	2
#define STAR	4
#define MARK    8

/*
 ** Softening types!
 */
#define PLUMMER		1
#define SPLINE		2

#define KD_SUBSUMED 0
#define KD_IGNORED 1

typedef struct pInitial {
	float r[3];
	float v[3];
	float fMass;
	float fPhi;
	float fTemp;
	float fBall2;
	float fDensity;
	int iOrder;
        int iGrp;
        int nSubsumed;
        int nIgnored;
	} PINIT;

typedef struct pMoved {
	float r[3];
	float rOld[3];
	float a[3];
	int iOrder;
	} PMOVE;

typedef struct pGroup {
	float rel[3];
	float rCenter[3];
	float vcm[3];
	float fMass;
	float fRadius;
	int nMembers;
	int pStart;
	int pCurr;
	} PGROUP;

typedef struct bndBound {
	float fMin[3];
	float fMax[3];
	} BND;

typedef struct kdNode {
	float fSplit;
	BND bnd;
	int iDim;
	int pLower;
	int pUpper;
	} KDN;

typedef struct grpNode {
    int index;
    float pos[3];
    float vcm[3];
    float fRgtp;    /* Radius of input GTP group */
    float fGTPMass; /* Mass of GTP input group (used for mass sort) */
    float fMvir;
    float fRvir;
    float fVcirc[NVCIRC];
    float fRmass[2];  /* quarter and half mass radii */
    float fRmax;
    float fVmax;
    float fDark[NMASSPROFILE];
    float fGas[NMASSPROFILE];
    float fStar[NMASSPROFILE];
    float fMark[NMASSPROFILE];
} GRPNODE;

typedef struct kdContext {
	int nBucket;
        int bPeriodic;
	float fPeriod[3];
	float fCenter[3];
	float G;
	CSM csm;
	float z;
        float fMassUnit;
        float fMpcUnit;
	int nParticles;
	int nDark;
	int nGas;
	int nStar;
	int inType;
	float fTime;
	int nLevels;
	int nNodes;
	int nSplit;
	int nMove;
	int nActive;
	int nInitActive;
	PINIT *pInit;
	PMOVE *pMove;
	PGROUP *pGroup;
	KDN *kdNodes;
	int nGroup;
	int *piGroup;
	int uSecond;
	int uMicro;
	int bOutDiag;
        GRPNODE *grps;
        int nGrps;
        int nMembers;        
        int bDark;
        int bGas;
        int bStar;
        int bMark;
        char *bMarkList;
        int bPot;
        int nInGTP;
    int iGroupsRemoved;
    int iGroupsSlurped;
    int iParticlesRemoved;
    int iParticlesIgnored;
    float fMassRemoved;
    float fMassIgnored;
	} * KD;


#define INTERSECT(pkdn,fBall2,lx,ly,lz,x,y,z,sx,sy,sz,bPeriodic,label)\
{\
	FLOAT INTRSCT_dx,INTRSCT_dy,INTRSCT_dz;\
	FLOAT INTRSCT_dx1,INTRSCT_dy1,INTRSCT_dz1;\
    FLOAT INTRSCT_fDist2;\
	INTRSCT_dx = (pkdn)->bnd.fMin[0]-x;\
	INTRSCT_dx1 = x-(pkdn)->bnd.fMax[0];\
	INTRSCT_dy = (pkdn)->bnd.fMin[1]-y;\
	INTRSCT_dy1 = y-(pkdn)->bnd.fMax[1];\
	INTRSCT_dz = (pkdn)->bnd.fMin[2]-z;\
	INTRSCT_dz1 = z-(pkdn)->bnd.fMax[2];\
	if (INTRSCT_dx > 0.0) {\
		INTRSCT_dx1 += lx;\
		if (INTRSCT_dx1 < INTRSCT_dx) {\
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
			sx = x+lx;\
			bPeriodic = 1;\
			}\
		else {\
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;\
			sx = x;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto label;\
		}\
	else if (INTRSCT_dx1 > 0.0) {\
		INTRSCT_dx += lx;\
		if (INTRSCT_dx < INTRSCT_dx1) {\
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;\
			sx = x-lx;\
			bPeriodic = 1;\
			}\
		else {\
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
			sx = x;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto label;\
		}\
	else {\
		INTRSCT_fDist2 = 0.0;\
		sx = x;\
		}\
	if (INTRSCT_dy > 0.0) {\
		INTRSCT_dy1 += ly;\
		if (INTRSCT_dy1 < INTRSCT_dy) {\
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
			sy = y+ly;\
			bPeriodic = 1;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;\
			sy = y;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto label;\
		}\
	else if (INTRSCT_dy1 > 0.0) {\
		INTRSCT_dy += ly;\
		if (INTRSCT_dy < INTRSCT_dy1) {\
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;\
			sy = y-ly;\
			bPeriodic = 1;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
			sy = y;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto label;\
		}\
	else {\
		sy = y;\
		}\
	if (INTRSCT_dz > 0.0) {\
		INTRSCT_dz1 += lz;\
		if (INTRSCT_dz1 < INTRSCT_dz) {\
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
			sz = z+lz;\
			bPeriodic = 1;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;\
			sz = z;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto label;\
		}\
	else if (INTRSCT_dz1 > 0.0) {\
		INTRSCT_dz += lz;\
		if (INTRSCT_dz < INTRSCT_dz1) {\
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;\
			sz = z-lz;\
			bPeriodic = 1;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
			sz = z;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto label;\
		}\
	else {\
		sz = z;\
		}\
	}




void kdTime(KD,int *,int *);
int kdInit(KD *,int,float *,float *,int,int,int,int,int,int,int,int);
void kdSetUniverse(KD,float,float,float,float,float,float,float);
int kdParticleType(KD,int);
int kdReadMark(KD, char *);
int kdReadStat(KD, char *);
int kdReadGTPList(KD, char *, char *, float, int);
int kdReadTipsy(KD,FILE *,int);
void kdSO(KD, float, int);
void kdWriteProfile(KD, char *, time_t, FILE *, int);
void kdWriteOut(KD, FILE *);
int kdBuildTree(KD);
void kdFinish(KD);
void kdWriteConflict(KD kd, char *achOutFileBase, int iOpt);
void kdOutStats(KD kd, FILE *fpOutFile);
void kdWriteArray(KD, char *);
void kdWriteGTP(KD, char *, int);

void indexx(int n, float arr[], int indx[]);

#endif



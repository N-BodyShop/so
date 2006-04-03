#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <sys/resource.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <assert.h>
#include <limits.h>
#include "smooth2.h"
#include "tipsydefs.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

int xdrHeader(XDR *pxdrs,struct dump *ph)
{
	int pad = 0;
	
	if (!xdr_double(pxdrs,&ph->time)) return 0;
	if (!xdr_int(pxdrs,&ph->nbodies)) return 0;
	if (!xdr_int(pxdrs,&ph->ndim)) return 0;
	if (!xdr_int(pxdrs,&ph->nsph)) return 0;
	if (!xdr_int(pxdrs,&ph->ndark)) return 0;
	if (!xdr_int(pxdrs,&ph->nstar)) return 0;
	if (!xdr_int(pxdrs,&pad)) return 0;
	return 1;
	}

void kdTime(KD kd,int *puSecond,int *puMicro)
{
	struct rusage ru;

	getrusage(0,&ru);
	*puMicro = ru.ru_utime.tv_usec - kd->uMicro;
	*puSecond = ru.ru_utime.tv_sec - kd->uSecond;
	if (*puMicro < 0) {
		*puMicro += 1000000;
		*puSecond -= 1;
		}
	kd->uSecond = ru.ru_utime.tv_sec;
	kd->uMicro = ru.ru_utime.tv_usec;
	}


int kdInit(KD *pkd,int nBucket,float *fPeriod,float *fCenter,int bOutDiag,
	   int nMembers, int bPeriodic, int bDark, int bGas, int bStar,
	   int bMark, int bPot)
{
	KD kd;
	int j;

	kd = (KD)malloc(sizeof(struct kdContext));
	assert(kd != NULL);
	kd->nBucket = nBucket;
	for (j=0;j<3;++j) {
		kd->fPeriod[j] = fPeriod[j];
		kd->fCenter[j] = fCenter[j];
		}
	kd->bOutDiag = bOutDiag;
	kd->G = 1.0;
	kd->csm = NULL;
	kd->z = 0.0;
	kd->nParticles = 0;
	kd->nDark = 0;
	kd->nGas = 0;
	kd->nStar = 0;
	kd->inType = 0;
	kd->fTime = 0.0;
	kd->nLevels = 0;
	kd->nNodes = 0;
	kd->nSplit = 0;
	kd->nMove = 0;
	kd->nActive = 0;
	kd->nInitActive = 0;
	kd->pMove = NULL;
	kd->pInit = NULL;
	/*	kd->pGroup = NULL; */
	kd->kdNodes = NULL;
	kd->piGroup = NULL;
	kd->nGroup = 0;

	kd->nMembers= nMembers;
	kd->grps = NULL;
	kd->nGrps = 0;
	*pkd = kd;
	kd->bPeriodic = bPeriodic;
        kd->bDark = bDark;
        kd->bGas = bGas;
        kd->bStar = bStar;
        kd->bMark = bMark;
	kd->bMarkList = NULL;
	kd->bPot = bPot;
	return(1);
	}


void kdSetUniverse(KD kd, float G, float Omega0, float Lambda, float H0, float z,
		   float fMassUnit, float fMpcUnit)
{
	kd->G = G;
	csmInitialize(&kd->csm);
	kd->csm->dOmega0 = Omega0;
	kd->csm->dLambda = Lambda;
	kd->csm->dHubble0 = H0;
	/*
	 * XXX bComove is never used below, nor is there a command line
	 * flag to set it.  Caveat Utor.
	 */
	if(H0 != 0.0) kd->csm->bComove = 1;
	kd->z = z;
	kd->fMassUnit=fMassUnit;
	kd->fMpcUnit=fMpcUnit;
	}


int kdParticleType(KD kd,int iOrder)
{
	if (iOrder < kd->nGas) return(GAS);
	else if (iOrder < (kd->nGas + kd->nDark)) return(DARK);
	else if (iOrder < kd->nParticles) return(STAR);
	else return(0);
	}


int kdReadMark(KD kd, char *achMarkFile)
{
    int i,j,k,nmark;
    FILE *in;

    in = fopen(achMarkFile,"r");
    assert(in != NULL);

    assert(kd->bMarkList == NULL);
    kd->bMarkList = (char *) malloc(kd->nParticles*sizeof(char));
    assert(kd->bMarkList != NULL);
    for(i=0;i<kd->nParticles;++i)
	kd->bMarkList[i]=0;

    fscanf(in,"%d %d %d",&i,&j,&k);  /* header */
    nmark=0;
    while(fscanf(in,"%d",&i) != EOF) {
	--i;  /* Mark file indexing begins at 1 */
	assert(i < kd->nParticles);
	kd->bMarkList[i]=1;
	++nmark;
    }
    
    fclose(in);
    return(nmark);
}
    

int kdReadGTPList(KD kd, char *achGTPFile, char *achListFile,
		     float fMinMass, int bStandard)
{
	int i,j,k;
	int fofgrpnum;
	struct dump h;
	struct star_particle pp,*GTP_particles;
	int *fofList,fofAllocSize=1028,nFof=0;
	    
	FILE *fp;
	XDR xdrs;


	/* 
	 * Read in file with list of GTP group numbers
	 */
	if (achListFile != NULL) {
	    fofList=(int *) malloc(fofAllocSize*sizeof(int));
	    assert(fofList != NULL);
	    fp = fopen(achListFile,"r");
	    assert(fp != NULL);
	    while(fscanf(fp,"%d",&fofgrpnum) != EOF) {
		++nFof;
		if (nFof > fofAllocSize) {
		    fofAllocSize *= 2;
		    fofList = (int *) realloc(fofList, fofAllocSize*sizeof(int));
		    assert(fofList != NULL);
		}
		fofList[nFof-1]=fofgrpnum;
	    }
	    fclose(fp);
	}
		
	/*
	 * Read in GTP file
	 */
	fp = fopen(achGTPFile,"rb");
	if (bStandard) {
	    assert(sizeof(Real)==sizeof(float)); /* Otherwise, this XDR stuff
						    ain't gonna work */
	    xdrstdio_create(&xdrs, fp, XDR_DECODE);
	    xdrHeader(&xdrs,&h);
	}
	else {
	    fread(&h,sizeof(struct dump),1,fp);
	}
	if ( h.ndark > 0 || h.nsph > 0 ) {
	    fprintf(stderr," FILE TYPE MISMATCH: GTP file contains non-star particles!\n");
	    exit(1);
	}
	GTP_particles = (struct star_particle *) 
	                       malloc(h.nstar*sizeof(struct star_particle));
	assert(GTP_particles != NULL);
	if (bStandard) {
	    for (i=0;i<h.nstar;++i) {
		xdr_vector(&xdrs, (char *) &GTP_particles[i], 11,
			   sizeof(Real), (xdrproc_t) xdr_float);
	    }
	    xdr_destroy(&xdrs);
	}
	else {
	    for (i=0;i<h.nstar;++i) {
		fread(&GTP_particles[i],sizeof(struct star_particle),1,fp);
	    }
	}
	fclose(fp);

	/* 
	 * Find centers of all relevant groups
	 */
	if (nFof) {
	    kd->grps = (GRPNODE *) malloc(nFof*sizeof(GRPNODE));
	    assert(kd->grps!=NULL);
	    for (i=0,k=0;i<nFof;++i) {
		if (GTP_particles[fofList[i]-1].mass >= fMinMass) {
		    for (j=0;j<3;++j) {
			kd->grps[k].pos[j] = GTP_particles[fofList[i]-1].pos[j];
		    }
		    kd->grps[k].fRgtp = GTP_particles[fofList[i]-1].eps;
		    kd->grps[k].index = fofList[i];
		    ++k;
		}
	    }
	    assert(k<=nFof);
	    free(fofList);
	} else {
	    kd->grps = (GRPNODE *) malloc(h.nstar*sizeof(GRPNODE));
	    assert(kd->grps!=NULL);
	    for (i=0,k=0;i<h.nstar;++i) {
		if (GTP_particles[i].mass >= fMinMass) {
		    for (j=0;j<3;++j) {
			kd->grps[k].pos[j] = GTP_particles[i].pos[j];
		    }
		    kd->grps[k].fRgtp = GTP_particles[i].eps;
		    kd->grps[k].index = i+1;
		    ++k;
		}
	    }
	    assert(k<=h.nstar);
	}    
	kd->nGrps = k;
	kd->nInGTP = h.nstar;
	free(GTP_particles);
	return(k);
}	


int kdReadStat(KD kd, char *achStatFile)
{
    int i,j,k,itemp,grpnum;
    float ftemp, r[3];
    FILE *fp;

    assert(achStatFile != NULL);
    fp = fopen(achStatFile,"r");
    assert(fp != NULL);
    k=0;
    while(fscanf(fp,"%d %d",&grpnum,&itemp) != EOF) {
	for (i=0;i<16;++i) {
	    fscanf(fp,"%f",&ftemp);
	}
	fscanf(fp,"%f %f %f",&r[0],&r[1],&r[2]);
	if (grpnum == kd->grps[k].index) {   /* Replace center */
	    /*fprintf(stderr,"Replaced grp %d pos: %g %g %g  to  pos: %g %g %g\n",
		    kd->grps[k].index,kd->grps[k].pos[0],kd->grps[k].pos[1],
		    kd->grps[k].pos[2],r[0],r[1],r[2]);*/
	    for (j=0;j<3;++j) {
		kd->grps[k].pos[j] = r[j];
	    }
	    ++k;
	}
    }
    fclose(fp);
    return(k);
}


int kdReadTipsy(KD kd,FILE *fp, int bStandard)
{
	PINIT *p;
	int i,j;
	struct dump h;
	struct gas_particle gp;
	struct dark_particle dp;
	struct star_particle sp;
	XDR xdrs;

	if (kd->bOutDiag) puts(">> kdReadTipsy()");
	fflush(stdout);
	if (bStandard) {
	    assert(sizeof(Real)==sizeof(float)); /* Otherwise, this XDR stuff
						    ain't gonna work */
	    xdrstdio_create(&xdrs, fp, XDR_DECODE);
	    xdrHeader(&xdrs,&h);
	    }
	else {
	    fread(&h,sizeof(struct dump),1,fp);
	    }
	kd->inType = 0;
	kd->nDark = h.ndark;
	if (kd->nDark) kd->inType |= DARK;
	kd->nGas = h.nsph;
	if (kd->nGas) kd->inType |= GAS;
	kd->nStar = h.nstar;
	if (kd->nStar) kd->inType |= STAR;
	kd->fTime = h.time;
	kd->nParticles = kd->nDark + kd->nGas + kd->nStar;
	kd->nInitActive = kd->nParticles;
	/*
	 ** Allocate arrays.
	 */
	kd->pInit = (PINIT *)malloc(kd->nParticles*sizeof(PINIT));
	assert(kd->pInit != NULL);
	fprintf(stderr,"nDark:%d nGas:%d nStar:%d\n",kd->nDark,kd->nGas,kd->nStar);
	fflush(stdout);
	p = kd->pInit;
	/*
	 ** Read Stuff!
	 */
	for (i=0;i<kd->nParticles;++i) {
		p[i].iOrder = i;
		p[i].iGrp = 0;
		p[i].fDensity = 0.0;
		switch (kdParticleType(kd,i)) {
		case (GAS):
			if (bStandard) {
			    xdr_vector(&xdrs, (char *) &gp, 12,
				       sizeof(Real), (xdrproc_t) xdr_float);
			    }
			else {
			    fread(&gp,sizeof(struct gas_particle),1,fp);
			    }
			p[i].fMass = gp.mass;
			p[i].fPhi = gp.phi;
			p[i].fTemp = gp.temp;
			for (j=0;j<3;++j) {
				p[i].r[j] = gp.pos[j];
				p[i].v[j] = gp.vel[j];
				}
			break;
		case (DARK):
			if (bStandard) {
			    xdr_vector(&xdrs, (char *) &dp, 9,
				       sizeof(Real), (xdrproc_t) xdr_float);
			    }
			else {
			    fread(&dp,sizeof(struct dark_particle),1,fp);
			    }
			p[i].fMass = dp.mass;
			p[i].fPhi = dp.phi;
			p[i].fTemp = 0.0;
			for (j=0;j<3;++j) {
				p[i].r[j] = dp.pos[j];
				p[i].v[j] = dp.vel[j];
				}
			break;
		case (STAR):
			if (bStandard) {
			    xdr_vector(&xdrs, (char *) &sp, 11,
				       sizeof(Real), (xdrproc_t) xdr_float);
			    }
			else {
			    fread(&sp,sizeof(struct star_particle),1,fp);
			    }
			p[i].fMass = sp.mass;
			p[i].fPhi = sp.phi;
			p[i].fTemp = 0.0;
			for (j=0;j<3;++j) {
				p[i].r[j] = sp.pos[j];
				p[i].v[j] = sp.vel[j];
				}
			break;
			}
		}
	if (bStandard) xdr_destroy(&xdrs);
	if (kd->bOutDiag) puts("<< kdReadTipsy()");
	fflush(stdout);
	return(kd->nParticles);
	}



int CmpList(const void *p1,const void *p2)			
{
    NN *a = (NN *)p1;
    NN *b = (NN *)p2;
    
    if (a->fDist2 < b->fDist2)
	return(-1);
    if (a->fDist2 > b->fDist2)
	return(1);
    return(0);
}

void addProfileMass(GRPNODE *grp, int i, float mass, int ptype)
{
	/* Now we have all masses */
	switch (ptype) {
	case (DARK):
	    grp->fDark[i]=mass;
	    break;
	case (GAS):
	    grp->fGas[i]=mass;
	    break;
	case (STAR):
	    grp->fStar[i]=mass;
	    break;
	case (MARK):
	    grp->fMark[i]=mass;
	    break;
	    assert(0);
	}
}


void kdMassProfile(KD kd, SMX smx, GRPNODE *grp, float rvir, int nParticles, int ptype)
{
    int i,j;
    float f,fmin,r,r2,mass;

    fmin=2.0/NVCIRC;  /* NVCIRC is number if Vcirc bins...this sets increments */
    j=0;
    mass=0.0;
    for(f=fmin,i=0;i<NVCIRC-1;++i,f+=fmin) {  /* f is fraction of Rvir */
	r=f*rvir;
	r2=r*r;
    	while(smx->nnList[j].fDist2 < r2 && j < nParticles) {
	    if (ptype == MARK) {
		if (kd->bMarkList[smx->nnList[j].pInit->iOrder])
		    mass += smx->nnList[j].pInit->fMass;
	    } else {
		if(kdParticleType(kd,smx->nnList[j].pInit->iOrder) == ptype)
		    mass += smx->nnList[j].pInit->fMass;
	    }
	    ++j;
	}
	
	addProfileMass(grp,i,mass,ptype);
    }

    /* For the last bin, just add up the remaining particles */
    while(j<nParticles) {
	if (ptype == MARK) {
	    if (kd->bMarkList[smx->nnList[j].pInit->iOrder])
		mass += smx->nnList[j].pInit->fMass;
	} else {
	    if(kdParticleType(kd,smx->nnList[j].pInit->iOrder) == ptype)
		mass += smx->nnList[j].pInit->fMass;
	}
	++j;
    }
    addProfileMass(grp,NVCIRC-1,mass,ptype);

}

float kdVcirc(KD kd, SMX smx, GRPNODE *grp)
{
    int i,j,jlast,nParticles;
    float fBall,fBall2,rvir,mvir,mass,m,r,r2,vm,rm,vc;
    float fmin,f;

    /*
     * Find circular velocities at (0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0)Rvir
     */
    
    fmin=2.0/NVCIRC;  /* NVCIRC is number if Vcirc bins...this sets increments */
    rvir = grp->fRvir;
    mvir = grp->fMvir;
    fBall=2.*rvir;
    fBall2=fBall*fBall;
    nParticles = smBallGather(smx,fBall2,grp->pos);  /* Find particles within 2 Rvir */
    qsort(smx->nnList,nParticles,sizeof(NN),CmpList);
    j=0;
    mass=0.0;
    for(f=fmin,i=0;i<NVCIRC-1;++i,f+=fmin) {  /* f is fraction of Rvir */
	r=f*rvir;
	r2=r*r;
	while(smx->nnList[j].fDist2 < r2 && j < nParticles) {
	    mass += smx->nnList[j].pInit->fMass;
	    ++j;
	}
	/* Now we have mass enclosed */
	grp->fVcirc[i]=sqrt(kd->G * mass / r);
    }
    /* For the last bin, just add up the remaining particles */
    while(j<nParticles) {
	mass += smx->nnList[j].pInit->fMass;
	++j;
    }
    grp->fVcirc[NVCIRC-1]=sqrt(kd->G * mass / fBall);
    
    /*
     * Find 0.25 Mvir and 0.5 Mvir radii
     */
    for(f=0.25,i=0;i<2;++i,f+=0.25) {
	m=f*mvir;
	j=0;
	mass=smx->nnList[0].pInit->fMass;
	while(mass < m) {
	    ++j;
	    mass += smx->nnList[j].pInit->fMass;
	}
	grp->fRmass[i]=sqrt(smx->nnList[j].fDist2);
    }	
    
    /*
     * Find Vcirc_max
     */
    vm=0;
    rm=0;
    mass=0.;
    for(j=0;j<kd->nMembers;++j) {
	mass += smx->nnList[j].pInit->fMass;
    }
    rm = sqrt(smx->nnList[kd->nMembers-1].fDist2);
    vm = sqrt(kd->G * mass / rm);
    for(j=kd->nMembers;j<nParticles;++j) {
	mass += smx->nnList[j].pInit->fMass;
	r = sqrt(smx->nnList[j].fDist2);
	vc = sqrt(kd->G * mass / r);
	if (vc > vm) {
	    vm = vc;
	    rm = r;
	}
    }
    grp->fRmax = rm;
    grp->fVmax = vm;

    /*
     * Find radial profiles of other things specified by user
     */

    if (kd->bDark)
	kdMassProfile(kd,smx,grp,rvir,nParticles,DARK);
    if (kd->bGas)
	kdMassProfile(kd,smx,grp,rvir,nParticles,GAS);
    if (kd->bStar)
	kdMassProfile(kd,smx,grp,rvir,nParticles,STAR);
    if (kd->bMark)
	kdMassProfile(kd,smx,grp,rvir,nParticles,MARK);


    return(vm);
}
    
float rhoEnclosed(float mass, float r2)
{
    float r3;
    r3 = r2*sqrt(r2);
    return mass / (1.33333333*M_PI*r3);  /* (1.333333=4/3) */
}

    

float kdRvir(KD kd, SMX smx, float fRhoVir, GRPNODE *grp)
{
    int i,j,k,l,jlast,nParticles;
    float mass,fBall,fBall2,r,r3,rho,minPot,ballCtr[3];
    
    /*
     *  Beginning at 1.2 Rgtp, gather all particles and sort them in order
     *  of ascending radius.  Find the rho enclosed by kd->nMembers particles,
     *  then incrementally add additional particles until rho is less than
     *  fRhoVir for two consecutive particles.  If this condition is not met
     *  by 1.2 Rgtp, increase radius by 20% and redo BallSearch.
     *  During these additional iterations, start with the accumulations from
     *  the particles at the last tested radius (oohhh...aaahhhh....efficiency!).
     *  The virial radius is defined as the radius at which the mass enclosed
     *  has density fRhoVir.
     *  Errors are returned if nParticles within 1.2 Rgtp is less than nMembers,
     *  if rho < rhovir with only nMembers particles, or if rho never < rhovir within
     *  3 Rgtp.
     */
    jlast=0;
    mass = 0.0;
    fBall = grp->fRgtp;                  /* Start at 1.2 Rvir */
    for (k=0;k<3;++k) ballCtr[k] = grp->pos[k];  /* Always do first search on GRP center */
    /*fprintf(stderr,"In kdRvir: index: %d  fRgtp: %g  pos: %g %g %g\n",grp->index,
      grp->fRgtp,grp->pos[0],grp->pos[1],grp->pos[2]);*/

    /* Find most bound particle if applicable */
    if (kd->bPot) {              
	fBall2 = fBall*fBall;
	nParticles = smBallGather(smx,fBall2,ballCtr);
	/*fprintf(stderr,"Old center: %g %g %g\n",ballCtr[0],ballCtr[1],ballCtr[2]);*/
	minPot = smx->nnList[0].pInit->fPhi;
	for (k=0;k<3;k++) ballCtr[k] = smx->nnList[0].pInit->r[k];
	for (i=1;i<nParticles;++i) {
	    if (smx->nnList[i].pInit->fPhi < minPot) {
		/*if (grp->index==220) {
		    fprintf(stderr,"minPot0: %g  minPot1: %g  pos: %g %g %g\n",minPot,
			    smx->nnList[i].pInit->fPhi,smx->nnList[i].pInit->r[0],
			    smx->nnList[i].pInit->r[1],smx->nnList[i].pInit->r[2]);
		    minPot = smx->nnList[i].pInit->fPhi;
		    }*/
		for (k=0;k<3;k++) ballCtr[k] = smx->nnList[i].pInit->r[k];
	    }		    
	}
	/*fprintf(stderr,"%d: %g %g %g\n",grp->index,ballCtr[0],ballCtr[1],ballCtr[2]);*/
    }
    

    while(fBall < (3.0 * grp->fRgtp)) {  /* Declare defeat by 3 Rgtp */
	fBall *= 1.2;                           /* Increase by 20% each iteration */
	fBall2 = fBall*fBall;
	nParticles = smBallGather(smx,fBall2,ballCtr);
	
	if (!jlast) {                    /* First time through */
	    if (nParticles < kd->nMembers) {   /* Declare failure if not enough particles */
		grp->fRvir=-1.0;
		grp->fMvir=-1.0;
		return(-1.0);   
	    }
	}
		
	/* Sort according to distance from center */
	/*fprintf(stderr,"fBall: %g  fBall2: %g  nParticles: %d\n",fBall,fBall2,nParticles);*/
	qsort(smx->nnList,nParticles,sizeof(NN),CmpList);
	/*for (i=0;i<8;++i) {
	    fprintf(stderr,"%g  ",smx->nnList[i-1].fDist2);
	    }*/

	if (!jlast) {                          /* First time through */
	    for(j=0;j<kd->nMembers;++j) {
		mass += smx->nnList[j].pInit->fMass;
	    }
	    /*r3 = smx->nnList[j-1].fDist2*sqrt(smx->nnList[j-1].fDist2);
	      rho = mass / ((4./3.)*M_PI*r3); */
	    /* fprintf(stderr,"\n%d: %g  %g\n,",j-1,r3,rho); */
            /* Declare failure if already below density for this particle *and* the next one */
	    if (rhoEnclosed(mass,smx->nnList[j-1].fDist2) < fRhoVir &&
		rhoEnclosed(mass+smx->nnList[j].pInit->fMass,smx->nnList[j].fDist2) < fRhoVir ) {
		grp->fRvir=-2.0;
		grp->fMvir=-2.0;
		return(-2.0);   
	    }
	    jlast=j;
	}
	for(j=jlast;j<nParticles-1;j++) {
	    mass += smx->nnList[j].pInit->fMass;
	    /*r = sqrt(smx->nnList[j].fDist2);
	    r3 = r*smx->nnList[j].fDist2;
	    rho = mass / ((4./3.)*M_PI*r3);*/
	    /* Again, virial radius criterion is that density for this and next
	     * Next particle is below fRhoVir */
	    if(rhoEnclosed(mass,smx->nnList[j].fDist2) < fRhoVir &&
	       rhoEnclosed(mass+smx->nnList[j+1].pInit->fMass,smx->nnList[j+1].fDist2) < fRhoVir) {
		mass -= smx->nnList[j].pInit->fMass;  /* now "mass" is mass *within* Rvir */
		r3=mass / ((4./3.)*M_PI*fRhoVir); /* Now we know Mvir and RhoVir -> find Rvir  */
		r=pow(r3,0.3333333333);
		grp->fRvir=r;
		grp->fMvir=mass;
		/* Tag all particles that are within fRvir and find grp mean velocity */
		for (k=0;k<j;++k) {
		    smx->nnList[k].pInit->iGrp=grp->index;
		    /*fprintf(stderr,"  iGrp1: %d   iOrder: %d  pos: %g %g %g\n",
			    smx->nnList[k].pInit->iGrp,smx->nnList[k].pInit->iOrder,
			    smx->nnList[k].pInit->r[0],smx->nnList[k].pInit->r[1],
			    smx->nnList[k].pInit->r[2]);*/
		    for (l=0;l<3;l++)
			grp->vcm[l] += smx->nnList[k].pInit->fMass*smx->nnList[k].pInit->v[l];
		}
		for (l=0;l<3;l++)
		    grp->vcm[l] /= mass;

		return(r);
	    }
	}
	jlast=j;     /* Rvir not found...redo loop with larger ball size */
    }

    grp->fRvir=-r;
    grp->fMvir=-mass;
    return(-3.0);
}
	
void kdSO(KD kd, float rhovir, int nSmooth)
{
    int i;
    float fBall2, rvir;
    SMX smx;

    smInit(&smx,kd,nSmooth);
    for (i=0;i<kd->nGrps;++i) {
	/*
	 * Find virial radius
	 */
	rvir=kdRvir(kd,smx,rhovir,&kd->grps[i]);
	/*
	 * Find circular velocities
	 */ 
	if (rvir > 0.0) {
	    kdVcirc(kd,smx,&kd->grps[i]);
	}
    }
    smFinish(smx);
}

	    

#define GRAV       6.6726e-8   /* G in cgs */

void kdWriteProfile(KD kd, char *achOutFileBase,time_t timeRun, FILE *fpOutFile, int ptype)
{
    int i,j;
    float massunit;
    char pstring[5],longword[128];
    GRPNODE *gg;
    FILE *outfile;

    switch (ptype) {
    case (DARK):
	assert(kd->bDark);
	sprintf(longword,"%s.sodark",achOutFileBase);
	strcpy(pstring,"dark");
	break;
    case (GAS):
	assert(kd->bGas);
	sprintf(longword,"%s.sogas",achOutFileBase);
	strcpy(pstring,"gas");
	break;
    case (STAR):
	assert(kd->bStar);
	sprintf(longword,"%s.sostar",achOutFileBase);
	strcpy(pstring,"star");
	break;
    case (MARK):
	assert(kd->bMark);
	sprintf(longword,"%s.somark",achOutFileBase);
	strcpy(pstring,"marked");
	break;
	assert(0);
    }
	
    outfile = fopen(longword,"w");
    assert(outfile != NULL);
    fprintf(fpOutFile,"# Radial mass profile for %s particles written to %s\n",pstring,longword);

    if (kd->fMassUnit < 0.0) {  /* If units not set, do not convert */
	massunit = 1.0;
    } else {
	massunit = kd-> fMassUnit;
    }

    fprintf(outfile,"# Radial mass profile for %s particles\n",pstring);
    fprintf(outfile,"# Run on %s",ctime(&timeRun));
    fprintf(outfile,"# grp# Mass(R = %4.2f ... 2 Rvir)\n",2.0/NVCIRC);
    for(gg=kd->grps,i=0;i<kd->nGrps;i++,gg++) {
	fprintf(outfile,"%d ",gg->index);
	for(j=0;j<NVCIRC;++j) {
	    switch (ptype) {
	    case (DARK):
		fprintf(outfile,"%g ",gg->fDark[j]*massunit);
		break;
	    case (GAS):
		fprintf(outfile,"%g ",gg->fGas[j]*massunit);
		break;
	    case (STAR):
		fprintf(outfile,"%g ",gg->fStar[j]*massunit);
		break;
	    case (MARK):
		fprintf(outfile,"%g ",gg->fMark[j]*massunit);
		break;
		assert(0);
	    }
	}
	fprintf(outfile,"\n");
    }
    fclose(outfile);
}

void kdWriteOut(KD kd, FILE *fpOutFile)
{
    int i,j;
    float kmsecunit,massunit,kpcunit;
    double dtemp;
    GRPNODE *gg;
    

    /* 
     * Find conversion for km/sec
     */ 
    if (kd->fMassUnit < 0.0) {  /* If units not set, do not convert */
	kmsecunit = 1.0;
	kpcunit = 1.0;
	massunit = 1.0;
    } else {
	dtemp = GRAV * kd->fMassUnit * (1.0+kd->z) / kd->fMpcUnit;
	dtemp = 25388.8 * sqrt(dtemp) / 100000.0;
	kmsecunit = (float) dtemp;
	kpcunit = kd->fMpcUnit*1000.0;
	massunit = kd-> fMassUnit;
    }

    fprintf(fpOutFile,"#\n# grp# Mvir Rvir R(0.25Mvir) R(0.5Mvir)  R(Vc_max)  Vc_max  Vc(R = %4.2f ... 2 Rvir)\n",
	   2.0/NVCIRC);
    for (i=0,gg=kd->grps;i<kd->nGrps;++i,++gg) {
	if (gg->fMvir < 0.0) {  /* Error condition */
	    fprintf(fpOutFile,"%i %g %g ",gg->index,gg->fMvir,gg->fRvir);
	} else {
	    fprintf(fpOutFile,"%i %g %g ",gg->index,gg->fMvir*massunit,gg->fRvir*kpcunit);
	}
	fprintf(fpOutFile,"%g %g %g %g ",gg->fRmass[0]*kpcunit,gg->fRmass[1]*kpcunit,gg->fRmax*kpcunit,
	       gg->fVmax*kmsecunit);
	for(j=0;j<NVCIRC;++j) {
	    fprintf(fpOutFile,"%g ",gg->fVcirc[j]*kmsecunit);
	}
	fprintf(fpOutFile,"\n");
    }
}




void kdSelectInit(KD kd,int d,int k,int l,int r)
{
	PINIT *p,t;
	double v;
	int i,j;

	p = kd->pInit;
	while (r > l) {
		v = p[k].r[d];
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (p[++i].r[d] >= v) break;
			while (i < j) if (p[--j].r[d] <= v) break;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			if (j <= i) break;
			}
		p[j] = p[i];
		p[i] = p[r];
		p[r] = t;
		if (i >= k) r = i - 1;
		if (i <= k) l = i + 1;
		}
	}


void Combine(KDN *p1,KDN *p2,KDN *pOut)
{
	int j;

	/*
	 ** Combine the bounds.
	 */
	for (j=0;j<3;++j) {
		if (p2->bnd.fMin[j] < p1->bnd.fMin[j])
			pOut->bnd.fMin[j] = p2->bnd.fMin[j];
		else
			pOut->bnd.fMin[j] = p1->bnd.fMin[j];
		if (p2->bnd.fMax[j] > p1->bnd.fMax[j])
			pOut->bnd.fMax[j] = p2->bnd.fMax[j];
		else
			pOut->bnd.fMax[j] = p1->bnd.fMax[j];
		}
	}


void UpPassInit(KD kd,int iCell)
{
	KDN *c;
	int l,u,pj,j;

	c = kd->kdNodes;
	if (c[iCell].iDim != -1) {
		l = LOWER(iCell);
		u = UPPER(iCell);
		UpPassInit(kd,l);
		UpPassInit(kd,u);
		Combine(&c[l],&c[u],&c[iCell]);
		}
	else {
		l = c[iCell].pLower;
		u = c[iCell].pUpper;
		for (j=0;j<3;++j) {
			c[iCell].bnd.fMin[j] = kd->pInit[u].r[j];
			c[iCell].bnd.fMax[j] = kd->pInit[u].r[j];
			}
		for (pj=l;pj<u;++pj) {
			for (j=0;j<3;++j) {
				if (kd->pInit[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = kd->pInit[pj].r[j];
				if (kd->pInit[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = kd->pInit[pj].r[j];
				}
			}
		}
	}


int kdBuildTree(KD kd)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	BND bnd;

	if (kd->bOutDiag) puts(">> kdBuildTree()");
	fflush(stdout);
	if (kd->nInitActive == 0) {
		if (kd->kdNodes) free(kd->kdNodes);
		kd->kdNodes = NULL;
		return(1);
		}
	assert(kd->nInitActive > 0);
	n = kd->nInitActive;
	kd->nLevels = 1;
	l = 1;
	while (n > kd->nBucket) {
		n = n>>1;
		l = l<<1;
		++kd->nLevels;
		}
	kd->nSplit = l;
	kd->nNodes = l<<1;
	if (kd->kdNodes) {
		free(kd->kdNodes);
		kd->kdNodes = NULL;
		}
	kd->kdNodes = (KDN *)malloc(kd->nNodes*sizeof(KDN));
	assert(kd->kdNodes != NULL);
	/*
	 ** Calculate Bounds.
	 */
	for (j=0;j<3;++j) {
		bnd.fMin[j] = kd->pInit[0].r[j];
		bnd.fMax[j] = kd->pInit[0].r[j];
		}
	for (i=1;i<kd->nInitActive;++i) {
		for (j=0;j<3;++j) {
			if (bnd.fMin[j] > kd->pInit[i].r[j]) 
				bnd.fMin[j] = kd->pInit[i].r[j];
			else if (bnd.fMax[j] < kd->pInit[i].r[j])
				bnd.fMax[j] = kd->pInit[i].r[j];
			}
		}
	/*
	 ** Set up ROOT node
	 */
	c = kd->kdNodes;
	c[ROOT].pLower = 0;
	c[ROOT].pUpper = kd->nInitActive-1;
	c[ROOT].bnd = bnd;
	i = ROOT;
	while (1) {
		assert(c[i].pUpper - c[i].pLower + 1 > 0);
		if (i < kd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;

			m = (c[i].pLower + c[i].pUpper)/2;
			kdSelectInit(kd,d,m,c[i].pLower,c[i].pUpper);

			c[i].fSplit = kd->pInit[m].r[d];
			c[LOWER(i)].bnd = c[i].bnd;
			c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m-1;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			c[UPPER(i)].pLower = m;
			c[UPPER(i)].pUpper = c[i].pUpper;
			diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
			assert(diff == 0 || diff == 1);
			i = LOWER(i);
			}
		else {
			c[i].iDim = -1;
			SETNEXT(i);
			if (i == ROOT) break;
			}
		}
	UpPassInit(kd,ROOT);
	if (kd->bOutDiag) puts("<< kdBuildTree()");
	fflush(stdout);
	return(1);
	}


void kdFinish(KD kd)
{
	if (kd->pMove) free(kd->pMove);
	if (kd->pInit) free(kd->pInit);
	/*	if (kd->pGroup) free(kd->pGroup); */
	if (kd->kdNodes) free(kd->kdNodes);
	if (kd->piGroup) free(kd->piGroup);
	if (kd->grps) free(kd->grps);
	free(kd);
	}


int CmpInit(const void *p1,const void *p2)
{
	PINIT *a = (PINIT *)p1;
	PINIT *b = (PINIT *)p2;

	return(a->iOrder - b->iOrder);
	}

void Order(KD kd)
{
    if (kd->pInit != NULL) {
	qsort(kd->pInit,kd->nParticles,sizeof(PINIT),CmpInit);
    }
}


void kdWriteArray(KD kd, char *achOutFileBase)
{
    FILE *fp;
    int i,j;
    char longword[128];
    
    /* Order particles */
    Order(kd);
    /* Open grp file */
    sprintf(longword,"%s.sogrp",achOutFileBase);
    fp = fopen(longword,"w");
    assert(fp != NULL);
    fprintf(fp,"%d\n",kd->nParticles);
    for (i=0;i<kd->nParticles;++i) {
	fprintf(fp,"%d\n",kd->pInit[i].iGrp);
	/*if (kd->pInit[i].iGrp != 0)
	    fprintf(stderr,"i: %d  iOrder: %d  iGrp: %d\n",i,kd->pInit[i].iOrder,
	    kd->pInit[i].iGrp);*/
    }
    fclose(fp);
}


void kdWriteGTP(KD kd, char *achOutFileBase, int bStandard)
{
	FILE *fp;
	int i,j;
	char longword[128];
	struct dump h;
	struct star_particle sp;
	GRPNODE *gg;
	XDR xdrs;

	/*
	 ** Now output the GTP file.
	 */

	sprintf(longword,"%s.sogtp",achOutFileBase);
	fp = fopen(longword,"wb");
	assert(fp != NULL);
	h.nbodies = kd->nInGTP;
	h.nstar = kd->nInGTP;
	h.ndark = 0;
	h.nsph = 0;
	h.ndim = 3;
	h.time = kd->fTime;
	if (bStandard) {
	    assert(sizeof(Real)==sizeof(float)); /* Otherwise, this XDR stuff
						    ain't gonna work */
	    xdrstdio_create(&xdrs, fp, XDR_ENCODE);
	    xdrHeader(&xdrs,&h);
	}
	else {
	    fwrite(&h,sizeof(struct dump),1,fp);
	}
	for (i=0,gg=kd->grps;i<kd->nInGTP;++i) {
	    if (gg->index == i+1) {
		sp.mass = max(gg->fMvir,0.0);  /* Set mass to zero for errorcode groups */
		for (j=0;j<3;++j) {
			sp.pos[j] = gg->pos[j];
			sp.vel[j] = gg->vcm[j];
			}
		sp.eps = gg->fRvir;  /* Preserve errorcodes from eps field */
		sp.metals = 0.0;
		sp.tform = (float)gg->index;
		sp.phi = 0.0;
		++gg;
	    } else {
		sp.mass = 0.0;
		for (j=0;j<3;++j) {
			sp.pos[j] = 0.0;
			sp.vel[j] = 0.0;
			}
		sp.eps = 0.0;
		sp.metals = 0.0;
		sp.tform = (float)(i+1);
		sp.phi = 0.0;
	    }		
	    if (bStandard) {
		xdr_vector(&xdrs, (char *) &sp, 11,
			   sizeof(Real), (xdrproc_t) xdr_float);
	    }
	    else {
		fwrite(&sp,sizeof(struct star_particle),1,fp);
	    }
	}
	if (bStandard) xdr_destroy(&xdrs);
	fclose(fp);
}


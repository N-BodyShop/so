#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include "smooth2.h"


int smInit(SMX *psmx,KD kd,int nSmooth)
{
	SMX smx;

	assert(nSmooth <= kd->nInitActive);
	smx = (SMX)malloc(sizeof(struct smContext));
	assert(smx != NULL);
	smx->kd = kd;
	smx->nSmooth = nSmooth;
	smx->pq = (PQ *)malloc(nSmooth*sizeof(PQ));
	assert(smx->pq != NULL);
	PQ_INIT(smx->pq,nSmooth);
	smx->iMark = (char *)malloc(kd->nInitActive*sizeof(char));
	assert(smx->iMark);
	smx->pp = NULL;
	/*
	 ** Initialize arrays for calculated quantities.
	 */
	*psmx = smx;	
	/* 
	 * Initialize added quantities
	 */ 
	smx->bPeriodic = kd->bPeriodic;
	smx->nListSize = smx->nSmooth;
	smx->nnList = (NN *)malloc(smx->nListSize*sizeof(NN));	
	assert(smx->nnList != NULL);
	return(1);
	}


void smFinish(SMX smx)
{
	if (smx->pp) free(smx->pp);
	free(smx->iMark);
	free(smx->pq);
	free(smx->nnList);
	free(smx);
	}



void smGrowList(SMX smx)
{
    smx->nListSize *= 1.5;
    
    smx->nnList = (NN *) realloc(smx->nnList, smx->nListSize*sizeof(NN));
    assert(smx->nnList != NULL);
}


int smBallGather(SMX smx,FLOAT fBall2,FLOAT *ri)
{
	KDN *c = smx->kd->kdNodes;
	PINIT *p = smx->kd->pInit;
	int pj,nCnt,cp,pUpper;
	FLOAT dx,dy,dz,x,y,z,lx,ly,lz,sx,sy,sz,fDist2;
	int iDum;

	x = ri[0];
	y = ri[1];
	z = ri[2];
	assert(smx->bPeriodic);
	lx = smx->kd->fPeriod[0];
	ly = smx->kd->fPeriod[1];
	lz = smx->kd->fPeriod[2];
	nCnt = 0;
	cp = ROOT;
	while (1) {
		INTERSECT(&c[cp],fBall2,lx,ly,lz,x,y,z,sx,sy,sz,iDum,GetNextCell);
		/*
		 ** We have an intersection to test.
		 */
		if (c[cp].iDim >= 0) {
			cp = LOWER(cp);
			continue;
			}
		else {
			pUpper = c[cp].pUpper;
			/* fprintf(stderr,"fBall2: %g  pos: %g,%g,%g  (sx,sy,sz): %g,%g,%g\n",
			   fBall2,x,y,z,sx,sy,sz); */
			for (pj=c[cp].pLower;pj<=pUpper;++pj) {
				dx = sx - p[pj].r[0];
				dy = sy - p[pj].r[1];
				dz = sz - p[pj].r[2];
				fDist2 = dx*dx + dy*dy + dz*dz;
				/* fprintf(stderr,"  (px,py,pz): %g %g %g  fDist2: %g\n",
				   p[pj].r[0],p[pj].r[1],p[pj].r[2],fDist2); */
				if (fDist2 <= fBall2) {
					if(nCnt >= smx->nListSize)
					    smGrowList(smx);
					smx->nnList[nCnt].fDist2 = fDist2;
					smx->nnList[nCnt].dx = dx;
					smx->nnList[nCnt].dy = dy;
					smx->nnList[nCnt].dz = dz;
					smx->nnList[nCnt].pInit = &p[pj];
					smx->nnList[nCnt].iIndex = pj;
					nCnt++;
					}
				}
			}
	GetNextCell:
		SETNEXT(cp);
		if (cp == ROOT) break;
		}
	assert(nCnt <= smx->nListSize);
	return(nCnt);
	}






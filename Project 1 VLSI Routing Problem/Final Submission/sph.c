/* Matlab interface for a one to many shortest path routine.
   The shortest path routine is LQUEUE from Annals of Operations 
   Research, 1988. The routine for path backtracking is from 
   the disagregated simplicial decomposition code by tolar, mipat.

   Syntax:  [I,P]=SP(START,END,LENGTH,SOURCE_NODE,DESTINATION_NODES,
                     NO_OF_NODES)

   where the input parameters are,

    START = Integer vector of start nodes for the 
            arcs (in forward-star order).
    END = Integer vector of end nodes for the arcs.
    LENGHT = Vector of link costs (real valued).
    SOURCE_NODE = The root node, integer.
    DESTINATION_NODES = Vector of destination nodes, integral.
    NO_OF_NODES = The total number of nodes, integer.

   and the output,

    I = Vector of indeces giving start poisiton in the P vector 
        of the next path.  
    P = Vector of links in the paths.

   (mex) Clas Rydergren, 1997.
*/

#include <math.h>
#include "mex.h"

const INF=999999999;


void l2queue(int A[],int ND[], double LNGT[], 
            double D[], int P[], int Q[], int N, int R){

int i,j;
int U,V, NN, LAST, INIT, IFIN, PNTR; 
double DV;

/* Init */
  for(i=0; i<=N; i++)
    Q[i] = 0;          
    D[i] = INF;
  Q[R] = -1;       
  D[R] = 0;
  P[R] = 0;
  NN = N+1;
  Q[NN] = NN;
  LAST = NN;
  PNTR = NN;   
  U = R;
  for(i=1; i<=N; i++)
/* Explore forward star of U */
l20:INIT=A[U];
    IFIN=A[U+1]-1;
    if (INIT > IFIN) goto l60;   
    for(j=INIT; j<=IFIN; j++){
      V=ND[j];
      DV=D[U]+LNGT[j];
      if (D[V] <= DV) goto l50;  
      D[V]=DV;
      P[V]=U;
      if (Q[V] < 0) 
	goto l30;
      else if (Q[V] > 0) 
	goto l50; 
      else 
	goto l40; 
l30:  Q[V] = Q[PNTR];  
      Q[PNTR] = V;     
      if ( LAST == PNTR ) LAST=V;  
      PNTR = V;
      goto l50; 
l40:  Q[LAST]=V;
      Q[V]=NN;
      LAST=V;
l50:  continue;
    }
/* Remove the new current node U */    
l60:U=Q[NN];
    Q[NN]=Q[U];
    Q[U]=-1;
    if (LAST == U) LAST = NN;
    if (PNTR == U) PNTR = NN;   
    /* Check whether the list is empty */
    if (U <= N) goto l20; 
    return;
}

int backtrack(int PRED[], int  SOURCE, int SINK, int NODER, int ROAD[],
               int BAAGAR, int START[], int SLUT[], int nodlista[]){

int i, INDEX, START_INDEX, SLUT_INDEX, LENGTH;
int *SLASK;

     SLASK = (int*) mxCalloc(BAAGAR+1, sizeof(int));
    
     LENGTH=0;

     SLUT_INDEX=SINK;
l10: START_INDEX=PRED[SLUT_INDEX];
     if (START_INDEX == 0 ) goto l100;
     for(i=nodlista[START_INDEX]; i<=BAAGAR; i++)
       if(SLUT_INDEX==SLUT[i]) goto l50;
l50: INDEX=i;
     ROAD[LENGTH+1]=INDEX;
     LENGTH=LENGTH+1;
     SLUT_INDEX=START_INDEX;
     goto l10;
l100:for(i=1; i<=LENGTH; i++)
       SLASK[i]=ROAD[i];
     for(i=1; i<=LENGTH; i++)
       ROAD[i]=SLASK[LENGTH-i+1];
     mxFree(SLASK);
     return LENGTH;
}


void mexFunction(
  int nlhs,       mxArray *plhs[],
  int nrhs, const mxArray *prhs[]){

  /* Declaration of input variables */

  double *LNGT, *ND, *ST, *N, *R, *DEST, *TH;

  /* Declaration of local copies */

  int *lND, *lP, *lQ, *lST; 
  double *lLNGT, *lD;

  /* Declaration of local variables */

  int *lA, *lNODELIST, *lROAD;
  int lN, lR;
  int lB, lSINK, lLENGTH;
  int i, j, cp, NODEST, RINDEX;
  int cntZL;
  double lTH, lRC;
  
  /* Declaration of output variables */
  double *utLNGT, *utROAD, *utPR, *utZL, *utNZL;

  /* Catch input parameters */
 
  ST=mxGetPr(prhs[0]);    /* Start nodes of arcs */
  ND=mxGetPr(prhs[1]);    /* End nodes of arcs */
  LNGT=mxGetPr(prhs[2]);  /* Arcs lenghts */
  R=mxGetPr(prhs[3]);     /* Source node*/
  DEST=mxGetPr(prhs[4]);  /* Destinations */
  N=mxGetPr(prhs[5]);     /* Number of nodes */
  TH=mxGetPr(prhs[6]);    /* Treshold for the zero costs */

  /* Allocation of local variables */

  lN=(int) *N;         /* Number of nodes */
  lR=(int) *R;         /* Source node */
  lTH=(double) *TH;    /* Treshold */
  lB=mxGetM(prhs[0]);  /* Number of arcs */

  lA = (int*) mxCalloc(lN+2, sizeof(int));
  lND = (int*) mxCalloc(lB+1, sizeof(int));
  lD = (double*) mxCalloc(lN+1, sizeof(double));
  lP = (int*) mxCalloc(lN+1, sizeof(int));
  lQ = (int*) mxCalloc(lN+1, sizeof(int));
  lLNGT = (double*) mxCalloc(lB+1, sizeof(double));
  lNODELIST = (int*) mxCalloc(lN+2, sizeof(int));
  lST = (int*) mxCalloc(lB+1, sizeof(int));
  lROAD = (int*) mxCalloc(lB+1, sizeof(int));



  for(i=0; i<=lB; i++)
    lND[i+1]=(int) ND[i];
  for(i=0; i<=lB; i++)
    lLNGT[i+1]=LNGT[i];
  for(i=0; i<=lB; i++)
    lST[i+1]=(int)ST[i];


  for(i=0; i<=lN; i++){
    lD[i]=INF;
  }


  for(i=1;i<=lN+1;i++)
    lA[i]=lB+1;
  for(i=lB;i>=1;i--)
    if(i < lA[lST[i]])
      lA[lST[i]]=i;


  l2queue(lA, lND, lLNGT, lD, lP, lQ, lN, lR);


  for(i=1; i<=lN;i++)
    lNODELIST[i]=0;
  cp=lST[1];
  lNODELIST[cp]=1;
  for(i=2; i<=lB;i++)
    if(lST[i] != cp) {
      cp=lST[i];
      lNODELIST[cp]=i;
  }

  
  NODEST=mxGetN(prhs[4]);  
  plhs[0]=mxCreateDoubleMatrix(NODEST, 1, mxREAL);
  utLNGT=mxGetPr(plhs[0]);

  plhs[1]=mxCreateDoubleMatrix(1, NODEST*lB, mxREAL);
  utROAD=mxGetPr(plhs[1]);
  
  plhs[2]=mxCreateDoubleMatrix(1, lN, mxREAL);
  utPR=mxGetPr(plhs[2]);
  for(i=0;i<=lN-1;i++)
    utPR[i]=0;

  plhs[3]=mxCreateDoubleMatrix(1, lB+1, mxREAL);
  utZL=mxGetPr(plhs[3]);
  for(i=0;i<=lB-1;i++)
    utZL[i]=0;

  plhs[4]=mxCreateDoubleMatrix(1, 1, mxREAL);
  utNZL=mxGetPr(plhs[4]);
    utNZL[0]=0;

  for(i=0;i<=lN-1;i++)
    utPR[i]=lP[i+1];
      
  cntZL=0;
  for(i=0;i<=lB-1;i++){
    lRC=lD[lST[i+1]]-lD[lND[i+1]]+lLNGT[i+1];
    if (lRC < lTH){
      if (lP[lND[i+1]] != lST[i+1]){
	    utZL[cntZL]=i+1;
	    cntZL++;
	    }
    }
  }

  utNZL[0]=cntZL;
  mxSetN(plhs[3],cntZL);    
    
  RINDEX=0;
  for(j=0;j<NODEST;j++){
    lSINK=(int) DEST[j];

    lLENGTH=backtrack(lP, lR, lSINK, lN, lROAD, lB, lST, 
                     lND, lNODELIST);

                     
    utLNGT[j]=lLENGTH;
    for(i=1;i<=lLENGTH;i++){
      utROAD[RINDEX]=lROAD[i];
      RINDEX=RINDEX+1;
    }
  }
  mxSetN(plhs[1],RINDEX);    
  mxSetM(plhs[1],1);    
  
}
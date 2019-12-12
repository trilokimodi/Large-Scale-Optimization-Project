/* Matlab interface for a many to many shortest path routine.
   The code includes a grid construction route for a Manhattan wiring grid 
   a la Feo and Hochbaum "Lagrangian relaxation for testing infeasibility
   in VLSI routing".
   The shortest path routine is LQUEUE from Annals of Operations 
   Research, 1988. 
   
   Syntax:  [routenodes]=GSP(DimX,DimY,pi,k,com)

   where the input parameters are,

    DimX = Number of horizontal nodes
    DimY = Number of nodes vertically
    pi = Link costs for the. Vector of size DimX*DimY*2 
    k = Number of commodities
    com = Matrix of commodity pairs. Size k times 2.  
 
    and the output,

    routenodes = A vector of node numbers corresponding to the used routes for all commodities. Size unknown. 

   (mex) Clas Rydergren, 2002.
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

     int INDEX, START_INDEX, SLUT_INDEX, LENGTH;
    
     LENGTH=1;
     ROAD[LENGTH++]=SINK;
     SLUT_INDEX=SINK;
     START_INDEX=PRED[SLUT_INDEX];
     while (START_INDEX != 0 ){
        ROAD[LENGTH++]=START_INDEX;
        SLUT_INDEX=START_INDEX;
        START_INDEX=PRED[SLUT_INDEX];
     }
     return LENGTH-1;
}


void mexFunction(
    int nlhs,       mxArray *plhs[],
    int nrhs, const mxArray *prhs[]){

    /* Declaration of input variables */
    double *dimX, *dimY, *pi, *k, *com;

    /* Declaration of local copies */

    int *lND, *lP, *lQ, *lST; 
    double *lLNGT, *lcom, *lD;

    /* Declaration of local variables */

    int *lA, *lNODELIST, *lROAD;
    int lN, lR, lDEST, ldimX, ldimY, lk;
    int lB, lSINK, lLENGTH;
    int i, j, count, cp, NODEST, RINDEX, commod;
  
    /* Declaration of output variables */
    double *utROUTENODES;

    /* Catch input parameters */
 
    dimX=mxGetPr(prhs[0]);    /*  */
    dimY=mxGetPr(prhs[1]);    /*  */
    pi=mxGetPr(prhs[2]);  /* Arcs lenghts */
    k=mxGetPr(prhs[3]);     /* commodities */
    com=mxGetPr(prhs[4]);  /* Matrix of commodity start and end nodes */

    
    
    ldimX=(int)*dimX;
    ldimY=(int)*dimY;
    lk=(int)*k;

    if(mxGetM(prhs[2])!=ldimX*ldimY*2){
        printf("Wrong size of the pi vector! Size: %d, required size: %d\n",mxGetM(prhs[2]),ldimX*ldimY*2);
        plhs[0]=mxCreateDoubleMatrix(1, 1, mxREAL);
        return;
    }

    if(mxGetN(prhs[4])!=2 && mxGetM(prhs[4])!=lk){
        printf("Wrong number of rows in the commodity matrix! Rows: %d, required rows: %d\n",mxGetM(prhs[4]),lk);
        plhs[0]=mxCreateDoubleMatrix(1, 1, mxREAL);
        return;
    }
    j=0;
    for (i=0;i<ldimY*ldimY*2;i++){
      if (pi[i]<0) j=1;
    }
    if(j==1){
        printf("At least one negativ node cost detected!");
        plhs[0]=mxCreateDoubleMatrix(1, 1, mxREAL);
        return;
    }
      

    utROUTENODES = (double*) mxCalloc(ldimX*ldimY, sizeof(double));    

    lST = (int*) mxCalloc(2*ldimX*ldimY*3+1, sizeof(int));
    lND = (int*) mxCalloc(2*ldimX*ldimY*3+1, sizeof(int));

    count=1;
    for (i=1;i<=ldimY;i++){
        for (j=1;j<=ldimX;j++){
            if (j>1){
                lST[count]=j+ldimX*(i-1);
                lND[count++]=j+ldimX*(i-1)-1;
            }
            if(j<ldimX){
                lST[count]=j+ldimX*(i-1);
                lND[count++]=j+ldimX*(i-1)+1;
            }
            lST[count]=j+((ldimX)*(i-1));
            lND[count++]=ldimX*ldimY+i+((ldimY)*(j-1));
        }
    }
    for(i=1;i<=ldimX;i++){
        for(j=1;j<=ldimY;j++){
            lST[count]=ldimX*ldimY+j+((ldimY)*(i-1));
            lND[count++]=i+((ldimX)*(j-1));
            if(j>1){
                lST[count]=ldimX*ldimY+j+ldimY*(i-1);
                lND[count++]=ldimX*ldimY+j+ldimY*(i-1)-1;
            }
            if(j<ldimY){
                lST[count]=ldimX*ldimY+j+ldimY*(i-1);
                lND[count++]=ldimX*ldimY+j+ldimY*(i-1)+1;
            }   
        }
    }

    lB=count-1;               /* number of arcs */
    lN=ldimX*ldimY*2;         /* Number of nodes */
    
    /* Alloc mem for local variables */

    lROAD = (int*) mxCalloc(lB+1, sizeof(int));
    lA = (int*) mxCalloc(lN+2, sizeof(int));
    lD = (double*) mxCalloc(lN+2, sizeof(double));
    lP = (int*) mxCalloc(lN+2, sizeof(int));
    lQ = (int*) mxCalloc(lN+2, sizeof(int));
    lLNGT = (double*) mxCalloc(lB+1, sizeof(double));
    lNODELIST = (int*) mxCalloc(lN+2, sizeof(int));


    for(i=1; i<=lN;i++)
        lNODELIST[i]=0;
    cp=lST[1];
    lNODELIST[cp]=1;
    for(i=2; i<=lB;i++)
        if(lST[i] != cp) {
            cp=lST[i];
            lNODELIST[cp]=i;
        }
    for(i=1;i<=lN+1;i++)
        lA[i]=lB+1;
    for(i=lB;i>=1;i--)
        if(i < lA[lST[i]])
            lA[lST[i]]=i;

    plhs[0]=mxCreateDoubleMatrix(ldimX*ldimY*lk+500, 1, mxREAL);
    utROUTENODES=mxGetPr(plhs[0]);
    RINDEX=0;
    for(commod=1;commod<=lk;commod++){
            
        for(i=1;i<=lB;i++){
            lLNGT[i]=pi[lST[i]-1]+(1+rand())*1e-12;
        }
        for(i=0; i<=lN; i++){
            lD[i]=INF;
            lP[i]=0;
            lQ[i]=0;
        }
        lR=(int)com[commod-1];
        lDEST=(int)com[lk+commod-1];
        l2queue(lA, lND, lLNGT, lD, lP, lQ, lN, lR);
        lLENGTH=backtrack(lP, lR, lDEST, lN, lROAD, lB, lST, lND, lNODELIST);
        for(i=1;i<=lLENGTH;i++){
            utROUTENODES[RINDEX]=lROAD[i];
            RINDEX=RINDEX+1;
        }
    }
    mxSetN(plhs[0],1);    
    mxSetM(plhs[0],RINDEX);    

}
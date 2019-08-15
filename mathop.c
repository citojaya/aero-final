#include "common.h"

/* Return particle volume*/
real partVol(int p){
    return (4.0/3.0)*PI*pow(demPart[p].dia*0.5,3);
}

/*Find solid fraction within cell radius*/
real solidFraction(int ip){
    //Find cell center
    int iIndex = floor((demPart[ip].pos[0] - xmin)/domainDx);
    int jIndex = floor((demPart[ip].pos[1] - ymin)/domainDy);
    int kIndex = floor((demPart[ip].pos[2] - zmin)/domainDz);
    int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;     
    
    real cX = iIndex*domainDx + 0.5*domainDx;
    real cY = jIndex*domainDy + 0.5*domainDy;
    real cZ = kIndex*domainDz + 0.5*domainDz;
    real vol = 0.0;
    real totVol = 0.0;
    for(int i=0; i<bdBox[cellIndex].noOfParticles; i++){
        int jp = bdBox[cellIndex].parts[i];
        real r = 0.5*demPart[jp].dia;
        real jpX = demPart[jp].pos[0];
        real jpY = demPart[jp].pos[1];
        real jpZ = demPart[jp].pos[2];
        
        real rR = sqrt((jpX-cX)*(jpX-cX) + (jpY-cY)*(jpY-cY) + (jpZ-cZ)*(jpZ-cZ)); 
        
        if( rR <= (0.5*domainDx-r)){
            vol = 0.8*partVol(jp);
        }
        else if(rR > (cellRadius-r)){
            vol = partVol(jp);
        }
        else if(rR <= (cellRadius-r) && rR  > (domainDx*0.5-r)){
            vol = 0.5*partVol(jp); 
        }
        totVol += vol;       
    }
    //writeLogNum("logfile3.log","VOL ",totVol/(lengthFactor*lengthFactor*lengthFactor));
    //bdBox[cellIndex].fluidVolF = totVol/((4.0/3.0)*PI*pow(cellRadius,3));
    return totVol/((4.0/3.0)*PI*pow(cellRadius,3));
}


void getNearFaceRot(Particle *p, struct wallFace *fc){
    //Find cell center
    int ip = p->part_id;
    int iIndex = floor((demPart[ip].pos[0] - xmin)/wBoxDx);
    int jIndex = floor((demPart[ip].pos[1] - ymin)/wBoxDy);
    int kIndex = floor((demPart[ip].pos[2] - zmin)/wBoxDz);
    real prevDist = 1.e3;
    //real prevDist2 = 1.e3;
    real dist = 0.;
    nFace = -1;
    //writeLog3Num("wallmesh.log","getNearF ",nFace);
    for(int k=kIndex-1; k<kIndex+2; k++){
        for(int j=jIndex-1; j<jIndex+2; j++){
            for(int i=iIndex-1; i<iIndex+2; i++){
                if(i>=0 && j>=0 && k>=0){
                    if(i<wBXDiv && j<wBYDiv && k<wBZDiv){
                int cI = i+ j*wBXDiv + k*wBXDiv*wBYDiv;
                for(int n=0; n<wBox[cI].totalFaces; n++){
                    int fId = wBox[cI].face[n];
                    //writeLogNum("wallmesh.log","fID ",fId);
                    real ipX = demPart[ip].pos[0];
                    real ipY = demPart[ip].pos[1];
                    real ipZ = demPart[ip].pos[2];
                    
                    // Rotor face
                    //if(wBox[cI].faceType[n] == rotZoneID){
                    real jpX = fc[fId].centroid[0];
                    real jpY = fc[fId].centroid[1];
                    real jpZ = fc[fId].centroid[2];
                    //}
                    
                    dist = sqrt((ipX-jpX)*(ipX-jpX) + (ipY-jpY)*(ipY-jpY) + (ipZ-jpZ)*(ipZ-jpZ));
                    
                    if(dist < prevDist){
                        nFace = fId;
                        faceType = wBox[cI].faceType[n];
                        prevDist = dist;
                    }
                }}
                }
                    
            }
        }
    }
    if(nFace > 0){
        real *dVec = allocateDoubleArray(DIM);
        getProjectOnFace(demPart[ip].pos, fc[nFace].node1,fc[nFace].node2,fc[nFace].node3, dVec);
        findContactFromMesh(p, dVec);
        free(dVec);
    }
    //return face;
}

/*
  Search through all the faces in neighbour cells and return the nearest face
*/
void getNearFace(Particle *p, struct wallFace *fc){
    //Find cell center
    int ip = p->part_id;
    int iIndex = floor((demPart[ip].pos[0] - xmin)/domainDx);
    int jIndex = floor((demPart[ip].pos[1] - ymin)/domainDy);
    int kIndex = floor((demPart[ip].pos[2] - zmin)/domainDz);
    real prevDist = 1.e3;
    real dist = 0.;
    nFace = -1;

       
    //writeLog3Num("wallmesh.log","getNearF ",nFace);
    //for(int k=kIndex-1; k<kIndex+2; k++){
        //for(int j=jIndex-1; j<jIndex+2; j++){
            //for(int i=iIndex-1; i<iIndex+2; i++){
                //int cI = i+ j*xDiv + k*yDiv*zDiv;
                int cI = iIndex+ jIndex*xDiv + kIndex*xDiv*yDiv;
                for(int n=0; n<bdBox[cI].totalFaces; n++){
                    int fId = bdBox[cI].face[n];
                    //distance vector (from proejction to particle center)
                    real *dVec = allocateDoubleArray(DIM);
                    // struct wallFace *fc;
                    // // Rotor face
                    // if(bdBox[cI].faceType[n] == rotZoneID){
                    //     fc = wFaceRot;
                    // }
                    // // Stator face
                    // else if(bdBox[cI].faceType[n] == statZoneID){
                    //     fc = wFace;
                    // }
                    //real *temp = allocateDoubleArray(DIM);
                    //vecSub(demPart[ip].pos,fc[fId].centroid,temp);
                    //if(vecMag(temp) < demPart[ip].dia*3.){
                        if(getProjectOnFace(demPart[ip].pos, fc[fId].node1,fc[fId].node2,fc[fId].node3, dVec) == 1){
                            nFace = fId;
                            faceType = bdBox[cI].faceType[n];
                            findContactFromMesh(p, dVec);
                        }
                    //}
                    //free(temp);
                    free(dVec);
                 }
            //}
        //}
    //}
}

real triArea(real *v1, real *v2){
    real a = 0;
    //Cross product
    real *v1v2 = allocateDoubleArray(DIM);
    crossProd(v1,v2,v1v2);
    a = 0.5*vecMag(v1v2);
    free(v1v2);
    return a;
}

/* If porjection is inside the triangle return 1 else return */
int getProjectOnFace(real *p, real *n1, real *n2, real *n3, real *dVec){
    real *n1n2 = allocateDoubleArray(DIM);
    real *n1n3 = allocateDoubleArray(DIM);
    real *pn1 = allocateDoubleArray(DIM);
    vecSub(n2,n1,n1n2);
    vecSub(n3,n1,n1n3);
    vecSub(p, n1, pn1);
    
    //Temp unit vector
    real *uVec = allocateDoubleArray(DIM);
    getUnitVector(n1,n2, n3,uVec);
    real *pj = allocateDoubleArray(DIM);
    projVec(pn1 ,uVec, pj, 0);

    //Projection to particle vector
    vecSub(pn1,pj,dVec);
    
    //Projection to n2 vector
    real *pn2 = allocateDoubleArray(DIM);
    vecSub(n1n2,pj, pn2);
    //Projection to n1 vector
    real *pn3 = allocateDoubleArray(DIM); 
    vecSub(n1n3,pj,pn3);
    //writeLog3Num("wallmesh.log","proj face ",pj[0],pj[1],pj[2]);

    real a1 = triArea(n1n2, pj);
    real a2 = triArea(n1n3, pj);
    real a3 = triArea(pn2,pn3);
    real totalArea = triArea(n1n2,n1n3);

    free(pn2);
    free(pn3);
    free(pn1);
    free(pj);
    free(uVec);
    free(n1n2);
    free(n1n3);

    // writeLogNum("wallmesh.log","totalArea ",totalArea);
    // writeLogNum("wallmesh.log","a1+a2+a3 ",a1+a2+a3);
    //if(fabs(a1+a2+a3 - totalArea) < 1.e-10*pow(lengthFactor,3)){
    //     writeLogNum("wallmesh.log","totalArea ",totalArea);
    //     writeLogNum("wallmesh.log","a1+a2+a3 ",a1+a2+a3);
    if(a1+a2+a3 == totalArea){
        return 1;
    }
    return 0;
}

/* Return center distance of two particles
param:
p1 - particle 1
p2 - particle 2 */
real getCenterDist(int ip, int jp){


    real ipX = demPart[ip].pos[0];
    real ipY = demPart[ip].pos[1];
    real ipZ = demPart[ip].pos[2];
    real jpX = demPart[jp].pos[0];
    real jpY = demPart[jp].pos[1];
    real jpZ = demPart[jp].pos[2];

    real val = sqrt((ipX-jpX)*(ipX-jpX) + (ipY-jpY)*(ipY-jpY) + (ipZ-jpZ)*(ipZ-jpZ));
    return val;
 }

/* Add two vecotrs
param:
v1 - input vector 1
v2 - input vector 2
vec - resultant vector */
void vecAdd(real *v1, real *v2, real *vec){
   vec[0] = v1[0]+v2[0];
   vec[1] = v1[1]+v2[1];
   vec[2] = v1[2]+v2[2];
}

/* Substract two vecotrs
param:
v1 - input vector 1
v2 - input vector 2
vec - resultant vector */
void vecSub(real *v1, real *v2, real *vec){
   vec[0] = v1[0]-v2[0];
   vec[1] = v1[1]-v2[1];
   vec[2] = v1[2]-v2[2];
}

/* Find relative velocity magnitude between two particles */
real relVel(int ip, int jp){
    real *relVelVec = allocateDoubleArray(DIM);
    vecSub(demPart[ip].vel, demPart[jp].vel, relVelVec);
    real vel = vecMag(relVelVec); 
    free(relVelVec);
    return vel;
}

/* Find vector cross product
param:
v1 - input vector 1
v2 - input vector 2
vec - resultant vector */
void crossProd(real *v1, real *v2, real *vec){
    real temp1 = v1[1]*v2[2];
    real temp2 = v2[1]*v1[2];
    vec[0] = temp1 - temp2;

    temp1 = v2[0]*v1[2];
    temp2 = v1[0]*v2[2];
    vec[1] = temp1 - temp2;

    temp1 = v1[0]*v2[1];
    temp2 = v2[0]*v1[1];
    vec[2] = temp1 - temp2;   
}

/*
Dot product of two vectors
param:
v1 - vector 1
v2 - vector 2
*/
real dotProduct(real *v1, real *v2){
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

/*
Find normal unit vector to the surface defined by vector v1v2 and v2v3
*/
void getUnitVector(real *v1, real *v2, real *v3, real *uVec){
    real *v1v2 = allocateDoubleArray(DIM);
    real *v1v3 = allocateDoubleArray(DIM);

    vecSub(v1, v2, v1v2);
    vecSub(v1, v3, v1v3);
    crossProd(v1v2, v1v3, uVec);
    unitVec(uVec, uVec);

    free(v1v2);
    free(v1v3);
}

/* Find unit vector
param: 
v - input vector
vec - unit vector */
void unitVec(real *v, real *vec){
    if(vecMag(v) > 0){
        real temp1 = v[0]/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        real temp2 = v[1]/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        real temp3 = v[2]/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        vec[0] = temp1;
        vec[1] = temp2;
        vec[2] = temp3;
    }
    else{
        vec[0] = 0.0;
        vec[1] = 0.0;
        vec[2] = 0.0;
    }
}

/* Find magnitude of a vector*/
real vecMag(real *vec){
    return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}

/* Multiply input vector by a scaler
param: 
scl - scaler
vec - vector to be multiplied*/
void sclMult(real scl, real *vec){
    vec[0] = scl*vec[0];
    vec[1] = scl*vec[1];
    vec[2] = scl*vec[2];
}

/* Multiply input vector by a scaler and assign to a vector
param: 
scl - scaler
inVec - input vector to be multiplied
outVec - reusltant output vector*/
void sclVecMult(real scl, real *inVec, real *outVec){
    outVec[0] = scl*inVec[0];
    outVec[1] = scl*inVec[1];
    outVec[2] = scl*inVec[2];
}


/* Returns the project vector on the plane defined by normal unit vector
param:
v - input vector
n - unit vector
vec - resultant project vector
type - either 0 or 1, 0 gives project vector, 1 gives project vector scaled by input vector */
void projVec(real *v, real *n, real *vec, int type){
    real *tV1 = malloc(DIM*sizeof(real));
    real *tV2 = malloc(DIM*sizeof(real));
    if(type == 0){
        crossProd(n,v, tV1);
        tV2[0] = -n[0];
        tV2[1] = -n[1];
        tV2[2] = -n[2];
        //printf("tV1 %lf,%lf,%lf\n ",tV1[0],tV1[1],tV1[2]);
        //printf("tV2 %lf,%lf,%lf\n ",tV2[0],tV2[1],tV2[2]);
        crossProd(tV2, tV1, vec);
        //printf("Type 0\n");
    }
    else{
        real *tV3 = malloc(DIM*sizeof(real));
        crossProd(n, v, tV1);
        tV2[0] = -n[0];
        tV2[1] = -n[1];
        tV2[2] = -n[2];
        crossProd(tV2, tV1, tV3);
        real temp;
        if((v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) != 0.0){
            temp = sqrt((tV3[0]*tV3[0]+tV3[1]*tV3[1]+tV3[2]*tV3[2])/
                            (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));  
        }
        else{
            temp = 0.0;
        }
        vec[0] = tV3[0]*temp;
        vec[1] = tV3[1]*temp;
        vec[2] = tV3[2]*temp;      
        free(tV3);
    }
    free(tV1);
    free(tV2);
}




#include "common.h"

//allocated = 0;
/*Read particle information and domain information
*/
void demInit(){
    xmin = 1.e6;
    ymin = 1.e6;
    zmin = 1.e6;
    xmax = -1.e6;
    ymax = -1.e6;
    zmax = -1.e6;  
    statZoneID = 3;
    maxFaces = 0;
    maxRotFaces = 0;
    noOfRotorFaces = 1;
    noOfStatorFaces = 1;
    maxCnt = 0;
    maxCharge = 0.0;
    maxES = 0.0;
    maxVF = 0.0;
    saveDEM = 0;
    updateDPM = 0;
    //dump_count = 0;
    cfdcycles = 0;
    demcycles = 0;
    //time_count = 0;// time counter used for saving output file
    particle_counter = 0; //counter for keeping the track of number of particles
    demTime = 0.0;

    readDomain("infile");

    //Read material data
    readInput("infile", &parArraySize, &dens, &ymod, &pois, &sfc, &rec, &dmpn, &rf, &cyldia, &timeStep, 
            &noOfWalls, &updateDPM, &maxVel, &rotZoneID);
    
    //Read particle-wall contact surfaces
    readWalls("infile", walls);
    readInputVelocity("infile");
    readGeom("infile", &ductxmin, &ductxmax, &ductxedge1, &ductxedge2, &ductymin, 
            &ductymax, &ductzmin, &ductzmax, &ductzedge);

    writeLog3Num("logfile2.log", "Initial min ",xmin,ymin,zmin);
    writeLog3Num("logfile2.log", "Initial max ",xmax,ymax,zmax);
    writeLogNum("logfile2.log","ductxmin ",ductxmin);
    writeLogNum("logfile2.log","ductxmax ",ductxmax);
    writeLogNum("logfile2.log","ductxedge1 ",ductxedge1);
    writeLogNum("logfile2.log","ductxedge2 ",ductxedge2);
    writeLogNum("logfile2.log","ductymin ",ductymin);
    writeLogNum("logfile2.log","ductymax ",ductymax);
    writeLogNum("logfile2.log","ductzmin ",ductzmin);
    writeLogNum("logfile2.log","ductzmax ",ductzmax);
    writeLogNum("logfile2.log","ductzedge ",ductzedge);
    writeLogNum("logfile2.log","maxVel ",maxVel);
    writeLogNum("logfile2.log","PP haConst (*1e20) ",haPP*1e20);
    writeLogNum("logfile2.log","PW haConst (*1e20) ",haPW*1e20);
    writeLogNum("logfile2.log","lamda1 ",lamda1*1e9);
    writeLogNum("logfile2.log","lamda2 ",lamda2*1e9);
    writeLogNum("logfile2.log","rms1 ",rms1*1e9);
    writeLogNum("logfile2.log","rms2 ",rms2*1e9);
    writeLogNum("logfile2.log","Update DPM ",updateDPM);
    writeLogNum("logfile2.log","Particle Array Size ",parArraySize);
    writeLogNum("logfile2.log","density ",dens);
    writeLogNum("logfile2.log","Youngs Modulus ",ymod);
    writeLogNum("logfile2.log","largestParDia ",largestParDia);
    writeLogNum("logfile2.log","V1 ",V1);
    writeLogNum("logfile2.log","V2 ",V2);
    writeLogNum("logfile2.log","imageConst ",imageConst);
    writeLogNum("logfile2.log","Zs (x1e9) ",Zs*1.e9);
    writeLogNum("logfile2.log","alpha ",alpha);
    writeLogNum("logfile2.log","permitivity (x1e-12)",permitivity*1e12);
    writeLogNum("logfile2.log","K0 ",imageConst*Zs/permitivity/(4.0*PI*pow(largestParDia*0.5,2)));
    writeLogNum("logfile2.log","esfTrue ",esfTrue);
    real qMax = (V1-V2)/(imageConst*Zs/permitivity/(4.0*PI*pow(largestParDia*0.5,2)));
    real maxESForce = qMax*qMax/(4.*PI*permitivity*1.0e-6*1.0e-6);
    writeLogNum("logfile2.log","Qmax (x1e12)",qMax*1e12);
    writeLogNum("logfile2.log","maxESForce (x1e6) ",maxESForce*1e6);
    writeLogNum("logfile2.log","Timestep (x1e6)",timeStep*1.e6);
    writeLogNum("logfile2.log","noOfCFDCells ",noOfCFDCells);
    writeLogNum("logfile2.log","maxBound ",maxBound);
    writeLogNum("logfile2.log","rotor zone id ",rotZoneID);

    noOfRotorFaces = getFaceNumber(rotZoneID);
    // noOfStatorFaces = getFaceNumber(statZoneID);

    readWallMesh("stator-wall-mesh", 0); //read number of wall faces
    noOfStatorFaces = noOfFaces;
    readWallMesh("stator-wall-mesh", 3); //set boundary
        
    writeLogNum("wallmesh.log","no of rotor faces ",noOfRotorFaces);
    writeLogNum("wallmesh.log","no of stator faces ",noOfStatorFaces);
    //Initial cell divisions for particles
    xDiv = floor((xmax-xmin)/(largestParDia*multif3));
    yDiv = floor((ymax-ymin)/(largestParDia*multif3));
    zDiv = floor((zmax-zmin)/(largestParDia*multif3));

    writeLogNum("wallmesh.log","FINALE xDiv ",xDiv);
    writeLogNum("wallmesh.log","FINALE yDiv ",yDiv);
    writeLogNum("wallmesh.log","FINALE zDiv ",zDiv);
    domainDx = (xmax-xmin)/xDiv;
    domainDy = (ymax-ymin)/yDiv;
    domainDz = (zmax-zmin)/zDiv;

    //Adjust boundary
    xmin += -domainDx*3.;
    ymin += -domainDy*3.;
    zmin += -domainDz*3.;

    xmax += domainDx*3.;
    ymax += domainDy*3.;
    zmax += domainDz*3.;

    //Set grid size for wall mesh
    setWallGridSize();
    writeLogLine("logfile2.log"," setWallGridSize\n ");

    //New cell divisions
    //Initial cell divisions for particles
    xDiv = floor((xmax-xmin)/(largestParDia*multif3));
    yDiv = floor((ymax-ymin)/(largestParDia*multif3));
    zDiv = floor((zmax-zmin)/(largestParDia*multif3));

    writeLogNum("logfile2.log","FINALE xDiv ",xDiv);
    writeLogNum("logfile2.log","FINALE yDiv ",yDiv);
    writeLogNum("logfile2.log","FINALE zDiv ",zDiv);
    domainDx = (xmax-xmin)/xDiv;
    domainDy = (ymax-ymin)/yDiv;
    domainDz = (zmax-zmin)/zDiv;

    writeLogNum("wallmesh.log","FINAL XMIN ",xmin);
    writeLogNum("wallmesh.log","YMIN ",ymin);
    writeLogNum("wallmesh.log","ZMIN ",zmin);
    writeLogLine("wallmesh.log","----------------\n");
    writeLogNum("wallmesh.log","FINAL XMAX ",xmax);
    writeLogNum("wallmesh.log","YMAX ",ymax);
    writeLogNum("wallmesh.log","ZMAX ",zmax);
    //readWallMesh("capsule-wall-mesh", 0);
    //noOfRotorFaces = noOfFaces;
    //readWallMesh("stator-wall-mesh", 0); //read number of wall faces
    //noOfStatorFaces = noOfFaces;
    //Get number of rotor faces

    //Allocate memory
    allocate();
    writeLogLine("logfile2.log"," allocate\n ");
   
    //Read wall face nodes
    //readWallMesh("capsule-wall-mesh", 1);
    //Read stator wall mesh
    readWallMesh("stator-wall-mesh", 2);  
}

/*Get number of rotor faces*/
int getFaceNumber(int zID){
  // Extract wall faces on rotor 
  cell_t c;
  face_t f;
  int noOfFaces = 0;
  real x[ND_ND];
  Domain *d = Get_Domain(1);
  Thread *th = Lookup_Thread(d,zID);
    begin_c_loop(c,th)
      int n;
      c_face_loop(c,th,n) 
      { 
        Thread *tf = C_FACE_THREAD(c,th,n); 
        if(THREAD_TYPE(tf) == THREAD_F_WALL){
            f = C_FACE(c,th,n);
            int m;
            int j=0;
            real *n1 = allocateDoubleArray(DIM);
            real *n2 = allocateDoubleArray(DIM);
            real *n3 = allocateDoubleArray(DIM);
            f_node_loop(f,tf,m)//Loop over face nodes
            {
                Node *node = F_NODE(f,tf,m);
                //real tLength = sqrt(NODE_X[node]); 
                //minEdgeLength = fmin(minEdgeLength,tLength);
                //maxEdgeLength = fmax(maxEdgeLength,tLength);
                switch(j){
                case 0:
                    n1[0] = NODE_X(node);
                    n1[1] = NODE_Y(node);
                    n1[2] = NODE_Z(node);
                case 1:
                    n2[0] = NODE_X(node);
                    n2[1] = NODE_Y(node);
                    n2[2] = NODE_Z(node);
                case 2:
                    n3[0] = NODE_X(node);
                    n3[1] = NODE_Y(node);
                    n3[2] = NODE_Z(node);
                }
                j++;
                xmin = fmin(xmin,n1[0]);
                xmin = fmin(xmin,n2[0]);
                xmin = fmin(xmin,n3[0]);
                xmax = fmax(xmax,n1[0]);
                xmax = fmax(xmax,n2[0]);
                xmax = fmax(xmax,n3[0]);

                ymin = fmin(ymin,n1[1]);
                ymin = fmin(ymin,n2[1]);
                ymin = fmin(ymin,n3[1]);
                ymax = fmax(ymax,n1[1]);
                ymax = fmax(ymax,n2[1]);
                ymax = fmax(ymax,n3[1]);

                zmin = fmin(zmin,n1[2]);
                zmin = fmin(zmin,n2[2]);
                zmin = fmin(zmin,n3[2]);
                zmax = fmax(zmax,n1[2]);
                zmax = fmax(zmax,n2[2]);
                zmax = fmax(zmax,n3[2]);
            }

            minMaxEdgeLength(n1, n2, n3, zID);
            free(n1);
            free(n2);
            free(n3);
            noOfFaces++;
        }
      }
    end_c_loop(c,rot_t)

    return noOfFaces;
}

/* Calculate divisions for wall mesh grid*/
void setWallGridSize(){
    wBXDiv = ceil((xmax-xmin)/(maxRotEdgeLength*2.0));
    wBYDiv = ceil((ymax-ymin)/(maxRotEdgeLength*2.0));
    wBZDiv = ceil((zmax-zmin)/(maxRotEdgeLength*2.0));
    wBoxDx = (xmax-xmin)/wBXDiv;
    wBoxDy = (ymax-ymin)/wBYDiv;
    wBoxDz = (zmax-zmin)/wBZDiv;
    writeLog3Num("wallmesh.log","Rotor mesh grid wBXDiv,wBYDiv,wBZDiv ",wBXDiv,wBYDiv,wBZDiv);
    writeLog3Num("wallmesh.log","Rotor mesh grid wBoxDw,wBoxDy,wBoxDz ",wBoxDx,wBoxDy,wBoxDz);
}
/*
Add wall face to BD box
fc - face index
type - face type (rotor=1, stator=2)
cI - BD box index
*/
void addFaceToBDBox(int fc, int type, int cI){
    bdBox[cI].face[bdBox[cI].totalFaces] = fc;
    bdBox[cI].faceType[bdBox[cI].totalFaces] = type;
    bdBox[cI].totalFaces++;
    //wFaceRot[fc].bdBoxIndex = cI; //assign new cell index
}

/*
Delete wall face from BD box
fc - face index
type - face type (rotor=1, stator=2)
cI - BD box index
*/
void deleteFace(int fc, int type, int cI){
    for (int i=0; i<bdBox[cI].totalFaces; i++){
        //int wallF = bdBox[cI].face[i];
        //int fType = bdBox[cI].faceType[i];
        //if(bdBox[cI].face[i] == fc && bdBox[cI].faceType[i] == type){
            //writeLogNum("logfile5.log","DELETE ", cnt);
            bdBox[cI].face[i] = bdBox[cI].face[bdBox[cI].totalFaces-1];
            bdBox[cI].faceType[i] = bdBox[cI].faceType[bdBox[cI].totalFaces-1];
            bdBox[cI].totalFaces--;
            //break;
        //}
    }
}

/*
Add wall face to BD box
fc - face index
type - face type (rotor=1, stator=2)
cI - BD box index
*/
void addRotFaceToBDBox(int fc, int type, int cI){
    wBox[cI].face[wBox[cI].totalFaces] = fc;
    wBox[cI].faceType[wBox[cI].totalFaces] = type;
    wBox[cI].totalFaces++;
    //wFaceRot[fc].wBoxIndex = cI; //assign new cell index
}

/*
Delete wall face from BD box
fc - face index
type - face type (rotor=1, stator=2)
cI - BD box index
*/
void deleteRotFace(int fc, int type, int cI){
    for (int i=0; i<wBox[cI].totalFaces; i++){
        //int wallF = wBox[cI].face[i];
        //int fType = wBox[cI].faceType[i];
        if(wBox[cI].face[i] == fc){// && wBox[cI].faceType[i] == type){
            //writeLogNum("logfile5.log","DELETE ", cnt);
            wBox[cI].face[i] = wBox[cI].face[wBox[cI].totalFaces-1];
            wBox[cI].faceType[i] = wBox[cI].faceType[wBox[cI].totalFaces-1];
            wBox[cI].totalFaces--;
            //break;
        }
    }
}

/* 
* Assign wall faces to BDBox
*/
void assignFaceToBDBox(struct wallFace *wF, int noOfF, int type){
    for(int i=0; i<noOfF; i++){   
        int iIndex = floor((wF[i].centroid[0]-xmin)/domainDx);
        int jIndex = floor((wF[i].centroid[1]-ymin)/domainDy);
        int kIndex = floor((wF[i].centroid[2]-zmin)/domainDz);
        int cI = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;
//        if(cI < 0){
        //writeLogNum("wallmesh.log","cell ID ",cI);
//        }
        bdBox[cI].face[bdBox[cI].totalFaces] = i;
        bdBox[cI].faceType[bdBox[cI].totalFaces] = type;
        bdBox[cI].totalFaces += 1;
        
        maxFaces = fmax(maxFaces,bdBox[cI].totalFaces);

        //Assign as a boundary cell
        // for(int k=kIndex-1; k<kIndex+2; k++){
        //     for(int j=jIndex-1; j<jIndex+2; j++){
        //         for(int i=iIndex-1; i<iIndex+2; i++){
        //             cI = i + j*xDiv + k*xDiv*yDiv;
        //             bdBox[cI].isWall = 1;
        //         }
        //     }
        // }
        // real cx = wFace[i].centroid[0]/lengthFactor;
        // if(cx > 0.0048 && cx < 0.0052){
        //     writeLogNum("wallmesh.log","ccent ",cx);
        // }
    }
    // writeLogNum("wallmesh.log","max faces in a BD box ",maxFaces);
    // writeLogNum("wallmesh.log","no of rotor faces ",noOfRotorFaces);
    // writeLogNum("wallmesh.log","no of stator faces ",noOfStatorFaces);
}

/* Deallocate memory*/
void deallocate(){
    if(sizeof(walls) != 0){
        free(walls);
        free(uVec);
        free(ipRVec);
        free(jpRVec);
        free(bdBox);
        free(ijVec);
        free(rotVel);
        free(ipCntPntVel);
        free(jpCntPntVel);
        free(cntPntVel);
        free(demPart);
        free(wFace);
        free(wFaceRot);
        //free(tbdBox);

    } 

}

/* Allocate arrays */
void allocate(){
    //If allocated then dealllocate
    deallocate();

    uVec = allocateDoubleArray(DIM);
    ipRVec = allocateDoubleArray(DIM);
    jpRVec = allocateDoubleArray(DIM);
    ijVec = allocateDoubleArray(DIM);
    //tempVec = allocateDoubleArray(dim);
    rotVel = allocateDoubleArray(DIM);
    ipCntPntVel = allocateDoubleArray(DIM);
    jpCntPntVel = allocateDoubleArray(DIM);
    cntPntVel = allocateDoubleArray(DIM);
    wFace = allocateWallFaceArray(noOfStatorFaces);
    wFaceRot = allocateWallFaceArray(noOfRotorFaces);
    // cfdcell = allocateCFDCell(noOfCFDCells);
    
    for(int i=0; i<noOfStatorFaces; i++){
        wFace[i].bdBoxIndex = -1;
    }
    for(int i=0; i<noOfRotorFaces; i++){
        wFaceRot[i].bdBoxIndex = -1;
    }  

    bdBox = allocateBdBoxArray(xDiv*yDiv*zDiv); //bounding box array
    wBox = allocateWallBoxArray(wBXDiv*wBYDiv*wBZDiv);

    for(int i=0; i<wBXDiv*wBYDiv*wBZDiv; i++){
        wBox[i].totalFaces = 0;
    }
    writeLogLine("logfile2.log"," allocate stage 1\n ");

    for(int i=0; i<xDiv*yDiv*zDiv; i++){
        bdBox[i].noOfParticles = 0;
        bdBox[i].totalFaces = 0;
        //bdBox[i].isWall = 0;
    }
    
    demPart = allocatePar(parArraySize);
    for(int i=0; i<parArraySize; i++){
        demPart[i].dt = 0.0;
        demPart[i].currentTime = 0.0;
        demPart[i].nrmDisp = 0.0;
        demPart[i].noOfNeigh = 0;
        demPart[i].cordNo = 0;
        demPart[i].pos = allocateDoubleArray(DIM);
        demPart[i].angVel = allocateDoubleArray(DIM);
        demPart[i].vel = allocateDoubleArray(DIM);
        demPart[i].hisDisp = allocateDoubleArray(DIM);

        demPart[i].hisDisp[0]  = 0.0;
        demPart[i].hisDisp[1]  = 0.0;
        demPart[i].hisDisp[2]  = 0.0;   

        demPart[i].angVel[0] = 0.0;
        demPart[i].angVel[1] = 0.0;
        demPart[i].angVel[2] = 0.0;

        demPart[i].force = allocateDoubleArray(DIM);
        demPart[i].dragForce = allocateDoubleArray(DIM);
        demPart[i].momentum = allocateDoubleArray(DIM);
        demPart[i].momentum[0] = 0.0;
        demPart[i].momentum[1] = 0.0;
        demPart[i].momentum[2] = 0.0;
        demPart[i].dragForce[0] =0.0;
        demPart[i].dragForce[1] =0.0;
        demPart[i].dragForce[2] =0.0;
        demPart[i].displacement = 0.0;
        demPart[i].insertable = 1;
        demPart[i].active = 1;
        
        demPart[i].noOfCnt = 0;
        demPart[i].noOfWallCnt = 0;
        demPart[i].incontact = 0;
        demPart[i].eCharge = 0.0;//1.03e-4*pow(0.5*largestParDia,1.7);
        demPart[i].noOfPartColl = 0;
        demPart[i].noOfWallColl = 0;
	    demPart[i].maxPartCollE = 0.0;
        demPart[i].maxWallCollE = 0.0;
        //demPart[i].preFluidCell = -1;
        //demPart[i].wallType = 0;
    }

    walls = allocateIntArray(noOfWalls);
    writeLogLine("logfile2.log"," allocate completed\n ");
}

/* Generate cells for the problem domain which is used by particle and fluent cells. 
These cells are used to reduce the number of searches when particle and fluid cells exchange 
information 
*/
void buildDEMCFDCellMap(){
    printf("buildDEMCFDCellMap\n");
}


/*
Add contact faces to bounding box
*/
// void addFaceToBdBox(){
//     int noOfF=0;
//     Domain *d = Get_Domain(1);
    
//     for(int i=0; i<noOfWalls; i++){
//         Thread *tf;
//         face_t f;
//         Node *node;
//         tf = Lookup_Thread(d, walls[i]);
//         int n;
//         begin_f_loop(f, tf)//Loop over faces
//         {
//             real f_cent[ND_ND];
//             F_CENTROID(f_cent,f,tf);
//             int iIndex = floor((f_cent[0]*lengthFactor-xmin)/domainDx);
//             int jIndex = floor((f_cent[1]*lengthFactor-ymin)/domainDy);
//             int kIndex = floor((f_cent[2]*lengthFactor-zmin)/domainDz);
//             int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;
//             noOfF++;
//         }
//         end_f_loop(f,tf)
        
//     }
//     printf("FACES %d\n",noOfF);
// }


/*
Initially neighbourlist array is empty. Fill neighbourlist. 
*/
void initialize(real *nbList, int *parIndex, int *cellSE, int np,
    real *pos, real *parDia)
{
 
}

/*
Assign scale factors
*/
void setReduceUnits()
{
    //Scale factors for reduced units
	refLength = largestParDia;
	refDensity = largestParDensity;
	lengthFactor = 1.0/refLength;
	volumeFactor = pow(lengthFactor,3);
	massFactor = 6.0/(PI*pow(refLength,3)*refDensity);
	timeFactor = sqrt(gravity/refLength);
	densityFactor = 6.0/(PI*refDensity);
	forceFactor = 6.0/(gravity*PI*pow(refLength,3)*refDensity);
	pressureFactor = 6.0/(gravity*PI*refLength*refDensity);
	StressFactor = pressureFactor;
	energyFactor = 6.0/(gravity*PI*pow(refLength,4)*refDensity);
	momentFactor = energyFactor;
	powerFactor = 6.0/(pow(gravity,1.5)*PI*pow((real)refLength,3.5)*refDensity);
	velocityFactor = 1.0/sqrt(refLength*gravity);
	accFactor = 1.0/gravity;
	angVelFactor = sqrt(refLength/gravity);
	angAccFactor = refLength/gravity;
	freqFactor = sqrt(refLength/gravity);
	inertiaFactor = 6.0/(PI*pow(refLength,5)*refDensity);

    cutGap = 1.2*largestParDia*lengthFactor;

    dsmaxCff = sfc*(2.0-pois)/(2.0*(1.0-pois));
    writeLogNum("logfile2.log"," mg in reduced units",(4.0/3.0)*PI*pow(largestParDia*0.5,3)*largestParDensity*massFactor);

    
    dti = 0.0;
    dd = 0.0;
    dsmax = 0.0;

    // scale particle properties
    ymod = ymod*pressureFactor;
    cyldia = cyldia*lengthFactor;
	timeStep = timeStep*timeFactor;
	maxTime	= maxTime*timeFactor;
    //haConst = 1.6E-19;
    haPP = haPP*forceFactor*lengthFactor;
    haPW = haPW*forceFactor*lengthFactor;
    prevCPUTime = 0.;

    elasticMod = ymod/(1.0-pow(pois,2));
   //Find allowed displacment for neighbourlist update
    rIn = largestParDia*lengthFactor;
    rOut = 1.55*rIn; //By definition
    allowedDisp = 0.5*(rOut-rIn);

    for (int i=0; i<np; i++){
        demPart[i].haPp = haPP;
        //demPart[i].haPw = haPW;
        demPart[i].dt = timeStep;
        demPart[i].displacement = 2.0*allowedDisp;
    }

    //Multiply face nodes by lengthFactor
    initializeFaces(wFace, noOfStatorFaces);
    //initializeFaces(wFaceRot, noOfRotorFaces);    

    //Adjust boundary limits
    xmin = xmin*lengthFactor;
    ymin = ymin*lengthFactor;
    zmin = zmin*lengthFactor;
    xmax = xmax*lengthFactor;
    ymax = ymax*lengthFactor;
    zmax = zmax*lengthFactor;

    domainDx = domainDx*lengthFactor;
    domainDy = domainDy*lengthFactor;
    domainDz = domainDz*lengthFactor;

    wBoxDx = wBoxDx*lengthFactor;
    wBoxDy = wBoxDy*lengthFactor;
    wBoxDz = wBoxDz*lengthFactor;

    cellRadius = 0.5*sqrt(domainDx*domainDx+domainDy*domainDy);

    iter = CURRENT_TIMESTEP/(timeStep/timeFactor);
    lamda1 = lamda1*lengthFactor;
    lamda2 = lamda2*lengthFactor;
    rms1 = rms1*lengthFactor;
    rms2 = rms2*lengthFactor;
}

/* Multiply wall face coordinate by lenghFactors*/
void initializeFaces(struct wallFace *wF, int nF){
    for(int i=0; i<nF; i++){
        minMaxEdgeLength(wF[i].node1,wF[i].node2,wF[i].node3, 2);
        wF[i].centroid[0] = wF[i].centroid[0]*lengthFactor;
        wF[i].centroid[1] = wF[i].centroid[1]*lengthFactor;
        wF[i].centroid[2] = wF[i].centroid[2]*lengthFactor;

        wF[i].node1[0] = wF[i].node1[0]*lengthFactor;
        wF[i].node1[1] = wF[i].node1[1]*lengthFactor;
        wF[i].node1[2] = wF[i].node1[2]*lengthFactor;
        wF[i].node2[0] = wF[i].node2[0]*lengthFactor;
        wF[i].node2[1] = wF[i].node2[1]*lengthFactor;
        wF[i].node2[2] = wF[i].node2[2]*lengthFactor;
        wF[i].node3[0] = wF[i].node3[0]*lengthFactor;
        wF[i].node3[1] = wF[i].node3[1]*lengthFactor;
        wF[i].node3[2] = wF[i].node3[2]*lengthFactor;
        

    }
}

void minMaxEdgeLength(real *n1, real *n2, real *n3, int zId){
    // Find min and mx of edge length
    real *tempLength = allocateDoubleArray(DIM);
    vecSub(n1,n2,tempLength);
    real tLength = vecMag(tempLength);
    minEdgeLength = fmin(minEdgeLength,tLength);
    maxEdgeLength = fmax(maxEdgeLength,tLength);

    vecSub(n2,n3,tempLength);
    tLength = vecMag(tempLength);
    minEdgeLength = fmin(minEdgeLength,tLength);
    maxEdgeLength = fmax(maxEdgeLength,tLength);

    vecSub(n3,n1,tempLength);
    tLength = vecMag(tempLength);
    minEdgeLength = fmin(minEdgeLength,tLength);
    maxEdgeLength = fmax(maxEdgeLength,tLength);

    if(zId == rotZoneID){
        minRotEdgeLength = fmin(minRotEdgeLength,tLength);
        maxRotEdgeLength = fmax(maxRotEdgeLength,tLength);
    }       

    free(tempLength); 
}

/* Go through all DEM cells and filter faces within the cell.
    Filtering is done by checking the distance from cell DEM center 
    to projection on face
*/
void refineBdBoxFace(struct wallFace *wF){
    //When searching near faces each cell checks neighbouring cells
    //extent of neighbouring cells is decided by offDiv index
    int offDiv = floor(maxEdgeLength*lengthFactor/domainDx);


    // if(tbdBox !=NULL){
    //     free(tbdBox);
    // }
    tbdBox = allocateTBdBoxArray(xDiv*yDiv*zDiv);
    //Make a temporary copy of from bdBox to tbdBox
    for(int k=0; k<zDiv; k++){
       for(int j=0; j<yDiv; j++){
           for(int i=0; i<xDiv; i++){
               int cI = i + j*xDiv + k*xDiv*yDiv;
               tbdBox[cI].totalFaces = bdBox[cI].totalFaces;
               for(int m=0; m<bdBox[cI].totalFaces; m++){
                   tbdBox[cI].face[m] = bdBox[cI].face[m];
                    tbdBox[cI].faceType[m] = bdBox[cI].faceType[m];
               }
               bdBox[cI].totalFaces = 0;
           }
       }
    }

    //
    //writeLogNum("wallmesh.log","maxEdgeLength ",maxEdgeLength);
    //writeLogNum("wallmesh.log","domainDx ",domainDx);
    writeLogNum("wallmesh.log","div ",offDiv);
    for(int k=0; k<zDiv; k++){
       for(int j=0; j<yDiv; j++){
           for(int i=0; i<xDiv; i++){
                if(i>=0 && j>=0 && k>=0){
                    if(i<xDiv && j<yDiv && k<zDiv){
                        int cI = i + j*xDiv + k*xDiv*yDiv;

                        real *cellCent = allocateDoubleArray(DIM);
                        real *dVec = allocateDoubleArray(DIM);
                        cellCent[0] = i*domainDx + 0.5*domainDx + xmin;
                        cellCent[1] = j*domainDy + 0.5*domainDy + ymin;
                        cellCent[2] = k*domainDz + 0.5*domainDz + zmin;

                        //Go through neighbouring cells for face hunting
                        for(int kk=k-offDiv; kk<k+offDiv+1; kk++){
                            for(int jj=j-offDiv; jj<j+offDiv+1; jj++){
                                for(int ii=i-offDiv; ii<i+offDiv+1; ii++){
                                    if(ii>=0 && jj>=0 && kk>=0){
                                        if(ii<xDiv && jj<yDiv && kk<zDiv){
                                            int cII = ii + jj*xDiv + kk*xDiv*yDiv;
                                            ///if(tbdBox[cII].totalFaces>0){
                                                            //}
                                            for(int m=0; m<tbdBox[cII].totalFaces; m++){  
                                                //writeLogNum("wallmesh.log","tbdBox[cII].totalFaces ",wF[tbdBox[cII].face[m]].node1[0]/lengthFactor);
    
                                                getProjectOnFace(cellCent, wF[tbdBox[cII].face[m]].node1, 
                                                    wF[tbdBox[cII].face[m]].node2,
                                                    wF[tbdBox[cII].face[m]].node3,dVec);
                                                if(vecMag(dVec) < 0.5*sqrt(domainDx*domainDx+
                                                    domainDy*domainDy+domainDz*domainDz)+
                                                    0.5*largestParDia*lengthFactor){
                                                    addFaceToBDBox(tbdBox[cII].face[m], tbdBox[cII].faceType[m], cI);
                                                    // writeLogNum("wallmesh.log","dVec ",vecMag(dVec)/lengthFactor);
                                                    // writeLogNum("wallmesh.log","limit ",(0.5*sqrt(domainDx*domainDx+
                                                    //domainDy*domainDy+domainDz*domainDz)+
                                                    //0.6*largestParDia*lengthFactor)/lengthFactor);

                                                }
                                    
                                            }
                                            
                                            
                                        }
                                    }
                                }
                            }
                        }
                        free(dVec);
                        free(cellCent);
                        // free(fc);
                        // free(fcType);
                        maxFaces = fmax(maxFaces,bdBox[cI].totalFaces); 
                    }//end of i<iDiv && j<yDiv && k<zDiv
                }// end of i>=0 && j>=0 && k>=0 
                          
           }
       }
    }
    //if(tbdBox != NULL){
    free(tbdBox);
    writeLogNum("wallmesh.log","MAXIMUM FACES IN BD BOX (stator)",maxFaces);
    //}
}

/* Assign rotor wall faces to wbdBox 
    Go through all CFD cells and assign face ids to DEM cells
*/
void assignRotFace(int step, struct wallFace *wF, int zID){
  Domain *d = Get_Domain(1);
  cell_t c;
  face_t f;
  int faceCount = 0;
  //zoneID = 12;
  //Get rotor thread
  Thread *rot_t = Lookup_Thread(d,zID);
  //thread_loop_c (rot_t,d)
  //Update node and centroid positions of rotor 
    begin_c_loop(c,rot_t)
      int n;
      c_face_loop(c,rot_t,n) 
      { 
        Thread *tf = C_FACE_THREAD(c,rot_t,n); 
        if(THREAD_TYPE(tf) == THREAD_F_WALL){
          f = C_FACE(c,rot_t,n);
          int m;
          int j=0;
          f_node_loop(f,tf,m)//Loop over face nodes
          {
            Node *node = F_NODE(f,tf,m);
            switch(j){
              case 0:
                wF[faceCount].node1[0] = NODE_X(node)*lengthFactor;
                wF[faceCount].node1[1] = NODE_Y(node)*lengthFactor;
                wF[faceCount].node1[2] = NODE_Z(node)*lengthFactor;
              case 1:
                wF[faceCount].node2[0] = NODE_X(node)*lengthFactor;
                wF[faceCount].node2[1] = NODE_Y(node)*lengthFactor;
                wF[faceCount].node2[2] = NODE_Z(node)*lengthFactor;
              case 2:
                wF[faceCount].node3[0] = NODE_X(node)*lengthFactor;
                wF[faceCount].node3[1] = NODE_Y(node)*lengthFactor;
                wF[faceCount].node3[2] = NODE_Z(node)*lengthFactor;
              }
              j++; 
              //writeLogLine("logfile4.log","face nodes");
          }//end of loop over face nodes

          real x[ND_ND];
          F_CENTROID(x,f,tf);
          wF[faceCount].centroid[0] = x[0]*lengthFactor;
          wF[faceCount].centroid[1] = x[1]*lengthFactor;
          wF[faceCount].centroid[2] = x[2]*lengthFactor;

        //   if(faceCount == 1){
        //       writeLog3Num("wallmesh.log","cen ",wF[faceCount].centroid[0]/lengthFactor,
        //       wF[faceCount].centroid[1]/lengthFactor,wF[faceCount].centroid[2]/lengthFactor);
        //   }
          
          int iIndex = floor((wF[faceCount].centroid[0]-xmin)/wBoxDx);
          int jIndex = floor((wF[faceCount].centroid[1]-ymin)/wBoxDy);
          int kIndex = floor((wF[faceCount].centroid[2]-zmin)/wBoxDz);
          int cI = iIndex + jIndex*wBXDiv + kIndex*wBXDiv*wBYDiv;
          
          
          //Only at the begining assign face to bdBox
          if(step == 1){
            wBox[cI].face[wBox[cI].totalFaces] = faceCount;
            //if(zID == rotZoneID){
                wBox[cI].faceType[wBox[cI].totalFaces] = zID; 
            //}

            wBox[cI].totalFaces += 1;
            wF[faceCount].bdBoxIndex = cI;
            //writeLogNum("wallmesh.log","cI ",cI);
            
            maxRotFaces = fmax(maxRotFaces,wBox[cI].totalFaces);
            //minMaxEdgeLength(wF, faceCount);
          }
          else{
            //delete from previous bdBox
            int prevCI = wF[faceCount].bdBoxIndex; //previous BD box index
            if(prevCI > 0){
                deleteRotFace(faceCount, zID , prevCI); //delete face of type 1 from bdBox
            
                //Add to new bdBox
                addRotFaceToBDBox(faceCount, zID , cI);
                wF[faceCount].bdBoxIndex = cI;
                //writeLogNum("wallmesh.log","delete ",cI);
            }
          }

          faceCount++;
          //writeLogNum("wallmesh.log","Cent ",x[1]);
        }
        
      }
    end_c_loop(c,rot_t)
}

/* Assign wall faces to bdBox 
    Go through all CFD cells and assign face ids to DEM cells
*/
void assignFace(int step, struct wallFace *wF, int zID){
  Domain *d = Get_Domain(1);
  cell_t c;
  face_t f;
  int faceCount = 0;
  //zoneID = 12;
  //Get rotor thread
  Thread *rot_t = Lookup_Thread(d,zID);
  //thread_loop_c (rot_t,d)
  //Update node and centroid positions of rotor 
    begin_c_loop(c,rot_t)
      int n;
      c_face_loop(c,rot_t,n) 
      { 
        Thread *tf = C_FACE_THREAD(c,rot_t,n); 
        if(THREAD_TYPE(tf) == THREAD_F_WALL){
          f = C_FACE(c,rot_t,n);
          int m;
          int j=0;
          f_node_loop(f,tf,m)//Loop over face nodes
          {
            Node *node = F_NODE(f,tf,m);
            switch(j){
              case 0:
                wF[faceCount].node1[0] = NODE_X(node)*lengthFactor;
                wF[faceCount].node1[1] = NODE_Y(node)*lengthFactor;
                wF[faceCount].node1[2] = NODE_Z(node)*lengthFactor;
              case 1:
                wF[faceCount].node2[0] = NODE_X(node)*lengthFactor;
                wF[faceCount].node2[1] = NODE_Y(node)*lengthFactor;
                wF[faceCount].node2[2] = NODE_Z(node)*lengthFactor;
              case 2:
                wF[faceCount].node3[0] = NODE_X(node)*lengthFactor;
                wF[faceCount].node3[1] = NODE_Y(node)*lengthFactor;
                wF[faceCount].node3[2] = NODE_Z(node)*lengthFactor;
              }
              j++; 
              //writeLogLine("logfile4.log","face nodes");
          }//end of loop over face nodes

          real x[ND_ND];
          F_CENTROID(x,f,tf);
          wF[faceCount].centroid[0] = x[0]*lengthFactor;
          wF[faceCount].centroid[1] = x[1]*lengthFactor;
          wF[faceCount].centroid[2] = x[2]*lengthFactor;

        //   if(faceCount == 1){
        //       writeLog3Num("wallmesh.log","cen ",wF[faceCount].centroid[0]/lengthFactor,
        //       wF[faceCount].centroid[1]/lengthFactor,wF[faceCount].centroid[2]/lengthFactor);
        //   }
          
          int iIndex = floor((wF[faceCount].centroid[0]-xmin)/domainDx);
          int jIndex = floor((wF[faceCount].centroid[1]-ymin)/domainDy);
          int kIndex = floor((wF[faceCount].centroid[2]-zmin)/domainDz);
          int cI = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;
          
          
          //Only at the begining assign face to bdBox
          if(step == 1){
            bdBox[cI].face[bdBox[cI].totalFaces] = faceCount;
            bdBox[cI].faceType[bdBox[cI].totalFaces] = zID; 
            bdBox[cI].totalFaces += 1;
            wF[faceCount].bdBoxIndex = cI;
            maxFaces = fmax(maxFaces,bdBox[cI].totalFaces);
          }
          else{
            //delete from previous bdBox
            int prevCI = wF[faceCount].bdBoxIndex; //previous BD box index
            if(prevCI > 0){
                deleteFace(faceCount, zID , prevCI); //delete face of type 1 from bdBox
            
                //Add to new bdBox
                addFaceToBDBox(faceCount, zID , cI);
                wF[faceCount].bdBoxIndex = cI;
                //writeLogNum("wallmesh.log","delete ",cI);
            }
          }

          faceCount++;
          //writeLogNum("wallmesh.log","Cent ",x[1]);
        }
        
      }
    end_c_loop(c,rot_t)
}

/*
Allocate an integer type array
return: int* 
*/
int *allocateIntArray(int size)
{
    int *val = (int*)malloc(size*sizeof(int));
    memset(val,0,size*sizeof(int));
    return val;
}

/*
Allocate a real* type array
return: real* 
*/
real *allocateDoubleArray(int size)
{
    real *val = (real*)malloc(size*sizeof(real));
    memset(val,0.0,size*sizeof(real));
    return val;
}

/*
Allocate a char* type array
return: char* 
*/
char *allocateCharArray(int size)
{
    char *val = (char*)malloc(size*sizeof(char));
    return val;
}

/*
Allocate bounding box type array
return: BdBox*
*/
struct BdBox *allocateBdBoxArray(int size)
{
    struct BdBox *bdB = (struct BdBox*)malloc(size*sizeof(struct BdBox));
    return  bdB;
}

/*
Allocate bounding box type array
return: tempBdBox*
*/
struct tempBdBox *allocateTBdBoxArray(int size)
{
    struct tempBdBox *bdB = (struct tempBdBox*)malloc(size*sizeof(struct tempBdBox));
    return  bdB;
}

/*
Allocate bounding box type array
return: BdBox*
*/
struct wallBox *allocateWallBoxArray(int size)
{
    struct wallBox *wB = (struct wallBox*)malloc(size*sizeof(struct wallBox));
    return  wB;
}

/*
Allocate bounding box type array
return: BdBox*
*/
struct wallFace *allocateWallFaceArray(int size)
{
    struct wallFace *wf = (struct wallFace*)malloc(size*sizeof(struct wallFace));
    return  wf;
}

/*
Allocate particle array
*/
struct demParticle *allocatePar(int np)
{
    struct demParticle *par = (struct demParticle*)malloc(np*sizeof(struct demParticle));
    return par;
}

/*
Allocate cfd eell array
*/
struct cfdCell *allocateCFDCell(int nc)
{
    struct cfdCell *cfdcell = (struct cfdCell*)malloc(nc*sizeof(struct cfdCell));
    return cfdcell;
}


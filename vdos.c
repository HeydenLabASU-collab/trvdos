#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double vector[3];

void vecCopy(double *a, vector *res) {
    (*res)[0]=a[0];
    (*res)[1]=a[1];
    (*res)[2]=a[2];
}

void vecSub(vector a, vector b, vector *res) {
    (*res)[0]=a[0]-b[0];
    (*res)[1]=a[1]-b[1];
    (*res)[2]=a[2]-b[2];
}

double vecDot(vector a, vector b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double vecNorm(vector a) {
    return sqrt(vecDot(a,a));
}

void vecNormalize(vector a, vector *res) {
    double norm;
    norm=vecNorm(a);
    (*res)[0]=a[0]/norm;
    (*res)[1]=a[1]/norm;
    (*res)[2]=a[2]/norm;
}

void vecScale(double s, vector a, vector *res) {
    (*res)[0]=a[0]*s;
    (*res)[1]=a[1]*s;
    (*res)[2]=a[2]*s;
}

void vecCross(vector a, vector b, vector *res) {
    (*res)[0] = a[1] * b[2] - a[2] * b[1];
    (*res)[1] = a[2] * b[0] - a[0] * b[2];
    (*res)[2] = a[0] * b[1] - a[1] * b[0];
}

void vecCenter(vector a, vector b, vector *res) {
    (*res)[0] = (a[0] + b[0]) / 2.0;
    (*res)[1] = (a[1] + b[1]) / 2.0;
    (*res)[2] = (a[2] + b[2]) / 2.0;
}

typedef struct {
    int bondAtomIndices[2];
    int sat0AtomIndices[4];
    int nSat0;
    int sat3AtomIndices[4];
    int nSat3;
    double inertia;
    double logInertia;
    double *wOmegaBuffer;
} t_rotBond;

typedef struct {
    /* constant data */
    int nAtoms;
    double resMass;
    int nRotBonds;
    double *atomMasses;
    double *atomMassesX3;
    t_rotBond *rotBonds;
    /* transient data */
    double inertia[3];
    double inertiaTensor[9];
    double rotAxes[9];
    double angMomLab[3];
    double angMomMol[3];
    double omegaMol[3];
    double omegaLab[3];
    double *COMposBuffer;
    double *COMvelBuffer;
    double *wOmegaBuffer;
    /* accumulating data */
    double sumInertia[3];
    double sumLogInertia[3];
    double *totCorr;
    double *trCorr;
    double *rotCorr;
    double *rotBondCorr;
    int corrCnt;
    /* offset in MD info arrays */
    int offset;
} t_residue;

typedef struct {
    t_residue *residues;
    int nResidues;
    double *totCorr;
} t_residueList;

typedef struct {
    int frameSize;
    double *atomCrd;
    int bufferSize;
    double *atomVelBuffer;
    double *atomVel;
    int nCorr;
    int nAtomsSel;
} t_MDinfo;

typedef struct {
    int *indices;
    long long rank;
} t_dih;

int allocResidueList(t_residueList *residueList,int nRes,int nCorr) {
    residueList->residues=(t_residue*)malloc(nRes*sizeof(t_residue));
    residueList->nResidues=nRes;
    residueList->totCorr=(double*)calloc(nRes*nCorr,sizeof(double));
    return 0;
}

int allocMDinfo(t_MDinfo *MDinfo, int nCorr, int nAtomsSel) {
    int i;

    MDinfo->frameSize=nAtomsSel*3;
    MDinfo->atomCrd=(double*)malloc(MDinfo->frameSize*sizeof(double));
    MDinfo->bufferSize=nCorr*MDinfo->frameSize;
    MDinfo->atomVelBuffer=(double*)malloc(MDinfo->bufferSize*sizeof(double));
    MDinfo->atomVel=MDinfo->atomVelBuffer;
    MDinfo->nCorr=nCorr;
    MDinfo->nAtomsSel=nAtomsSel;
    return 0;
}

int allocResidue(t_residue *res,int nAtoms,double *atomMasses,double resMass,int nAtomsSel,int nCorr) {
    int i,j;
    
    res->nAtoms=nAtoms;
    res->resMass=resMass;
    res->atomMasses=(double*)malloc(nAtoms*sizeof(double));
    res->atomMassesX3=(double*)malloc(nAtoms*3*sizeof(double));
    for(i=0;i<nAtoms;i++) {
        res->atomMasses[i]=atomMasses[i];
        for(j=0;j<3;j++) {
            res->atomMassesX3[i*3+j]=atomMasses[i];
        }
    }
    for(i=0;i<3;i++) {
        res->sumInertia[i]=0.0;
        res->sumLogInertia[i]=0.0;
    }
    res->totCorr=(double*)calloc(nCorr, sizeof(double));
    res->corrCnt=0;
    return 0;
}

int setArrayIndexOffsets(t_residueList *residueList,t_MDinfo *MDinfo) {
    int i;
    int offset=0;

    for(i=0;i<residueList->nResidues;i++) {
        residueList->residues[i].offset=offset;
        offset+=residueList->residues[i].nAtoms*3;
    }
    return 0;
}

int computeTotalVDoS(t_residue *res,t_MDinfo *MDinfo) {
    int i,j,k;
    double mass;
    double *m,*v1,*v2;

    m=res->atomMassesX3;
    v1=MDinfo->atomVel+res->offset;
    for(i=0;i<MDinfo->nCorr;i++) {
        k = ((v1+(i*MDinfo->frameSize))-MDinfo->atomVelBuffer) % MDinfo->bufferSize;
        v2=MDinfo->atomVelBuffer+k;
        for(j=0;j<res->nAtoms*3;j++) {
            res->totCorr[i]+=m[j]*v1[j]*v2[j];
        }
    }
    res->corrCnt++;
 }

int processStep(int tStep,t_MDinfo *MDinfo,double *crds,double *vels,t_residueList *residueList,int nCorr) {
    int i;
    int atomVelBufferShift;

    for(i=0;i<3*MDinfo->nAtomsSel;i++) {
        MDinfo->atomCrd[i]=crds[i];
        MDinfo->atomVel[i]=vels[i];
    }
    /* prep MDinfo->atomVel for next time step and analysis functions*/
    atomVelBufferShift=((tStep+1)%MDinfo->nCorr)*MDinfo->nAtomsSel*3;
    MDinfo->atomVel=MDinfo->atomVelBuffer+atomVelBufferShift;
    
    if(tStep>=nCorr-1) {
        #pragma omp parallel for
        for(i=0;i<residueList->nResidues;i++) {
            computeTotalVDoS(&residueList->residues[i],MDinfo);
        }
    }
    return 0;
}

int postProcess(t_residueList *residueList,int nCorr) {
    int i,j;
    int normFactor;

    for(i=0;i<residueList->nResidues;i++) {
        normFactor=residueList->residues[i].corrCnt;
        for(j=0;j<nCorr;j++) {
            residueList->residues[i].totCorr[j]/=normFactor;
            residueList->totCorr[j]+=residueList->residues[i].totCorr[j];
        }
    }
    normFactor=residueList->nResidues;
    for(j=0;j<nCorr;j++) {
        residueList->totCorr[j]/=normFactor;
    }
    return 0;
}

// int inertiaTensorAngularMomentum(t_residue *res,int tStep,int nCorr) {
//     int i,j;
//     int idx;
//     double m;
//     double *crd;
//     double *vel;

//     idx = tStep % nCorr;
//     if(res->nAtoms==1) {
//         for(i=0;i<9;i++) {
//             res->inertiaTensor[i]=0.0;
//         }
//         for(i=0;i<3;i++) {
//             res->angMomLab[i]=0.0;
//         }
//     } else {
//         for(i=0;i<res->nAtoms;i++) {
//             m=res->atomMasses[i];
//             crd=&res->atomCrdList[i*3];
//             vel=&res->atomVelListBuffer[idx*res->nAtoms*3+i*3];
//             res->inertiaTensor[0]+=m*(crd[1]*crd[1]+crd[2]*crd[2]);
//             res->inertiaTensor[4]+=m*(crd[0]*crd[0]+crd[2]*crd[2]);
//             res->inertiaTensor[8]+=m*(crd[0]*crd[0]+crd[1]*crd[1]);
//             res->inertiaTensor[1]-=m*crd[0]*crd[1];
//             res->inertiaTensor[2]-=m*crd[0]*crd[2];
//             res->inertiaTensor[5]-=m*crd[1]*crd[2];
//             res->inertiaTensor[3]=res->inertiaTensor[1];
//             res->inertiaTensor[6]=res->inertiaTensor[2];
//             res->inertiaTensor[7]=res->inertiaTensor[5];
//             res->angMomLab[0]+=m*(crd[1]*vel[2]-crd[2]*vel[1]);
//             res->angMomLab[1]+=m*(crd[2]*vel[0]-crd[0]*vel[2]);
//             res->angMomLab[2]+=m*(crd[0]*vel[1]-crd[1]*vel[0]);
//         }
//     }
//     return 0;
// }

int applyAxes(double *inertia, double *rotAxesTrans, double *angMomLab, double *angMomMol, double *wOmegaMol,int nAtomsInRes, double *pos, double *vel) {
    int i,j;
    vector omegaMol;
    vector omegaLab;
    vector tmp;
    vector radVel;
    
    /*rotate angular momentum from lab frame into molecular frame*/
    for(i=0;i<3;i++) {
        angMomMol[i]=0.0;
        for(j=0;j<3;j++) {
            angMomMol[i]+=rotAxesTrans[j*3+i]*angMomLab[j];
        }
    }
    /*convert angular momentum into weighted rotational velocity*/
    for(i=0;i<3;i++) {
        wOmegaMol[i]=angMomMol[i]/sqrt(inertia[i]);
    }
    /*convert angular momentum into rotational velocity*/
    for(i=0;i<3;i++) {
        omegaMol[i]=angMomMol[i]/inertia[i];
    }
    /*rotate rotational velocity into lab frame*/
    for(i=0;i<3;i++) {
        omegaLab[i]=0.0;
        for(j=0;j<3;j++) {
            omegaLab[i]+=rotAxesTrans[i*3+j]*omegaMol[j];
        }
    }
    for(i=0;i<nAtomsInRes;i++) {
        vecCopy(&pos[i*3],&tmp);
        vecCross(omegaLab,tmp,&radVel);
        for(j=0;j<3;j++) {
            vel[i*3+j]-=radVel[j];
        }
    }
    return 0;
}

/* qsort int comparison function */
int int_cmp(const void *a, const void *b)
{
        const int *ia = (const int *)a; /* casting pointer types  */
        const int *ib = (const int *)b;
        return *ia  - *ib;
        /* integer comparison: returns negative if b > a
        and positive if a > b */
}

/* qsort comparison function for dihedral index list */
int dih_cmp(const void *a, const void *b)
{
        const t_dih *ia = (const t_dih *)a; /* casting pointer types  */
        const t_dih *ib = (const t_dih *)b;

        long long diff = ia[0].rank  - ib[0].rank;

        if(diff<0) {
            return -1;
        } else if(diff>0) {
            return 1;
        } else {
            return 0;
        }
        /* integer comparison: returns negative if b > a
        and positive if a > b */
}

int getRotBonds(t_residueList *sets,int nSets,int *dihedAtomIndices,int nDih,int *resAtomIdxRange, int nCorr) {
    int i,j,k,l;
    int dih[4];
    int *tmp;
    int min,max;
    long long range;
    int cnt,idx;
    t_dih *dihList;
    t_rotBond *rotBonds;
    int flag;
    int nRotBonds;

    if(nDih>0) {
        dihList=(t_dih*)malloc(nDih*sizeof(t_dih));
        /*ensure order for dihedral list: [X a b Y] with a<b*/
        for(i=0;i<nDih;i++) {
            if(dihedAtomIndices[i*4+1] > dihedAtomIndices[i*4+2]) {
                for(k=0;k<4;k++) {
                    dih[k]=dihedAtomIndices[i*4+k];
                }
                for(k=0;k<4;k++) {
                    dihedAtomIndices[i*4+k]=dih[3-k];
                }
            }
        }
        
        /*find largest index in dihedral list*/
        min=dihedAtomIndices[0];
        max=-1;
        for(i=0;i<nDih*4;i++) {
            if(dihedAtomIndices[i]<min) {
                min=dihedAtomIndices[i];
            }
            if(dihedAtomIndices[i]>max) {
                max=dihedAtomIndices[i];
            }
        }
        range=(long long)(max-min+1);
        /*sort dihedral lists*/
        for(i=0;i<nDih;i++) {
            dihList[i].indices=&dihedAtomIndices[4*i];
            /*computing rank for sorting*/
            /*priority is on atom indices for rotatable bond*/
            dihList[i].rank =((long long)(dihList[i].indices[1]-min+1))*range*range*range;
            dihList[i].rank+=((long long)(dihList[i].indices[2]-min+1))*range*range;
            dihList[i].rank+=((long long)(dihList[i].indices[0]-min+1))*range;
            dihList[i].rank+=((long long)(dihList[i].indices[3]-min+1));
        }
        qsort(dihList,nDih,sizeof(t_dih),dih_cmp);
        tmp=(int*)malloc(nDih*sizeof(int));

        /*counting unique rotatable bonds: easy for sorted dihedrals*/
        i=0;
        j=dihList[i].indices[1];
        k=dihList[i].indices[2];
        cnt=1;
        for(i=1;i<nDih;i++) {
            if(dihList[i].indices[1]!=j || dihList[i].indices[2]!=k) {
                /*new unique rotBond detected*/
                j=dihList[i].indices[1];
                k=dihList[i].indices[2];
                cnt++;
            }
        }
        /*temporarily store all unique rotatable bonds in dihedral list in array 'rotBonds'*/
        rotBonds=(t_rotBond*)malloc(cnt*sizeof(t_rotBond));
        idx=0;
        i=0;
        j=dihList[i].indices[1];
        k=dihList[i].indices[2];
        rotBonds[idx].bondAtomIndices[0]=j;
        rotBonds[idx].bondAtomIndices[1]=k;
        rotBonds[idx].sat0AtomIndices[0]=dihList[i].indices[0];
        rotBonds[idx].nSat0=1;
        rotBonds[idx].sat3AtomIndices[0]=dihList[i].indices[3];
        rotBonds[idx].nSat3=1;
        for(i=1;i<nDih;i++) {
            if(dihList[i].indices[1]!=j || dihList[i].indices[2]!=k) {
                /*new unique rotBond detected*/
                j=dihList[i].indices[1];
                k=dihList[i].indices[2];
                idx++;
                rotBonds[idx].bondAtomIndices[0]=j;
                rotBonds[idx].bondAtomIndices[1]=k;
                rotBonds[idx].sat0AtomIndices[0]=dihList[i].indices[0];
                rotBonds[idx].nSat0=1;
                rotBonds[idx].sat3AtomIndices[0]=dihList[i].indices[3];
                rotBonds[idx].nSat3=1;
            } else {
                /*rotatable bond is the same*/
                /*test if the satellites are the same, otherwise add to list of satellites*/
                /*check satellite index 0*/
                flag=0;
                for(l=0;l<rotBonds[idx].nSat0;l++) {
                    if(dihList[i].indices[0]==rotBonds[idx].sat0AtomIndices[l]) {
                        flag=1;
                        break;
                    }
                }
                if(flag==0) {
                    if(rotBonds[idx].nSat0<4) {
                        rotBonds[idx].sat0AtomIndices[rotBonds[idx].nSat0]=dihList[i].indices[0];
                        rotBonds[idx].nSat0++;
                    } else {
                        printf("ERROR: found too many satellites for rotatable bond\n");
                        printf("       rotatable bond atom indices: %d %d\n",rotBonds[idx].bondAtomIndices[0],rotBonds[idx].bondAtomIndices[1]);
                        return 1;
                    }
                }
                /*check satellite index 3*/
                flag=0;
                for(l=0;l<rotBonds[idx].nSat3;l++) {
                    if(dihList[i].indices[3]==rotBonds[idx].sat3AtomIndices[l]) {
                        flag=1;
                        break;
                    }
                }
                if(flag==0) {
                    if(rotBonds[idx].nSat3<4) {
                        rotBonds[idx].sat3AtomIndices[rotBonds[idx].nSat3]=dihList[i].indices[3];
                        rotBonds[idx].nSat3++;
                    } else {
                        printf("ERROR: found too many satellites for rotatable bond\n");
                        printf("       rotatable bond atom indices: %d %d\n",rotBonds[idx].bondAtomIndices[0],rotBonds[idx].bondAtomIndices[1]);
                        return 1;
                    }
                }
            }
        }
        nRotBonds=idx+1;

        /*just for style, let's sort the satellite indices for each unique rotatable bond*/
        for(i=0;i<nRotBonds;i++) {
            qsort(rotBonds[i].sat0AtomIndices,rotBonds[i].nSat0,sizeof(int),int_cmp);
            qsort(rotBonds[i].sat3AtomIndices,rotBonds[i].nSat3,sizeof(int),int_cmp);
        }

        /*assign unique rotatable bonds to residues (one set of unique rotatable bonds per residue)*/
        /* loops over residues, i.e., sets*/
        // FILE *debug;
        // debug = fopen("/Users/mheyden/Documents/GitHub/HeydenLabASU/MDA-3D-2PT/debug.txt", "w");
        // fprintf(debug, "unique rotatable bonds\n");
        // for(i=0;i<nRotBonds;i++) {
        //     fprintf(debug, "rotBond %d\n", i+1);
        //     fprintf(debug, "bondAtomIndices: %d %d\n", rotBonds[i].bondAtomIndices[0], rotBonds[i].bondAtomIndices[1]);
        //     fprintf(debug, "sat0AtomIndices:");
        //     for(j=0;j<rotBonds[i].nSat0;j++) {
        //         fprintf(debug, " %d", rotBonds[i].sat0AtomIndices[j]);
        //     }
        //     fprintf(debug, "\n");
        //     fprintf(debug, "sat3AtomIndices:");
        //     for(j=0;j<rotBonds[i].nSat3;j++) {
        //         fprintf(debug, " %d", rotBonds[i].sat3AtomIndices[j]);
        //     }
        //     fprintf(debug, "\n");
        // }
        for(i=0;i<nSets;i++) {
            /*min max atom indices of residue i*/
            // fprintf(debug, "Residue %d\n", i+1);
            min=resAtomIdxRange[i*2+0];
            max=resAtomIdxRange[i*2+1];
            // fprintf(debug, "min: %d, max: %d\n", min, max);
            cnt=0;
            /*loop over unique rotatable bonds*/
            for(j=0;j<nRotBonds;j++) {
                flag=0;
                for(k=0;k<2;k++) {
                    l=rotBonds[j].bondAtomIndices[k];
                    if(l<min || l>max) {
                        flag=1;
                        break;
                    }
                }
                for(k=0;k<rotBonds[j].nSat0;k++) {
                    l=rotBonds[j].sat0AtomIndices[k];
                    if(l<min || l>max) {
                        flag=1;
                        break;
                    }
                }
                for(k=0;k<rotBonds[j].nSat3;k++) {
                    l=rotBonds[j].sat3AtomIndices[k];
                    if(l<min || l>max) {
                        flag=1;
                        break;
                    }
                }
                /*true if dihedral is part of residue i*/
                if(flag==0) {
                    tmp[cnt]=j;
                    // fprintf(debug, "rotBond %d %d\n",rotBonds[j].bondAtomIndices[0],rotBonds[j].bondAtomIndices[1]);
                    // fprintf(debug, "sat0AtomIndices:");
                    // for(k=0;k<rotBonds[j].nSat0;k++) {
                    //     fprintf(debug, " %d", rotBonds[j].sat0AtomIndices[k]);
                    // }
                    // fprintf(debug, "\n");
                    // fprintf(debug, "sat3AtomIndices:");
                    // for(k=0;k<rotBonds[j].nSat3;k++) {
                    //     fprintf(debug, " %d", rotBonds[j].sat3AtomIndices[k]);
                    // }
                    // fprintf(debug, "\n");
                    cnt++;
                }
            }
            /*allocate memory for rotatable bonds of residue i in sets->rotBondSets[i]*/
            sets->residues[i].nRotBonds=cnt;
            sets->residues[i].rotBonds=(t_rotBond*)malloc(cnt*sizeof(t_rotBond));
            /*copy rotatable bonds into for residue 'i' into sets->rotBondSets[i].rotBonds */
            for(j=0;j<cnt;j++) {
                l=tmp[j];
                for(k=0;k<2;k++) {
                    sets->residues[i].rotBonds[j].bondAtomIndices[k]=rotBonds[l].bondAtomIndices[k]-min;
                }
                for(k=0;k<rotBonds[l].nSat0;k++) {
                    sets->residues[i].rotBonds[j].sat0AtomIndices[k]=rotBonds[l].sat0AtomIndices[k]-min;
                }
                sets->residues[i].rotBonds[j].nSat0=rotBonds[tmp[j]].nSat0;
                for(k=0;k<rotBonds[l].nSat3;k++) {
                    sets->residues[i].rotBonds[j].sat3AtomIndices[k]=rotBonds[l].sat3AtomIndices[k]-min;
                }
                sets->residues[i].rotBonds[j].nSat3=rotBonds[tmp[j]].nSat3;
                sets->residues[i].rotBonds[j].inertia=0.0;
                sets->residues[i].rotBonds[j].logInertia=0.0;
                sets->residues[i].rotBonds[j].wOmegaBuffer=(double*)malloc(nCorr*sizeof(double));
            }
        }
        // fclose(debug);
    }
    return 0;
}

int analyzeRotBondsInResidue(t_residue *set,double *masses, double *pos, double *vel, int tStep, int nCorr) {
    int i,j,k,l;
    int idx;
    vector b0, b1, sat, satVel, center, axis, proj, perp, mom, angMom;
    int nRotBonds=set->nRotBonds;
    t_rotBond *rotBond;
    double inertia1,inertia2;
    double reducedInertia;
    double angMom1,angMom2;
 
    idx = tStep % nCorr;

    for(i=0;i<nRotBonds;i++) {
        inertia1=0.0;
        inertia2=0.0;
        angMom1=0.0;
        angMom2=0.0;
        rotBond=&set->rotBonds[i];

        /*compute center of rotatable bond*/
        vecCopy(&pos[rotBond->bondAtomIndices[0]*3],&b0);
        vecCopy(&pos[rotBond->bondAtomIndices[1]*3],&b1);
        vecCenter(b0,b1,&center);

        /*compute axis of rotation*/
        vecSub(b1,b0,&axis);
        vecNormalize(axis,&axis);

        /*compute inertia*/
        /*add inertia contributions from satellites with index 0 in dihedral*/
        for(j=0;j<rotBond->nSat0;j++) {
            vecCopy(&pos[rotBond->sat0AtomIndices[j]*3],&sat);
	        vecCopy(&vel[rotBond->sat0AtomIndices[j]*3],&satVel);
            /*coordinate of satellite relative to bond center*/
            vecSub(sat,center,&sat);
            /*project on rotational axis*/
            vecScale(vecDot(sat,axis),axis,&proj);
            /*isolate perpendicular component*/
            vecSub(sat,proj,&perp);
            /*moment of inertia*/
            inertia1+=masses[rotBond->sat0AtomIndices[j]]*vecDot(perp,perp);
            /*momentum*/
            vecScale(masses[rotBond->sat0AtomIndices[j]],satVel,&mom);
            /*angular momentum with respect to axis*/
            vecCross(perp,mom,&angMom);
            /*isolate angular momentum around axis*/
            angMom1+=vecDot(angMom,axis);
        }
        /*add inertia contributions from satellites with index 3 in dihedral*/
        for(j=0;j<rotBond->nSat3;j++) {
            vecCopy(&pos[rotBond->sat3AtomIndices[j]*3],&sat);
	        vecCopy(&vel[rotBond->sat3AtomIndices[j]*3],&satVel);
            /*coordinate of satellite relative to bond center*/
            vecSub(sat,center,&sat);
            /*project on rotational axis*/
            vecScale(vecDot(sat,axis),axis,&proj);
            /*isolate perpendicular component*/
            vecSub(sat,proj,&perp);
            /*moment of inertia*/
            inertia2+=masses[rotBond->sat3AtomIndices[j]]*vecDot(perp,perp);
            /*momentum*/
            vecScale(masses[rotBond->sat3AtomIndices[j]],satVel,&mom);
            /*angular momentum with respect to axis*/
            vecCross(perp,mom,&angMom);
            /*isolate angular momentum around axis*/
            angMom2+=vecDot(angMom,axis);
        }
        /*equivalent to reduced mass of vibrating bond in diatomic molecule*/
        reducedInertia=inertia1*inertia2/(inertia1+inertia2);
        /*accumulate reduced inertia for averaging*/
        rotBond->inertia+=reducedInertia;
        rotBond->logInertia+=log(reducedInertia);
	/*1) convert angular momenta into angular velocities*/
	/*2) compute difference (twist velocity of dihedral angle)*/
	/*3) multiply with square root of "reduced" inertia for weighted rotational velocity*/
        rotBond->wOmegaBuffer[idx]=sqrt(reducedInertia)*(angMom1/inertia1-angMom2/inertia2);
    }
    return 0;
}

int corrRotBonds(t_residueList *sets, double *corr,int start, int nCorr) {
    int i,j,k,l,m;
    int nSets;
    t_residue *set;
    t_rotBond *rotBond;

    nSets=sets->nResidues;
    for(i=0;i<nSets;i++) {
        set=&sets->residues[i];
        for(j=0;j<nCorr;j++) {
            k = start % nCorr;
            l = (k + j) % nCorr;
            for(m=0;m<set->nRotBonds;m++) {
                rotBond=&set->rotBonds[m];
                corr[j*nSets+i]+=rotBond->wOmegaBuffer[k]*rotBond->wOmegaBuffer[l];
            }
        }
    }
    return 0;
}

int corrAtomsRes(int nAtomsRes, int residueIdx, int nRes, double *masses, double *atVelBuffer, double *corr,int start, int nCorr) {
    int i,i2,j,k,l,m,n;
    double tmp;
    
    for(i=0;i<nCorr;i++) {
        i2 = i*nRes+residueIdx;
        k = start % nCorr;
        l = (k + i) % nCorr;
        for(j=0;j<nAtomsRes;j++) {
            m=0+3*(j+nAtomsRes*k);
            n=0+3*(j+nAtomsRes*l);
            tmp =atVelBuffer[m]*atVelBuffer[n];
            m++; n++;
            tmp+=atVelBuffer[m]*atVelBuffer[n];
            m++; n++;
            tmp+=atVelBuffer[m]*atVelBuffer[n];
            corr[i2]+=masses[j]*tmp;
        }
    }
    return 0;
}

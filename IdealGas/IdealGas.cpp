#include <iostream>
#include <math.h>
#include <chrono>
#include <random>
#include <fstream>

using namespace std;


//User-defined constants
const int NUMBER_OF_ATOMS = 100;
const int NUMBER_OF_ITERATIONS = 1000;
const float MIN_GEN_DIST = 0.85;                    //Minimal allowed distance between atoms when they are spawned
const float BOX_SIZE = 10;                       //Size of the box for periodic boundary conditions
const float TIMESTEP = 0.01;                     //dt that we use for integration
const float R_CUTOFF = 3;                        //Cutoff range for Lennard-Jones potential
const char* SAVE_FILE_NAME = "Trajectories.txt"; //Text file in which we save coordinates of all atoms every frame
const char* ATOM_TYPE = "H";                     //Only affects Ovito representation


//Auxilary constants
const float R_MIN = pow(2, 1.0 / 6.0);     //Distance at which there is a potential minimum
const float CUTOFF_STEP = 0.0109438300771; //MUST BE PRECALCULATED//Value of LJ potential derivative at cutoff range


//Random number generators
unsigned seed = chrono::system_clock::now().time_since_epoch() / chrono::microseconds(1) % 1000000;
default_random_engine rnd(seed);
uniform_real_distribution<float> rndPos(0, BOX_SIZE);


//Functions
float LJPotDeriv(float r) {
    if (r > R_CUTOFF) { return 0; }
    return 24 * (pow(r, -7) - 2 * pow(r, -13)) - CUTOFF_STEP; //Subtract the value of function at cutoff range to avoid abrupt force change
}


//Classes
struct Atom {
    float m;
    float rPrev[3];
    float rCurr[3];
    float rNext[3];
    float rAbs[3];
    float a[3];

    Atom(float m, float x, float y, float z) : m(m) {
        rCurr[0] = x;
        rCurr[1] = y;
        rCurr[2] = z;
        for (int i = 0; i != 3; i++) {
            rPrev[i] = rCurr[i];
            rAbs[i] = rCurr[i];
            rNext[i] = 0;
            a[i] = 0;
        }
    }

    void Interact(Atom* other) {
        if (other == this) { return; }
        float imageCoord[3];
        for (int i = 0; i != 3; i++) {
            imageCoord[i] = other->rCurr[i];
        }

        findClosestImageCoord(imageCoord); //Function modifies the array
        float dist = getDistTo(imageCoord[0], imageCoord[1], imageCoord[2]);
        if (dist > R_CUTOFF) { return; }

        float aMod = LJPotDeriv(dist) / m;                    //aMod > 0 means atoms are pulled together
        for (int i = 0; i != 3; i++) {                        //aMod < 0 means atoms are pushed apart
            a[i] += aMod * (imageCoord[i] - rCurr[i]) / dist; //The signs of acceleration projections are all accounted for
        }
    }

    void CalcSummaryAcc(Atom** allAtoms, int size) {
        for (int i = 0; i != size; i++) {
            Interact(allAtoms[i]);
        }
    }

    void VerletPosUpd(float dt) {
        for (int i = 0; i != 3; i++) {
            rNext[i] = 2 * rCurr[i] - rPrev[i] + a[i] * dt * dt; //Use Verlet integration scheme

            rAbs[i] += rNext[i] - rCurr[i]; //Save absolute coordinates, then apply PBC
            while (rNext[i] < 0) {          
                rNext[i] += BOX_SIZE;       
                rCurr[i] += BOX_SIZE;       
            }                              
            while (rNext[i] > BOX_SIZE) {   
                rNext[i] -= BOX_SIZE;       
                rCurr[i] -= BOX_SIZE;       
            }                               

            rPrev[i] = rCurr[i]; //rPrev may be out of box, but it allows to avoid jumps in Verlet scheme
            rCurr[i] = rNext[i]; 
            a[i] = 0;            //Reset summary acceleration to zero
        }

    }

    float* findClosestImageCoord(float* coordArr) {
        float minD[3];
        for (int i = 0; i != 3; i++) {
            minD[i] = abs(rCurr[i] - coordArr[i]);
            if (abs(rCurr[i] - coordArr[i] + BOX_SIZE) < minD[i]) {
                coordArr[i] -= BOX_SIZE;
            }
            if (abs(rCurr[i] - coordArr[i] - BOX_SIZE) < minD[i]) {
                coordArr[i] += BOX_SIZE;
            }
        }
        return minD;
    }

    bool isTooClose(Atom** allAtoms, int size) { 
        for (int i = 0; i != size; i++) { 
            if (getSqrDistTo(allAtoms[i]) < MIN_GEN_DIST * MIN_GEN_DIST) { //Be careful with size you pass, so that the atom
                return true;                                               //doesn't try to calculate the distance to itself
            }
        }
        return false;
    }

    float getDistTo(float x, float y, float z) {
        return sqrt((rCurr[0] - x) * (rCurr[0] - x) + (rCurr[1] - y) * (rCurr[1] - y) + (rCurr[2] - z) * (rCurr[2] - z));
    }
    float getDistTo(Atom* other) {
        return getDistTo(other->rCurr[0], other->rCurr[1], other->rCurr[2]);
    }

    float getSqrDistTo(float x, float y, float z) {
        return (rCurr[0] - x) * (rCurr[0] - x) + (rCurr[1] - y) * (rCurr[1] - y) + (rCurr[2] - z) * (rCurr[2] - z);
    }
    float getSqrDistTo(Atom* other) {
        return getSqrDistTo(other->rCurr[0], other->rCurr[1], other->rCurr[2]);
    }

};


int main() {

    Atom* allAtoms[NUMBER_OF_ATOMS];
    for (int i = 0; i != NUMBER_OF_ATOMS; i++) {                           //TODO: Initialization process of atoms should be optimized later
        allAtoms[i] = new Atom(1, rndPos(rnd), rndPos(rnd), rndPos(rnd));  //
        if (allAtoms[i]->isTooClose(allAtoms, i)) {                        //
            delete allAtoms[i];                                            //
            i--;                                                           //
            cout << "Retry" << "\n";
        }
    }
    cout << "All atoms successfully generated" << "\n";

    ofstream SaveFile(SAVE_FILE_NAME);                              //Open file in which trajectories will be saved
    for (int iter = 0; iter != NUMBER_OF_ITERATIONS; iter++) {      //
                                                                    //
        if (iter % 100 == 0) {                                      //Display iteration counter in console just for convenience
            cout << "Iteration " << iter << "\n";                   //
        }                                                           //
                                                                    //
        SaveFile << NUMBER_OF_ATOMS << "\n\n";                      //This is for Ovito to work propperly
                                                                    //
        for (int i = 0; i != NUMBER_OF_ATOMS; i++) {                //          
                                                                    //
            SaveFile << ATOM_TYPE << " "                            //Save coordinates of every atom to the text file
                     << allAtoms[i]->rCurr[0] << " "                //
                     << allAtoms[i]->rCurr[1] << " "                //
                     << allAtoms[i]->rCurr[2] << "\n";              //
                                                                    //
            allAtoms[i]->CalcSummaryAcc(allAtoms, NUMBER_OF_ATOMS); //Make each atom interact with all other atoms
        }                                                           //
                                                                    //
        for (int i = 0; i != NUMBER_OF_ATOMS; i++) {                //Update every atom's position
            allAtoms[i]->VerletPosUpd(TIMESTEP);
        }

    }
    SaveFile.close();                  //Close the trajectories file
    cout << "Finished saving" << "\n"; //

    cout << "\n";
    for (int i = 0; i != NUMBER_OF_ATOMS; i++) {
        cout << allAtoms[i]->rAbs[0] << " " << allAtoms[i]->rAbs[1] << " " << allAtoms[i]->rAbs[2] << "\n";
    }

}
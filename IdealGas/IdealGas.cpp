#include <iostream>
#include <math.h>
#include <chrono>
#include <random>
#include <fstream>

using namespace std;


//User-defined constants
const int NUMBER_OF_ATOMS = 150;
const int NUMBER_OF_ITERATIONS = 1000;
const float MIN_GEN_DIST = 0.85;                 //Minimal allowed distance between atoms when they are spawned
const float BOX_SIZE = 10;                       //Size of the box for periodic boundary conditions
const float TIMESTEP = 0.01;                     //dt that we use for integration
const float R_CUTOFF = 3;                        //Cutoff range for Lennard-Jones potential

const char* SAVE_FILE_NAME = "Trajectories.txt"; //Text file in which we save coordinates of all atoms every frame
const char* ATOM_TYPE = "H";                     //Only affects Ovito representation


//Auxiliary functions
float LJPot(float r) {
    return 4 * (pow(r, -12) - pow(r, -6));
}

float LJPotDeriv(float r) {
    return 24 * (pow(r, -7) - 2 * pow(r, -13));
}


//Auxiliary constants
const float R_MIN = pow(2, 1.0 / 6.0);                //Distance at which there is a potential minimum
const float ENERGY_CUTOFF_STEP = LJPot(R_CUTOFF);     //Value of LJ potential at cutoff range
const float FORCE_CUTOFF_STEP = LJPotDeriv(R_CUTOFF); //Value of LJ potential derivative at cutoff range


//Random number generators
unsigned seed = chrono::system_clock::now().time_since_epoch() / chrono::microseconds(1) % 1000000;
default_random_engine rnd(seed);
uniform_real_distribution<float> rndPos(0, BOX_SIZE);


//Functions
float getPotEnergy(float r) {
    if (r > R_CUTOFF) { return 0; }
    return 4 * (pow(r, -12) - pow(r, -6)) - ENERGY_CUTOFF_STEP; //Subtract step to prevent abrupt energy change
}

float getForce(float r) {
    if (r > R_CUTOFF) { return 0; }
    return 24 * (pow(r, -7) - 2 * pow(r, -13)) - FORCE_CUTOFF_STEP; //Subtract step to prevent abrupt force change
}


//Classes
struct Vec3D {
    float x, y, z;
    
    Vec3D(float x=0, float y=0, float z=0): x(x), y(y), z(z) {}
    Vec3D(const Vec3D& other): Vec3D(other.x, other.y, other.z) {}

    Vec3D& operator=(const Vec3D& other) {
        x = other.x;
        y = other.y;
        z = other.z;
        return *this;
    }

    Vec3D& operator+=(const Vec3D& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
    Vec3D operator+(const Vec3D& other) {
        Vec3D buff(*this);
        buff += other;
        return buff;
    }Vec3D& operator-=(const Vec3D& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }
    Vec3D operator-(const Vec3D& other) {
        Vec3D buff(*this);
        buff -= other;
        return buff;
    }
    Vec3D& operator*=(float a) {
        x *= a;
        y *= a;
        z *= a;
        return *this;
    }
    Vec3D& operator/=(float a) {
        x /= a;
        y /= a;
        z /= a;
        return *this;
    }
    Vec3D operator*(float a) {
        Vec3D buff(*this);
        buff *= a;
        return buff;
    }Vec3D operator/(float a) {
        Vec3D buff(*this);
        buff /= a;
        return buff;
    }
};
Vec3D operator*(float a, Vec3D vec) {
    return vec * a;
}

struct Atom {
    float m;
    float rPrev[3];
    float rCurr[3];
    float rNext[3];
    float rAbs[3];
    float aPrev[3];
    float aCurr[3];

    Atom(float m, float x, float y, float z) : m(m) {
        rCurr[0] = x;
        rCurr[1] = y;
        rCurr[2] = z;
        for (int i = 0; i != 3; i++) {
            rPrev[i] = rCurr[i];
            rAbs[i] = rCurr[i];
            rNext[i] = 0;
            aCurr[i] = 0;
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

        float aMod = getForce(dist) / m;                          //aMod > 0 means atoms are pulled together
        for (int i = 0; i != 3; i++) {                            //aMod < 0 means atoms are pushed apart
            aCurr[i] += aMod * (imageCoord[i] - rCurr[i]) / dist; //The signs of acceleration projections are all accounted for
        }
    }

    void CalcSummaryAcc(Atom** allAtoms, int size) {
        for (int i = 0; i != size; i++) {
            Interact(allAtoms[i]);
        }
    }

    void VerletPosUpd(float dt) {
        for (int i = 0; i != 3; i++) {
            rNext[i] = 2 * rCurr[i] - rPrev[i] + aCurr[i] * dt * dt; //Use Verlet integration scheme

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
            aCurr[i] = 0;            //Reset summary acceleration to zero
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
                return true;                                               //doesn't try to calculate distance to itself
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
    float tmp = R_MIN - 0.01;

    cout << tmp << "\n";

    cout << getPotEnergy(tmp) << "\n" << getForce(tmp);

    /*
    Atom* allAtoms[NUMBER_OF_ATOMS];
    for (int i = 0; i != NUMBER_OF_ATOMS; i++) {                           //TODO: Initialization process of atoms should be optimized
        allAtoms[i] = new Atom(1, rndPos(rnd), rndPos(rnd), rndPos(rnd));  //
        if (allAtoms[i]->isTooClose(allAtoms, i)) {                        //
            delete allAtoms[i];                                            //
            i--;                                                           //
            cout << "Retry" << "\n";
        }
    }
    cout << "\n" << "All atoms successfully generated" << "\n\n";

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
    cout << "\n" << "Finished saving" << "\n"; //

    cout << "\n";
    for (int i = 0; i != NUMBER_OF_ATOMS; i++) {
        cout << allAtoms[i]->rAbs[0] << " " << allAtoms[i]->rAbs[1] << " " << allAtoms[i]->rAbs[2] << "\n";
    }
    */
}
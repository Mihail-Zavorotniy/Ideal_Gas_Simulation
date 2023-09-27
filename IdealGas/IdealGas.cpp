#include <iostream>
#include <math.h>
#include <chrono>
#include <random>
#include <fstream>

using namespace std;


//User-defined constants
const int NUMBER_OF_ATOMS = 150;
const int NUMBER_OF_ITERATIONS = 1000;
const float MIN_GEN_DIST = 1;                 //Minimal allowed distance between atoms when they are spawned
const float BOX_SIZE = 12;                       //Size of the box for periodic boundary conditions
const float R_CUTOFF = 3;                        //Cutoff range for Lennard-Jones potential
const float TIMESTEP = 0.001;                     //dt that we use for integration

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
    
    explicit Vec3D(float x=0, float y=0, float z=0): x(x), y(y), z(z) {}
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
    Vec3D operator+(const Vec3D& other) const {
        Vec3D buff(*this);
        buff += other;
        return buff;
    }
    Vec3D& operator-=(const Vec3D& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }
    Vec3D operator-(const Vec3D& other) const {
        Vec3D buff(*this);
        buff -= other;
        return buff;
    }
    Vec3D& operator*=(const float a) {
        x *= a;
        y *= a;
        z *= a;
        return *this;
    }
    Vec3D operator*(const float a) const {
        Vec3D buff(*this);
        buff *= a;
        return buff;
    }
    Vec3D& operator/=(const float a) {
        x /= a;
        y /= a;
        z /= a;
        return *this;
    }
    Vec3D operator/(const float a) const {
        Vec3D buff(*this);
        buff /= a;
        return buff;
    }
};
Vec3D operator*(float a, const Vec3D& vec) {
    return vec * a;
}
ostream& operator<<(ostream& os, const Vec3D& vec) {
    os << vec.x << " " << vec.y << " " << vec.z;
    return os;
}
float modSqr(const Vec3D& vec) { //Returns squared module of vector
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z; 
}
float mod(const Vec3D& vec) { //Returns module of vector
    return sqrt(modSqr(vec));
}


struct Atom {
    float m;
    Vec3D rCurr, rPrev, rNext, rDispl;
    Vec3D aCurr, aPrev;

    Atom(float m, float x, float y, float z) : m(m), rCurr(x, y, z), rPrev(rCurr), rNext(), rDispl(), aCurr(), aPrev() {}

    void Interact(Atom* other) {
        if (other == this) { return; } //Do not interact if passed atom is this atom itself

        Vec3D rImage = findClosestImageCoord(other->rCurr); //Function returns radius vector of the closest image
        float dist = getDistTo(rImage);  

        if (dist > R_CUTOFF) { return; } //Do not interact if beyond cutoff range

        float aMod = getForce(dist) / m;         //aMod > 0 means atoms are pulled together
        aCurr += aMod * (rImage - rCurr) / dist; //aMod < 0 means atoms are pushed apart
    }

    void CalcSummaryAcc(Atom** allAtoms, int size) {
        for (int i = 0; i != size; i++) {
            Interact(allAtoms[i]);
        }
    }

    void VerletPosUpd(float dt) {
        rNext = 2 * rCurr - rPrev + aCurr * dt * dt; //Use Verlet integration scheme
        rDispl += rNext - rCurr; //Add dr to displacement vector, then apply PBC

        float stepN;
        stepN = floor(rNext.x / BOX_SIZE); //Applying PBC
        rNext.x -= BOX_SIZE * stepN;       //
        rCurr.x -= BOX_SIZE * stepN;       //
        stepN = floor(rNext.y / BOX_SIZE); //
        rNext.y -= BOX_SIZE * stepN;       //
        rCurr.y -= BOX_SIZE * stepN;       //
        stepN = floor(rNext.z / BOX_SIZE); //
        rNext.z -= BOX_SIZE * stepN;       //
        rCurr.z -= BOX_SIZE * stepN;       //                     

        rPrev = rCurr; //rPrev may be out of box, but it allows to avoid jumps in Verlet scheme
        rCurr = rNext; 

        aPrev = aCurr;   //Save current acceleration
        aCurr = Vec3D(); //Reset current acceleration to zero
    }

    Vec3D findClosestImageCoord(const Vec3D& rOther) {
        Vec3D rImage(rOther);

        if (abs(rCurr.x - (rImage.x - BOX_SIZE)) < abs(rCurr.x - rImage.x)) {
            rImage.x -= BOX_SIZE;
        }
        else if (abs(rCurr.x - (rImage.x + BOX_SIZE)) < abs(rCurr.x - rImage.x)) {
            rImage.x += BOX_SIZE;
        }

        if (abs(rCurr.y - (rImage.y - BOX_SIZE)) < abs(rCurr.y - rImage.y)) {
            rImage.y -= BOX_SIZE;
        }
        else if (abs(rCurr.y - (rImage.y + BOX_SIZE)) < abs(rCurr.y - rImage.y)) {
            rImage.y += BOX_SIZE;
        }

        if (abs(rCurr.z - (rImage.z - BOX_SIZE)) < abs(rCurr.z - rImage.z)) {
            rImage.z -= BOX_SIZE;
        }
        else if (abs(rCurr.z - (rImage.z + BOX_SIZE)) < abs(rCurr.z - rImage.z)) {
            rImage.z += BOX_SIZE;
        }

        return rImage;
    }

    bool isTooClose(Atom** allAtoms, int size) { 
        for (int i = 0; i != size; i++) { 
            if (getSqrDistTo(allAtoms[i]) < MIN_GEN_DIST * MIN_GEN_DIST) { //Be careful with size you pass, so that the atom
                return true;                                               //doesn't try to calculate distance to itself
            }
        }
        return false;
    }

    float getSqrDistTo(const Vec3D& rOther) {
        return modSqr(rCurr - rOther);
    }
    float getSqrDistTo(float x, float y, float z) {
        Vec3D tmp(x, y, z);
        return getSqrDistTo(tmp);
    }
    float getSqrDistTo(Atom* other) {
        return getSqrDistTo(other->rCurr);
    }

    float getDistTo(const Vec3D& rOther) {
        return sqrt(getSqrDistTo(rOther));
    }
    float getDistTo(float x, float y, float z) {
        return sqrt(getSqrDistTo(x, y, z));
    }
    float getDistTo(Atom* other) {
        return sqrt(getSqrDistTo(other));
    }
};


int main() {
    
    Atom* allAtoms[NUMBER_OF_ATOMS];                                       
    for (int i = 0; i != NUMBER_OF_ATOMS; i++) {                           //TODO: Initialization process of atoms should be optimized
        allAtoms[i] = new Atom(1, rndPos(rnd), rndPos(rnd), rndPos(rnd));  //
        if (allAtoms[i]->isTooClose(allAtoms, i)) {                        //
            delete allAtoms[i];                                            //
            i--;                                                           //
            cout << "Retry" << "\n";                                       //
        }
    }
    
    cout << "\n" << "All atoms successfully generated" << "\n\n";

    ofstream SaveFile(SAVE_FILE_NAME);                              //Open file in which trajectories will be saved
    for (int iter = 0; iter != NUMBER_OF_ITERATIONS; iter++) {      
                                                                    
        if (iter % 100 == 0) {                                      //Display iteration counter in console just for convenience
            cout << "Iteration " << iter << "\n";                   //
        }                                                           //
                                                                       
        SaveFile << NUMBER_OF_ATOMS << "\n\n";                      //This is for Ovito to work propperly
                                                                    
        for (int i = 0; i != NUMBER_OF_ATOMS; i++) {                //Iterate through all atoms
            SaveFile << ATOM_TYPE << " "                            //Save coordinates of every atom to the text file
                     << allAtoms[i]->rCurr << "\n";                 //
                                                                    
            allAtoms[i]->CalcSummaryAcc(allAtoms, NUMBER_OF_ATOMS); //Make each atom interact with all other atoms
        }                                                           
                                                                    
        for (int i = 0; i != NUMBER_OF_ATOMS; i++) {                //Update every atom's position
            allAtoms[i]->VerletPosUpd(TIMESTEP);                    //
        }                                                           //
    }                                                               
    SaveFile.close();                                               //Close the trajectories file
    cout << "\n" << "Finished saving" << "\n";                      //

    cout << "\n";
    for (int i = 0; i != NUMBER_OF_ATOMS; i++) {
        cout << allAtoms[i]->rDispl << "\n";
    }
    
}
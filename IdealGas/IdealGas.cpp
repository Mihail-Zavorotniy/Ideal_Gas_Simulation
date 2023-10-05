#include <iostream>
#include <math.h>
#include <chrono>
#include <random>
#include <fstream>

using namespace std;


//User-defined constants
const int NUMBER_OF_ATOMS = 150;                           //TODO: Figure out why system explodes so often
const int NUMBER_OF_ITERATIONS = 1000;
const float MIN_GEN_DIST = 1;                              //Minimal allowed distance between atoms when they are spawned
const float BOX_SIZE = 12;                                 //Size of the box for periodic boundary conditions
const float R_CUTOFF = 3;                                  //Cutoff range for Lennard-Jones potential
const float TIMESTEP = 0.01;                               //dt that we use for integration

const char* TRAJECTORY_FILE = "Trajectories.txt";          //Text file in which we save coordinates of all atoms every frame
const char* VELOCITY_FILE = "Velocities.txt";              //Text file in which we save velocities of all atoms every frame
const char* POT_ENERGY_FILE = "PotentialEnergies.txt";     //Text file in which we save summary energy of system every frame
const char* ATOM_TYPE = "H";                               //Only affects Ovito representation


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
    
    explicit Vec3D(float x=0, float y=0, float z=0) : x(x), y(y), z(z) {}
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
float modSqr(const Vec3D& vec) { //Returns squared module of a vector
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z; 
}
float mod(const Vec3D& vec) { //Returns module of a vector
    return sqrt(modSqr(vec));
}


struct Atom {
    float m;
    Vec3D rCurr, rNext, rDispl;
    Vec3D vCurr;
    Vec3D aCurr, aNext;

    Atom(float m, const Vec3D& r, const Vec3D& v) : m(m), rCurr(r), rNext(), rDispl(), vCurr(v), aCurr(), aNext() {}

    Vec3D getAcc(const Vec3D& r, Atom* other) {     //Returns acceleration that would be caused by other atom if this atom was at position r
        Vec3D rImage = getClosestImageCoord(r, other->rCurr);    //Function returns radius vector of the closest image
        float dist = mod(r - rImage);
        float aMod = getForce(dist) / m;     //aMod > 0 means atoms are pulled together, aMod < 0 means atoms are pushed apart
        return aMod * (rImage - r) / dist;
    }

    Vec3D getSumAcc(Vec3D r, Atom** allAtoms, int size) { //Returns acceleration this atom would experience if it was at position r
        Vec3D sumAcc{};
        for (int i = 0; i != size; i++) {
            if (allAtoms[i] != this) {
                sumAcc += getAcc(r, allAtoms[i]);
            }
        }
        return sumAcc;
    }

    float getSumPotEnergy(Atom** allAtoms, int size) {
        float sumPotEnergy = 0;
        for (int i = 0; i != size; i++) {
            if (allAtoms[i] != this) {
                Vec3D rImage = getClosestImageCoord(rCurr, allAtoms[i]->rCurr);
                sumPotEnergy += getPotEnergy(mod(rCurr - rImage));
            }
        }
        return sumPotEnergy;
    }

    void aCurrUpd(Atom** allAtoms, int size) {
        aCurr = getSumAcc(rCurr, allAtoms, size);
    }

    void verletPosUpd(float dt, Atom** allAtoms, int size) {
        rNext = rCurr + vCurr * dt + 0.5 * aCurr * dt * dt; //Use velocity Verlet scheme
        aNext = getSumAcc(rNext, allAtoms, size);           //
        vCurr += 0.5 * (aCurr + aNext) * dt;                //

        rDispl += rNext - rCurr; //Update displacement vector

        rNext.x -= BOX_SIZE * floor(rNext.x / BOX_SIZE); //Apply PBC
        rNext.y -= BOX_SIZE * floor(rNext.y / BOX_SIZE); //
        rNext.z -= BOX_SIZE * floor(rNext.z / BOX_SIZE); //                   

        rCurr = rNext; //Move position and acceleration a step forward
        aCurr = aNext; //
    }

    Vec3D getClosestImageCoord(const Vec3D& r, const Vec3D& rOther) {
        Vec3D rImage(rOther);

        if (abs(r.x - (rImage.x - BOX_SIZE)) < abs(r.x - rImage.x)) {
            rImage.x -= BOX_SIZE;
        }
        else if (abs(r.x - (rImage.x + BOX_SIZE)) < abs(r.x - rImage.x)) {
            rImage.x += BOX_SIZE;
        }

        if (abs(r.y - (rImage.y - BOX_SIZE)) < abs(r.y - rImage.y)) {
            rImage.y -= BOX_SIZE;
        }
        else if (abs(r.y - (rImage.y + BOX_SIZE)) < abs(r.y - rImage.y)) {
            rImage.y += BOX_SIZE;
        }

        if (abs(r.z - (rImage.z - BOX_SIZE)) < abs(r.z - rImage.z)) {
            rImage.z -= BOX_SIZE;
        }
        else if (abs(r.z - (rImage.z + BOX_SIZE)) < abs(r.z - rImage.z)) {
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
        return getSqrDistTo(Vec3D(x, y, z));
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
    float energyBuff(0);
    
    Atom* allAtoms[NUMBER_OF_ATOMS];                                       
    for (int i = 0; i != NUMBER_OF_ATOMS; i++) {                                  //TODO: Initialization process of atoms should be optimized
        Vec3D r(rndPos(rnd), rndPos(rnd), rndPos(rnd));                           //
        Vec3D v(0, 0, 0);                                                         //TODO: Add gamma distribution for velocities
        allAtoms[i] = new Atom(1, r, v);                                          //
        if (allAtoms[i]->isTooClose(allAtoms, i)) {                               //
            delete allAtoms[i];                                                   //
            i--;                                                                  //
            cout << "Retry" << "\n";                                              //
        }                                                                         //
    }                                                                             //
    for (int i = 0; i != NUMBER_OF_ATOMS; i++) {                                  //
        allAtoms[i]->aCurrUpd(allAtoms, NUMBER_OF_ATOMS);                         //
    }
    
    cout << "\n" << "All atoms successfully generated" << "\n\n";

    ofstream trajFile(TRAJECTORY_FILE);                                           //Open file in which trajectories will be saved
    ofstream velFile(VELOCITY_FILE);                                              //Open file in which velocities will be saved
    ofstream potEnergyFile(POT_ENERGY_FILE);                                      //Open file in which energies will be saved
    for (int iter = 0; iter != NUMBER_OF_ITERATIONS; iter++) {                                               
        if (iter % 100 == 0) {                                                    //Display iteration counter in console just for convenience
            cout << "Iteration " << iter << "\n";                                 //
        }                                                                         //
                                                                       
        trajFile << NUMBER_OF_ATOMS << "\n\n";                                    //This is for Ovito to work propperly
        for (int i = 0; i != NUMBER_OF_ATOMS; i++) {                             
            trajFile << "H " << allAtoms[i]->rCurr << "\n";                                   //Save positions
            velFile << allAtoms[i]->vCurr << "\n";                                            //Save velocities
            potEnergyFile << allAtoms[i]->getSumPotEnergy(allAtoms, NUMBER_OF_ATOMS) << "\n"; //Save potential energies
        }                                                                                     //
        velFile << "\n";                                                                      //
        potEnergyFile << "\n";                                                                //
                                                                             
        for (int i = 0; i != NUMBER_OF_ATOMS; i++) {                              //Iterate through all atoms again       
            allAtoms[i]->verletPosUpd(TIMESTEP, allAtoms, NUMBER_OF_ATOMS);       //Update every atom's position
        }
    }                                                               
    trajFile.close();                                                             //Close the trajectory file
    velFile.close();                                                              //Close the velocity file
    potEnergyFile.close();                                                        //Close the energy file
    cout << "\n" << "Finished saving" << "\n";                                    //

    cout << "\n";
    for (int i = 0; i != NUMBER_OF_ATOMS; i++) {
        cout << allAtoms[i]->rDispl << "\n";
    }
}
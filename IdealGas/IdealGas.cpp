#include <iostream>
#include <math.h>
#include <chrono>
#include <random>
#include <fstream>

using namespace std;


//User-defined constants
const int ATOM_NUM = 250;                               //TODO: Figure out why system explodes so often
const int ITER_NUM = 1000;
const double MIN_GEN_DIST = 1;                          //Minimal allowed distance between atoms when they first created
const double BOX_SIZE = 12;                             //Size of the box for periodic boundary conditions
const double R_CUTOFF = 3;                              //Cutoff range for Lennard-Jones potential
const double TIMESTEP = 0.001;                          //dt that we use for integration
const double ATOM_MASS = 1;                             //For now all atoms have the same mass

const char* TRAJECTORY_FILE = "Trajectories.txt";       //Text file in which we save coordinates of all atoms every frame
const char* VELOCITY_FILE = "Velocities.txt";           //Text file in which we save velocities of all atoms every frame
const char* POT_ENERGY_FILE = "PotentialEnergies.txt";  //Text file in which we save summary energy of system every frame
const char* ATOM_TYPE = "H";                            //Only affects Ovito representation


//Auxiliary functions
double LJPot(double r) {
    return 4 * (pow(r, -12) - pow(r, -6));
}
double LJPotDeriv(double r) {
    return 24 * (pow(r, -7) - 2 * pow(r, -13));
}


//Auxiliary constants
const double R_MIN = pow(2, 1.0 / 6.0);                //Distance at which there is a potential minimum
const double ENERGY_CUTOFF_STEP = LJPot(R_CUTOFF);     //Value of LJ potential at cutoff range
const double FORCE_CUTOFF_STEP = LJPotDeriv(R_CUTOFF); //Value of LJ potential derivative at cutoff range


//Random number generators
unsigned seed = chrono::system_clock::now().time_since_epoch() / chrono::microseconds(1) % 1000000;
default_random_engine rnd(seed);
uniform_real_distribution<double> rndPos(0, BOX_SIZE);


//Functions
double getPotEnergy(double r) {
    if (r > R_CUTOFF) { return 0; }
    return 4 * (pow(r, -12) - pow(r, -6)) - ENERGY_CUTOFF_STEP; //Subtract step to prevent abrupt energy change
}
double getForce(double r) {
    if (r > R_CUTOFF) { return 0; }
    return 24 * (pow(r, -7) - 2 * pow(r, -13)) - FORCE_CUTOFF_STEP; //Subtract step to prevent abrupt force change
}


//Classes
struct Vec3D {
    double x, y, z;
    
    explicit Vec3D(double x=0, double y=0, double z=0) : x(x), y(y), z(z) {}
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
    Vec3D& operator*=(const double a) {
        x *= a;
        y *= a;
        z *= a;
        return *this;
    }
    Vec3D operator*(const double a) const {
        Vec3D buff(*this);
        buff *= a;
        return buff;
    }
    Vec3D& operator/=(const double a) {
        x /= a;
        y /= a;
        z /= a;
        return *this;
    }
    Vec3D operator/(const double a) const {
        Vec3D buff(*this);
        buff /= a;
        return buff;
    }
};
Vec3D operator*(double a, const Vec3D& vec) {
    return vec * a;
}
ostream& operator<<(ostream& os, const Vec3D& vec) {
    os << vec.x << " " << vec.y << " " << vec.z;
    return os;
}

double modSqr(const Vec3D& vec) { //Returns squared module of a vector
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z; 
}
double mod(const Vec3D& vec) { //Returns module of a vector
    return sqrt(modSqr(vec));
}

Vec3D getClosestImage(const Vec3D& r, const Vec3D& rOther) {
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


struct Box {
    int atomN = ATOM_NUM;
    double m[ATOM_NUM];
    Vec3D rCurr[ATOM_NUM], rNext[ATOM_NUM], rDispl[ATOM_NUM];
    Vec3D vCurr[ATOM_NUM];
    Vec3D aCurr[ATOM_NUM], aNext[ATOM_NUM];

    Box() {
        for (int i = 0; i != atomN; i++) {
            m[i] = ATOM_MASS;

            Vec3D r(rndPos(rnd), rndPos(rnd), rndPos(rnd));
            for (int j = 0; j != i; j++) {
                if (mod(r - getClosestImage(r, rCurr[j])) < MIN_GEN_DIST) {
                    cout << "Retry" << "\n";
                    r = Vec3D(rndPos(rnd), rndPos(rnd), rndPos(rnd));
                    j = 0;
                }
            }

            rCurr[i] = r;
            rNext[i] = Vec3D();
            rDispl[i] = Vec3D();
            vCurr[i] = Vec3D();   //TODO: add gamma distribution
            aNext[i] = Vec3D();
        }
        
        calcAcc(rCurr, aCurr);

        cout << "\n" << "All atoms successfully generated" << "\n\n";
    }

    void calcAcc(Vec3D* rArr, Vec3D* aArr) { //Used to calculate both aCurr and aNext
        Vec3D r, rImage, force;
        double dist;

        for (int i = 0; i != atomN; i++) {     //Set acceleration to zero for all atoms
            aArr[i] = Vec3D();                 // 
        }                                      //

        for (int i = 0; i != atomN; i++) {
            r = rArr[i];
            for (int j = i + 1; j != atomN; j++) {
                rImage = getClosestImage(r, rArr[j]);
                dist = mod(r - rImage);
                force = getForce(dist) * (rImage - r) / dist;
                aArr[i] += force / m[i];
                aArr[j] -= force / m[j];
            }
        }
    }

    void verletUpdate(double dt=TIMESTEP) {
        for (int i = 0; i != atomN; i++) {
            rNext[i] = rCurr[i] + vCurr[i] * dt + 0.5 * aCurr[i] * dt * dt;
        }
        calcAcc(rNext, aNext);

        for (int i = 0; i != atomN; i++) {
            vCurr[i] += 0.5 * (aCurr[i] + aNext[i]) * dt;
            rDispl[i] += rNext[i] - rCurr[i];
            rNext[i].x -= BOX_SIZE * floor(rNext[i].x / BOX_SIZE); //Apply PBC
            rNext[i].y -= BOX_SIZE * floor(rNext[i].y / BOX_SIZE); //
            rNext[i].z -= BOX_SIZE * floor(rNext[i].z / BOX_SIZE); //
            rCurr[i] = rNext[i];
            aCurr[i] = aNext[i];
        }
    }

    double getSystemPotEnergy() {
        Vec3D r, rImage;
        double sumE{0};

        for (int i = 0; i != atomN; i++) {
            r = rCurr[i];
            for (int j = i + 1; j != atomN; j++) {
                rImage = getClosestImage(r, rCurr[j]);
                sumE += getPotEnergy(mod(r - rImage));
            }
        }

        return sumE;
    }
};


int main() {
    Box box = Box(); //Create box with all atoms                                                             

    ofstream trajFile(TRAJECTORY_FILE);      //Open file in which trajectories will be saved
    ofstream velFile(VELOCITY_FILE);         //Open file in which velocities will be saved
    ofstream potEnergyFile(POT_ENERGY_FILE); //Open file in which energies will be saved

    for (int iter = 0; iter != ITER_NUM; iter++) {                                               
        if (iter % 100 == 0) {                     //Display iteration counter in console just for convenience
            cout << "Iteration " << iter << "\n";  //
        }                                          //
                                                                       
        trajFile << box.atomN << "\n\n";                       //Save number of atoms for Ovito to work propperly
        for (int i = 0; i != box.atomN; i++) {                 //                               
            trajFile << "H " << box.rCurr[i] << "\n";          //Save positions
            velFile << box.vCurr[i] << "\n";                   //Save velocities
        }                                                      //
        potEnergyFile << box.getSystemPotEnergy() << "\n";     //Save potential energy
        velFile << "\n";                                       //                               
        potEnergyFile << "\n";                                 //                               
                                                                             
        box.verletUpdate(); //Use verlet integration scheme
    }      

    trajFile.close();      //Close the trajectory file
    velFile.close();       //Close the velocity file
    potEnergyFile.close(); //Close the energy file

    cout << "\n" << "Finished saving" << "\n\n";

    for (int i = 0; i != box.atomN; i++) {
        cout << box.rDispl[i] << "\n";
    }
}
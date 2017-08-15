#include <math.h>
#include <string>
#include <vector>
#include "Pythia8/Pythia.h"
#include "HepMC/GenEvent.h"

/*
Method to run Python script to get PDG data and save it in ~/.pdg_data/
This Python script should already be located in that directory. See
get_pdg_data.py for more info.
*/
void save_pdg_data(void);

/*
Particle class to contain relevant properties of particles.

Intended only to be used internally for the BackgroundGenerator.

Members (with accessors defined):
  pid :
    standard particle ID as defined by Particle Data Group
  mass :
    mass of particle (should be up-to-date if local copy of PDG listing
    is recent)
  name :
    non-standard name of particle from PDG listing (should not be used as
    an identifier)
  charge :
    string representation of charged used to compute charge to total particle
    ratio
  px, py, pz :
    components of the 3-momentum of the particle
  charge :
    NOT IMPLEMENTED YET

Methods:
  Particle(int, double, std::string) :
    standard constructor
  Particle(const Particle&) :
    copy constructor
  Particle& operator=(const Particle& that) :
    copy operator
  double e() :
    computes energy of the particle
  void print() :
    pretty prints the particle's defining properties

Additional Notes:


*/
class Particle {

private:

  int pid_;
  double mass_;
  std::string name_;
  std::string charge_;
  double px_ = 0.0;
  double py_ = 0.0;
  double pz_ = 0.0;

public:

  // Constructors
  Particle(int pid, double mass, std::string name, std::string charge);
  Particle(const Particle& that);
  Particle& operator=(const Particle& that);

  // Methods
  double e() { return sqrt(mass_ * mass_ + px_ * px_ + py_ * py_ + pz_ * pz_); }

  // Debug
  void print();

  // Accessors
  int pid() { return pid_; }
  double mass() { return mass_; }
  std::string name() { return name_; }
  std::string charge() { return charge_; }
  double px() { return px_; }
  double py() { return py_; }
  double pz() { return pz_; }
  void pid(int x) { pid_ = x; }
  void mass(double x) { mass_ = x; }
  void name(std::string x) { name_ = std::string(x); }
  void charge(std::string x) { charge_ = std::string(x); }
  void px(double x) { px_ = x; }
  void py(double x) { py_ = x; }
  void pz(double x) { pz_ = x; }

};

/*
Class for a Cumulative Mass Function allowing for easy sampling.

Intended only to be used internally for the BackgroundGenerator.

Members:
  v :
    ordered vector of values in range [0, 1)
  p :
    vector of particles with positions corresponding to values in v
  temp :
    temperature of freezeout used to sample momenta
  charge_ratio :
    ratio of total particles to charged particles

Methods:
  CMF() :
    empty constructor
  CMF(std::vector<double>, std::vector<Particle>, double) :
    standard constructor
  Particle sample(double) :
    returns a particle with sampled momenta based on a value on [0, 1]. For a
    value x, returns a particle corresponding to the index
    where v[index] <= x < v[index + 1] (edge case handle appropriately)
  void print() :
    prints the values of the CMF and their associated particles

Additional Notes:


*/
class CMF {

private:

  std::vector<double> v;
  std::vector<Particle> p;
  double temp;
  double charge_ratio_;

public:

  // Constructors
  CMF();
  CMF(std::vector<double> values, std::vector<Particle> particles,
      double temperature);

  // Methods
  Particle sample(double x);
  double charge_ratio() { return charge_ratio_; }

  // Debug
  void print();

};

/*
Class to initialize and quickly generate heavy-ion backgrounds using the
TRENTO model and a thermal model for the resulting particle multiplicities
and momenta.

Initial construction sets up the TRENTO command.
Successive runs of generate run one instance of TRENTO and take the
entropy to generate a set of particles.

The particle types are sampled from a thermal model.
The particle momenta are sampled from a 3-D Boltzmann distribution.

Members:
  cmd :
    TRENTO command to be run from which entropy is parsed
  f :
    CMF object that will be initialized by the constructor
  pythia :
    Pythia8 object that is used to handle decays of unstable particles
  directory :
    directory where event data will be saved if it is written to file
  degen_factor :
    factor to ensure that two events generated in the same second do not
    overwrite eachother when written to file
  energy_ :
    energy of the event used to scale multiplicity
  proj :
    notation for projectile in event used to scale multiplicity

Methods:
  BackgroundGenerator(std::string, std::string, double, double, std::string,
                      std::string, std::string) :
    most general constructor (all allowable overloads of it are supported)
    non-optional arguments are p1, p2, and energy
  std::string generate() :
    method to generate event and save in HepMC format to file
    filename returned as string
  void generate(HepMC::GenEvent*) :
    method to generate event and store in provided GenEvent object

Additional Notes:
- Current bottleneck is generating TRENTO
- Possible performance improvement obtainable by just running a huge set of
  TRENTO events first and saving the file and sampling from that
- C++ implementation of TRENTO to save IO time?

*/
class BackgroundGenerator {

private:

  std::string cmd;
  CMF f;
  Pythia8::Pythia pythia;
  std::string directory;
  int degen_factor = 0;
  double energy_;
  std::string proj;

  // Private constructor
  void init(std::string p1, std::string p2, double energy, double temp,
            std::string preset, std::string save_dir, std::string options);

  // Private internal methods
  void gen_event(HepMC::GenEvent* evt, int mult);
  int get_event_mult();
  std::string get_filename();

public:

  // Overloaded constructors
  BackgroundGenerator(std::string p1, std::string p2, double energy);
  BackgroundGenerator(std::string p1, std::string p2, double energy,
                      double temp);
  BackgroundGenerator(std::string p1, std::string p2, double energy,
                      std::string preset);
  BackgroundGenerator(std::string p1, std::string p2, double energy,
                      double temp, std::string preset);
  BackgroundGenerator(std::string p1, std::string p2, double energy, 
                      std::string preset, std::string save_dir);
  BackgroundGenerator(std::string p1, std::string p2, double energy,
                      double temp, std::string preset,
                      std::string save_dir, std::string options);

  // Methods
  std::string generate();
  void generate(HepMC::GenEvent* bg_evt);

};

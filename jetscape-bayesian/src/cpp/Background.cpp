/*
TODO:
- Follow style guide more closely
*/


#include "Background.h"
#include <math.h>
#include <unistd.h>
#include <pwd.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/algorithm/string.hpp>
#include "json.hpp"
#include "HepMC/GenEvent.h"
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"


/* --------------------- Helper Function Prototypes ----------------------- */


void sample_momenta(Particle* p, double temp);
std::string exec(const char* cmd);
double particle_density(double temp, double mass, int degen,
                        double chem_pot=0.0, double eps=1.0 * pow(10, -8));
CMF generate_cmf(double temp, bool get_new_data=false);
int compute_mult(double entropy, double scale, double ch_factor);
int bisect(const std::vector<double>& values, double n);


/* ------------------------- BackgroundGenerator -------------------------- */


BackgroundGenerator::BackgroundGenerator(std::string p1, std::string p2,
                                         double energy) {
  init(p1, p2, energy, 0.15, "LHC", "", "");
}

BackgroundGenerator::BackgroundGenerator(std::string p1, std::string p2,
                                         double energy, double temp) {
  init(p1, p2, energy, temp, "LHC", "", "");
}
BackgroundGenerator::BackgroundGenerator(std::string p1, std::string p2,
                                         double energy, std::string preset) {
  init(p1, p2, energy, 0.15, preset, "", "");
}
BackgroundGenerator::BackgroundGenerator(std::string p1, std::string p2,
                                         double energy, double temp,
                                         std::string preset) {
  init(p1, p2, energy, temp, preset, "", "");
}
BackgroundGenerator::BackgroundGenerator(std::string p1, std::string p2,
                                         double energy, std::string preset,
                                         std::string save_dir) {
  init(p1, p2, energy, 0.15, preset, save_dir, "");
}

BackgroundGenerator::BackgroundGenerator(std::string p1, std::string p2,
                                         double energy, double temp,
                                         std::string preset,
                                         std::string save_dir,
                                         std::string options) {
  init(p1, p2, energy, temp, preset, save_dir, options);
}

void BackgroundGenerator::init(std::string p1, std::string p2, double energy,
                          double temp, std::string preset,
                          std::string save_dir, std::string options) {

  // Save directory
  directory = save_dir;

  // Save TRENTO command
  cmd = "trento " + p1 + " " + p2 + " " + options;

  // Generate CMF and save for later use
  f = generate_cmf(temp);

  // Save energy
  energy_ = energy;

  // Seed random number
  srand(time (NULL));

  // Call rand() once to start off randomness
  int a = rand();

  // Set up PYTHIA to do decays
  pythia.readString("Print:quiet = on");
  pythia.readString("ProcessLevel:all = off");
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Next:numberCount = 0");

  // Initialize PYTHIA with given settings
  pythia.init();
}

std::string BackgroundGenerator::generate() {
  
  // Get multiplicity of background event
  int mult = get_event_mult();

  // Create HepMC event
  HepMC::GenEvent* bg_evt = new HepMC::GenEvent();
  bg_evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);

  // Generate the event based on computed multiplicity
  gen_event(bg_evt, mult);

  // Get filename
  std::string filename = get_filename();

  // Save event to file (ascii_out closes as it goes out of scope)
  HepMC::IO_GenEvent ascii_out(directory + filename, std::ios::out);
  ascii_out << bg_evt;

  // Delete GenEvent object
  delete bg_evt;

  return directory + filename;
}

void BackgroundGenerator::generate(HepMC::GenEvent* bg_evt) {

  // Get multiplicity of background event
  int mult = get_event_mult();

  // Generate the event based on computed multiplicity
  gen_event(bg_evt, mult);
}

void BackgroundGenerator::gen_event(HepMC::GenEvent* evt, int mult) {

  // Reset PYTHIA object
  pythia.event.reset();


  // Loop over multiplicity
  for (int i = 0; i < mult; i++) {
    // Sample particle from previously computed CMF
    double n = ((double) rand() / (RAND_MAX));
    Particle a = f.sample(n);

    // Add particle to PYTHIA
    pythia.event.append(a.pid(), 1, 0, 0, a.px(), a.py(),
                        a.pz(), a.e(), a.mass());
  }

  // Run PYTHIA to do the decays
  pythia.next();

  // Convert PYTHIA event to HepMC::GenEvent object
  HepMC::Pythia8ToHepMC converter;
  converter.fill_next_event(pythia, evt);

  // For debugging
  // pythia.event.list();
  // bg_evt->print();

  // Reset PYTHIA object
  pythia.event.reset();
}

int BackgroundGenerator::get_event_mult() {
 
  // Run TRENTO and save output
  std::string trento = exec(cmd.c_str());

  // Parse TRENTO output
  int number;
  double impact_param;
  int npart;
  double entropy;
  double e2;
  sscanf(trento.c_str(), "%d %lf %d %lf %lf", &number, &impact_param,
         &npart, &entropy, &e2);

  // TODO: Get correct scale factor and charged to neutral particle ratio
  double ch_to_n_ratio = f.charge_ratio();

  // Compute multiplicity
  int mult = compute_mult(entropy, 4.65905256, ch_to_n_ratio);

  return mult;
}

std::string BackgroundGenerator::get_filename() {

  // Get time
  std::time_t time = std::time(nullptr);

  // Convert time to string format
  char date_buffer[80];
  struct tm * timeinfo = gmtime(&time);
  strftime(date_buffer, 80, "%Y.%m.%d.%T", timeinfo);
  std::string date_str(date_buffer);

  // Assemble filename
  std::string filename("EventBackground_" + date_str + "-"
                       + std::to_string(degen_factor) + ".hepmc");

  // Increment degeneracy factor to prevent overwriting
  degen_factor += 1;

  return filename;
}

/* ----------------------------- Particle ---------------------------- */

Particle::Particle(int pid, double mass, std::string name,
                   std::string charge) {
  pid_ = pid;
  mass_ = mass;
  name_ = name;
  charge_ = charge;
}

Particle::Particle(const Particle& that) {
  pid_ = that.pid_;
  mass_ = that.mass_;
  name_ = std::string(that.name_);
  charge_ = std::string(that.charge_);
  px_ = that.px_;
  py_ = that.py_;
  pz_ = that.pz_;
}

Particle& Particle::operator=(const Particle& that) {
  if (this != &that) {
    name_ = std::string(that.name_);
    charge_ = std::string(that.charge_);
    mass_ = that.mass_;
    pid_ = that.pid_;
    px_ = that.px_;
    py_ = that.py_;
    pz_ = that.pz_;
  }
  return *this;
}

void Particle::print() {
  printf("pid: %d -- mass: %lf, p: %lf %lf %lf\n", pid_, mass_, px_, py_, pz_);
}

/* ----------------------- Cumulative Mass Function ----------------------- */


CMF::CMF() {

}

CMF::CMF(std::vector<double> values, std::vector<Particle> particles,
    double temperature) {

  v = std::vector<double>(values);
  p = std::vector<Particle>(particles);
  temp = temperature;

  // Compute total to charged particle ratio (before decays)
  charge_ratio_ = 0.0;
  for (int i = 0; i < v.size(); i++) {
    if (p[i].charge() != "0") {
      if (i == v.size() - 1) {
        charge_ratio_ += 1.0 - v[i];
      } else {
        charge_ratio_ += v[i + 1] - v[i];
      }
    }
  }
  if (charge_ratio_ != 0.0) {
    charge_ratio_ = 1.0 / charge_ratio_;
  }
}

Particle CMF::sample(double x) {

  // Get index corresponding to x in v
  int index = bisect(v, x);

  // Create new particle from particle at index in p
  Particle part = Particle(p[index]);

  // Sample particle 3-momentum
  sample_momenta(&part, temp);

  // Return particle with 3-momentum already sampled
  return part;
}

void CMF::print() {

  // Iterate over particles
  for (int i = 0; i < v.size(); i++) {
    Particle a = p[i];
    double val = v[i];

    // Print value and particle properties
    printf("%lf - %d %lf %s\n", val, a.pid(), a.mass(), a.name().c_str());
  }
}

/* ------------------------- Helper Functions -------------------------- */

/*
Function to call Python script (2 or 3) to save data in JSON format
*/
void save_pdg_data(void) {
  FILE* p = popen("python ~/.pdg_data/get_pdg_data.py", "r");
  pclose(p);
}

/*
Function to sample particle 3-momentum from a 3-D Boltzmann distribution
*/
void sample_momenta(Particle* p, double temp) {

  // Get particle mass
  double mass = (*p).mass();

  // Loop until Boltzman distribution is satisfied
  double r1, r2, r3, b_sample, pmag, e, ctheta, stheta, phi;
  do {
    // Sample four random numbers
    r1 = ((double) rand() / (RAND_MAX));
    r2 = ((double) rand() / (RAND_MAX));
    r3 = ((double) rand() / (RAND_MAX));
    b_sample = ((double) rand() / (RAND_MAX));

    // Compute 3-momentum magnitude
    pmag = -1 * temp * log(r1 * r2 * r3);

    // Compute energy
    e = sqrt(pmag * pmag + mass * mass);
  } while (b_sample >= exp((pmag - e) / temp));

  // Compute angles based on random numbers
  ctheta = log(r1 / r2) / log(r1 * r2);
  stheta = sqrt(1 - ctheta * ctheta);
  phi = temp * temp * pow(log(r1 * r2), 2) / (pmag * pmag);
  if (phi > 1.0) {
    printf("phi = %lf, out of range\n", phi);
    // TODO raise error here
  }
  phi *= 2.0 * M_PI;

  // Set particle 3-momentum
  (*p).pz(pmag * ctheta);
  (*p).px(pmag * stheta * cos(phi));
  (*p).py(pmag * stheta * sin(phi));
}

/*
Function to run cmd on shell and save and return output
*/
std::string exec(const char* cmd) {
  std::array<char, 128> buffer;
  std::string result;
  std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
  if (!pipe) throw std::runtime_error("popen() failed!");
  while (!feof(pipe.get())) {
    if (fgets(buffer.data(), 128, pipe.get()) != nullptr) {
      result += buffer.data();
    }
  }
  return result;
}

/*
Function to compute particle density based on properties of 
particle and system
*/
double particle_density(double temp, double mass, int degen, double chem_pot,
                        double eps) {

  // Compute scale factor for sum
  double scale = temp * degen / (2 * pow(M_PI, 2));

  // Set flag based on fermion or boson
  int flag = 1;
  if (degen % 2 == 0) {
    flag = -1;
  }

  // Compute lambda
  // Note: this will be different if chem_pot != 0
  double lam = 1.0;

  // Compute series until sufficiently converged
  double series_sum = 0;
  double term = 0;
  int k = 1;
  do {
    term = pow(flag, k + 1) * pow(lam, 2) * pow(mass, 2) / k;
    // Use boost/math modified Bessel function of second kind
    term *= boost::math::cyl_bessel_k(2, k * mass / temp);
    series_sum += term;
    k += 1;
  } while (std::abs(term / series_sum) > eps);

  // Scale sum and return
  return scale * series_sum;
}

/*
Function to compute CMF for particle types based on thermal model

get_new_data runs the Python script to fetch most recent data from
PDG website. This is a heavy bottleneck and causes poor performance.
Simply, loading the already existing JSON should be sufficient in
most cases.

If you are using the same CMF a lot of times, it may be
good to get_new_data to guarantee that it is up-to-date. For many sampling
runs, the time to fetch the data will be dwarfed by the actual calculation.
*/
CMF generate_cmf(double temp, bool get_new_data) {

  // Initialize empty vectors for particles and values
  std::vector<Particle> particles;
  std::vector<double> values;

  // Density sum for normalization
  double density_sum = 0.0;

  if (get_new_data){
    save_pdg_data();
  }

  // Get system home directory
  const char *homedir;
  if ((homedir = getenv("HOME")) == NULL) {
      homedir = getpwuid(getuid())->pw_dir;
  }
  std::string data_dir = std::string(homedir) + "/.pdg_data/mass_width.json";

  // Open json file
  std::ifstream i(data_dir);

  // Check that it was actually opened
  if (!i.good()) {
    printf("Failed to find mass_width data to load.\n");
    printf("Exiting...\n");
    exit(0);
  }

  // Load file into JSON data object
  nlohmann::json data;
  i >> data;

  // Close file
  i.close();

  // Iterate over all particles in data
  for (int i = 0; i < data.size(); i++) {

    // Get pid and mass
    int id = data[i]["pid"];
    double mass = data[i]["mass"];

    // Skip over partons, some bosons, and K0
    // Draw mass cutoff at 2.5
    if ((std::abs(id) >= 23) && (std::abs(id) != 311) && (mass <= 2.5)) {

      // Get degeneracy based on pid meaning
      int degen = std::abs(id) % 10;

      // Handle K(L) and K(S) separately
      if ((id == 130) || (id == 310)) {
        degen = 1;
      }

      // Create particle
      Particle particle = Particle(id, mass, data[i]["name"],
                                   data[i]["charge"]);

      // Compute density
      double density = particle_density(temp, mass, degen);

      // Push particle and cumulative density to vectors and accumulate
      particles.push_back(particle);
      values.push_back(density_sum);
      density_sum += density;
    }
  }

  // Normalize values
  for (int i = 0; i < values.size(); i++) {
    values[i] = values[i] / density_sum;
  }

  // Create CMF object and return
  return CMF(values, particles, temp);
}

// Compute multiplicity from entropy
// Sample from POISSON ???
int compute_mult(double entropy, double scale, double ch_factor) {
  return round(entropy * scale * ch_factor);
}

// Binary search bisect for performance
int bisect(const std::vector<double>& values, double n) {
  int low = 0;
  int high = values.size();
  int index = 0;
  while (high - low > 1) {
    index = floor((high + low) / 2.0);
    if (values[index] > n) {
      high = index;
    } else if (values[index] < n) {
      low = index;
    } else {
      high = index + 1;
      low = index;
    }
  }
  index = low;
  return index;
}

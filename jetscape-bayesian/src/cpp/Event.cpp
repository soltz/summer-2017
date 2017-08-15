#include "Event.h"
#include <vector>
#include <string>
#include <iostream>
#include <ctime>
#include "HepMC/GenEvent.h"
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"


/* ------------------------- EventGenerator -------------------------- */


EventGenerator::EventGenerator() {
  std::vector<std::string> a;
  init(2760.0, "", a);
}

EventGenerator::EventGenerator(double energy) {
  std::vector<std::string> a;
  init(energy, "", a);
}

EventGenerator::EventGenerator(double energy, double quench) {
  std::vector<std::string> a;
  q = quench;
  init(energy, "", a);
}

EventGenerator::EventGenerator(std::string save_dir) {
  std::vector<std::string> a;
  init(2760.0, save_dir, a);
}

EventGenerator::EventGenerator(const std::vector<std::string>& options) {
  init(2760.0, "", options);
}

EventGenerator::EventGenerator(double energy, std::string save_dir) {
  std::vector<std::string> a;
  init(energy, save_dir, a);
}

EventGenerator::EventGenerator(double energy, double quench,
                               std::string save_dir) {
  std::vector<std::string> a;
  q = quench;
  init(energy, save_dir, a);
}

EventGenerator::EventGenerator(double energy,
                               const std::vector<std::string>& options) {
  init(energy, "", options);
}

EventGenerator::EventGenerator(double energy, double quench,
                               const std::vector<std::string>& options) {
  q = quench;
  init(energy, "", options);
}

EventGenerator::EventGenerator(std::string save_dir,
                               const std::vector<std::string>& options) {
  init(2760.0, save_dir, options);
}

EventGenerator::EventGenerator(double energy, std::string save_dir,
                               const std::vector<std::string>& options) {
  init(energy, save_dir, options);
}

EventGenerator::EventGenerator(double energy, double quench,
                               std::string save_dir,
                               const std::vector<std::string>& options) {
  q = quench;
  init(energy, save_dir, options);
}

void EventGenerator::init(double energy, std::string save_dir,
                          const std::vector<std::string>& options) {

  // Save directory
  directory = save_dir;

  // Set up PYTHIA to generate jet events
  pythia.readString("Print:quiet = on");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PromptPhoton:all = on");
  pythia.readString("Random:setseed = on");
  pythia.readString("Random:seed = 0");
  pythia.readString("Beams:eCM = " + std::to_string(energy));

  // Add custom options
  for (int i = 0; i < options.size(); i++) {
    pythia.readString(options[i]);
  }

  // Initialize PYTHIA with given settings
  pythia.init();

}

std::string EventGenerator::generate() {

  // Create HepMC event
  HepMC::GenEvent* evt = new HepMC::GenEvent();
  evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);

  // Generate jet event
  gen_event(evt, q);  

  // Get filename
  std::string filename = get_filename();

  // Save event to file (ascii_out closes as it goes out of scope)
  HepMC::IO_GenEvent ascii_out(directory + filename, std::ios::out);
  ascii_out << evt;

  // Delete GenEvent object
  delete evt;

  return directory + filename;
}

void EventGenerator::generate(HepMC::GenEvent* evt) {
  gen_event(evt, q);
}

void EventGenerator::gen_event(HepMC::GenEvent* evt, double quench) {
  // Reset PYTHIA object
  pythia.event.reset();

  // Generate jet event
  pythia.next();

  // If jets need to be quenched
  if (quench != 1.0) {

    // Select jet to quench (if photon jet, choose opposite)
    int index = 5 + rand() % 2;
    if (pythia.event[5].id() == 22) {
      index = 6;
    } else if (pythia.event[6].id() == 22) {
      index = 5;
    }

    // Go through quenched jet daughters and quench all of them
    std::vector<int> d = pythia.event[index].daughterListRecursive();
    for (int i = 0; i < d.size(); i++) {
      double px = pythia.event[d[i]].px() * quench;
      double py = pythia.event[d[i]].py() * quench;
      double pz = pythia.event[d[i]].pz() * quench;
      double m = pythia.event[d[i]].m();
      double e = sqrt(m * m + px * px + py * py + pz * pz);
      pythia.event[d[i]].px(px);
      pythia.event[d[i]].py(py);
      pythia.event[d[i]].pz(pz);
      pythia.event[d[i]].e(e);
    }
  }

  // Convert PYTHIA event to HepMC::GenEvent object
  HepMC::Pythia8ToHepMC converter;
  converter.fill_next_event(pythia, evt);

  // For debugging
  // pythia.event.list();
  // evt->print();

  // Reset PYTHIA object
  pythia.event.reset();
}

std::string EventGenerator::get_filename() {

  // Get time
  std::time_t time = std::time(nullptr);

  // Convert time to string format
  char date_buffer[80];
  struct tm * timeinfo = gmtime(&time);
  strftime(date_buffer, 80, "%Y.%m.%d.%T", timeinfo);
  std::string date_str(date_buffer);

  // Assemble filename
  std::string filename("EventJets_" + date_str + "-"
                       + std::to_string(degen_factor) + ".hepmc");

  // Increment degeneracy factor to prevent overwriting
  degen_factor += 1;

  return filename;
}

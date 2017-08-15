#include <vector>
#include <string>
#include "Pythia8/Pythia.h"
#include "HepMC/GenEvent.h"
#include <fastjet/PseudoJet.hh>


/*
Class to initialize and quickly generate heavy-ion jet events using PYTHIA.

Initial construction initializes the PYTHIA object.
Successive runs of generate create an event and populate a HepMC::GenEvent
object.

Members:
  pythia :
    PYTHIA object to generate jet events
  degen_factor :
    factor to ensure that two events generated in the same second do not
    overwrite eachother when written to file
  directory :
    directory where event data will be saved if it is written to file
  q :
    quench factor for quenched jets (non-photon or randomly chosen parton jet)

Methods:
  EventGenerator(double, double, std::string,
                 const std::vector<std::string>&) :
    most general constructor (all allowable overloads of it are supported)
  std::string generate() :
    method to generate event and save in HepMC format to file
    filename returned as string
  void generate(HepMC::GenEvent*) :
    method to generate event and store in provided GenEvent object
  void set_quench(double) :
    method to set the quench factor for future events

Additional Notes:

*/
class EventGenerator {

private:
  
  Pythia8::Pythia pythia;
  int degen_factor = 0;
  std::string directory;
  double q = 1.0;

  // Private constructor
  void init(double energy, std::string save_dir,
            const std::vector<std::string>& options);

  // Private internal methods
  void gen_event(HepMC::GenEvent* evt, double quench);
  std::string get_filename();

public:

  // Overloaded constructors
  EventGenerator();
  EventGenerator(double energy);
  EventGenerator(double energy, double quench);
  EventGenerator(std::string save_dir);
  EventGenerator(const std::vector<std::string>& options);
  EventGenerator(double energy, std::string save_dir);
  EventGenerator(double energy, double quench, std::string save_dir);
  EventGenerator(double energy, const std::vector<std::string>& options);
  EventGenerator(double energy, double quench,
                 const std::vector<std::string>& options);
  EventGenerator(std::string save_dir,
                 const std::vector<std::string>& options);
  EventGenerator(double energy, std::string save_dir,
                 const std::vector<std::string>& options);
  EventGenerator(double energy, double quench, std::string save_dir,
                 const std::vector<std::string>& options);

  // Methods
  std::string generate();
  void generate(HepMC::GenEvent* evt);
  void set_quench(double quench) { q = quench; }

};

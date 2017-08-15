# Dependencies

- (HepMC)[http://hepmc.web.cern.ch/hepmc/]
- (Pythia8)[http://home.thep.lu.se/~torbjorn/Pythia.html]
- (fastjet)[http://fastjet.fr/]
- (Boost)[http://www.boost.org/]
- (JSON)[https://github.com/nlohmann/json]

# Setup

After building HepMC, Pythi8, and fastjet to /usr/local/ and installing Boost and JSON via homebrew (`brew install boost` and `brew install nlohmann/json/nlohmann_json`):

1. Run `mkdir ~/.pdg_data`
2. Move `get_pdg_data.py` and `whitelist.txt` into the directory you just made
3. Run `python get_pdg_data.py` once to create the JSON file

# Compiling

- Must link explicitly against HepMC, Pythia8, and fastjet
- Must set `std=c++11`
- Must compile .cpp file used

# To-do

- Stricter following of style-guide
- Set up some meaningful namespace for the library
- Actual form of library should not change, just some `const` in places where it makes explicit assumptions about things not changing

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include <args.hxx>

double delta = 0.5;
double divisor = 2 / (256.0 * 2 * delta * delta);

#define EXACT 1
#define ENSEMBLE_N_EVENTS 2000
#define CUT 0.4

#if EXACT
double my_exp(double x) { return std::exp(-x / (delta * delta)); }
#else
double my_exp(double x) {
  x = 1 - x * divisor;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  return x;
}
#endif

struct Event {
  double s12;
  double s13;
  double s24;
  double s34;
  double s134;
  double half_mag_squared;
};

const std::vector<Event> read_file(const std::string filename, const size_t n_events, const int subvector_begin, const int subvector_end) {
  std::fstream file(filename, std::ios_base::in);
  if (!file)
    throw std::runtime_error("Error opening file " + filename);

  std::vector<Event> events;
  events.reserve(std::min((size_t)5000000, n_events));

  std::string line;
  while (std::getline(file, line) && events.size() < n_events) {
    std::istringstream iss(line);
    Event e;
    iss >> e.s12 >> e.s13 >> e.s24 >> e.s34 >> e.s134;
    if (iss.fail())
      throw std::runtime_error("Error reading line " + std::to_string(events.size()+1) + " in " + filename);
    e.half_mag_squared = 0.5 * (e.s12*e.s12 + e.s13*e.s13 + e.s24*e.s24 + e.s34*e.s34 + e.s134*e.s134);
    events.push_back(e);
  }
    
  std::vector<Event> subevents;
  subevents.reserve(subvector_end - subvector_begin + 1); 
  
  std::copy(events.begin()+subvector_begin, events.begin()+subvector_end, std::back_inserter(subevents)); 

  return subevents;
}

// WARNING: THIS HAS BEEN MODIFIED TO RETURN EUCLIDEAN DISTANCE
// MAKE SURE TO CHANGE THIS IF RETURNING ANYTHING OTHER THAN THIS
double compute_distance(const std::vector<Event> &data_1, const std::vector<Event> &data_2, const bool upper_only) {
    double total = 0.0;
#pragma omp parallel for reduction(+ : total) schedule(static, 1)
    for (size_t i=0; i < data_1.size(); ++i) {
        auto event_1 = data_1[i];
        for (size_t j=(upper_only ? i+1 : 0); j < data_2.size(); ++j) {
            auto event_2 = data_2[j];
            double x1 = event_1.s13 - event_2.s13; 
            double x2 = event_1.s12 - event_2.s12; 
            double x3 = event_1.s24 - event_2.s24;
            double x4 = event_1.s34 * event_2.s34;
            double x5 = event_1.s134 * event_2.s134;
            double distance =
                sqrt(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5);
            //if (distance < 0.4387862)
            //if (distance > 0.1793)
            if (distance > CUT) 
            //if (distance > 0.6)
            total++; 
            //total = total + distance;
        }
    }
    return total;
}

double compute_statistic(const std::vector<Event> &data_1, const std::vector<Event> &data_2, const bool debug = false) {
  double dist_11 = compute_distance(data_1, data_1, true);
  dist_11 /= data_1.size() * (data_1.size() - 1);
  if (debug) 
      std::cout << "    dist_11 = " << dist_11 << std::endl;

  double dist_22 = compute_distance(data_2, data_2, true);
  dist_22 /= data_2.size() * (data_2.size() - 1);
  if (debug) 
      std::cout << "    dist_22 = " << dist_22 << std::endl;

  double dist_12 = compute_distance(data_1, data_2, false);
  dist_12 /= data_2.size() * data_1.size();
  if (debug)
      std::cout << "    dist_12 = " << dist_12 << std::endl; 
    

  double difference; 
  difference = dist_11+dist_22-dist_12; 
  difference = difference*difference; 
  if (debug)
      std::cout << "    difference in counts (T) is " << difference << std::endl; 

  return difference;
}


double run_energy_test(int argc, char *argv[], int subvector_begin, int subvector_end) {
#if EXACT
  std::cout << "Running in exact mode" << std::endl;
#else
  std::cout << "Running in approximate mode" << std::endl;
#endif

  // Parse the command line arguments
  args::ArgumentParser parser("CPU based energy test");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::Flag permutations_only(parser, "permutations only", "Only calculate permutations", {"permutations-only"});
  args::Flag output_write(parser, "output write", "write output Tvalues", {"output-write"});
  args::ValueFlag<size_t> n_permutations(parser, "n_permutations", "Number of permutations to run", {"n-permutations"});
  args::ValueFlag<size_t> max_events_1(parser, "max events 1", "Maximum number of events to use from dataset 1", {"max-events-1"});
  args::ValueFlag<size_t> max_events_2(parser, "max events 2", "Maximum number of events to use from dataset 2", {"max-events-2"});
  args::ValueFlag<size_t> max_events(parser, "max events", "Max number of events in each dataset", {"max-events"});
  args::ValueFlag<size_t> seed(parser, "seed", "seed for permutations", {"seed"});
  args::ValueFlag<size_t> max_permutation_events_1(parser, "max permutation events 1", "Max number of events in dataset 1 for permutations",
                                                   {"max-permutation-events-1"});

  args::ValueFlag<double> delta_value(parser, "delta value", "delta_value", {"delta-value"});
  args::Positional<std::string> filename_1(parser, "dataset 1", "Filename for the first dataset");
  args::Positional<std::string> filename_2(parser, "dataset 2", "Filename for the second dataset");
  args::Positional<std::string> output_fn(parser, "output filename", "Output filename for the permutation test statistics", {"output-fn"});

  try {
    parser.ParseCLI(argc, argv);
    if (!filename_1 || !filename_2)
      throw args::ParseError("Two dataset filenames must be given");
    if ((max_events_1 || max_events_2) && max_events)
      throw args::ParseError("--max-events cannot be used with --max-events-1 or --max-events-2");
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // set delta
  if (delta_value) {
    delta = args::get(delta_value);
  }
  std::cout << "Distance parameter set to " << delta << std::endl;
  divisor = 2 / (256.0 * 2 * delta * delta);

  // Parse the maximum number of events to use
  size_t data_1_limit = std::numeric_limits<size_t>::max();
  size_t data_2_limit = std::numeric_limits<size_t>::max();

  if (max_events) {
    data_1_limit = args::get(max_events);
    data_2_limit = args::get(max_events);
  } else {
    if (max_events_1)
      data_1_limit = args::get(max_events_1);
    if (max_events_2)
      data_2_limit = args::get(max_events_2);
  }

  // Load the data
  const auto dataset_1 = read_file(args::get(filename_1), data_1_limit, subvector_begin, subvector_end);
  std::cout << "Dataset 1 size is " << dataset_1.size() << std::endl;
  const auto dataset_2 = read_file(args::get(filename_2), data_2_limit, subvector_begin, subvector_end);
  std::cout << "Dataset 2 size is " << dataset_2.size() << std::endl;
  std::cout << std::endl;

  double real_test_statistic = -999;

  if (!permutations_only) {
    // Compute the test statistic for the current dataset
    std::cout << "Calculating test statistic for nominal dataset:" << std::endl;
    real_test_statistic = compute_statistic((std::vector<Event>)dataset_1, (std::vector<Event>)dataset_2, true);
    std::cout << "    T = " << real_test_statistic << std::endl;
  }
  std::cout << std::endl;
  if (n_permutations) {

    // Merge the vectors of events so we can shuffle them
    std::vector<Event> all_events;
    all_events.insert(all_events.end(), dataset_1.begin(), dataset_1.end());
    all_events.insert(all_events.end(), dataset_2.begin(), dataset_2.end());

    size_t N = args::get(n_permutations);

    // Scale the number of events used in permutations
    size_t n_events_1 = dataset_1.size();
    size_t n_events_2 = dataset_2.size();
    if (max_permutation_events_1) {
      n_events_1 = std::min(args::get(max_permutation_events_1), n_events_1);
      n_events_2 = std::round(n_events_1 * ((double)dataset_2.size() / (double)dataset_1.size()));
    }

    double factor = 1.0 * (n_events_1 + n_events_2) / (1.0 * (dataset_1.size() + dataset_2.size()));
    // Set up the random number generator
    int random_seed = std::mt19937::default_seed;
    if (seed)
      random_seed = args::get(seed);
    std::mt19937 random_generator(random_seed);

    // Open the output file
    std::string output_filename;

    if (output_fn) {
      output_filename = args::get(output_fn);
    } else {
      output_filename = "Ts." + std::to_string(dataset_1.size()) + "_";
      output_filename += std::to_string(dataset_2.size()) + "_";
      output_filename += std::to_string(N) + "_";
      output_filename += std::to_string(random_seed) + ".txt";
    }
    
    std::string pvalues_filename;
    pvalues_filename = "pvalues_c" + std::to_string(CUT) + ".txt"; 

    std::ofstream output_file;
    if (output_write)
      output_file.open(output_filename);
    std::ofstream p_output_file;
    if (!permutations_only)
      p_output_file.open(pvalues_filename, std::iostream::out | std::iostream::app);
    int nsig = 0;

    std::cout << "Running " << N << " permutations of " << n_events_1 << " and "
              << n_events_2 << " events using seed " << random_seed << std::endl;
    if (output_write)
      std::cout << "Output filename is " << output_filename << std::endl;

    // Counter to avoid shuffling every time, start large to do a first shuffle
    size_t events_to_skip = all_events.size() + 1;

    for (size_t i = 0; i < N; ++i) {
      if ((i + 1) % std::max((size_t)100, N / 100) == 0)
        std::cout << "Calculating permutation " << i + 1 << " of " << N << std::endl;

      // Reshuffle if we've ran out of events
      if (events_to_skip + n_events_1 + n_events_2 > all_events.size()) {
        std::shuffle(all_events.begin(), all_events.end(), random_generator);
        events_to_skip = 0;
      }

      const std::vector<Event> data_1(all_events.begin() + events_to_skip, all_events.begin() + events_to_skip + n_events_1);
      events_to_skip += n_events_1;
      const std::vector<Event> data_2(all_events.begin() + events_to_skip, all_events.begin() + events_to_skip + n_events_2);
      events_to_skip += n_events_2;

      const double test_statistic = compute_statistic(data_1, data_2);
      if (output_write) {
        output_file << test_statistic << std::endl;
      }
      if (!permutations_only) {
        if (test_statistic > real_test_statistic) {
          nsig++;
          //std::cout << "BOOM!" << std::endl; 
        }
      }
    }
    if (!permutations_only) {
      p_output_file << (1.0 * nsig) / (1.0 * N) << std::endl;
      std::cout << "p value is: " << (1.0 * nsig) / (1.0 * N) << std::endl;
    }
    p_output_file.close();
    output_file.close();
  }

  return 0;
}



int main(int argc, char *argv[]) {

  for (int i = 0; i < 100; i++) {

    std::cout << "\n \nRunning on ensemble " << i << " of 100" << std::endl;  

    //try {
    run_energy_test(argc, argv, i*ENSEMBLE_N_EVENTS, (i+1)*ENSEMBLE_N_EVENTS);
    //} 
  
  
    //catch (std::runtime_error &e) {
     //   std::cerr << e.what() << std::endl;
     
    }
    return 0; 
}


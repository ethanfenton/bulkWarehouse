#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
CharacterMatrix ampSearch(CharacterVector long_strings, CharacterVector short_strings, bool show_progress=false){
  // Initialize a vector of strings to store the results
  vector<string> long_results;
  vector<string> short_results;
  
  // Total iterations for progress calculation
  size_t total_iterations = long_strings.size();
  size_t iteration_count = 0;
  
  // Iterate over each pair of long and short strings
  for (size_t i = 0; i < long_strings.size(); i++) {
    string long_str = as<string>(long_strings[i]);
    for (size_t j = 0; j < short_strings.size(); j++) {
      string short_str = as<string>(short_strings[j]);
      
      // Check if the short string is found in the long string
      if (long_str.find(short_str) != string::npos) {
        // If found, add the pair to the results
        long_results.push_back(long_str);
        short_results.push_back(short_str);
      }
    }
    // Update and print progress
    iteration_count++;
    if (show_progress) {
      double progress = (static_cast<double>(iteration_count) / total_iterations) * 100.0;
      Rcpp::Rcout << "\rProgress: " << std::fixed << std::setprecision(2) << progress << "%" << std::flush;    }
  }
  Rcpp::Rcout << endl; // Print a new line at the end of the process
  
  // Create a CharacterMatrix from the results
  int n = long_results.size();
  CharacterMatrix mat(n, 2);
  for (int i = 0; i < n; i++) {
    mat(i, 0) = long_results[i];
    mat(i, 1) = short_results[i];
  }
  
  return mat;
}

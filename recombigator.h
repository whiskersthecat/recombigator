#ifndef recombigator
#define recombigator

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <set>
#include <map>
#include <sstream>
#include <cstring>
#include <chrono>
#include <experimental/filesystem>
#include <sys/types.h>
#include <sys/stat.h>
#include <iomanip>
#include <chrono>
using namespace std;

void printProgramInfo();
void printCategorizationResults(double time);
void makeDirectories(string base_dir);
void readSNPs(string polymorphism_file);
void readChromosomes(string alignment_genome_Afile, string alignment_genome_Bfile);
void analyzeReads();
void memoryCleanup();
//string parseMD(ifstream* sam_file, vector<int>* mm_loc_ptr, int* len_ptr, int* loc_ptr, int* chrom_ptr, int* num_aa_ptr, string read_name);
void parseCIGAR(stringstream& read, vector<int>& mm_locs, int& len, unsigned long& loc, string& chrom_name);
void categorizeRead(int categorization, string& read_name, string& sam_line_1, string& sam_line_2, string& textoutput, bool silenced);
void createSNPfile(string polymorphism_file);

#endif

// testing

// control
// ./recombigator ./minimap_alignments_variant_calling_UxS/eqx_Alignment_CONTROL_toSalinas_14.SAM ./minimap_alignments_variant_calling_UxS/eqx_Alignment_CONTROL_toUS96UC23_14.SAM ./syri_US_to_SA/syri.out control

// downsampled
// ./recombigator ./minimap_alignments_nanopore_UxS/eqx_Alignment_UxS_downsampled1000_toSalinas14.SAM  ./minimap_alignments_nanopore_UxS/eqx_Alignment_UxS_downsampled1000_toUS96UC23_14.SAM ./syri_US_to_SA/syri.out UxS_downsampled1000

// final
// ./recombigator ./minimap_alignments_nanopore_UxS/eqx_Alignment_UxS_all_toSalinas14.SAM  ./minimap_alignments_nanopore_UxS/eqx_Alignment_UxS_all_toUS96UC23_14.SAM ./syri_US_to_SA/syri.out UxS_all

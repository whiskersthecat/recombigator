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
void parseCIGAR(stringstream& read, vector<int>& mm_locs, int& len, unsigned long& loc, string& chrom_name, bool& dir);
void categorizeRead(int categorization, string& read_name, string& sam_line_1, string& sam_line_2, string& textoutput, bool silenced, int len);
void createSNPfile(string polymorphism_file);
void generatePlots(string base_dir);

#endif

// testing

// control
// ./recombigator ./minimap_alignments_nanopore_UxS/eqx_Alignment_CONTROL_toSalinas_14.SAM ./minimap_alignments_nanopore_UxS/eqx_Alignment_CONTROL_toUS96UC23_14.SAM ./syri_US_to_SA/syri.out control

// UxS downsampled
// ./recombigator ./minimap_alignments_nanopore_UxS/eqx_Alignment_UxS_downsampled1000_toSalinas14.SAM  ./minimap_alignments_nanopore_UxS/eqx_Alignment_UxS_downsampled1000_toUS96UC23_14.SAM ./syri_US_to_SA/syri.out UxS_downsampled1000

// UxS
// ./recombigator ./minimap_alignments_nanopore_UxS/eqx_Alignment_UxS_all_toSalinas14.SAM  ./minimap_alignments_nanopore_UxS/eqx_Alignment_UxS_all_toUS96UC23_14.SAM ./syri_US_to_SA/syri.out UxS_all

// UxS retest recombinants
// ./recombigator ./recombigator_output_on_UxS_all/sorted_alignments/toGenomeA/toGenomeA_Recombinant.SAM  ./recombigator_output_on_UxS_all/sorted_alignments/toGenomeB/toGenomeB_Recombinant.SAM ./syri_US_to_SA/syri.out UxS_recombinants

// PxU_22G86T downsampled
// ./recombigator ./minimap_alignments_nanopore_PxU_22G86T/eqx_Alignment_PxU_22G86T_downsampled1000_toPI251246_14.SAM  ./minimap_alignments_nanopore_PxU_22G86T/eqx_Alignment_PxU_22G86T_downsampled1000_toUS96UC23_14.SAM ./syri_US_to_PI/syri.out PxU_22G86T_downsampled1000

// PxU_22G86T
// ./recombigator ./minimap_alignments_nanopore_PxU_22G86T/eqx_Alignment_PxU_22G86T_all_toPI251246_14.SAM  ./minimap_alignments_nanopore_PxU_22G86T/eqx_Alignment_PxU_22G86T_all_toUS96UC23_14.SAM ./syri_US_to_PI/syri.out PxU_22G86T_all
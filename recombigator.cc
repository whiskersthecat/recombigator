#include "recombigator.h"
#define NUM_CATEGORIES 8
#define NUM_RECOMBINANT_CATEGORIES 3
#define CHOMP_CATEGORIES_INDEX 4

bool*** SNP_coordinates;
int num_chromosomes;

ofstream*** output_SAM_files;
ofstream output_TXT_files[NUM_CATEGORIES + 1];
ofstream summary_file;
ofstream log_file;
ofstream recombinantinfo_files[NUM_RECOMBINANT_CATEGORIES];
ofstream readlength_files[NUM_CATEGORIES + 1];
ofstream error_file;

ifstream sam_files[2];

unsigned long long totalErrors = 0;
unsigned long long totalErrorBP = 0;

map<string, int> chromosome_numbers;
map<int, string> chromosome_names[2];

static string SNP_symbols[] = {"A","B","!","?"};
static string SNP_meanings[] = {"GenotypeA","GenotypeB","BothMismatch","NeitherMismatch"};
static char ALIGN_DIRECTIONS[] = {'f','r'};

static string CATEGORY_NAMES[] = {"GenotypeA", "GenotypeB", "LowQualityRecombinant", "MediumQualityRecombinant", "HighQualityRecombinant","ChompBadAlign", "ChompNoSNPs", "ChompBadGT"};
int num_reads_in_category[NUM_CATEGORIES];

int main(int argc, char* argv[]) {
    if(argc < 4) {
        printProgramInfo();
    }
    string dir_string = "recombigator_output";
    if(argc > 4) {
        string append = argv[4];
        dir_string += "_on_" + append;
    }

    makeDirectories(dir_string);
    readChromosomes(argv[1], argv[2]);
    readSNPs(argv[3]);

    auto begin = chrono::high_resolution_clock::now();
    analyzeReads();
    cout << "[recombigator] Done analyzing reads" << endl;

    auto end = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::microseconds>(end - begin).count();

    printCategorizationResults((double)elapsed / (double) 1000000);

    memoryCleanup();

    generatePlots(dir_string);
}

void categorizeRead(int categorization, string& read_name, string& sam_line_A, string& sam_line_B, string& textoutput, bool silenced, int len) {

    // 1. print textoutput to standard output, if not silenced
    textoutput.append(" ***** Read is " + CATEGORY_NAMES[categorization] + "\n");
    if(!silenced) cout << textoutput << endl;

    // 2. increment counter based on categorization
    num_reads_in_category[categorization]++;

    // 3. write the SAM files with the read
    *output_SAM_files[0][categorization] << read_name << sam_line_A << endl;
    *output_SAM_files[1][categorization] << read_name << sam_line_B << endl;

    // 4. write the text file with the read and the all reads text file
    output_TXT_files[categorization] << textoutput << endl;
    output_TXT_files[NUM_CATEGORIES] << textoutput << endl;

    // 5. write the length of the read to its corresponding file
    readlength_files[categorization] << len << endl;
    readlength_files[NUM_CATEGORIES] << len << endl;

    return;
}

void analyzeReads() {
    
    // the header fields are already consumed when finding the chromosome sizes

    cout << "[recombigator] Analyzing reads" << endl;

    string current_read_name = "";
    string next_read_name = "";

    sam_files[0] >> next_read_name;
    sam_files[1] >> next_read_name;

    int read_number = 0;
    
    while(sam_files[0]) {

        // 1. get mismatch locations for both alignments
        vector<int> mmlocations[2];
        int align_len[2];
        unsigned long align_loc[2];
        unsigned long align_end_loc[2];
        int align_chrom[2];
        string align_chrom_name[2];
        int n_alternate_alignments[2];
        string SAM_line[2];
        bool align_dir[2];

        vector<int> expectedSNPlocations[2];
        vector<bool> genotypeMatrix[2];
        vector<int> concensusMatrix;
        int num_gt_switches = -1;

        string textoutput;

        bool silenced = true;
        if(read_number % 100 == 0)
            silenced = false;
        read_number++;

        current_read_name = next_read_name;
        textoutput.append("#### [" + to_string(read_number) + "] Processing Read " + current_read_name + '\n');

        string line;
        for(int i = 0; i < 2; i++) {
            getline(sam_files[i], SAM_line[i]);
            stringstream read(SAM_line[i]);

            parseCIGAR(read, (mmlocations[i]), align_len[i], align_loc[i], align_chrom_name[i], align_dir[i]);
            align_end_loc[i] = align_loc[i] + align_len[i];
            align_chrom[i] = chromosome_numbers[align_chrom_name[i]];

            n_alternate_alignments[i] = 0;
            sam_files[i] >> next_read_name;
            while(next_read_name == current_read_name && sam_files[i]) { // stop looking for alternates if a new read name comes up or end of file reached
                // alternate alignment
                n_alternate_alignments[i]++;
                getline(sam_files[i], line);
                sam_files[i] >> next_read_name;
            }
        }

        // 2. create SNP matrix for both alignments
        for(int i = 0; i < 2; ++i) {
            auto mismatch_location_it = mmlocations[i].begin();
            // stores the last checked index for a mismatch
            for(unsigned long pos = align_loc[i]; pos < align_end_loc[i]; ++pos) {
                if(SNP_coordinates[i][align_chrom[i]][pos] == 1) {
                    int relative_align_pos = pos - align_loc[i];
                    // Add this to the list of SNP locations
                    expectedSNPlocations[i].push_back(relative_align_pos);
                    bool SNP_exists = false;
                    // Check if Mismatch position exists
                    while(true) {
                        if(mismatch_location_it == mmlocations[i].end() || *mismatch_location_it > relative_align_pos) {
                            break;
                        }        
                        if(*mismatch_location_it == relative_align_pos) {
                            SNP_exists = true;
                            break;
                        }
                        mismatch_location_it++;
                    }
                    if(SNP_exists)
                        genotypeMatrix[i].push_back(i ^ true);
                    else
                        genotypeMatrix[i].push_back(i);
                }
            }
        }

        // 3. create consensus SNP matrix
        vector<bool>::iterator matrix_it[2];
        int numSNPs[2] = {0};
        int totalSNPs = 0;
        int numIgnoredSNPs = 0;
        for(int i = 0; i < 2; i++) matrix_it[i] = genotypeMatrix[i].begin();
        int previous_gt = -1, breakpoint = -1, index = 0, start_gt = -1;
        while(matrix_it[0] != genotypeMatrix[0].end() && matrix_it[1] != genotypeMatrix[1].end()) {

            if(*matrix_it[0] == *matrix_it[1]) {
                // A) both agree on categorization
                int SNP_categorization = *matrix_it[0];
                concensusMatrix.push_back(SNP_categorization);
                if(start_gt == -1)
                    start_gt = SNP_categorization;
                if(previous_gt != SNP_categorization) {
                    previous_gt = SNP_categorization;
                    num_gt_switches++;
                    breakpoint = totalSNPs;
                }
                ++numSNPs[SNP_categorization];
                ++totalSNPs;
            } else if(*matrix_it[0] == 1 && *matrix_it[1] == 0) {
                // B) both have a mismatch, ignore this SNP
                ++numIgnoredSNPs;
                concensusMatrix.push_back(2);
            } else {
                // C) neither have a mismatch, ignore this SNP
                ++numIgnoredSNPs;
                concensusMatrix.push_back(3);
            }
            ++matrix_it[0]; ++ matrix_it[1]; ++ index;
        }

        // 4. Write to file
        for(int i = 0; i < 2; i++) {
            textoutput.append("  >> Aligned to chromosome " + to_string(align_chrom[i] + 1)  + " (called " + align_chrom_name[i] + ") from position " + to_string(align_loc[i]) + " to " + to_string(align_end_loc[i]) + " (Total Length " + to_string(align_len[i]) + ")");
            if(n_alternate_alignments[i] > 0)
                textoutput.append(" Note: " + to_string(n_alternate_alignments[i]) + " alternate alignment(s)");
            textoutput.append("\n");

            textoutput.append("Mismatch Locations     (Total " + to_string(mmlocations[i].size()) + "):");
            for(auto it = mmlocations[i].begin(); it != mmlocations[i].end(); ++it) {
                textoutput.append('\t' + to_string(*it));
            } textoutput.append("\n");

            textoutput.append("Expected SNP Locations (Total " + to_string(expectedSNPlocations[i].size()) + "):");
            for(auto it = expectedSNPlocations[i].begin(); it != expectedSNPlocations[i].end(); ++it) {
               textoutput.append('\t' + to_string(*it));
            } textoutput.append("\n");

            textoutput.append("Genotype Matrix:");
            for(auto it = genotypeMatrix[i].begin(); it != genotypeMatrix[i].end(); ++it) {
                textoutput.append("\t[" + (SNP_symbols[*it]) + "]");
            } textoutput.append("\n");
        }
        textoutput.append(" **** Concensus Matrix ");
        if(genotypeMatrix[0].size() != genotypeMatrix[1].size()) {
            // log_file << "[anomolous case] read " << current_read_name << " has genotype matrixes with different number of SNPs" << endl;
            textoutput.append("(Warning: different GT Matrix sizes)");
        }
        textoutput.append(":");
        for(auto it = concensusMatrix.begin(); it != concensusMatrix.end(); ++it) {
            textoutput.append("\t" + (SNP_symbols[*it]) + "");
        } textoutput.append("\n");

        // 5. categorization
        // #CHOMP CATEGORY 1
        if(align_chrom[0] != align_chrom[1]) {
            categorizeRead(CHOMP_CATEGORIES_INDEX + 1, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced, align_len[1]);
            continue;
        }
        // #CHOMP CATEGORY 2
        if(totalSNPs < 1) {
            categorizeRead(CHOMP_CATEGORIES_INDEX + 2, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced, align_len[1]);
            continue;
        }
        // #CHOMP CATEGORY 3
        if(num_gt_switches > 1) {
            categorizeRead(CHOMP_CATEGORIES_INDEX + 3, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced, align_len[1]);
            continue;
        }

        // Check for alternate alignments
        if(n_alternate_alignments[0] > 0 || n_alternate_alignments[1] > 0) log_file << "[anomolous case] Read " << current_read_name << " (which is not thrown out) has alternate alignment(s)" << endl;

        // #GENOTYPE A or GENOTYPE B
        if(num_gt_switches == 0) {
            // 1. Count Errors
            totalErrors += mmlocations[1].size() - expectedSNPlocations[1].size();
            totalErrorBP += align_len[1];

            // 2. Categorize
            categorizeRead(start_gt, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced, align_len[1]);
            continue;
        }
        
        // ### RECOMBINANT
        int quality = 0;

        // #LOW QUALITY RECOMBINANT
        if(breakpoint == 1 || breakpoint == totalSNPs - 1) quality = 0;

        // #MEDIUM QUALITY RECOMBINANT
        else if(breakpoint == 2 || breakpoint == totalSNPs - 2) quality = 1;

        // #HIGH QUALITY RECOMBINANT
        else quality = 2;

        // write to recombinant info
        double relativeBreakpoint = (double) breakpoint / (double) totalSNPs;
        
        // name alignmentSNPs(A:B)  numSNPs(A:B)    relativeBreakpoint  numIgnoredSNPs  numAlternateAlignments(A:B)   firstSNP    alignDir(A:B)   chromosome(A:B)   len(A:B)    alignloc(A:B)
        recombinantinfo_files[quality] << current_read_name << '\t' << genotypeMatrix[0].size() << ':' << genotypeMatrix[1].size() << '\t' << numSNPs[0] << ':' << numSNPs[1] << '\t' <<  relativeBreakpoint << '\t'
            << numIgnoredSNPs << '\t' << n_alternate_alignments[0] << ':' << n_alternate_alignments[1] << '\t' << SNP_symbols[start_gt] << '\t' << ALIGN_DIRECTIONS[align_dir[0]] << ':' << ALIGN_DIRECTIONS[align_dir[1]] << '\t'
            << align_chrom_name[0] << ':' << align_chrom_name[1] << '\t' << align_len[0] << ':' << align_len[1] << '\t' << align_loc[0] << ':' << align_loc[1] << endl;

        if(align_dir[0] != align_dir[1]) log_file << "[anomolous case] Recombinant read " << current_read_name << " has alignments in different directions" << endl;
        categorizeRead(2 + quality, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced, align_len[1]);
    }

    error_file << "Total Error Rate: " << 100 * ((double) totalErrors / (double) totalErrorBP) << "%" << endl;

    return;
}

void parseCIGAR(stringstream& read, vector<int>& mm_locs, int& len, unsigned long& loc, string& chrom_name, bool& dir) {
    // parse the CIGAR string of the next read in the sam file
    
    // Lsat_Sal__Chrom2_(7016910..7020531)	2064	SERH_U_Chr_2	35850229	1	4=1X82=1X79=1X60=1X111=1X20=1X3=1X97=1X194=1X12=1X28=1X12=1X21=1X15=1X3=1I49=1X10=1X1=1X144=1D76=1X9=1X1=1X5=2565H
    // readname                             flag    chrom_aligned   position    qual    CIGARstring
    
    int num, flag;
    char token;

    read >> flag;
    read >> chrom_name;
    read >> loc;
    read >> num;

    dir = flag & 16; // check if 16 present in flag, indicating aligned reverse direction

    int pos = 0; // current position on reference

    while(read >> num) {
        read >> token;
        if(token == '=' || token == 'D') { // = match // D deletion
            pos += num;
        } else if(token == 'X') {
            for(int i = 0; i < num; i++) {
                mm_locs.push_back(pos);
                pos++;
            }
        } else {
            // ignore
        }
        
        // X mismatch
        // H hard clip (ignore, does not consume reference)
        // S soft clip (ignore, does not consume reference)
        // D deletion in reference (consume the reference)
        // I insertion (ignore, does not consume reference)
    }
    len = pos;
    return;
}

void printProgramInfo() {
    const char* info = "\nrecombigator version 2.0\nPeter Reifenstein June 2024\n            _____        ____\n           / (^) \\______/ (^) \\ \n          /      ________      \\ \n/\\__/\\____|__o___\\______/__o___|____/\\/\\/\\ \n\nA tool to classify Nanopore or HiFi DNA reads as belonging to Genome A, Genome B, or Recombinant\n\ncompile:\ng++ recombigator.cc -std=c++11 -O3 -o recombigator\n\nusage: \n./recombigator alignment_A.sam alignment_B.sam variations.tsv [output_name]\n\n-- alignment_A.sam\nincludes all reads aligned to genome A, with match/mismatch information in the CIGAR string\nfor example, to align using the minimap2 aligner (https://github.com/lh3/minimap2):\n./minimap2 --eqx -a --secondary=no genomeA.fasta reads.fasta\n\n-- variations.tsv\na tab seperated value file with the format specified by SYRI (https://schneebergerlab.github.io/syri/fileformat.html):\n  genomeA_chrom  genomeA_pos  ~  ~  ~  genomeB_chrom  genomeB_pos  ~  ~  region_type  variation_type  ~ \nregion_type beginning with \"SYN\" and variation_type = \"SNP\" indicates a SNP in a syntenic region\ngenomeA_chrom and genomeA_pos give the coordinates of the SNP in genome A\nfields marked with ~ are not used and can be any string\n";
    cout << info << endl;
    exit(1);
}

void printCategorizationResults(double time) {
    stringstream summary;
    long total_reads = 0;
    long total_discarded_reads = 0;

    summary << left;
    for(int i = 0; i < NUM_CATEGORIES; i++) {
        if(i > CHOMP_CATEGORIES_INDEX) {
            if(i == CHOMP_CATEGORIES_INDEX + 1)
                summary << endl;
            total_discarded_reads += num_reads_in_category[i];
            summary << "Total " << setw(16) << (CATEGORY_NAMES[i] + ":") << num_reads_in_category[i] << endl;
        } else {
            summary << "Total " << setw(16) << (CATEGORY_NAMES[i] + ":") << num_reads_in_category[i] << endl;
        }
        total_reads += num_reads_in_category[i];
    }

    double attrition_rate = (double)total_discarded_reads / (double)total_reads;

    summary << endl << "Processed [" << total_reads << "] reads, " << 100 * attrition_rate << "% were thrown out (chomped)" << endl;

    
    summary << "Elapsed time: " << time << " seconds" << endl;

    string summary_str = summary.str();
    cout << summary_str << endl;
    summary_file << summary_str << endl;

    return;
}

void makeDirectories(string base_dir) {
    mkdir((base_dir).c_str(), 0777);

    string alignments_dir = base_dir + "/sorted_alignments";
	mkdir((alignments_dir).c_str(), 0777);

    // open ofstreams for writing SAM files
    output_SAM_files = (ofstream***) malloc(2 * sizeof(ofstream**));

    for(int i = 0; i < 2; i++) {
        string SAM_names_prefix = "toGenome";
        SAM_names_prefix += SNP_symbols[i];
        output_SAM_files[i] = (ofstream**) malloc(NUM_CATEGORIES * sizeof(ofstream*));
        string gt_dir = alignments_dir + "/" + SAM_names_prefix;
        mkdir((gt_dir).c_str(), 0777);
        for(int j = 0; j < NUM_CATEGORIES; j++) {
            output_SAM_files[i][j] = new ofstream(gt_dir + "/" + SAM_names_prefix + "_" + CATEGORY_NAMES[j] + ".SAM");
        }
    }

    // open ofstreams for writing to text files
    string txt_dir = base_dir + "/text_files";
	mkdir((txt_dir).c_str(), 0777);
    for(int i = 0; i < NUM_CATEGORIES; i++) {
        output_TXT_files[i].open(txt_dir + "/" + CATEGORY_NAMES[i] + ".txt");
    }
    output_TXT_files[NUM_CATEGORIES].open(txt_dir + "/" + "AllReads" + ".txt");

    // open ofstream for writing summary
    summary_file.open(base_dir + "/summary.txt");

    // open ofstream for writing log
    log_file.open(base_dir + "/log.txt");

    // open ofstreams for writing recombinant info
    string info_dir = base_dir + "/recombinant_info";
	mkdir((info_dir).c_str(), 0777);
    for(int i = 0; i < NUM_RECOMBINANT_CATEGORIES; i++) {
        recombinantinfo_files[i].open(info_dir + "/" + CATEGORY_NAMES[i + 2] + "Info.tsv");
        recombinantinfo_files[i] << "name alignmentSNPs(A:B)  numSNPs(A:B)    relativeBreakpoint  numIgnoredSNPs  numAlternateAlignments(A:B)   firstSNP    alignDir(A:B)   chromosome(A:B)   len(A:B)    alignloc(A:B)" << endl;
    }

    string misc_dir = base_dir + "/misc";
    mkdir((misc_dir).c_str(), 0777);

    // open ofstreams for writing lengths of reads
    string len_dir = misc_dir + "/read_lengths";
    mkdir((len_dir).c_str(), 0777);
    for(int i = 0; i < NUM_CATEGORIES; i++) {
        readlength_files[i].open(len_dir + "/" + CATEGORY_NAMES[i] + "_lengths" + ".tsv");
    }
    readlength_files[NUM_CATEGORIES].open(len_dir + "/" + "AllReads_lengths" + ".tsv");

    // open ofstream for writing error rate
    error_file.open(misc_dir + "/Error_Rate.txt");

    return;
}

void createSNPfile(string polymorphism_file) {
    ifstream p_file(polymorphism_file);
    ofstream p_SNP_file(polymorphism_file + ".recombigator.SNPs");
    ofstream p_info_file(polymorphism_file + ".recombigator.info");
    ofstream p_VCF_files[2];
    for(int i = 0; i < 2; i++) {
        p_VCF_files[i].open(polymorphism_file + ".recombigator.SNPs" + ".genome" + SNP_symbols[i] + ".VCF");
        p_VCF_files[i] << "##fileformat=VCFv4.2\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" << endl;
    }
    p_SNP_file << "chromosome_num\tgenome_A_location\tgenome_B_location" << endl;

    if(!p_file) { cout << "[ERROR] cannot open " << polymorphism_file << endl; exit(1); }
    string chrBloc_str, unique_id, parent_id, type, prev_parent_id, str;
    string seq[2];
    string chrName[2];
    //int chrAloc, chrBloc, chrAnum, chrBnum, chrAname, chrBname;
    int chrLoc[2];
    int chrNum[2];
    unsigned long totalcount = 0;
    int num_SNPs_on_chrom [num_chromosomes] = {0};
    int num_syntenic_regions_on_chrom [num_chromosomes] = {0};
    while(p_file >> chrName[0]) {
        p_file >> chrLoc[0]; p_file >> chrLoc[0]; p_file >> seq[0]; p_file >> seq[1];
        p_file >> chrName[1]; p_file >> chrBloc_str;  p_file >> chrBloc_str;
        p_file >> unique_id; p_file >> parent_id; p_file >> type; p_file >> str;

        // SALH_C_Chr_1	1	26851	-	-	-	-	-	NOTAL1	-	NOTAL	-
        if(parent_id.substr(0, 3) == "SYN" && type == "SNP") {
            // found SNP in syntenic region
            chrLoc[1] = stoi(chrBloc_str);
            totalcount++;
            chrNum[0] = chromosome_numbers[chrName[0]];
            chrNum[1] = chromosome_numbers[chrName[1]];
            num_SNPs_on_chrom[chrNum[0]] ++;
            if(chrNum[0] != chrNum[1]) {
                log_file << "[createSNPfile] [warning] Differnent chromosomes in 'syntenic region'" << endl;
            }

            // write to SNP file
            p_SNP_file << chrNum[0] + 1 << '\t' << chrLoc[0] << '\t' << chrLoc[1] << endl;

            // write to both VCF files
            for(int i = 0; i < 2; i++) {
                // CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
                p_VCF_files[i] << (chromosome_names[i][chrNum[i]] + "\t" + to_string(chrLoc[i]) + "\t" + unique_id + "\t" + seq[i] + "\t" + seq[i ^ true] + "\t.\tPASS\t." ) << endl;
            }

            if(parent_id != prev_parent_id) {
                num_syntenic_regions_on_chrom[chrNum[0]]++;
            }
            prev_parent_id = parent_id;
            //if(totalcount % 400000 == 0)
            //    cout << "[SNP parse] Found SNP in syntenic region on chromosome " << chrAnum << " A coordinate: " << chrAloc << " B coordinate: " << chrBloc << endl;
        }
    }
    for(int i = 0; i < num_chromosomes; i++) {
        p_info_file << "Chromosome " << i + 1 << " has total " << num_SNPs_on_chrom[i] << " SNPs among " << num_syntenic_regions_on_chrom[i] << " syntenic regions" << endl;
    }
    
}

void readSNPs(string polymorphism_file) {

    // Lsat_PI251246_v14_Chr1  3435181 3435181 C       T       Lsat_PI251246_v14_Chr4  360187218       360187218       SNP13644        INVTR7135       SNP     -
    // ref_chrom               ref_pos ref_pos ref_seq qry_seq qry_chrom               qry_pos         qry_pos         unique_id       parent_id       type    copy-status
    //genomeA_chrom  genomeA_pos  ~  ~  ~  genomeB_chrom  genomeB_pos  ~  ~  region_type  variation_type  ~ \nwhere the region_type beginning with "SYN" indicates a syntenic region and variation_type = "SNP" indicates a SNP\nand genomeA_chrom and genomeA_pos give the coordinates of the SNP in genome A\n(the other fields marked with ~ are not used and can be any string)
    //
    //
    // (the other fields marked with ~ are not used and can be any string)
    // look for parent_id = SNYx (i.e. the variation occurs in the xth syntenic alignment region) and type = SNP   

    cout << "[recombigator] Reading in Polymorphism file" << endl;

    ifstream p_SNP_file_test(polymorphism_file + ".recombigator.SNPs");
    if(!p_SNP_file_test) {
        cout << "[SNP parse] No SNP file created yet from polymorphisms file, creating SNP file and VCF files..." << endl;
        createSNPfile(polymorphism_file);
    }

    cout << "[recombigator] Reading in SNP file" << endl;
    ifstream p_SNP_file(polymorphism_file + ".recombigator.SNPs");
    string header;
    getline(p_SNP_file, header);
    unsigned long chrAloc, chrBloc;
    int chromosome;
    int count = 0;
    while(p_SNP_file >> chromosome >> chrAloc >> chrBloc) {
        SNP_coordinates[0][chromosome - 1][chrAloc] = 1;
        SNP_coordinates[1][chromosome - 1][chrBloc] = 1;
        count ++;
    }
}

void readChromosomes(string alignment_genome_Afile, string alignment_genome_Bfile) {
    sam_files[0].open(alignment_genome_Afile);
    sam_files[1].open(alignment_genome_Bfile);

    if(!sam_files[0] || !sam_files[1]) { cout << "[ERROR] cannot open one or more SAM files" << endl; exit(1); }

    cout << "[recombigator] Reading chromosomes from SAM files" << endl;
    vector<int> genome_sizes[2];
    string chrom_name;

    char c;
    string str;
    int len;
    int n = 0;
    for(int i = 0; i < 2; i++) {
        n = 0;
        while(true) {
            sam_files[i] >> c;
            if(c!= '@') {
                sam_files[i].unget();
                break;
            }
            
            getline(sam_files[i], str);
            str = "@" + str;
            stringstream header_line(str);
            header_line << str;
            for(int k = 0; k < NUM_CATEGORIES; k++) {
                *output_SAM_files[i][k] << str << endl;
            }

            header_line >> str;
            if(str == "@SQ") { // genomic alignment reference sequence (i.e. chromosome)
                header_line >> str;
                str = str.substr(3);
                header_line >> c;header_line >> c;header_line >> c;
                header_line >> len;

                log_file << "[chromosome] " << str << " len: " << len << endl;
                chromosome_numbers[str] = n;
                chromosome_names[i][n] = str;
                genome_sizes[i].push_back(len);
                n++;
            } else { // consume the rest of the field
                getline(header_line, str); 
            }
            
        }
    }
    num_chromosomes = genome_sizes[1].size();

    cout << "[recombigator] Allocating Memory for SNP locations" << endl;

    SNP_coordinates    = (bool***) malloc(2 * sizeof(bool**));
    SNP_coordinates[0] = (bool**)  malloc(num_chromosomes * sizeof(bool*));
    SNP_coordinates[1] = (bool**)  malloc(num_chromosomes * sizeof(bool*));

    for(int i = 0; i < num_chromosomes; i++) {
        // memory allocated with calloc will be zero initialized
        SNP_coordinates[0][i] = (bool*) calloc((genome_sizes[0][i] + 1), sizeof(bool));
        SNP_coordinates[1][i] = (bool*) calloc((genome_sizes[1][i] + 1), sizeof(bool));
        log_file << "[calloc] Allocated memory for chromosome " << i + 1 << " SNPs (total " << genome_sizes[0][i] + 1 << " bytes)" << endl;
    }
    return;
}

void memoryCleanup() {
    log_file << "[recombigator] Deallocating heap memory" << endl;
    for(int i = 0; i < num_chromosomes; i++) {
        free(SNP_coordinates[0][i]);
        free(SNP_coordinates[1][i]);
    }
    free(SNP_coordinates[0]);
    free(SNP_coordinates[1]);
    free(SNP_coordinates);

    for(int i = 0; i < NUM_CATEGORIES; i++) {
        free(output_SAM_files[0][i]);
        free(output_SAM_files[1][i]);
    }
    free(output_SAM_files[0]);
    free(output_SAM_files[1]);
    free(output_SAM_files);

    log_file << "[recombigator] Deallocated heap memory" << endl;
    return;
}

void generatePlots(string base_dir) {
    string command = "python3 generatePlots.py " + base_dir;
    cout << "[recombigator] Generating plots" << endl;
    int result = system(command.c_str());
    if(result == 1) {
        cout << "[recombigator] Error generating plots" << endl;
    } else {
        cout << "[recombigator] Plots generated" << endl;
    }
}
#pragma once

#ifdef COMPILING_GA_DLL
#define CTM_GA_API __declspec(dllexport)
#else
#define CTM_GA_API __declspec(dllimport)
#endif

#include <vector>
#include <string>
#include <unordered_map>

#include "GA_Structs.h"
#include "GA_InternalStructs.h"

namespace CTM_GA {

/**
 1 intersection gene including:
 offset, cycle, phases all in binary and decimal (eg."0100110" or int 45)
*/
struct intersec_gene {
	// only 1 offset among the whole simulation, it ignorant to the typeOfplan, no_Cycle_intotal
	// decimal form
	int dec_Offset;

	// if FGFC,only one dec and binary cycle inside the following vector, eg. dec_Offset.at(0)
	// if VGVC,mutiple dec and binary cycles inside the following vector
	// decimal form
	std::vector<int> dec_Cycle;

	//------------- 1 phase duration gene in 1 intersection (eg."0100110" or int 45) -------------//
   // the inner std::vec => numbers of cycle within that a certain phase	 ==  std::vector<int> dec_phase_duration;	
   // if FGFC,only one dec phase duration inside the following vector, eg. dec_phase_duration.at(0)
   // if VGFC,mutiple dec phase durations inside the following vector
   // decimal form  
   // collectinos of phase_gene, the order in g1, g2 , g3 ....
	std::vector<std::vector<int>> phase_duration_collections;

};

// 1 Chromosome among all population, it includes all intersections inf inside
class CTM_GA_API Chromosome {
private:
	// Gene in string binary form, created in create_wholegene()
	// it may contain offset, cycle and phase 
	std::vector<std::string> wholegene_ref;

	// index location for referring back to the wholegene_ref string operational set
	std::unordered_map <int, int> bin_offset_locations;
	std::unordered_map <int, std::vector<int>> bin_cycle_locations;
	std::unordered_map <int, std::vector<std::vector<int>> > bin_phase_locations;

	// real data struct in decimal integer
	std::unordered_map<int, intersec_gene> dec_data_structure;

	// leave a place for the input check
	const input_param& ref_sourceData;

	// leave a place for the input check
	std::unordered_map<int, intersec_inf>& intersecs_inf;

public:
	// Default constructor
	Chromosome(const input_param& sourceData, std::unordered_map<int, intersec_inf>& intersecs_inf);
	// Default constructor
	Chromosome(const input_param& sourceData, std::unordered_map<int, intersec_inf>& intersecs_inf,
		std::unordered_map<int, intersec_gene> dec_data_structure);
	// destructor
	~Chromosome();

	// get the data inside chromosome object class
	const input_param& get_ref_sourceData() const;

	std::unordered_map<int, intersec_gene> get_dec_data_structure() const;

	void set_dec_data_structure_by_intersec_id(int intersec_id, intersec_gene&& dec_data_structure);


	// --- --- --- --- --- initiate the starting chromosome pop,  from initial setting --- --- --- --- ---
	void create_initial_seed();

	// --- --- --- --- --- generate the a population just on the dec_data_structure for the 1st generation--- --- --- --- ---
	//ref_chromo is used for when there is limit to changable range and is using initial guess/ SGA
	bool gen_random_timeplan();

	// --- --- --- --- --- (re)create whole gene and the bin location --- --- --- --- ---
	void create_wholegene();

	// --- --- --- --- --- according to bin location --- --- --- --- ---
	bool wholegene_operation2dec_data();


	// --- --- --- --- --- GA processing --- --- --- --- ---

	////............................... crossover by rate and mutation by rate ...............................////
	// do crossover by rate and mutation by rate on two selected chromosome in the  jobs_ptr_collections pool
	// it could be repeatedly called if the chromosome is failed in check validity
	void crossoverByrate_mutationByrate(Chromosome& c);

	////............................... crossover by crossover point and mutation by rate ...............................////
	// do crossover by crossover point and mutation by rate on two selected chromosome in the  jobs_ptr_collections pool
	// it could be repeatedly called if the chromosome is failed in check validity
	void crossoverBypt_mutationByrate(Chromosome& c);

	void print_dec_data_structure();

private:

	// --- --- --- --- --- generate set of crossover point locations for one crossover and mutation on a pair of chromosome--- --- --- --- ---
	// usde for the crossoverBypt_mutationByrate() func
	std::set<int> generate_rand_crossover_point();

	char bin_flip(char& ch);  // can be repeatedly call used for mutation()

	int bin2dec(std::string bin);

	std::string dec2bin(int dec);

	bool Check_crossPhaseMinGreen(const int& intersec_id, const intersec_gene& intersec_obj);

};

}	//namespace CTM_GA

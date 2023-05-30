#pragma once

#ifdef COMPILING_GA_DLL
#define DISCO2_GA_API __declspec(dllexport)
#else
#define DISCO2_GA_API __declspec(dllimport)
#endif

#include <set> 
#include <vector>
#include <unordered_set>
#include <unordered_map>

#include "GA_Structs.h"
#include "GA_InternalStructs.h"

namespace DISCO2_GA {

//------------- classes that needs to be pre-defined -------------

class I_SigDataConvertor;
class I_RunCase;
class I_RunCase_Factory;
class I_Runner;
class Chromosome;

class DISCO2_GA_API GA_Module {
private:

	// --- --- --- --- --- RunTime Storage --- --- --- --- --- //

	// holds ownership of Job ptrs for the current generation pre-run
	std::vector<I_RunCase*> jobs_ptr_collections;

	// ptr copy from the prev ptr collection, populated according to their fitness
	std::vector<I_RunCase*> jobs_intermediate_parents_collections;

	// holds ownership of Job ptrs after run(generated eval values) between generations
	std::vector<I_RunCase*> job_previousgen_ptr_collections;

	// fitness values of the Job ptrs, in same sequence of prev ptrs vector
	std::vector<double> pops_fit;

	// history of I_RunCase, inserted when it is deletable
	std::unordered_set<I_RunCase*> history_jobs;

	// --- --- --- --- --- RunTime Settings --- --- --- --- --- //

	input_param m_param;

	I_SigDataConvertor* m_convertor = nullptr;
	I_RunCase_Factory* m_runCaseFactory = nullptr;
	I_Runner* m_runner = nullptr;

public:

	// --- --- --- --- --- Base Information --- --- --- --- --- //

	// gather all interseciton information together into map
	std::unordered_map<int, intersec_inf> intersecs_inf;

	// --- --- --- --- --- Final Results --- --- --- --- --- //

	// holds ownership of the best I_RunCase of each (sub)generation
	std::vector <I_RunCase*> best_jobs;

	// stores the best I_RunCase's fitness of each (sub)generation
	std::vector<double> best_pops_fit;

public:

	// --- --- --- --- --- Constructor Destructor --- --- --- --- --- //

	GA_Module();

	~GA_Module();

	// --- --- --- --- --- input_param Sets --- --- --- --- --- //

	void set_rand_seed(unsigned int seed);
	void set_Elite_rate(double rate);
	void set_Crossover_rate(double rate);
	void set_crossover_multiplier(double rate);
	void set_Mutation_rate(double rate);
	void set_Power_factor(double rate);
	void set_Populcation_size(int integer);
	void set_no_generation(int integer);
	void set_no_sub_generation(int integer);
	void set_sga_use_cyc_subgen(bool useCycSubGen);
	//0 -> nothing; 1 -> minimal; 2 -> includes all best jobs details; 3 -> top & bottom 8
	void set_DebugPrintOutLevel(int level);

	// --- --- --- --- --- I_SigDataConvertor --- --- --- --- --- //

	//transfers ownership to this class, will be deleted on dtor
	void setSigDataConvertor(I_SigDataConvertor* sigDataConvertor);
	//return and dereference the convertor ptr from this (no longer owns the convertor)
	I_SigDataConvertor* popSigDataConvertor();	
	//must be called before the other settings
	void initSigDataConvertor();		

	// --- --- --- --- --- I_RunCase_Factory --- --- --- --- --- //

	//transfers ownership to this class, will be deleted on dtor
	void setRunCaseFactory(I_RunCase_Factory* runCaseFactory);
	//return and dereference the convertor ptr from this (no longer owns the factory)
	I_RunCase_Factory* popRunCaseFactory();

	// --- --- --- --- --- I_Runner --- --- --- --- --- //

	//transfers ownership to this class, will be deleted on dtor
	void setRunner(I_Runner* runCaseFactory);
	//return and dereference the convertor ptr from this (no longer owns the factory)
	I_Runner* popRunner();

	// --- --- --- --- --- Sets Value into the intersecs_inf Table --- --- --- --- --- //
	// overwrites(non-repairable) intersecs_inf, re-init sigConvertor to restore data

	void set_use_initial_seed(bool use);
	void set_initial_phase(int intersec_id, int phase_id, int initial_phase);
	void set_initial_Cycle_time(int intersec_id, int initial_Cycle_time);
	void set_initial_offset(int intersec_id, int offset);
	
	void set_Min_green(int intersec_id, int phase_id, int Min_green);
	void set_inter_green(int intersec_id, int phase_id, int inter_green);
	// only limits the lower bound of the phase green time, iff initial guess is used
	void set_phase_changeable_range(int intersec_id, int phase_id, int range);	

	void set_Max_Cycle_length(int intersec_id, int Max_Cycle);
	void set_Min_Cycle(int intersec_id, int Min_Cycle);						
	// only limits the lower bound of the cycle time, iff initial guess is used
	void set_cycle_time_changeable_range(int intersec_id, int range);	

	// unused
	void set_Stop_criteria(int integer);

	// --- --- --- --- --- How Values are Changeable --- --- --- --- --- //

	// Permanently changes the chamgeable setting for offset
	void set_Offset_Changeable(int intersec_id, bool newChangeable);
	// Permanently changes the chamgeable setting for offset
	void set_CycleTime_Changeable(int intersec_id, bool newChangeable);
	// Permanently changes the chamgeable setting for offset
	void set_PhaseG_Changeable_by_value(int intersec_id, int phase_id, bool newChangeable);
	void set_PhaseG_Changeable_by_ratio(int intersec_id, int phase_id, bool newChangeable);

	// --- --- --- --- --- RunTime Settings --- --- --- --- --- //

	// limits the number of tries to generate a chromosome from random gen or crossover&mutation
	void set_trial_count_limit(int trial_count_limit);

	// only applicable for SGA, defines the grouping of which the junctions will be 
	// allowed to be optimised at a time
	void set_unfreeze_list(std::vector < std::set <int>> unfreeze_list);

	// defines the groups which will share the same cycle time
	void set_intersections_gp(std::vector<std::vector<int>> intersec_ids_list);

	// sets the ordering of the intersections in a chromosome
	void set_intersecs_sequence(std::vector<int> id_list);

	// injects a seed chromosome, can be repeatedly called
	void set_seedjob_by_chromosome(Chromosome&& chromosome_seed);

	// --- --- --- --- --- Execution of GA Optimization --- --- --- --- ---

	void Run_NGA();
	void Run_SGA();

	// --- --- --- --- --- Analyse Results --- --- --- --- ---

	std::vector<double> get_best_pops_avd();
	std::vector<double> get_best_pops_fit();

private:

	// --- --- --- --- --- Run GA --- --- --- --- ---
	void initRuns();
	void set_no_bits_in_chromosome();
	void GA_input_check_before_run();
	void create_initial_seedjob();

	void fill_random_pop();
	void crossover_mutation();

	//UNUSED, for when VGVC, VGFC exist
	//to calculate the number required for the number of cycle time conservatively 
	void cal_max_no_Cycle();

	//copies the top fitness prev-object ptrs into the intermediate pop
	void cal_jobs_intermediate_parents_and_jobs_ptr();

	// move all the ptrs over from jobs_ptr_collections and rebuild the vector
	// holds pointers for copy in the next gen elites
	void save_job_previousgen_ptr_collections();

	void inheritance_of_elites();
	void sort_jobs_by_avd();
	void fitness_cal();
	//selection from the intermediate parents pool
	std::pair<int, int> generate_Rand_crossover_selection(std::set<int>& usedGenes);

	//adds job ptr to history and awaits delete
	//function performs checks to ensure input is not resevered elsewhwere
	void add_to_history(I_RunCase* jobPtr);
	// deletes all ptr contained within
	void clean_hisotry();

	// --- --- --- --- --- functions to deal with the Job Gene --- --- --- --- ---

	static bool sort_by_avd(I_RunCase* a, I_RunCase* b);
	static double accumulate_avd(double result, I_RunCase* obj);
	static bool search_max_avd(I_RunCase* a, I_RunCase* b);
	static bool search_min_avd(I_RunCase* a, I_RunCase* b);

	// --- --- --- --- --- add results --- --- --- --- ---

	void save_best_results();

	// --- --- --- --- --- Printouts --- --- --- --- ---

	//print out debugging info
	void debugPrintOut();
	//print out the best results in when GA finished
	void Print_best_gene();
	//print out the frist numOfjobs results in each generation
	void Print_best_genes_in_each_gen(int numOfjobs);
	//print out the end numOfjobs results in each generation
	void Print_worst_genes_in_each_gen(int numOfjobs);

	// --- --- --- --- --- Export Results --- --- --- --- ---

	//puts the best chromo ptr into convertor for export
	void exportBestChromosomeToConvertor();

};

}	//namespace DISCO2_GA

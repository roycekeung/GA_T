#include "GA_Module.h"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <thread>
#include <limits>
#include <numeric>
#include <cassert>

//internal impls
#include "Chromosome.h"
#include "I_SigDataConvertor.h"
#include "I_RunCase.h"
#include "I_Runner.h"

namespace DISCO2_GA {

// Util to measure performance
class Timer {
private:
	std::chrono::time_point<std::chrono::steady_clock > starttime;

public:

	void start() {
		starttime = std::chrono::steady_clock::now();
	}

	void Stop() {
		auto endtime = std::chrono::steady_clock::now();
		auto start = std::chrono::time_point_cast<std::chrono::milliseconds>(starttime).time_since_epoch().count();
		auto end = std::chrono::time_point_cast<std::chrono::milliseconds>(endtime).time_since_epoch().count();

		auto duration = (end - start);
		long long sec = (duration / 1000);

		std::cout << "\n" << sec << " seconds run time " << std::endl;
	}
};

// --- --- --- --- --- Constructor Destructor --- --- --- --- --- //

GA_Module::GA_Module() {
		srand(time(NULL));
}

GA_Module::~GA_Module() {
	// Created I_RunCase
	for (auto ptr : this->best_jobs)
		this->history_jobs.insert(ptr);
	this->clean_hisotry();
	// Other owned stuff
	if (this->m_convertor)
		delete this->m_convertor;
	if (this->m_runCaseFactory)
		delete this->m_runCaseFactory;
	if (this->m_runner)
		delete this->m_runner;
};

// --- --- --- --- --- input_param Sets --- --- --- --- --- //

void GA_Module::set_rand_seed(unsigned int seed) {
	std::srand(seed);
}

void GA_Module::set_Elite_rate(double rate) {
	if (rate <= 0) {
		throw std::runtime_error("this cant be 0 or smaller than 0");
	}
	else if (rate >= 1) {
		throw std::runtime_error("this cant be greater than 1");

	}
	// fixing the offspring_size and Elite_size even and odd number problems
	else if (std::round((int)std::round(this->m_param.Populcation_size * rate) % 2) == 0) {
		this->m_param.Elite_rate = rate;

		this->m_param.Elite_size = (int)std::round((this->m_param.Populcation_size * this->m_param.Elite_rate)) ;
		this->m_param.offspring_size = this->m_param.Populcation_size - this->m_param.Elite_size;
	}
	// fixing the offspring_size and Elite_size even and odd number problems
	else if (std::round((int)std::round(this->m_param.Populcation_size * rate) % 2) != 0) {
		this->m_param.Elite_rate = rate;

		this->m_param.Elite_size = (int)std::round((this->m_param.Populcation_size * this->m_param.Elite_rate)) + 1;
		this->m_param.offspring_size = this->m_param.Populcation_size - this->m_param.Elite_size;
	}
}

void GA_Module::set_Crossover_rate(double rate) {
	if (rate <= 0) {
		throw std::runtime_error("this cant be 0 or smaller than 0");
	}
	else if (rate >= 1) {
		throw std::runtime_error("this cant be greater than 1");
	}
	// fixing the offspring_size and Elite_size even and odd number problems
	else if (std::round((int)std::round(this->m_param.Populcation_size * 0.08) % 2) == 0) {
		this->m_param.Crossover_rate = rate;
	}
	// fixing the offspring_size and Elite_size even and odd number problems
	else if (std::round((int)std::round(this->m_param.Populcation_size * 0.08) % 2) != 0) {
		this->m_param.Crossover_rate = rate;
	}
}

void GA_Module::set_crossover_multiplier(double rate) {
	if (rate < this->m_param.bits_in_chromosome)
		this->m_param.crossover_multiplier = rate;
	else
		throw std::runtime_error("crossover multiplier is too large");
}

void GA_Module::set_Mutation_rate(double rate) {
	if (rate <= 0)
		throw std::runtime_error("this cant be 0 or smaller than 0");
	else if (rate >= 1)
		throw std::runtime_error("this cant be greater than 1");
	else
		this->m_param.Mutation_rate = rate;
}

void GA_Module::set_Power_factor(double rate) {
	if (rate <= 0)
		throw std::runtime_error("this cant be 0 or smaller than 0");
	else if (rate > 10)
		throw std::runtime_error("this cant be greater than 10");
	else
		this->m_param.Power_factor = rate;
}

void GA_Module::set_Populcation_size(int integer) {
	if (integer <= 0) {
		throw std::runtime_error("this cant be 0 or smaller than 0");
	}
	else if (integer % 2 == 0) {
		this->m_param.Populcation_size = integer;
		this->m_param.Elite_size = (int)std::round(this->m_param.Populcation_size * this->m_param.Elite_rate);   // the size of the elite to be inherit for the next generation
		//Elite_size = parent_size;
		this->m_param.offspring_size = this->m_param.Populcation_size - this->m_param.Elite_size;

	}
	else if (integer % 2 != 0) {
		throw std::runtime_error("Populcation_size cant be odd number");
	}
}

void GA_Module::set_no_generation(int integer) {
	if (integer <= 0)
		throw std::runtime_error("this cant be 0 or smaller than 0");
	else
		this->m_param.no_generation = integer;
}

void GA_Module::set_no_sub_generation(int integer) {
	if (integer < 0)
		throw std::runtime_error("this cant be 0 or smaller than 0");
	else
		this->m_param.no_sub_generation = integer;
}

void GA_Module::set_sga_use_cyc_subgen(bool useCycSubGen) {
	this->m_param.sga_cyc_subgen = useCycSubGen;
}

void GA_Module::set_DebugPrintOutLevel(int level) {
	level = std::min(0, level);
	this->m_param.debugPrintOutLevel = level;
}

// --- --- --- --- --- I_SigDataConvertor --- --- --- --- ---

void GA_Module::setSigDataConvertor(I_SigDataConvertor* sigDataConvertor) {
	//delete the exsistig one first
	if (this->m_convertor != nullptr) {
		delete this->m_convertor;
	}
	//move the new one in
	this->m_convertor = sigDataConvertor;
}

I_SigDataConvertor* GA_Module::popSigDataConvertor() {
	I_SigDataConvertor* temp = this->m_convertor;
	this->m_convertor = nullptr;
	return temp;
}

void GA_Module::initSigDataConvertor() {
	//use the default if a convertor is not supplied
	if (this->m_convertor == nullptr)
		throw std::runtime_error("I_SigDataConvertor not set");
	//invoke the init function
	this->m_convertor->init(this->intersecs_inf);
}

// --- --- --- --- --- I_RunCase_Factory --- --- --- --- --- //

void GA_Module::setRunCaseFactory(I_RunCase_Factory* runCaseFactory) {
	//delete the exsistig one first
	if (this->m_runCaseFactory != nullptr) {
		delete this->m_runCaseFactory;
	}
	//move the new one in
	this->m_runCaseFactory = runCaseFactory;
}

I_RunCase_Factory* GA_Module::popRunCaseFactory() {
	I_RunCase_Factory* temp = this->m_runCaseFactory;
	this->m_runCaseFactory = nullptr;
	return temp;
}

// --- --- --- --- --- I_Runner --- --- --- --- --- //

void GA_Module::setRunner(I_Runner* runCaseFactory) {
	//delete the exsistig one first
	if (this->m_runner != nullptr) {
		delete this->m_runner;
	}
	//move the new one in
	this->m_runner = runCaseFactory;
}

I_Runner* GA_Module::popRunner() {
	I_Runner* temp = this->m_runner;
	this->m_runner = nullptr;
	return temp;
}

// --- --- --- --- --- set the new value into the inf table --- --- --- --- ---

void GA_Module::set_Min_green(int intersec_id, int phase_id, int Min_green) {
	if (this->intersecs_inf.empty()) {
		throw std::runtime_error("plz initialize the intersection inf table first");
	}

	if (Min_green <= 0) {
		throw std::runtime_error("this cant be 0 or smaller than 0");
	}
	else {
		this->intersecs_inf.at(intersec_id).phase_duration_inf_collections[phase_id].min_green = Min_green;
	}
}

void GA_Module::set_inter_green(int intersec_id, int phase_id, int inter_green) {
	if (this->intersecs_inf.empty()) {
		throw std::runtime_error("plz initialize the intersection inf table first");
	}

	if (inter_green <= 0) {
		throw std::runtime_error("this cant be 0 or smaller than 0");
	}
	else {
		this->intersecs_inf.at(intersec_id).phase_duration_inf_collections[phase_id].intergreen = inter_green;
	}
}

void GA_Module::set_phase_changeable_range(int intersec_id, int phase_id, int range) {
	if (this->intersecs_inf.empty()) {
		throw std::runtime_error("plz initialize the intersection inf table first");
	}
	else if (this->intersecs_inf.count(intersec_id) && this->intersecs_inf.at(intersec_id).phase_duration_inf_collections.size() > phase_id)
		this->intersecs_inf.at(intersec_id).phase_duration_inf_collections[phase_id].changableInterval = range;
	else
		throw std::invalid_argument("no intersec id or phase found");
}

void GA_Module::set_Max_Cycle_length(int intersec_id, int Max_Cycle) {
	if (this->intersecs_inf.empty()) {
		throw std::runtime_error("plz initialize the intersection inf table first");
	}

	if (Max_Cycle <= 0) {
		throw std::runtime_error("this cant be 0 or smaller than 0");
	}
	else {
		this->intersecs_inf.at(intersec_id).max_Cycle= Max_Cycle;
	}
}

void GA_Module::set_Min_Cycle(int intersec_id, int Min_Cycle) {
	if (this->intersecs_inf.empty()) {
		throw std::runtime_error("plz initialize the intersection inf table first");
	}

	if (Min_Cycle <= 0) {
		throw std::runtime_error("this cant be 0 or smaller than 0");
	}
	else {
		this->intersecs_inf.at(intersec_id).min_Cycle = Min_Cycle;
	}
}

void GA_Module::set_initial_offset(int intersec_id, int offset) {
	if (this->intersecs_inf.empty()) {
		throw std::runtime_error("plz initialize the intersection inf table first");
	}

	if (offset <= 0) {
		throw std::runtime_error("this cant be 0 or smaller than 0");
	}
	else {
		this->intersecs_inf.at(intersec_id).initial_Offset = offset;
	}
}

void GA_Module::set_initial_Cycle_time(int intersec_id, int initial_Cycle_time) {
	if (this->intersecs_inf.empty()) {
		throw std::runtime_error("plz initialize the intersection inf table first");
	}

	if (initial_Cycle_time <= this->intersecs_inf.at(intersec_id).sum_PhaseInterG ) {
		throw std::runtime_error("this cycle time is too small");
	}
	else {

		//// ``````````````````````````````` requires to rebuilt the initial phases````````````````````````````````////
		

		//// ```````````````````````````````phases````````````````````````````````////   only have 1 no_cycle , dec_phase_duration size = 1
		//  sum of ( phase  - its min green)
		int sum_PhaseMinusMinG = 0;
		// sum the min green of all phases
		int sum_MinG = this->intersecs_inf.at(intersec_id).sum_PhaseMinG;

		// sum the inter green of all phases
		int sum_InterG = this->intersecs_inf.at(intersec_id).sum_PhaseInterG;


		//int int_phase;
		int tmp_phase = 0;
		std::vector<int> tmp_phase_per_intersec = {};


		// refresh phase_gen efor each no_phase
		std::vector <int> phase_gene1 = {};
		std::vector<std::vector <int>> phase_duration_collections1 = {};


		for (int no_phase = 0; no_phase < this->intersecs_inf.at(intersec_id).phase_duration_inf_collections.size(); ++no_phase) {
			// create that certain amouont of empty phase_gene1 first inside phase_duration_collections1
			phase_duration_collections1.emplace_back(phase_gene1);

			// all phase green in this intersection is fixed by ratio


							// the previous newly randomly generated cycle time from the above session
			tmp_phase = 
				this->intersecs_inf.at(intersec_id).phase_duration_inf_collections.at(no_phase).initial_phase_duration
				- this->intersecs_inf.at(intersec_id).phase_duration_inf_collections.at(no_phase).min_green;

			tmp_phase_per_intersec.emplace_back(tmp_phase);


			// sum up all unchanged phase green, for following scaled equation
			sum_PhaseMinusMinG += tmp_phase;			
		}

		// cal the range of the phase durations is allowed for variation
									// newly assign cycle time   -   all min green for all phases(sum_PhaseMinG)      -    all intergreen for all phases (sum_PhaseInterG)
		int phase_variation_range = initial_Cycle_time - this->intersecs_inf.at(intersec_id).sum_PhaseMinG - this->intersecs_inf.at(intersec_id).sum_PhaseInterG;

		// newly assign cycle time -   all intergreen for all phases (sum_PhaseInterG) 
		int unused_range = initial_Cycle_time - this->intersecs_inf.at(intersec_id).sum_PhaseInterG ;


		int new_phase;

		// To scale phases green based on ratio
		for (int no_phase = 0; no_phase < this->intersecs_inf.at(intersec_id).phase_duration_inf_collections.size(); ++no_phase) {

			// all phase green in this intersection is fixed by ratio
			
			new_phase= static_cast<int>(std::round(phase_variation_range * (tmp_phase_per_intersec.at(no_phase) / (double)sum_PhaseMinusMinG))
				+ this->intersecs_inf.at(intersec_id).phase_duration_inf_collections.at(no_phase).min_green);

			//--------- intialize the dec_phase_duration into phase_duration_collections1 
			phase_duration_collections1.at(no_phase).emplace_back(new_phase);

			// update the room left for the unused range
			unused_range -= new_phase;

			// to ensure the sum of all phases (including their mingreen) and intergreen = new cycle time
			if (no_phase == this->intersecs_inf.at(intersec_id).phase_duration_inf_collections.size() - 1 && unused_range != 0) {
				phase_duration_collections1.at(no_phase)[0] += unused_range;
			}
		}
		// clear for the next round diff intersection
		tmp_phase_per_intersec.clear();


		intersec_gene intersec_gene1{};
		intersec_gene1.dec_Cycle.emplace_back(initial_Cycle_time);
		intersec_gene1.dec_Offset = this->intersecs_inf.at(intersec_id).initial_Offset;
		//--------- intialize and add dec_oc_phases_per_intersec to dec_intersections_per_pop
		intersec_gene1.phase_duration_collections = phase_duration_collections1;
		



		// check the crossphase here, 
		// it would only check the unfreeze intersection phases, so its intersection by intersection
		
		// firstly assume all phases in this intersection do reach the crossphase requirement whether it has the crossphase or not
		bool bool_Check_gene_requirements = true;

		// --- --- --- --- --- --- --- --- check of  PhaseidList_CrossPhasesMinGreen std::pair first : list of phases index, second : crossphase min green  --- --- --- --- --- --- --- ---  //
		if (this->intersecs_inf.at(intersec_id).PhaseidList_CrossPhasesMinGreen.empty()) {
			// this intersection have no crossphase with empty vector here
			// do nothing

		}
		else {
			// this phase duration in this intersection have crossphase

			// define the temp value for following crosspahse value checking 
			int current_check_phase = 0, previous_phase_durations, previous_intergreens;
			int crossphase_druation;

			for (std::size_t no_cycle = 0; no_cycle < intersec_gene1.dec_Cycle.size(); ++no_cycle) {
				previous_phase_durations = 0;
				previous_intergreens = 0;

				//  std::pair first : list of phases index, second : crossphase min green  =====> this is for the crossphase only, the vector size depends on how many crossphase u have
				for (auto Crossphase1 : this->intersecs_inf.at(intersec_id).PhaseidList_CrossPhasesMinGreen) {

					// list of crossphases index (eg. <2,3,4> or <4,0>), looping each crossphase index
					for (auto crossphase_index_it = Crossphase1.first.begin(); crossphase_index_it != Crossphase1.first.end(); ++crossphase_index_it) {

						if (*crossphase_index_it == Crossphase1.first.back()) {
							// crossphase_index loop to the last element in which the real crossphase value checking starts
							// if loop to the last element phase of the crossphase, dont need to add up,saved that of phase duration into current_check_phase
							current_check_phase = intersec_gene1.phase_duration_collections[*crossphase_index_it][no_cycle];

						}
						else {
							//list of crossphases index (eg. <2,3,4> or <4,0>) the crossphase_index was still havent reached the last element

							//  add up previous phase of the crossphase
							previous_phase_durations += intersec_gene1.phase_duration_collections[*crossphase_index_it][no_cycle];
							//  add up previous inter green of the crossphase
							previous_intergreens += this->intersecs_inf.at(intersec_id).phase_duration_inf_collections[*crossphase_index_it].intergreen;

						}

					}

					crossphase_druation = current_check_phase + previous_phase_durations + previous_intergreens;

					if (crossphase_druation < Crossphase1.second) {
						// true, pass to here means crossPhaseMinGreen is failed for this no_cycle or for this intersection; failed
						bool_Check_gene_requirements = false;
						break;

					}
				}
			}
		}


		if (bool_Check_gene_requirements == false) {
			// failed in Check_gene_requirements, so abandon this generation chromosome

			// coz it would stop gen_pop func immediately, so delete the tmp intersec_gene1
			// stack object no need to delete

			// this equation would jump into the end, coz bool_Check_gene_requirements r false (fail in generating this chromosome)
			throw std::runtime_error("when phases split, this cycle time is too small, not enough room to distribute phases");
		}
		else {
			// passed

			this->intersecs_inf.at(intersec_id).initial_Cycle = initial_Cycle_time;

			for (int no_phase = 0; no_phase < this->intersecs_inf.at(intersec_id).phase_duration_inf_collections.size(); ++no_phase) {
				this->intersecs_inf.at(intersec_id).phase_duration_inf_collections.at(no_phase).initial_phase_duration = intersec_gene1.phase_duration_collections.at(no_phase)[0];
			}
		}


	}
}

void GA_Module::set_cycle_time_changeable_range(int intersec_id, int range) {
	if (this->intersecs_inf.empty()) {
		throw std::runtime_error("plz initialize the intersection inf table first");
	}
	else if (this->intersecs_inf.count(intersec_id))
		this->intersecs_inf.at(intersec_id).changableInterval = range;
	else
		throw std::invalid_argument("no intersec id found");
}

void GA_Module::set_Stop_criteria(int integer) {
	if (this->intersecs_inf.empty()) {
		throw std::runtime_error("plz initialize the intersection inf table first");
	}

	if (integer <= 0) {
		throw std::runtime_error("this cant be 0 or smaller than 0");
	}
	else if (this->m_param.no_generation > integer) {
		throw std::runtime_error("Stop_criteria cant be greater than the total generation");
	}
	else {
		this->m_param.Stop_criteria = integer;
	}
}

void GA_Module::set_use_initial_seed(bool use) {
	this->m_param.use_initial_seed = use;
}


void GA_Module::set_initial_phase(int intersec_id, int phase_id, int initial_phase) {

	if (this->intersecs_inf.count(intersec_id)) {

		this->intersecs_inf.at(intersec_id).phase_duration_inf_collections.at(phase_id).initial_phase_duration = initial_phase;
	}
	else {

		throw std::runtime_error("this intersection doesnt exist");
	}
}

// --- --- --- --- --- boolean of fixing the certain value Sets --- --- --- --- --- //

void GA_Module::set_Offset_Changeable(int intersec_id, bool newChangeable) {
	if (this->intersecs_inf.empty()) {
		throw std::runtime_error("plz initialize the intersection inf table first");
	}

	if (intersec_id < 0) {
		throw std::runtime_error("this cant be smaller than 0");
	}
	else {
		this->intersecs_inf.at(intersec_id).Offset_changeable = newChangeable;
	}
}

void GA_Module::set_CycleTime_Changeable(int intersec_id, bool newChangeable) {
	if (this->intersecs_inf.empty()) {
		throw std::runtime_error("plz initialize the intersection inf table first");
	}

	if (intersec_id < 0) {
		throw std::runtime_error("this cant be smaller than 0");
	}
	else {
		this->intersecs_inf.at(intersec_id).Cycle_changeable = newChangeable;
	}
}

void GA_Module::set_PhaseG_Changeable_by_value(int intersec_id, int phase_id, bool newChangeable) {
	if (this->intersecs_inf.empty()) {
		throw std::runtime_error("plz initialize the intersection inf table first");
	}

	if (intersec_id < 0) {
		throw std::runtime_error("this cant be smaller than 0");
	}
	else {
		this->intersecs_inf.at(intersec_id).phase_duration_inf_collections[phase_id].phase_changeable_by_value = newChangeable;
	}
}

void GA_Module::set_PhaseG_Changeable_by_ratio(int intersec_id, int phase_id, bool newChangeable) {
	if (this->intersecs_inf.empty()) {
		throw std::runtime_error("plz initialize the intersection inf table first");
	}

	if (intersec_id < 0) {
		throw std::runtime_error("this cant be smaller than 0");
	}
	else {
		this->intersecs_inf.at(intersec_id).phase_duration_inf_collections[phase_id].phase_changeable_by_ratio = newChangeable;
	}
}

void GA_Module::set_trial_count_limit(int trial_count_limit) {
	this->m_param.trial_count_limit = trial_count_limit;
}

void GA_Module::set_unfreeze_list(std::vector < std::set <int>> unfreeze_list) {
	this->m_param.unfreeze_list.clear();
	this->m_param.unfreeze_list = unfreeze_list;
}

void GA_Module::set_intersections_gp(std::vector<std::vector<int>> intersec_ids_list) {
	std::unordered_set<int> added{};
	this->m_param.intersections_gp.clear();

	for (auto intersec_ids : intersec_ids_list) {
		int coordinated_intersec_id = intersec_ids.front();

		for (int intersec_id : intersec_ids) {

			//check for overlapping
			if (added.count(intersec_id)) {
				this->m_param.intersections_gp.clear();
				throw std::invalid_argument("overlapping groups, all inputs cleared");
			}
			added.insert(intersec_id);

			intersecs_gps_inf intersecs_gp{};
			intersecs_gp.coordinated_intersec_id = coordinated_intersec_id;
			intersecs_gp.relation_group_id = intersec_ids;
			intersecs_gp.minimum_maxCycle = 0;
			intersecs_gp.maximum_minCycle = 300;
			this->m_param.intersections_gp.emplace(intersec_id, std::move(intersecs_gp));
		}
	}
}

void GA_Module::cal_max_no_Cycle() {
	// tmp  cycle_min on all min_Cycle in intersections
	int tmp_cycle_min =0;
	for (auto it = this->intersecs_inf.begin(); it != this->intersecs_inf.end(); ++it)
		if (tmp_cycle_min > it->second.min_Cycle) {
			tmp_cycle_min = it->second.min_Cycle;
		}
		else {
			// nth
		}
}

void GA_Module::set_no_bits_in_chromosome() {
	// tmp  cycle_max on all max_Cycle in intersections
	int tmp_cycle_max= 0;
	for (auto it = this->intersecs_inf.begin(); it != this->intersecs_inf.end(); ++it)
		if (tmp_cycle_max < it->second.max_Cycle) {
			tmp_cycle_max = it->second.max_Cycle;
		}
		else {
			// nth
		}

	if (tmp_cycle_max <= 63) {
		this->m_param.bits_in_chromosome = 6;
		this->m_param.bits_capability = 63;
	}
	else if (tmp_cycle_max > 63 && tmp_cycle_max <= 127) {
		this->m_param.bits_in_chromosome = 7;
		this->m_param.bits_capability = 127;
	}
	else if (tmp_cycle_max > 127 && tmp_cycle_max <= 255) {
		this->m_param.bits_in_chromosome = 8;
		this->m_param.bits_capability = 255;
	}
	else if (tmp_cycle_max > 255 && tmp_cycle_max <= 511) {
		this->m_param.bits_in_chromosome = 9;
		this->m_param.bits_capability = 511;
	}
}

void GA_Module::initRuns() {
	// GA input param check before run, to avoid error
	GA_input_check_before_run();

	// clear up all history that might have left from the last time GA run, to avoid memory leaks
	//direct insert to bypass check against best_jobs
	for (auto& jobPtr : this->best_jobs)
		this->history_jobs.insert(jobPtr);
	clean_hisotry();

	// clear up all history that might have left from the last time GA run, to avoid memory leaks
	set_no_bits_in_chromosome();

	//setup intersec gps
	for (auto& intersecGpData : this->m_param.intersections_gp) {
		int min_maxCycle = std::numeric_limits<int>::max();
		int max_minCycle = 0;

		for (auto& sigId : intersecGpData.second.relation_group_id) {
			// get the intersec data
			auto& intersecInfo = this->intersecs_inf.at(sigId);

			if (min_maxCycle > intersecInfo.max_Cycle)
				min_maxCycle = intersecInfo.max_Cycle;
			if (max_minCycle < intersecInfo.min_Cycle)
				max_minCycle = intersecInfo.min_Cycle;
		}

		intersecGpData.second.minimum_maxCycle = min_maxCycle;
		intersecGpData.second.maximum_minCycle = max_minCycle;
	}

	//counters
	this->m_param.current_generation = 0;
	this->m_param.current_sub_generation = 0;
	this->m_param.current_unfreeze_loc = 0;
	
	//reset all the lists
	this->m_convertor->setTChromosome(nullptr);
	this->m_param.best_ref_opt_job = nullptr;
	this->best_jobs.clear();
	this->best_pops_fit.clear();
	this->jobs_intermediate_parents_collections.clear();
	this->job_previousgen_ptr_collections.clear();
	//not clearing jobs_ptr_collections to allow injection of seeds

	//gen the intersec seq
	this->m_param.intersecs_sequence.clear();
	std::unordered_set<int> added{};
	//put in the intersec group ref jt first as their genes need to be generated first
	for (auto& set : this->m_param.intersections_gp) {
		int refJctId = set.second.coordinated_intersec_id;
		if (!added.count(refJctId)) {
			this->m_param.intersecs_sequence.push_back(refJctId);
			added.insert(refJctId);
		}
	}
	//dump in the rest
	for (auto& inf : this->intersecs_inf) {
		if (!added.count(inf.first)) {
			this->m_param.intersecs_sequence.push_back(inf.first);
			added.insert(inf.first);
		}
	}
}

//-- - -- - -- - -- - -- - GA input param checking before run -- - -- - -- - -- - -- - //
void GA_Module::GA_input_check_before_run() {
	switch (this->m_param.GA_type) {

	case input_param::GAtype::NGA:
	{
		if (!this->m_param.unfreeze_list.empty()) {
			std::cout << "unfreeze list is being cleared due to running NGA" << std::endl;
			this->m_param.unfreeze_list.clear();
		}
		break;
	}
	
	case input_param::GAtype::SGA:

		// check of unfreeze_list, has to be sth in it before SGA run
		if (this->m_param.unfreeze_list.empty()) {
			throw std::runtime_error("plz use set_unfreeze_list() to input unfreeze_list first in order to activate SGA\
				or restore_intersecs_sequence to get the random unorder_map sequence saved in the intersecs_inf");
		}

		// check of no_sub_generation, cant be less than or equal to 0
		if (this->m_param.no_sub_generation <= 0) {
			throw std::runtime_error("plz input a valid number for the sub generation");
		}
		break;
		
	}
}

void GA_Module::set_intersecs_sequence(std::vector<int> id_list) {
	this->m_param.intersecs_sequence.clear();
	this->m_param.intersecs_sequence = id_list ;
}

void GA_Module::set_seedjob_by_chromosome(Chromosome&& chromosome_seed) {
	if (this->m_runCaseFactory && this->m_convertor) {
		I_RunCase* seed = this->m_runCaseFactory->genRunCase(std::move(chromosome_seed), this->m_convertor);
		this->jobs_ptr_collections.emplace_back(seed);
	}
	else
		throw std::runtime_error("I_RunCase_Factory and/or I_SigDataConvertor not set!");
}

void GA_Module::create_initial_seedjob() {
	Chromosome chromo(this->m_param, this->intersecs_inf);
	// invoke the chromosome built in func for making the initial set
	chromo.create_initial_seed();
	// use create_wholegene to construct its whole gene string
	chromo.create_wholegene();

	I_RunCase* seed = this->m_runCaseFactory->genRunCase(std::move(chromo), this->m_convertor);
	this->jobs_ptr_collections.emplace_back(seed);
}

void GA_Module::fill_random_pop() {
	// process of the generation of random timing plan
	int tryCount = 0;
	while ((int)this->jobs_ptr_collections.size() < this->m_param.Populcation_size && tryCount < this->m_param.trial_count_limit) {
		Chromosome tmp_chromo(this->m_param, this->intersecs_inf);

		bool chromo_validation_check = tmp_chromo.gen_random_timeplan();
		if (chromo_validation_check == true) {
			// true, the randomly generated timing plan is valid after all requirement checking, so add that chromosome into newly job object
			I_RunCase* runCase = this->m_runCaseFactory->genRunCase(std::move(tmp_chromo), this->m_convertor);

			// if chromosome is validated then added into current generation jobs_ptr_collections for operation
			this->jobs_ptr_collections.emplace_back(runCase);
			// if chromosome is validated then use create_wholegene to construct its whole gene string
			runCase->get_Chromosome_obj().create_wholegene();

			//reset try count
			tryCount = 0;
		}
		else {
			// false, the randomly generated timing plan is invalid after all requirement checking, 
			// so do the while loop for gen_random_timeplan() until chromosome reach the requirements
			++tryCount;
		}

		//avoid infinite looping
		if (tryCount > this->m_param.trial_count_limit)
			break;
	}
}

void GA_Module::crossover_mutation() {
	// for the check of overlapping pops selected in selection_methods { wheeel_chart = 0, avoid_reuse };
	std::set<int> usedGenes{};

	std::pair<int, int> pair_chromo_location;

	bool parent_validation1, parent_validation2;
	bool deReselect = true;
	int trial_count = 0;
	int reselectTryCount = 0;
	// prepare for next round, undergoing crossover and mutaiton
	while (usedGenes.size() < this->m_param.Populcation_size		//normal number of times
		&& usedGenes.size() < this->jobs_intermediate_parents_collections.size()	//if the pop is small
		&& usedGenes.size() != this->jobs_intermediate_parents_collections.size() - 1 //if the pop is odd number
		&& reselectTryCount < this->m_param.trial_count_limit) {			//avoid infinite looping

		if (deReselect == true)
			pair_chromo_location = generate_Rand_crossover_selection(usedGenes);

		//// --- --- --- --- --- generate set of crossover point locations for one crossover and mutation on a pair of chromosome--- --- --- --- 

		Chromosome chromo1 = this->jobs_intermediate_parents_collections.at(pair_chromo_location.first)->get_Chromosome_obj();
		Chromosome chromo2 = this->jobs_intermediate_parents_collections.at(pair_chromo_location.second)->get_Chromosome_obj();

		/// calling the crossover and mutation function inside the chromosome object
		chromo1.crossoverBypt_mutationByrate(chromo2);

		// parent validation
		parent_validation1 = chromo1.wholegene_operation2dec_data();
		parent_validation2 = chromo2.wholegene_operation2dec_data();

		// for chromosome validation check 
		if (parent_validation1 == true && parent_validation2 == true) {
			// only both chromosome is passed in validation check, usedGenes is then added to usedGenes set, which means the two offsprings r formed

			// after verification, both chromosome passed and can be converted back dec2string and stored in wholegene string vector in order to update new changes
			chromo1.create_wholegene();
			chromo2.create_wholegene();
			this->jobs_ptr_collections.push_back(this->m_runCaseFactory->genRunCase(std::move(chromo1), this->m_convertor));
			this->jobs_ptr_collections.push_back(this->m_runCaseFactory->genRunCase(std::move(chromo2), this->m_convertor));

			// add both locations to usedGenes set
			usedGenes.insert(pair_chromo_location.first);
			usedGenes.insert(pair_chromo_location.second);

			// can select the new chromosome pop location
			deReselect = true;

			// trial count to limit the while loop of redoing the crossover and mutation, reset to 0
			trial_count = 0;
			reselectTryCount = 0;
		}
		else if (trial_count >= this->m_param.trial_count_limit) {
			// when trial_count is reached or out of the limit of trial_count
			// can select the new chromosome pop location
			deReselect = true;
			// trial count to limit the while loop of redoing the crossover and mutation,  reset to 0
			trial_count = 0;
			++reselectTryCount;
		}
		else {
			// one of chromosome is failed in validation check, 
			//  still using the same chromosome pop location
			deReselect = false;
			// trial count to limit the while loop of redoing the crossover and mutation, ++1
			++trial_count;
		}
	}
}

void GA_Module::cal_jobs_intermediate_parents_and_jobs_ptr() {
	//copies the top fitness prev-object ptrs into the intermediate pop

	this->jobs_intermediate_parents_collections.clear();
	this->jobs_ptr_collections.clear();

	int parent_number_select = 0;
	int np;

	do {
		// duplicate according to the fitness ratio in order to distribute the gene pop into wheel chart
		if (this->pops_fit.at(parent_number_select) * (double)this->pops_fit.size() > 1) {
			for (np = 0; 
				np < std::min((int)(std::round(this->pops_fit.at(parent_number_select) * (double)this->m_param.Populcation_size)), this->m_param.Populcation_size); ++np) {
				// copy the org ptr to fill up the spaces waiting for copy
				this->jobs_intermediate_parents_collections.emplace_back(this->job_previousgen_ptr_collections.at(parent_number_select));   
			}
			++parent_number_select;
		}
		else {
			// copy the org ptr to fill up the spaces waiting for copy
			this->jobs_intermediate_parents_collections.emplace_back(this->job_previousgen_ptr_collections.at(parent_number_select));  
			++parent_number_select;
		}
	} while ((int)this->jobs_intermediate_parents_collections.size() < this->m_param.Populcation_size && parent_number_select < (int)this->pops_fit.size());

	//  ~~~~~~~~~~~~~~~~~ deal with the obj out of the range ~~~~~~~~~~~~~~ ///
	if (this->jobs_intermediate_parents_collections.size() > this->m_param.Populcation_size) {
		//erase out the out of boundary intermediate_parents (within the range of population size), the intermediate_bin_parents.size = Populcation_size  here
		this->jobs_intermediate_parents_collections.erase(this->jobs_intermediate_parents_collections.begin() + this->m_param.Populcation_size\
			, this->jobs_intermediate_parents_collections.end());
	}
}

void GA_Module::save_job_previousgen_ptr_collections() {
	// transfer the existing ptrs to history for delete
	for (auto& ptr : this->job_previousgen_ptr_collections)
		this->add_to_history(ptr);

	// move all the ptrs over from jobs_ptr_collections and rebuild the vector
	this->job_previousgen_ptr_collections = std::move(this->jobs_ptr_collections);
	this->jobs_ptr_collections = {};
}

void GA_Module::inheritance_of_elites() {

	// --- --- --- --- --- to differentiate out SGA and NGA
	if (this->m_param.GA_type == input_param::GAtype::NGA) {
		if (this->m_param.current_generation > 0) {
			// the following generation, would have the previous records, elites could have inherited 

			// update the elite population (size = Elite_size ) into the current just generation calulation pointer list jobs_ptr_collections(inside r all offsprings)
			//     for 1st generation job_previousgen_ptr_collections is empty, so nth change in jobs_ptr_collections
			int i = 0;
			for (auto itr = this->job_previousgen_ptr_collections.begin(); itr != this->job_previousgen_ptr_collections.end() && i < this->m_param.Elite_size; ++itr) {
				this->jobs_ptr_collections.insert(this->jobs_ptr_collections.begin(), *itr);
				++i;
			}

			//remove copied over ptrs to avoid delete
			this->job_previousgen_ptr_collections.erase(
				this->job_previousgen_ptr_collections.begin(), 
				this->job_previousgen_ptr_collections.begin() + std::min((int)(this->job_previousgen_ptr_collections.size()), this->m_param.Elite_size));
		}
	}
	else if (this->m_param.GA_type == input_param::GAtype::SGA){
		if (this->m_param.current_sub_generation == 0) {
			if (this->m_param.best_ref_opt_job != nullptr) {
				// the following Large (outer) generation, but 1st sub generation and first intersections in current_unfreeze_loc

				// still need to inherit 1 best_ref_opt_job int  the beginning of job_previousgen_ptr_collections in order to start the next generation and inner sub generation
				this->jobs_ptr_collections.emplace_back(this->m_param.best_ref_opt_job);

				// clear previous outer generation wholegene string locations
				(this->m_param.best_ref_opt_job)->get_Chromosome_obj().clear_all_locations();

				// refresh the wholegene_ref because different sub generation might have diff intersection been optimized so wholegene string is also diff
				// reset wholegene_ref for the best_ref_opt_job
				(this->m_param.best_ref_opt_job)->get_Chromosome_obj().create_wholegene();

				//remove copied over ptrs to avoid delete
				if(this->job_previousgen_ptr_collections.size())
					this->job_previousgen_ptr_collections.erase(this->job_previousgen_ptr_collections.begin());
			}
		}
		else {
			// the following generation, would have the previous records, elites could have inherited 

			// update the elite population (size = Elite_size ) into the current just generation calulation pointer list jobs_ptr_collections(inside r all offsprings)
			//     for 1st generation job_previousgen_ptr_collections is empty, so nth change in jobs_ptr_collections
			int i = 0;
			for (auto itr = this->job_previousgen_ptr_collections.begin(); itr != this->job_previousgen_ptr_collections.end() && i < this->m_param.Elite_size; ++itr) {
				this->jobs_ptr_collections.insert(this->jobs_ptr_collections.begin(), *itr);
				++i;
			}

			//remove copied over ptrs to avoid delete
			this->job_previousgen_ptr_collections.erase(
				this->job_previousgen_ptr_collections.begin(),
				this->job_previousgen_ptr_collections.begin() + std::min((int)(this->job_previousgen_ptr_collections.size()), this->m_param.Elite_size));
		}
	}
}


// -- -- -- -- -- -- sort the job thru its fitness 
void GA_Module::sort_jobs_by_avd() {

	sort(this->jobs_ptr_collections.begin(), this->jobs_ptr_collections.end(), sort_by_avd);

	//  ~~~~~~~~~~~~~~~~~ deal with the obj out of the range ~~~~~~~~~~~~~~ ///
	if ((int)this->jobs_ptr_collections.size() > this->m_param.Populcation_size) {
		// schedule delete of excess ptrs, as it might be from elites, and will need to be checked against prev job ptrs
		for (auto job_ptr_it = this->jobs_ptr_collections.begin() + this->m_param.Populcation_size; 
				job_ptr_it != this->jobs_ptr_collections.end(); ++job_ptr_it) {
			this->add_to_history(*job_ptr_it);
		}

		//erase out the out of boundary jobs_ptr_collections (within the range of population size), the jobs_ptr_collections.size = Populcation_size  here
		this->jobs_ptr_collections.erase(
			this->jobs_ptr_collections.begin() + this->m_param.Populcation_size, 
			this->jobs_ptr_collections.end());
	}
}

//-- - -- - -- - -- - -- - fitness calculation -- - -- - -- - -- - -- - //
void GA_Module::fitness_cal() {
	//clear the pops_fit before following update new fitness from jobs in this current generation
	this->pops_fit.clear();

	double max_avd = (*max_element(this->jobs_ptr_collections.begin(), \
		this->jobs_ptr_collections.end(), search_max_avd))->getEvalValue();
	double min_avd = (*min_element(this->jobs_ptr_collections.begin(), \
		this->jobs_ptr_collections.end(), search_min_avd))->getEvalValue();

	// predefine the properties of A_i value and fitness
	std::vector<double> temp_A_i;
	double A_i, fit;

	// calculate the each A_i according to all cases max and min average delay  and its power amplification factor
	for (I_RunCase* n : this->jobs_ptr_collections) {
		A_i = exp((max_avd - n->getEvalValue()) / std::max(std::numeric_limits<double>::epsilon(), (max_avd - min_avd)) * this->m_param.Power_factor);
		// this it temporary vector for the collection of A_i value on each job cases
		temp_A_i.push_back(A_i);
	}

	// get the ratio of A_i collections
	double sumOftemp_A_i = accumulate(temp_A_i.begin(), temp_A_i.end(), 0);
	assert(sumOftemp_A_i != 0);
	for (double n : temp_A_i) {
		fit = n / sumOftemp_A_i;
		this->pops_fit.push_back(fit);
	}

}

void GA_Module::clean_hisotry() {
	for (I_RunCase* job : this->history_jobs)
		delete job;
	this->history_jobs.clear();
}

std::pair<int, int> GA_Module::generate_Rand_crossover_selection(std::set<int>& usedGenes) {
	//selection must not repeat
	int no_1_pop = -1;
	int no_2_pop = -1;
	int size = std::min((int)this->jobs_intermediate_parents_collections.size(), this->m_param.Populcation_size);

	do {
		// pick up 1st pop randomly as parent
		no_1_pop = rand() % std::max(size, 1);
	} while (usedGenes.count(no_1_pop));

	//    population cuold be randomly selected except from reselecting the same location of the population
	do {
		no_2_pop = rand() % std::max(size, 1);
	} while (no_1_pop == no_2_pop || usedGenes.count(no_2_pop));

	// push up no_1_pop, no_2_pop selected parent into list as first element
	return std::make_pair(no_1_pop, no_2_pop);
}

// --- --- --- --- --- execution of the optimization of GA --- --- --- --- ---
void GA_Module::Run_NGA() {
	// ------------------ initial set up -----------------------//
	// reset the type of GA strategy for this run
	this->m_param.GA_type = input_param::GAtype::NGA;
	//reset everything
	this->initRuns();
			
	//reset the Runner
	this->m_runner->reset();

	// start the timer
	Timer GA_timer;
	GA_timer.start();

	//------------------ true runnning of NGA -----------------------//

	// 1st generation
	if (this->m_param.use_initial_seed)
		this->create_initial_seedjob();
	this->fill_random_pop();

	for (this->m_param.current_generation; this->m_param.current_generation < this->m_param.no_generation; ++this->m_param.current_generation) {
		// print out level on the details 
		if (this->m_param.debugPrintOutLevel) {
			std::cout << "\n //--- --- --- --- --- --- --- " << this->m_param.current_generation + 1 << "th generation --- --- --- --- --- --- --- // " << std::endl;
			std::cout << "Generated " << this->jobs_ptr_collections.size() << " population" << std::endl;
		}

		// --- --- --- --- --- almost every generation undergoes following process --- --- --- --- ---

		// do computation
		this->m_runner->setCaseEditors(this->jobs_ptr_collections);
		this->m_runner->runAll();
		if (this->m_param.debugPrintOutLevel)
			std::cout << "Completed runAll" << std::endl;

		// add in the elites
		inheritance_of_elites();
		// sort everthing
		sort_jobs_by_avd();
		fitness_cal();

		save_job_previousgen_ptr_collections();

		// add to best result from jobs
		this->save_best_results();

		//print debug stuff
		debugPrintOut();

		//Exit iff last gen
		if (this->m_param.current_generation == this->m_param.no_generation - 1) {
			Print_best_gene();
			break;
		}

		// --- --- --- --- --- preparation for start of  crossover() --- --- --- --- ---
		cal_jobs_intermediate_parents_and_jobs_ptr();

		this->crossover_mutation();

		this->clean_hisotry();
	}

	//esc from last sub-gen
	//clear job_previousgen_ptr_collections
	for (auto& jobPtr : this->job_previousgen_ptr_collections)
		this->add_to_history(jobPtr);
	this->job_previousgen_ptr_collections.clear();
	//dump history
	this->clean_hisotry();
	//puts the best chromo ptr into convertor for export
	this->exportBestChromosomeToConvertor();
	// in order to count up the running time
	GA_timer.Stop();
}

void GA_Module::Run_SGA() {
	//// ------------------ initial set up -----------------------////
	// reset the type of GA strategy for this run
	this->m_param.GA_type = input_param::GAtype::SGA;
	//reset everything
	this->initRuns();

	//reset the Runner
	this->m_runner->reset();

	// start the timer
	Timer GA_timer;
	GA_timer.start();

	//make sure if there is a initial guess, the gene is setup properly
	if (this->m_param.use_initial_seed) {
		if (this->m_param.sga_cyc_subgen)
			this->m_param.current_isCycOnlyOptGen = true;

		this->create_initial_seedjob();
	}

	///------------------ true runnning of SGA -----------------------//
	for (this->m_param.current_generation; this->m_param.current_generation < this->m_param.no_generation; ++this->m_param.current_generation) {

		// print out level of the details on Large generation
		if (this->m_param.debugPrintOutLevel) {
			std::cout << "\n //--- --- --- --- --- --- --- " << this->m_param.current_generation + 1 << "th Large generation --- --- --- --- --- --- --- // " << std::endl;
		}

		///------------------ Gen to opt cycle only -----------------------//
		// this is a shit impl but merging sub-gen and normal gen to make a function is too complicated -JLo
		if (this->m_param.sga_cyc_subgen) {
			// print out level of the details on intersections in this session of unfreeze list 
			if (this->m_param.debugPrintOutLevel) {
				std::cout << "\n   * * *  * * *  * * *  * * * intersections' cycles are being optimized * * *  * * *  * * *  * * *    " << std::endl;
			}

			this->m_param.current_isCycOnlyOptGen = true;
			//------------------ loop through the cycle opt only sub gens -----------------------//

			// 1st generation 
			this->fill_random_pop();
			
			for (this->m_param.current_sub_generation = 0; this->m_param.current_sub_generation < this->m_param.no_sub_generation; ++this->m_param.current_sub_generation) {
				// print out level of the details on sub generation
				if (this->m_param.debugPrintOutLevel) {
					std::cout << "\n //~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ " << this->m_param.current_sub_generation + 1 << "th sub generation ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ // " << std::endl;
					std::cout << "Generated " << this->jobs_ptr_collections.size() << " population" << std::endl;
				}

				// --- --- --- --- --- almost every generation undergoes following  --- --- --- --- ---
					
				// do the computation
				this->m_runner->setCaseEditors(this->jobs_ptr_collections);
				this->m_runner->runAll();
				if (this->m_param.debugPrintOutLevel)
					std::cout << "Completed runAll" << std::endl;

				inheritance_of_elites();
				sort_jobs_by_avd();
				fitness_cal();

				save_job_previousgen_ptr_collections();

				// add to best result from jobs
				this->save_best_results();

				//print debug stuff
				debugPrintOut();

				//Exit iff end generation 
				if (this->m_param.current_sub_generation == this->m_param.no_sub_generation - 1) {
					//save the last best result as seed
					this->m_param.best_ref_opt_job = this->job_previousgen_ptr_collections.front();
					if (this->m_param.debugPrintOutLevel)
						Print_best_gene();
					break;
				}

				this->cal_jobs_intermediate_parents_and_jobs_ptr();
				this->crossover_mutation();
			}

			//esc from last sub-gen
			//clear job_previousgen_ptr_collections
			for (auto& jobPtr : this->job_previousgen_ptr_collections)
				this->add_to_history(jobPtr);
			this->job_previousgen_ptr_collections.clear();
			//dump history
			this->clean_hisotry();
		}
		this->m_param.current_isCycOnlyOptGen = false;

		///------------------ loop through the unfreeze list -----------------------//
		for (this->m_param.current_unfreeze_loc = 0; this->m_param.current_unfreeze_loc < (int)this->m_param.unfreeze_list.size(); ++this->m_param.current_unfreeze_loc) {

			// print out level of the details on intersections in this session of unfreeze list 
			if (this->m_param.debugPrintOutLevel) {
				std::cout << "\n   * * *  * * *  * * *  * * * intersections ";
				for (auto unfreezed_intersecs : this->m_param.unfreeze_list.at(this->m_param.current_unfreeze_loc)) {
					std::cout << unfreezed_intersecs << " ";
				}
				std::cout << " are being optimized * * *  * * *  * * *  * * *    " << std::endl;
			}

			this->fill_random_pop();

			///------------------ loop through the sub gens -----------------------//
			for (this->m_param.current_sub_generation = 0; this->m_param.current_sub_generation < this->m_param.no_sub_generation; ++this->m_param.current_sub_generation) {

				// print out level of the details on sub generation
				if (this->m_param.debugPrintOutLevel) {
					std::cout << "\n //~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ " << this->m_param.current_sub_generation + 1 << "th sub generation ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ // " << std::endl;
					std::cout << "Generated " << this->jobs_ptr_collections.size() << " population" << std::endl;
				}

				// --- --- --- --- --- almost every generation undergoes following  --- --- --- --- ---
					
				this->m_runner->setCaseEditors(this->jobs_ptr_collections);
				this->m_runner->runAll();
				if (this->m_param.debugPrintOutLevel)
					std::cout << "Completed runAll" << std::endl;

				inheritance_of_elites();
				sort_jobs_by_avd();
				fitness_cal();

				save_job_previousgen_ptr_collections();

				// add to best result from jobs
				this->save_best_results();

				//print debug stuff
				debugPrintOut();

				//Exit iff end generation 
				if (this->m_param.current_sub_generation == this->m_param.no_sub_generation - 1) {
					//save the last best result as seed
					this->m_param.best_ref_opt_job = this->job_previousgen_ptr_collections.front();
					if (this->m_param.debugPrintOutLevel)
						Print_best_gene();
					break;
				}

				// --- --- --- --- --- preparation for start of  crossover() --- --- --- --- ---
				this->cal_jobs_intermediate_parents_and_jobs_ptr();

				this->crossover_mutation();
			}

			//esc from last sub-gen
			//clear job_previousgen_ptr_collections
			for (auto& jobPtr : this->job_previousgen_ptr_collections)
				this->add_to_history(jobPtr);
			this->job_previousgen_ptr_collections.clear();
			//dump history
			this->clean_hisotry();
		}
		//reset counter as chromo depends on this value
		this->m_param.current_unfreeze_loc = 0;

	} // end of the Large (outer) generation

	//puts the best chromo ptr into convertor for export
	this->exportBestChromosomeToConvertor();
	// in order to count up the running time
	GA_timer.Stop();
}

void GA_Module::add_to_history(I_RunCase* jobPtr) {
	bool found = !jobPtr;
	// search in best_jobs
	if(!found)
		for(auto& bestPtr : this->best_jobs)
			if (bestPtr == jobPtr) {
				found = true;
				break;
			}
	// comp best_ref_opt_job
	if (!found)
		if (this->m_param.best_ref_opt_job == jobPtr)
			found = true;

	if(!found && !this->history_jobs.count(jobPtr))
		this->history_jobs.insert(jobPtr);
}

bool GA_Module::sort_by_avd(I_RunCase* a, I_RunCase* b) {
	return (a->getEvalValue() < b->getEvalValue());
}

double GA_Module::accumulate_avd(double result, I_RunCase* obj) {
	return  result + (obj->getEvalValue());
}

bool GA_Module::search_max_avd(I_RunCase* a, I_RunCase* b) {
	return (a->getEvalValue() < b->getEvalValue());
}

bool GA_Module::search_min_avd(I_RunCase* a, I_RunCase* b) {
	return (a->getEvalValue() < b->getEvalValue());
}

void GA_Module::save_best_results() {
	// copy directly to the best_jobs vector // call once per generation
	this->best_jobs.emplace_back(this->job_previousgen_ptr_collections.front());
	// copy directly to the best pops_fit vector // call once per generation
	this->best_pops_fit.emplace_back(this->pops_fit.front());
}

std::vector<double> GA_Module::get_best_pops_avd() {
	std::vector<double> best_pops_avd;
	for (auto best_job : this->best_jobs) {
		best_pops_avd.emplace_back(best_job->getEvalValue());
	};
	return best_pops_avd;
}
std::vector<double> GA_Module::get_best_pops_fit() {
	return this->best_pops_fit;
}

void GA_Module::exportBestChromosomeToConvertor() {
	if (this->m_convertor)
		this->m_convertor->setTChromosome(&(this->best_jobs.back()->get_Chromosome_obj()));
}

void GA_Module::debugPrintOut() {
	if (this->m_param.debugPrintOutLevel) {
		switch (this->m_param.debugPrintOutLevel) {
		case 1:
		{
			//quick print of best gene eval values
			std::cout << "< ";
			auto it = this->job_previousgen_ptr_collections.begin();
			for (int i = 0; i < this->m_param.Elite_size && it != this->job_previousgen_ptr_collections.end(); ++i) {
				std::cout << (*it)->getEvalValue() << " ";
				++it;
			}
			std::cout << ">" << std::endl;
		}
			break;

		case 2:
			//detailed all best gene
			Print_best_genes_in_each_gen(this->m_param.Elite_size);
			break;

		case 3:
			//detailed all best gene and worst
			Print_best_genes_in_each_gen(std::min(8, this->m_param.Populcation_size));
			Print_worst_genes_in_each_gen(std::min(8, this->m_param.Populcation_size));
			break;
		};
	}

}

void GA_Module::Print_best_gene() {
	std::cout << "->> the top chromosome are printed out." << std::endl;;
	I_RunCase* bestJob = this->job_previousgen_ptr_collections.front();
	std::cout << "average daley: " << bestJob->getEvalValue();
	bestJob->get_Chromosome_obj().print_dec_data_structure();
}

void GA_Module::Print_best_genes_in_each_gen(int numOfjobs) {

	int orderOfjob = 0;
	auto it = this->job_previousgen_ptr_collections.begin();

	std::cout << "->> the top " << numOfjobs << " chromosome are printed out.";

	while (orderOfjob < numOfjobs && it != this->job_previousgen_ptr_collections.end()) {

		std::cout << "\n"<< orderOfjob << "th average daley: " << (*it)->getEvalValue();

		(*it)->get_Chromosome_obj().print_dec_data_structure();

		++it;  
		++orderOfjob;
	}
}

void GA_Module::Print_worst_genes_in_each_gen(int numOfjobs) {

	int orderOfjob = 0;
	auto it = this->job_previousgen_ptr_collections.rbegin();

	std::cout << "->> the back " << numOfjobs << " chromosome are printed out.";

	while (orderOfjob < numOfjobs && it != this->job_previousgen_ptr_collections.rend()) {

		std::cout << "\n" << (this->m_param.Populcation_size - numOfjobs - orderOfjob) << "th average daley: " << (*it)->getEvalValue();

		(*it)->get_Chromosome_obj().print_dec_data_structure();

		++it;  
		++orderOfjob;
	}
}

}	//namespace DISCO2_GA



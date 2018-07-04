/*
* Generates square pulse current commands.
*/

#include <default_gui_model.h>
#include <data_recorder.h> /* including data recorder to lock the data acquisition on the protocol*/
#include <string>
#include <time.h>

const int Length_Randomization = 4000; 
const double Iout_Flag_for_Inactive = 0.0037; 

class cIstep : public DefaultGUIModel {
	
	public:
		cIstep(void);
		virtual ~cIstep(void);
		virtual void
		execute(void);
	
	protected:
		virtual void update(DefaultGUIModel::update_flags_t);
	
	private:
		int active;
		double V, Iout;
		double Start_Vector[Length_Randomization];
		double Stop_Vector[Length_Randomization];
		double Amplitude_Vector[Length_Randomization];
		bool step_on;
		double period;
		int steps;
		double systime;
		long long count;
		double dt;
		double delay;
		double random_delay;
		double Amin, Amax;
		int Nsteps;
		int step_counter;
		double duration;
		double deltaI;
		QString root_dir, filename;
		void initParameters();
		void initRandomization(void);
		void storeRandomization(void);
		void set_filename(void);
		double randZeroToOne();
};

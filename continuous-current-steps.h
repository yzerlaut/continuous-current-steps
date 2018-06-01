/*
Copyright (C) 2011 Georgia Institute of Technology

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

/*
* Generates square pulse current commands.
*/

#include <default_gui_model.h>
#include <string>

const int Length_Randomization = 2000; 

class cIstep : public DefaultGUIModel {
	
	public:
		cIstep(void);
		virtual ~cIstep(void);
		virtual void
		execute(void);
	
	protected:
		virtual void update(DefaultGUIModel::update_flags_t);
	
	private:
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
		void initParameters();
		void initRandomization(void);
		double randZeroToOne();
};

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

#include <continuous-current-steps.h>
#include <iostream>

extern "C" Plugin::Object *createRTXIPlugin(void) {
	return new cIstep();
}

static DefaultGUIModel::variable_t vars[] =
{
	{ "Vin", "", DefaultGUIModel::INPUT, },
	{ "Iout", "", DefaultGUIModel::OUTPUT, },
	{ "Duration (ms)", "Duration of one cycle", DefaultGUIModel::PARAMETER
		| DefaultGUIModel::DOUBLE, },
	{ "Delay (ms)", "Fixed Time until step starts from beginning of cycle",
		DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Random Delay (ms)", "Random Time until step starts from beginning of cycle",
		DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Min Amp (pA)", "Starting current of the steps",
		DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Max Amp (pA)", "Ending current of the steps", DefaultGUIModel::PARAMETER
		| DefaultGUIModel::DOUBLE, },
	{ "Increments", "How many steps to take between min and max",
		DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER, },
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

#include <stdlib.h>

double cIstep::randZeroToOne()
{
  return rand() / (RAND_MAX + 1.);
}

cIstep::cIstep(void) 
  : DefaultGUIModel("Continuous Randomized Current Steps", ::vars, ::num_vars)
{
  initParameters();
  setWhatsThis("<p><b>cont-I-Steps:</b><br>This module generates a continuous randomized train of current injection pulses with amplitudes between a user-specified minimum and maximum.</p>");
	createGUI(vars,
		  num_vars);
	update(INIT);
	refresh();
	QTimer::singleShot(0, this, SLOT(resizeMe()));
}

cIstep::~cIstep(void) {}

void cIstep::execute(void) {

  // V = input(0);
  systime = count * period; // time in seconds

  if ((systime>Start_Vector[step_counter % Length_Randomization]) && !(step_on)) {
      step_on = true; // turn on step
  } 

  if ((systime>Stop_Vector[step_counter % Length_Randomization]) && (step_on)) {
      step_on = false; // turn off step
      step_counter++;
  }

  // std::cout << ' ' << step_on << ' ' << step_counter << ' ' << Start_Vector[step_counter] << ' ' << systime << '\n'; 

  Iout = 0;
  if (step_on) {
    Iout = Amplitude_Vector[step_counter];
    }
  output(0) = Iout * 1e-12;
  count++;
}

void
cIstep::initParameters(void)
{
  systime = 0.0;
  count = 0;
  period = RT::System::getInstance()->getPeriod() * 1e-9; // s
  delay = 100; // ms
  duration = 200; // ms
  random_delay = 600; // ms
  Amin = -50. ; // pA
  Amax = 200. ; // pA
  Nsteps = 5 ;
  step_on = false;
  step_counter = 0;
}

void
cIstep::initRandomization(void)
{
  // prepare a vector for the random shuffling
  std::vector<int> myvector;
  myvector.push_back(Amin);
  for (int i=1; i<Nsteps; ++i) myvector.push_back(Amin+i*(Amax-Amin)/(Nsteps-1));
  // random shuffling
  int j=0;
  for (int i=0; i<Length_Randomization; ++i) {
    if ((i%Nsteps)==0) std::random_shuffle( myvector.begin(), myvector.end() );
    Amplitude_Vector[i] = myvector[i%Nsteps];
  }
  // now defining the temporal randomization
  Start_Vector[0] = 1e-3*delay ; // second
  Stop_Vector[0] = Start_Vector[0]+1e-3*duration ; // second
  for (int i=1;i<Length_Randomization;i++) {
    Start_Vector[i] = Stop_Vector[i-1]+1e-3*delay+
      1e-3*randZeroToOne()*random_delay; // second
    Stop_Vector[i] = Start_Vector[i]+1e-3*duration; // second
  }


  step_on = false;
}

void cIstep::update(DefaultGUIModel::update_flags_t flag) {
	switch (flag) {
		case INIT:
			setParameter("Duration (ms)", duration);
			setParameter("Delay (ms)", delay);
			setParameter("Random Delay (ms)", random_delay);
			setParameter("Min Amp (pA)", Amin);
			setParameter("Max Amp (pA)", Amax);
			setParameter("Increments", Nsteps);
			break;

		case MODIFY:
			duration = getParameter("Duration (ms)").toDouble();
			delay = getParameter("Delay (ms)").toDouble();
			random_delay = getParameter("Random Delay (ms)").toDouble();
			Amin = getParameter("Min Amp (pA)").toDouble();
			Amax = getParameter("Max Amp (pA)").toDouble();
			Nsteps = getParameter("Increments").toInt();
			break;

		case PAUSE:
			output(0) = 0;
	
		case PERIOD:
		  period = RT::System::getInstance()->getPeriod() * 1e-9; // s
	
		default:
			break;
	}
	
	// Some Error Checking for fun
	
	if (duration <= 0) {
		duration = 1;
		setParameter("Duration (ms)", duration);
	}
	
	if (Amin > Amax) {
		Amax = Amin;
		setParameter("Min Amp (pA)", Amin);
		setParameter("Max Amp (pA)", Amax);
	}
	
	if (Nsteps < 1) {
		Nsteps = 1;
		setParameter("Increments", Nsteps);
	}

	// if (delay <= 0 || delay > duration * duty / 100) {
	if (delay <= 1e-3) {
		delay = 1e-3;
		setParameter("Delay (ms)", delay);
	}
	if (random_delay <= 1e-3) {
		random_delay = 1e-3;
		setParameter("Random Delay (ms)", random_delay);
	}
	
	//Define deltaI based on params
	if (Nsteps > 1) {
		deltaI = (Amax - Amin) / (Nsteps - 1);
	}
	else {
		deltaI = 0;
	}
	
	// Initialize counters
	step_counter = 0;
	// Randomize
	initRandomization();
}


#include <TSystem.h>
#include <iostream>

using namespace std;

void calib()
{
//	gSystem->Load("/home/muzalevsky/soft/AculUti/libAculData.so");
//	gSystem->Load("/home/muzalevsky/soft/AculUti/libTELoss.so");

	AculCalibration cal;
	cal.SetWorkDirectory("/home/muzalevsky/work/exp1803/data/siCal/newCal/");
	cal.SetParFileName("/home/muzalevsky/soft/AculUti/macroexp18/parforcal.par");	
	cal.SetInputRootFile("/home/muzalevsky/work/temp/output.root");
	cal.Init();	//takes parameters from .par
	cal.PrintInputParameters();

	cal.FindPedestals(0, 150);
	cal.Mycalc(400, 1300);
	cal.FindEnergyPedestals();
	cal.ShowFullCalibratedSpectra(0, 4095);

}



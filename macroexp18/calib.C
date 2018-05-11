#include <TSystem.h>
#include <iostream>

using namespace std;

void calib()
{
	gSystem->Load("/home/muzalevsky/soft/AculUti/libAculData.so");
	gSystem->Load("/home/muzalevsky/soft/AculUti/libTELoss.so");

	AculCalibration cal;
	cal.SetWorkDirectory("/home/muzalevsky/work/exp1803/data/dataCal/");
	cal.SetParFileName("/home/muzalevsky/soft/AculUti/macroexp18/parforcal.par");	
	cal.SetInputRootFile("/media/analysis_nas/exp201804/calib/si_1000_LR_02_0001.root");
	cal.Init();	//takes parameters from .par
	cal.PrintInputParameters();


	cal.FindPedestals(0, 50);
	cal.Mycalc(150, 500);
	cal.FindEnergyPedestals();
	cal.ShowFullCalibratedSpectra(0, 4095);

}



#include <TSystem.h>
#include <iostream>

using namespace std;

void calib()
{
	gSystem->Load("/home/muzalevsky/AculUti/libAculData.so");
	gSystem->Load("/home/muzalevsky/AculUti/libTELoss.so");

	AculCalibration cal;
	cal.SetWorkDirectory("/home/muzalevsky/AculUti/exp1803/cali1mm/calDL/");
	cal.SetParFileName("/home/muzalevsky/AculUti/macroexp18/parforcal.par");	
	cal.SetInputRootFile("/media/users_NAS/Muzalevsky/si_1000_LR_02_0001.root");
	cal.Init();	//takes parameters from .par
	cal.PrintInputParameters();


	cal.FindPedestals(0, 50);
	cal.Mycalc(150, 400);
	cal.FindEnergyPedestals();
	cal.ShowFullCalibratedSpectra(0, 4095);

}



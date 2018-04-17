#include <TSystem.h>
#include <iostream>

using namespace std;

void calib()
{
	gSystem->Load("/home/muzalevsky/AculUtils/libAculData.so");
	gSystem->Load("/home/muzalevsky/AculUtils/libTELoss.so");

	AculCalibration cal;
	cal.SetWorkDirectory("/home/muzalevsky/AculUtils/exp1804/");
	cal.SetParFileName("/home/muzalevsky/AculUtils/macro/parforcal.par");	
	cal.SetInputRootFile("/home/muzalevsky/work/calib/si_20_04.root");
	cal.Init();	//takes parameters from .par


	cal.FindPedestals(0, 200);
	cal.Mycalc(500, 1000);
	cal.FindEnergyPedestals();
	cal.ShowFullCalibratedSpectra(0, 4095);

}



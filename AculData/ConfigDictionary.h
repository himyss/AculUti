#ifndef CONFIGDICTIONARY_H
#define CONFIGDICTIONARY_H


#include "TObject.h"
#include "ReturnCodes.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include <string>
#include <map>
#include <sstream>
#include "TLorentzVector.h"
#include <iostream>
#include "TMath.h"
#include "TString.h"
#include "TTree.h"


class ConfigDictionary{
public:
	typedef std::map<std::string,std::string>::iterator CDIter;
	ConfigDictionary();
	ConfigDictionary(std::string);
	virtual ~ConfigDictionary(){};//empty virtual destructor

	ClassDef(ConfigDictionary,1);
	std::string ToString();
	void FromString(std::string);

	//These throw errors if couldn't find key:
	std::string GetString(std::string)throw(std::string);
	int	GetInt(std::string)throw(std::string);
	double GetDouble(std::string)throw(std::string);
	bool GetBool(std::string)throw(std::string);

	//These will always set 'something' into map:
	void SetString(std::string,std::string);
	void SetDouble(std::string,double);
	void SetInt(std::string,int);
	void SetBool(std::string,bool);
	
	CDIter Begin(){return configMap.begin();};
	CDIter End(){return configMap.end();};

private:
	std::map<std::string,std::string> configMap;
};

#endif

#include "ConfigDictionary.h"

ClassImp(ConfigDictionary);

//////////////////////////////////////////////////////////////////////////////
//	BEGIN_HTML
//	<p><font size="4"><b>Config(uration)Dictionary class</b></font></p>
//  <br>
// 	<i>Author: Bartlomiej Hnatio 2012-08-06</i>
//	<br><br>
//	This is very useful class to convert strings containing pairs "key"="value"
//	into fast dictionaries. This strings are often read from external files.
//	From dictionary created in such way one can easily
//	extract values in few supported formats (integer, double, boolean, string)
//	It can also be used in opposite direction: when dictionary is created with
// 	given values and keys, one can create config string, which is most
//	convenient to write in external files.
//	<br><br><br>
//	<u>Most simple example of usage:</u>
//	<br><br>Suppose you have two variables fD and fI:
//	<pre>-------------------
//    Double_t fD;
//    Int_t fI;
//
//fD = 3.14;
//fI = 2012;
//-----------------------</pre>
//	To save its parameters into file you can create ConfigDictionary class
//	instance and use <b>SetDouble</b> and <b>SetInt</b> functions to
//	insert parameters values with arbitrarly defined keys (let them be
//	really long and descriptive in this example):
//	<pre>---------------------
//ConfigDictionary CD;
//CD.SetDouble("Most important variable D",fD);
//CD.SetInt("Equally important integer I",fI);
//---------------------</pre>
//	Now configuration string saved in <b>ConfigDictionary</b> CD
//	can be obtained with <b>ToString</b> method and should look like:
//<pre>
//"Most important variable D"="3.14" "Equally important integer I"="2012"
//</pre>
//	<br> It can be easily saved to a file using simplest <i>fstream</i>
//	methods. And the advantage is that as a key one can use any string,
//	so configuration file created in such way is very self-explanatory.
//	<br>
//	<br>
//	Now lets suppose the opposite action - loading config from file.
//	Imagine, that you have 1000 objects of A class, which config was saved
//	to file - one object per line:
//	<pre>-----------------
//"Most important variable D"="3.14" "Equally important integer I"="2012"
//"Equally important integer I"="1011" "Most important variable D"="8.15"
//"Most important variable D"="13.16" "Equally important integer I"="10"
//(...)
//----------------</pre>
//	Please notice that order in which pairs of keys and values are placed
//  in file doesn't make any difference to the ConfigDictionary class.
//	To recreate objects, you just read each	line containing single
//	config string and then create with it ConfigDictionary object
//  using special constructor with string as argument,
//	or using <b>FromString</b> method. After this, you can
// 	extract parameters (using <b>GetDouble</b> and <b>GetInt</b> methods)
//	with using same keys as were used for saving, and ConfigDictionary
//	will return their values:
//	<pre>----------------------
//string line;//This line you read from file
//ConfigDictionary CD(line);
//fD = CD.GetDouble("Most important variable D");
//fI = CD.GetInt("Equally important integer I");
//--------------------</pre>
//<br>
//	And last but not least: what happens, if key requested by user doesn't
//	exist in dictionary? (This can be caused by many reasons, mostly errors on
// 	user side, but not always).
//	When you try to extract non-existent key using <b>Get</b> functions,
//	an exception is risen. In normal program it should end execution and print
//	some information about where program stopped working.
//	But lets say that you don't want that, i.e. program may use default
//  configs instead of those from files, or not all keys were that important.
//	You can surround all uses of <b>Get</b> methods with try/catch clause.
//	So when exception is risen, you will catch it, and decide, if you want
//	the program to stop running, or anything else. Simple example is
//	shown below:
//<pre>-------------------
//ConfigDictionary CD(some_string);
//try{//try to read important variables:
//  double d = CD.GetDouble("crucial var");
//  bool b = CD.GetBool("most important b");
//  int i = CD.GetInt("unique ID");
//  //...and so on
//}catch(std::string & e){//catch any kind of exception
//  Error("Some crucial variable wasn't read, ending program!");
//  return SOME_REALLY_BAD_ERROR_CODE;
//}
//try{//now less important, or optional:
//  string name = CD.GetString("least important variable");
//  double p = CD.GetDouble("optional parameter");
//  //...and so on
//}catch(std::string & f){
//  Info("Some optional variables wasn't read!");
//}
//--------------------</pre>
//	END_HTML
//////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
ConfigDictionary::ConfigDictionary(){
	//just empty map...
}

//_____________________________________________________________________________
ConfigDictionary::ConfigDictionary(std::string params){
	//Just creates dictionary using FromString method
	FromString(params);
}

//_____________________________________________________________________________
std::string ConfigDictionary::ToString(){
	//Builds string that can be easily saved to file.
	//Same format that can be read from files or from QtGui
	//This should work same way for all uses

	std::map<std::string,std::string>::iterator it;//Iterator to map elements
	std::stringstream ss;//stream helpful with adding strings
	//iterate whole dictionary:
	for (it=configMap.begin();it!=configMap.end();++it){
		//insert pairs "key"="value":
		ss<<"\""<<it->first<<"\""<<"="<<"\""<<it->second<<"\""<<" ";
	}
	return ss.str();
}

//_____________________________________________________________________________
void ConfigDictionary::FromString(std::string params){
	//params - TString containing list of key=value pairs
	//
	//Changes string formatted:
	//"key1"="value1" "key2"="value with spaces 2" ...
	//into map with keys and values
	//Useful in lots of I/O functions,
	//when we want to have a nice, readable format of data

	std::stringstream loading(params);
	std::string k,v;

	while(!loading.fail()){
		getline(loading,k,'\"');//removes everything to first "
		getline(loading,k,'\"');//All chars between first two "" are the key
		getline(loading,v,'\"');//removes all until third "
		getline(loading,v,'\"');//All between another pair of "" is the value
		if (!loading.fail())
		    configMap[k]=v;
	}
}

//_____________________________________________________________________________
std::string ConfigDictionary::GetString(std::string key)throw(std::string){
	//Extracts string from given key
	//(if it exist, otherwise raises exception)

	if (configMap.find(key) == configMap.end()){
		Error("ConfigDictionary::GetString",
				"Couldn't find the key: %s!",key.c_str());
		throw(key);
	}
	return configMap[key];
}

//_____________________________________________________________________________
int	ConfigDictionary::GetInt(std::string key)throw(std::string){
	//Extracts integer from given key
	//(if it exist, otherwise raises exception)
	if (configMap.find(key) == configMap.end()){
		Error("ConfigDictionary::GetInt",
				"Couldn't find the key: %s!",key.c_str());
		throw(key);
	}
	int returned=0;
	//Convert string to int:
	std::stringstream ss(configMap[key]);
	ss>>returned;
	return returned;
}

//_____________________________________________________________________________
double ConfigDictionary::GetDouble(std::string key)throw(std::string){
	//Extracts integer from given key
	//(if it exist, otherwise raises exception)
	if (configMap.find(key) == configMap.end()){
		Error("ConfigDictionary::GetDouble",
				"Couldn't find the key: %s!",key.c_str());
		throw(key);
	}
	double returned=0.0;
	//Convert string to double:
	std::stringstream ss(configMap[key]);
	ss>>returned;
	return returned;
}

//_____________________________________________________________________________
bool ConfigDictionary::GetBool(std::string key)throw(std::string){
	//Extracts boolean from given key
	//(if it exist, otherwise raises exception)
	if (configMap.find(key) == configMap.end()){
		Error("ConfigDictionary::GetBool",
				"Couldn't find the key: %s!",key.c_str());
		throw(key);
	}
	//Convert string to bool:
	if (configMap[key].compare("true") == 0)
		return true;
	else return false;
}

//_____________________________________________________________________________
void ConfigDictionary::SetString(std::string key,std::string value){
	//Sets value to key, no comments needed here...
	configMap[key] = value;
}

//_____________________________________________________________________________
void ConfigDictionary::SetDouble(std::string key,double value){
	//Sets value to key, converts double to string first
	std::stringstream ss;
	ss<<value;
	configMap[key] = ss.str();
}

//_____________________________________________________________________________
void ConfigDictionary::SetInt(std::string key,int value){
	//Sets value to key, converts int to string first
	std::stringstream ss;
	ss<<value;
	configMap[key] = ss.str();
}

//_____________________________________________________________________________
void ConfigDictionary::SetBool(std::string key,bool value){
	//Sets value to key, converts bool to string first
	if (value == true)	configMap[key] = "true";
	else				configMap[key] = "false";
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//
//  StdArg.hpp: main's command line parser
//
//  Created: 04/29/2003, Andriy Zatserklyaniy
//  
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


#ifndef StdArg_hpp
#define StdArg_hpp

#include <map>
#include <sstream>
#include <iostream>
#include <cstring>

using std::cout;    using std::cerr;    using std::endl;
using std::string;  using std::map;


// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// -----------------------------------------------------------------------------------
// class StdArg

//    Note: 
//         Intended to parse command line arguments in form of key-value pairs and 
//         flags in any order like:
//
//            flag   key-value     key-value      flag
//
//         For example:  ./a.out -debug -in input.dat -out output.dat -mc
//
//
//  For insructions on how to compile and to use this class, see example commented at
//  the end of this file
//
// -----------------------------------------------------------------------------------
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

class StdArg 
{

public:
  

  // Create public structure to output message error when bad input are used
  // -----------------------------------------------------------------------

  struct BadInput  
  {
    BadInput(const string& mess, const string& key) 
    {
      cerr<< "*** ERROR StdArg: " << mess << key <<endl;
    }

    BadInput(const string& mess, const string& key, const string& value) 
    {
      cerr<< "*** ERROR StdArg: " << mess << ". Key-value pair: " << key << " " << value <<endl;
    }
  };


  // Create public structure to parse a flag
  // ---------------------------------------
  
  struct FlagMap: public std::map<string, int> 
  {
    FlagMap& operator << (const string& flag) 
    {
      insert( std::pair<string, int>(flag,0) );
      return *this;
    }
  };


  // Create public structure to parse a key and its value
  // ----------------------------------------------------

  struct KeyValueMap: public std::map<string, const char*> 
  {
    KeyValueMap& operator << (const string& key) 
    {
      insert( std::pair<string, const char*>(key,0) );
      return *this;
    }
  };
  

  // User interface for input flags and keys
  // ---------------------------------------

  FlagMap&     flags;  
  KeyValueMap& keys;   



private:

  int _argc;
  char** _argv;
  
  // input/output map for flags 
  // --------------------------

      // Note: '1' means 'has only one field: key' 

  FlagMap     _i1, _o1;


  // input, output map for pairs key-value 
  // -------------------------------------

      // Note: '2' means 'has two fields: key and value' 

  KeyValueMap _i2, _o2;
  bool _check_minus;



public:

  StdArg(int argc, char** argv, bool check_minus=true):
    _argc(argc), _argv(argv), flags(_i1), keys(_i2), _check_minus(check_minus)
  {}


  // Process the reading of the keys and the flags from command line
  // ---------------------------------------------------------------

  void Process() 
  {
    int iarg = 1;

    while (++iarg < _argc) 
      {
	//if (iarg == 0) continue;
	string key(_argv[iarg]);


	// First: look in defined flags
	// ...........................

	if (_i1.find(key) != _i1.end())
	  { 
	    
	    // If the flag is => throw a "defined twice" error

	    if (_o1[key] != 0) throw BadInput("defined twice flag: ", key);


	    // If the flag is not found => pass to the next argument

	    else 
	      {
	      _o1[key]++;
	      continue;  
	    }
	  }


	// Second: if not found in flags, look in defined keys for key-value
	// .................................................................

	if (_i2.find(key) != _i2.end())
	  { 

	    // If the key is found => throw a "defined twice" error

	    if (_o2.find(key) != _o2.end()) throw BadInput("defined twice key: ", key);


	    // Store value and point to next key, or throw and error if no value are defined in command line

	    if (iarg+1 < _argc) _o2[key] = _argv[++iarg];  
	    else throw BadInput("no value for key ", key);
	  }


	// Third: If not a key nor a flag, throw and error
	// ...............................................

	else throw BadInput("no such key: ", key);

      } // end while

  } // end Process



  // Function to show flags read from input
  // --------------------------------------

  void ShowFlags() const 
  {
    for (FlagMap::const_iterator p = _o1.begin(); p!=_o1.end(); p++) 
      {
	cout<< p->first <<endl;
      }
  }



  // Function to show key-value pairs read from input
  // ------------------------------------------------

  void ShowKeys() const 
  {
    for (KeyValueMap::const_iterator p = _o2.begin(); p!=_o2.end(); p++) 
      {
	cout<< p->first <<" "<< p->second <<endl;
      }
  }



  // Boolean function to tell if flags are all defined
  // -------------------------------------------------

  bool Flag(const string& flag) const { return _o1.find(flag) != _o1.end(); }


  // Boolean function to tell if key-value pairs are all defined
  // -----------------------------------------------------------

  bool Key(const string& key) const 
  {
    return _o2.find(key) != _o2.end();
  }


  // Doesn't throw exception
  // -----------------------

  const char* Value(const string& key) 
  {
    if (!Key(key)) return 0;
    else           return _o2[key];
  }



  // Template convertor for (almost) arbitrary type
  // ----------------------------------------------

      // Note:  usage: Type n = instance.Get<Type>("-key");
  
  template<class T> T Get(const string& key) 
  {
    T val;
    char ch;
    if (!Key(key)) throw BadInput("no such key: ", key);
    std::istringstream buf(_o2[key]);
    if (!(buf >> val))  throw BadInput("conversation error", key, Value(key));
    else if (buf >> ch) throw BadInput("conversation error", key, Value(key));  // ok if trailing '\0'
    return val;
  }
}; // end class






// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// -----------------------------------------------------------------------------------
// StdArg::Get

//    Note: 
//         Specializations of Get functions -> it should be defined out of class 
//         declaration
//
// -----------------------------------------------------------------------------------
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



// Specialization for const char*
// ------------------------------

template<>  
const char* StdArg::Get<const char*>(const string& key) 
{
  if (!Key(key)) throw BadInput("no such key: ", key);
  if (_check_minus && *Value(key) == '-')
    throw BadInput("value should not start from '-'", key, Value(key));
  return Value(key);
}


// Specialization for string: allows few words in string
// -----------------------------------------------------

template<>  
string StdArg::Get<string>(const string& key) 
{
  if (!Key(key)) throw BadInput("no such key: ", key);
  if (_check_minus && *Value(key) == '-')
    throw BadInput("value should not start from '-'", key, Value(key));
  return string(Value(key));
}


// Specialization for bool
// -----------------------

template<>  
bool StdArg::Get<bool>(const string& key) 
{
  bool val;
  char ch;
  if (!Key(key)) throw BadInput("no such key: ", key);
  if (strcmp(_o2[key],"true") ==0) return true;
  if (strcmp(_o2[key],"false")==0) return false;
  std::istringstream buf(_o2[key]);
  if (!(buf >> val))  throw BadInput("conversation error", key, Value(key));
  else if (buf >> ch) throw BadInput("conversation error", key, Value(key));  // ok if trailing '\0'
  return val;
}


#endif  // StdArg_hpp





// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// -----------------------------------------------------------------------------------
//
// EXAMPLES OF HOW TO COMPILE AND USE THE CLASS
//
//
//  To compile and run (after renaming to StdArg.cpp):
// 
//     KCC StdArg.cpp ; ./a.out -debug -in input.dat -out output.dat -num 100 -mc -somebool true
// or
//     g++ StdArg.cpp ; ./a.out -debug -in input.dat -out output.dat -num 100 -mc -somebool true
// 
//
//
//  To use in main code:
//
//
//    int main(int argc, char *argv[]) {
//      StdArg arg(argc,argv);
//      arg.flags << "-mc" << "-debug";
//      arg.keys << "-in" << "-out" << "-num" << "-somebool";
// 
//      const char* ifile  = 0;
//      const char* ofile;
//      int num = 0;
//      bool isMC  = false;
//      bool debug = false;
//   
//      bool somebool = false;
// 
//      try {
//        arg.Process();
//        // for required values
//        ifile = arg.Get<const char*>("-in");
//        ofile = arg.Get<const char*>("-out");
//        // for optional values
//        if (arg.Key("-num")) num = arg.Get<int>("-num");
//        isMC  = arg.Flag("-mc")? true: false;
//        debug = arg.Flag("-debug")? true: false;
//     
//        if (arg.Key("-somebool")) somebool = arg.Get<bool>("-somebool");
//      }
//      catch (StdArg::BadInput)
//      {
//        if (argc > 1) cout<< "Input error" <<endl;
//        // usage in case of error or no parameters
//        cout<< "Usage:" <<endl;
//        cout<< argv[0] << " -in input_file -out output_file [-num events_to_process] [-somebool somebool] [-mc]" <<endl;
//        return 0;
//      }
//      // ... usage of variables
//      cout<< "input file   = " << ifile <<endl;
//      cout<< "output file  = " << ofile <<endl;
//      if (num) cout<< "num = " << num <<endl;
//      if (isMC)  cout<< "MC mode"    <<endl;
//      if (debug) cout<< "debug mode" <<endl;
//   
//      if (arg.Key("-somebool")) cout<< "somebool = " << somebool <<endl;
//    }
//
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// -----------------------------------------------------------------------------------

// Time-stamp: </w/c/overprop/common2.h, Sun, 27 Apr 2014, 19:44:38 EDT, http://wrfranklin.org/>

// Common stuff for all my C++ programs.

// I keep adding things to this file.  Newer versions are mostly compatible.  If not, increment
// the file name.

// This obsoletes an earlier version called misc.h.

// Don't include this if including cuda.h

#ifndef __COMMON_H__
#define __COMMON_H__

//================================================================
// Includes
//================================================================

#include <algorithm>
#include <deque>
#include <errno.h>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <locale>
#include <math.h>
#include <queue>
#include <set>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <sys/time.h>
#include <sys/timeb.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#include <utility>
#include <vector>

#if 0
#include <boost/foreach.hpp>
#include <boost/geometry/arithmetic/dot_product.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/unordered_set.hpp>
#include <boost/multi_array.hpp>
#endif

//================================================================
// Usings
//================================================================
//using namespace std;  // might cause problems with array
//using boost::array;   // apparently is in std when not in cuda
using std::string;
using std::cout;
using std::cerr;
using std::ios_base;
using std::locale;
using std::endl;
using std::ofstream;
using std::ostream;
using std::setprecision;
using std::setw;
using std::fixed;


//================================================================
// Misc macros
//================================================================

// Expand and stringify a macro.
#define EXPANDANDSTRING(macro) EXPANDANDSTRING2(macro)
#define EXPANDANDSTRING2(macro)  #macro

//================================================================
// Printing
//================================================================

// http://en.wikipedia.org/wiki/ANSI_escape_code
const std::string boldblacktty("\033[1;30m");   // tell tty to switch to bold black
const std::string redtty("\033[1;31m");   // tell tty to switch to bold red
const std::string greentty("\033[1;32m");   // tell tty to switch to bright green
const std::string bluetty("\033[34m");   // tell tty to switch to blue
const std::string magentatty("\033[1;35m");   // tell tty to switch to bright magenta
const std::string yellowbgtty("\033[1;43m");   // tell tty to switch to bright yellow background
const std::string underlinetty("\033[4m");   // tell tty to switch to underline
const std::string deftty("\033[0m");      // tell tty to switch back to default color

#define PRINT(arg)  #arg "= " <<(arg)     // Print an expression's name then its value, possibly
// followed by a comma or std::endl.
// Ex: std::cout << PRINTC(x) << PRINTN(y);
#define PRINTC(arg)  #arg << "= " << bluetty << (arg) << deftty << ", "
#define PRINTN(arg)  #arg << "= " << bluetty << (arg) << deftty << std::endl

// Assuming that the expression is a qualified name, e.g., a.b, print the name after the dot etc.
#define PRINTT(arg)  (index(#arg,'.')+1) << "=" << (arg)
#define PRINTCT(arg)  (index(#arg,'.')+1) << "=" << (arg) << ", "


//================================================================
// Storage Allocation
//================================================================

// Create a new array of size var of type.  Report error.
#define NEWA(var,type,size) { try  { if(0==(var=new type [(size)])) throw;} catch(...)  { std::cerr << "NEWA failed on " #var "=new " #type "[" #size "=" <<(size) << "]" << std::endl; exit(1); }}

// Create a new vector.  Report error.
#define NEWV(var,type) { try  { var=new vector<type>; } catch(...)  { std::cerr << "NEWV failed on " #var "=new vector<" #type ">"  << std::endl; exit(1); }}

// Create a new vector of given size.  Report error.
#define NEWVS(var,type,size) { try  { var=new vector<type>(size); } catch(...)  { std::cerr << "ERROR: NEWVS failed on " #var "=new vector<" #type ">(" #size "), "  << PRINTC(sizeof(type)) << PRINTN(size); exit(1); }}

// Create a new square vector of vectors of given size.
#define NEWVVS(var,type,size) { try  { var=new vector<vector<type> >(size, vector<type>(size)); } catch(...)  { std::cerr << "NEWVVS failed on " #var "=new vector<vector<" #type "> >(" #size "), "  << PRINTC(sizeof(type)) << PRINTN(size); exit(1); }}

// Create a new rectangular vector of vectors of given size.
#define NEWVVSS(var,type,size1,size2) { try  { var=new vector<vector<type> >(size1, vector<type>(size2)); } catch(...)  { std::cerr << "NEWVVS failed on " #var "=new vector<vector<" #type "> >(" #size1 ", " #size2 "), "  << PRINTC(sizeof(type)) << PRINTC(size1) << PRINTN(size2); exit(1); }}

// Create a new vector of given capacity.  Report error.
#define NEWVR(var,type,size) { try  { var=new vector<type>; var->reserve(size); } catch(...)  { std::cerr << "NEWVR failed on " #var "=new vector<" #type ">(" #size "),"  << PRINTC(sizeof(type)) << PRINTN(size); exit(1); }}


//================================================================
// Compact 'FOR'
//================================================================

#define foreach BOOST_FOREACH

// Compact FOR loop.   Ex:   FOR(i,0,3,a[i]=0);
#define FOR(foriterator,formin,formaxplus1,forbody) {for(int foriterator=formin; foriterator<formaxplus1; foriterator++) {forbody;};}

// Compact FOR loop w large body.   Ex:   FORR(i,0,3) {foo; baz}
#define FORR(foriterator,formin,formaxplus1) for(unsigned int foriterator=formin; foriterator<formaxplus1; foriterator++)


//================================================================
// Misc
//================================================================

void system2(const char *s) {   // stop those 'ignoring return value' warnings
  (void)(system(s)); }          // However that doesn't always work. ?!


//================================================================
// Timing
//================================================================

// TIMING - 1 METHOD

// Exec a std::string print its time, e.g. TIME(init()), and return delta clock time.
#define TIME(arg) ( (arg), Print_Time(#arg) )

// TIMING - MORE DETAILED METHOD

double CPUTotTime(0.);
double ClockTotTime(0.);
struct timeval *watchdog_tv;
struct timezone *watchdog_tz;

void Print_Current_Process_Memory_Used() {
  std::ostringstream sout;
  pid_t p = getpid();
  sout << "cat /proc/" << p << "/status | grep Vm" << std::ends;
  cerr  << sout.str()  << endl;
  system2(sout.str().c_str());
}

double Process_CPU_Time()  {  // Return CPU(user+system) time since start of process.
  timespec tp;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tp);
  return tp.tv_sec + tp.tv_nsec*1e-9;
}

//  DELTA_CPU_TIME Returns time in seconds since last Delta_Time.  Automatically initializes
//  itself on 1st call and returns 0.  Also, set time from process start (or from 1st
//  call) into CPUTotTime.

double Delta_CPU_Time()  {
  static double old_time = 0.0;
  double  delta;
  CPUTotTime = Process_CPU_Time();
  delta = CPUTotTime - old_time;
  old_time = CPUTotTime;
  return delta;  }

double Delta_Clock_Time()  {
  timespec tp;
  static double Clockbasetime;
  static double old_time = 0.0;
  double  delta;
  static bool first = true;
  clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
  ClockTotTime = tp.tv_sec + tp.tv_nsec*1e-9;
  if (first) {
    Clockbasetime = ClockTotTime;
    first = false;
  }
  //  cout << PRINTC(first) << PRINTC(Clockbasetime) <<  PRINTN(ClockTotTime);
  ClockTotTime -= Clockbasetime;
  delta = ClockTotTime - old_time;
  old_time = ClockTotTime;
  return delta;  }

// PRINT_TIME If REALLY_PRINT_TIME is defined, then print the time since the last call, with a
// message.  In any case, return the delta clock time.

double Print_Time(const string &msg)  {
#ifdef CUDA
  cudaDeviceSynchronize();
#endif
  double  CPUincrtime = Delta_CPU_Time();
  double  Clockincrtime = Delta_Clock_Time();
  string msg2 = msg;
  const int lmsg = 30;
#ifdef REALLYPRINTTIME
  msg2.resize(lmsg, ' ');   // truncate long messages
#ifdef PRINTTIMESHORT
  cout << fixed << setprecision(3) << "TIME thru " << setw(lmsg) << msg2 << ": CPU: " << setw(6) << CPUTotTime << " (d="
       << setw(6) << CPUincrtime <<  ")" << endl;
#else
#ifdef PRINTTIMEMED
  cout << fixed << setprecision(3) << "TIME thru " << setw(lmsg) << msg2 << ": CPU: " << CPUTotTime << " (d="
       << CPUincrtime <<  "), Clock: " << setw(3) << ClockTotTime << " (d="
       << setw(3) << Clockincrtime << ")" << endl;
#else
  cout << fixed << setprecision(3) << "TIME thru " << setw(lmsg) << msg2 << ": CPU: " << CPUTotTime << " (d="
       << CPUincrtime <<  "), Clock: " << setw(3) << ClockTotTime << " (d="
       << setw(3) << Clockincrtime << "), R: " << setw(3) << (CPUTotTime/ClockTotTime)
       << " (d=" << setw(3) << (CPUincrtime/Clockincrtime) << ")" << endl;
#endif
#endif
#endif
  return Clockincrtime;
}

void Print_CPU_Time(string s)  {  // Print CPU (user+system) time since start of process.
  cout << "Elapsed CPU time at " << s << ": " << Process_CPU_Time() << endl;
}


// PTIME
 void ptime(const char *const msg) {         // Write cumulative CPU process time to std::cerr.
   float t= ((float)clock())/CLOCKS_PER_SEC;
   std::cerr << magentatty << "-Cumulative CPU time thru " << msg << "= " << t << deftty << std::endl;
 }

//================================================================
// Assertions, Exceptions, Die
//================================================================

// Test an assertion.  On failure, print assertion, and die.
#define ASSERTA(predicate) if(!(predicate)) die2(#predicate " failed at line ",(__LINE__));

// Test an assertion.  On failure, print assertion and context, and die.
#define ASSERTC(predicate,context) if(!(predicate)) die2(#predicate " failed in " #context " at line ",(__LINE__));

// Test an assertion.   On failure, exec an expression(perhaps to print some vars), print a message, and die.
#define ASSERTE(predicate,context,exec) if(!(predicate)) {(exec); die2(#predicate " failed in " #context " at line ",(__LINE__));};

// Test an assertion.   On failure, exec an expression(perhaps to print some vars), print a message, and throw an exception.
#define ASSERTT(predicate,context,exec) if(!(predicate)) {(exec); cerr << redtty << #predicate << " failed in " << #context << " at line " <<__LINE__ << deftty << endl; throw 1; };


// TRY

// Try an expression.   On failure, exec another expression and escalate.
#define TRY2(tryexpr, failexpr) try { tryexpr; } catch(...) { std::cerr << "ERROR: exception caught at line " << __LINE__ << std::endl;  failexpr; throw; }

// Try an expression.   On failure, print the expression and escalate.
#define TRY1(tryexpr) try {tryexpr;} catch(...) { die2(#tryexpr " failed at line ",(__LINE__)); throw; }



void die(const std::string msg) {                   // die
  std::cerr << "ERROR: " << msg << std::endl;
  exit(1);
}

void die2(const std::string msg, const int i) {
  std::cerr << "ERROR: " << msg << i << std::endl;
  exit(1);
}

#endif


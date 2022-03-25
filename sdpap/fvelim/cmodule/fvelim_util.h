/* ********************************************************************************
   fvelim_util.h
   A module of fvelim
   C interface to call free variable elimination

   December 2010, Kenta KATO
   ******************************************************************************** */

#ifndef __FVELIM_UTIL_H__
#define __FVELIM_UTIL_H__

#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;

/* ============================================================
   Message
   ============================================================ */
#define rMessage(message)                       \
    {cout << message << " :: line " << __LINE__ \
          << " in " << __FILE__ << endl; }

#define rError(message)                         \
    {cout << message << " :: line " << __LINE__ \
          << " in " << __FILE__ << endl;        \
        exit(false);}

/* ============================================================
   Allocate array
   ============================================================ */
#if 1
#define NewArray(val,type,number) \
  {val = NULL; \
    try{ val = new type[number]; } \
    catch(bad_alloc){ \
        rMessage("Memory Exhausted (bad_alloc)"); abort(); } \
    catch(...){ \
        rMessage("Fatal Error (related memory allocation"); abort(); } \
  }
#else
#define NewArray(val,type,number) \
  {rMessage("New Invoked"); \
   val = NULL; val = new type[number]; \
   if  (val==NULL) {rError("Over Memory");} \
  }
#endif

#define DeleteArray(val) \
  { if  (val!=NULL) { \
      delete[] val; \
      val = NULL; \
    } \
  }

#endif

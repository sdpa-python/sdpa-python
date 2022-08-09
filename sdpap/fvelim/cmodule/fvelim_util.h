/*
C Extension code to eliminate free variables
This file is a component of SDPAP
Copyright (C) 2010-2022 SDPA Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

December 2010: Originally written by Kenta Kato
*/

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

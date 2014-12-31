/*=============================================================================
#     FileName: config.h
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2013-02-19 19:11:15
#   LastChange: 2014-02-11 18:08:53
#      History:
=============================================================================*/

#ifndef  CONFIG_H
#define  CONFIG_H

#define VERSION "calcDescriptors.exe version 1.0"

#if !defined(PI)
  #define PI 3.1415926535897931
#endif

#if !defined(SEP)
  #if defined(WIN32)
    #define SEP "\\"
  #else
    #define SEP "/"
  #endif
#endif

// in order to build dynamic library, you need to define 'USING_DLL' and 'DESCRIPTORS_DLL_EXPORTS'
// when using the dynamic library, please define 'USING_DLL'
#if defined(USING_DLL)
 #if defined(DESCRIPTORS_DLL_EXPORTS)
  #define DESCRIPTOR_API __declspec(dllexport)
 #else
  #define DESCRIPTOR_API __declspec(dllimport)
 #endif
#else
 #define DESCRIPTOR_API
#endif

#endif   /* ----- #ifndef CONFIG_H  ----- */


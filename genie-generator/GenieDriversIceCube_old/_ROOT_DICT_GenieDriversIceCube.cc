// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME _ROOT_DICT_GenieDriversIceCube

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "GIceCubeDiffuseFlux.h"
#include "SimpleIceCubeGeomAnalyzer.h"

// Header files passed via #pragma extra_include

namespace genie {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *genie_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("genie", 0 /*version*/, "EVGDrivers/GFluxI.h", 33,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &genie_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *genie_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace genie {
   namespace flux {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *geniecLcLflux_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("genie::flux", 0 /*version*/, "GIceCubeDiffuseFlux.h", 47,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &geniecLcLflux_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *geniecLcLflux_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}
}

namespace genie {
   namespace geometry {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *geniecLcLgeometry_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("genie::geometry", 0 /*version*/, "SimpleIceCubeGeomAnalyzer.h", 27,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &geniecLcLgeometry_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *geniecLcLgeometry_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}
}

namespace ROOT {
   static TClass *geniecLcLfluxcLcLGIceCubeDiffuseFlux_Dictionary();
   static void geniecLcLfluxcLcLGIceCubeDiffuseFlux_TClassManip(TClass*);
   static void delete_geniecLcLfluxcLcLGIceCubeDiffuseFlux(void *p);
   static void deleteArray_geniecLcLfluxcLcLGIceCubeDiffuseFlux(void *p);
   static void destruct_geniecLcLfluxcLcLGIceCubeDiffuseFlux(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::genie::flux::GIceCubeDiffuseFlux*)
   {
      ::genie::flux::GIceCubeDiffuseFlux *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::genie::flux::GIceCubeDiffuseFlux));
      static ::ROOT::TGenericClassInfo 
         instance("genie::flux::GIceCubeDiffuseFlux", "GIceCubeDiffuseFlux.h", 49,
                  typeid(::genie::flux::GIceCubeDiffuseFlux), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &geniecLcLfluxcLcLGIceCubeDiffuseFlux_Dictionary, isa_proxy, 0,
                  sizeof(::genie::flux::GIceCubeDiffuseFlux) );
      instance.SetDelete(&delete_geniecLcLfluxcLcLGIceCubeDiffuseFlux);
      instance.SetDeleteArray(&deleteArray_geniecLcLfluxcLcLGIceCubeDiffuseFlux);
      instance.SetDestructor(&destruct_geniecLcLfluxcLcLGIceCubeDiffuseFlux);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::genie::flux::GIceCubeDiffuseFlux*)
   {
      return GenerateInitInstanceLocal((::genie::flux::GIceCubeDiffuseFlux*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::genie::flux::GIceCubeDiffuseFlux*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *geniecLcLfluxcLcLGIceCubeDiffuseFlux_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::genie::flux::GIceCubeDiffuseFlux*)0x0)->GetClass();
      geniecLcLfluxcLcLGIceCubeDiffuseFlux_TClassManip(theClass);
   return theClass;
   }

   static void geniecLcLfluxcLcLGIceCubeDiffuseFlux_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer_Dictionary();
   static void geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer_TClassManip(TClass*);
   static void delete_geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer(void *p);
   static void deleteArray_geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer(void *p);
   static void destruct_geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::genie::geometry::SimpleIceCubeGeomAnalyzer*)
   {
      ::genie::geometry::SimpleIceCubeGeomAnalyzer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::genie::geometry::SimpleIceCubeGeomAnalyzer));
      static ::ROOT::TGenericClassInfo 
         instance("genie::geometry::SimpleIceCubeGeomAnalyzer", "SimpleIceCubeGeomAnalyzer.h", 29,
                  typeid(::genie::geometry::SimpleIceCubeGeomAnalyzer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer_Dictionary, isa_proxy, 0,
                  sizeof(::genie::geometry::SimpleIceCubeGeomAnalyzer) );
      instance.SetDelete(&delete_geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer);
      instance.SetDeleteArray(&deleteArray_geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer);
      instance.SetDestructor(&destruct_geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::genie::geometry::SimpleIceCubeGeomAnalyzer*)
   {
      return GenerateInitInstanceLocal((::genie::geometry::SimpleIceCubeGeomAnalyzer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::genie::geometry::SimpleIceCubeGeomAnalyzer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::genie::geometry::SimpleIceCubeGeomAnalyzer*)0x0)->GetClass();
      geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer_TClassManip(theClass);
   return theClass;
   }

   static void geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_geniecLcLfluxcLcLGIceCubeDiffuseFlux(void *p) {
      delete ((::genie::flux::GIceCubeDiffuseFlux*)p);
   }
   static void deleteArray_geniecLcLfluxcLcLGIceCubeDiffuseFlux(void *p) {
      delete [] ((::genie::flux::GIceCubeDiffuseFlux*)p);
   }
   static void destruct_geniecLcLfluxcLcLGIceCubeDiffuseFlux(void *p) {
      typedef ::genie::flux::GIceCubeDiffuseFlux current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::genie::flux::GIceCubeDiffuseFlux

namespace ROOT {
   // Wrapper around operator delete
   static void delete_geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer(void *p) {
      delete ((::genie::geometry::SimpleIceCubeGeomAnalyzer*)p);
   }
   static void deleteArray_geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer(void *p) {
      delete [] ((::genie::geometry::SimpleIceCubeGeomAnalyzer*)p);
   }
   static void destruct_geniecLcLgeometrycLcLSimpleIceCubeGeomAnalyzer(void *p) {
      typedef ::genie::geometry::SimpleIceCubeGeomAnalyzer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::genie::geometry::SimpleIceCubeGeomAnalyzer

namespace {
  void TriggerDictionaryInitialization__ROOT_DICT_GenieDriversIceCube_Impl() {
    static const char* headers[] = {
"GIceCubeDiffuseFlux.h",
"SimpleIceCubeGeomAnalyzer.h",
0
    };
    static const char* includePaths[] = {
"/usr/include/libxml2",
"/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_7_x86_64/include/",
"/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/RHEL_7_x86_64/include",
"/srv/work/simulation/GENIE/include",
"/home/hignight/work/simulation/GENIE/R-2_12_8/src/",
"/home/hignight/work/simulation/GENIE/R-2_12_8/include/GENIE/",
"/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/RHEL_7_x86_64/include",
"/home/mliubar/Software/genie_workspace/genie-generator/GenieDriversIceCube/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "_ROOT_DICT_GenieDriversIceCube dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace genie{namespace flux{class __attribute__((annotate("$clingAutoload$GIceCubeDiffuseFlux.h")))  GIceCubeDiffuseFlux;}}
namespace genie{namespace geometry{class __attribute__((annotate("$clingAutoload$SimpleIceCubeGeomAnalyzer.h")))  SimpleIceCubeGeomAnalyzer;}}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "_ROOT_DICT_GenieDriversIceCube dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "GIceCubeDiffuseFlux.h"
#include "SimpleIceCubeGeomAnalyzer.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"genie::flux::GIceCubeDiffuseFlux", payloadCode, "@",
"genie::geometry::SimpleIceCubeGeomAnalyzer", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("_ROOT_DICT_GenieDriversIceCube",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization__ROOT_DICT_GenieDriversIceCube_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization__ROOT_DICT_GenieDriversIceCube_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization__ROOT_DICT_GenieDriversIceCube() {
  TriggerDictionaryInitialization__ROOT_DICT_GenieDriversIceCube_Impl();
}

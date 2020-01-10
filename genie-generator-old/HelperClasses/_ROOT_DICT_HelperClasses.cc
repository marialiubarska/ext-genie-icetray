// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME _ROOT_DICT_HelperClasses

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
#include "CrossSectionAccessor.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *CrossSectionAccessor_Dictionary();
   static void CrossSectionAccessor_TClassManip(TClass*);
   static void delete_CrossSectionAccessor(void *p);
   static void deleteArray_CrossSectionAccessor(void *p);
   static void destruct_CrossSectionAccessor(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CrossSectionAccessor*)
   {
      ::CrossSectionAccessor *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CrossSectionAccessor));
      static ::ROOT::TGenericClassInfo 
         instance("CrossSectionAccessor", "CrossSectionAccessor.h", 20,
                  typeid(::CrossSectionAccessor), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CrossSectionAccessor_Dictionary, isa_proxy, 0,
                  sizeof(::CrossSectionAccessor) );
      instance.SetDelete(&delete_CrossSectionAccessor);
      instance.SetDeleteArray(&deleteArray_CrossSectionAccessor);
      instance.SetDestructor(&destruct_CrossSectionAccessor);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CrossSectionAccessor*)
   {
      return GenerateInitInstanceLocal((::CrossSectionAccessor*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::CrossSectionAccessor*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CrossSectionAccessor_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::CrossSectionAccessor*)0x0)->GetClass();
      CrossSectionAccessor_TClassManip(theClass);
   return theClass;
   }

   static void CrossSectionAccessor_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_CrossSectionAccessor(void *p) {
      delete ((::CrossSectionAccessor*)p);
   }
   static void deleteArray_CrossSectionAccessor(void *p) {
      delete [] ((::CrossSectionAccessor*)p);
   }
   static void destruct_CrossSectionAccessor(void *p) {
      typedef ::CrossSectionAccessor current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::CrossSectionAccessor

namespace {
  void TriggerDictionaryInitialization__ROOT_DICT_HelperClasses_Impl() {
    static const char* headers[] = {
"CrossSectionAccessor.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/RHEL_7_x86_64/include",
"/home/hignight/work/simulation/GENIE/R-2_12_8/src",
"/home/hignight/work/simulation/GENIE/R-2_12_8/src/NuValidator",
"./",
"/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/RHEL_7_x86_64/include",
"/home/mliubar/Software/genie_workspace/genie-generator/HelperClasses/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "_ROOT_DICT_HelperClasses dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$CrossSectionAccessor.h")))  CrossSectionAccessor;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "_ROOT_DICT_HelperClasses dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "CrossSectionAccessor.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"CrossSectionAccessor", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("_ROOT_DICT_HelperClasses",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization__ROOT_DICT_HelperClasses_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization__ROOT_DICT_HelperClasses_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization__ROOT_DICT_HelperClasses() {
  TriggerDictionaryInitialization__ROOT_DICT_HelperClasses_Impl();
}

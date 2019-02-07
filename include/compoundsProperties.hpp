#ifndef COMPOUNDS_PROPERTIES
#define COMPOUNDS_PROPERTIES

#include "proteinManager.hpp"

#include <boost/preprocessor.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/tuple/elem.hpp>

namespace proteinManager {

#define MDL_PROPERTIES_LIST   ((Id    , id    ,int        ))

#define CHAIN_PROPERTIES_LIST ((Id    , id    ,std::string))

#define RES_PROPERTIES_LIST   ((Name       , name       ,std::string)) \
                              ((Seq        , seq        ,int        )) \
                              ((InsCode    , iCode      ,std::string))

#define ATOM_PROPERTIES_LIST  ((Serial     , serial     ,int        )) \
                              ((Name       , name       ,std::string)) \
                              ((AltLoc     , altLoc     ,std::string)) \
                              ((Coord      , coord      ,real3      )) \
                              ((Occupancy  , occupancy  ,real       )) \
                              ((TempFactor , tempFactor ,real       )) \
                              ((Element    , element    ,std::string)) \
                              ((Mass       , mass       ,real       )) \
                              ((Charge     , charge     ,real       )) \
                              ((Radius     , radius     ,real       )) \
                              ((C12        , c12        ,real       )) \
                              ((C6         , c6         ,real       )) \
                              ((SolvE      , solvE      ,real       )) \
                              ((SASA       , SASA       ,real       )) 
                              


//Get the Name (first letter capital) from a tuple in the property list
#define PROPNAME_CAPS(tuple) BOOST_PP_TUPLE_ELEM(3, 0 ,tuple)
//Get the name (no capital) from a tuple in the property list
#define PROPNAME(tuple) BOOST_PP_TUPLE_ELEM(3, 1 ,tuple)
//Get the type from a tuple in the property list
#define PROPTYPE(tuple) BOOST_PP_TUPLE_ELEM(3, 2 ,tuple)

//This macro iterates through all properties applying some macro
#define MDL_PROPERTY_LOOP(macro)    BOOST_PP_SEQ_FOR_EACH(macro, _, MDL_PROPERTIES_LIST)
#define CHAIN_PROPERTY_LOOP(macro)  BOOST_PP_SEQ_FOR_EACH(macro, _, CHAIN_PROPERTIES_LIST)
#define RES_PROPERTY_LOOP(macro)    BOOST_PP_SEQ_FOR_EACH(macro, _, RES_PROPERTIES_LIST)
#define ATOM_PROPERTY_LOOP(macro)   BOOST_PP_SEQ_FOR_EACH(macro, _, ATOM_PROPERTIES_LIST)
}

#endif

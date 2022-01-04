// This file has been generated by Py++.

#include "boost/python.hpp"
#include "wrappable_v3d.h"
#include "QPolygon.pypp.hpp"

namespace bp = boost::python;

void register_QPolygon_class(){

    { //::QPolygon
        typedef bp::class_< QPolygon, bp::bases< QVector< QPoint > > > QPolygon_exposer_t;
        QPolygon_exposer_t QPolygon_exposer = QPolygon_exposer_t( "QPolygon", bp::init< >() );
        bp::scope QPolygon_scope( QPolygon_exposer );
        QPolygon_exposer.def( bp::init< QPolygon const & >(( bp::arg("a") )) );
        QPolygon_exposer.def( bp::init< QVector< QPoint > const & >(( bp::arg("v") )) );
        bp::implicitly_convertible< QVector< QPoint > const &, QPolygon >();
        QPolygon_exposer.def( bp::init< QRect const &, bp::optional< bool > >(( bp::arg("r"), bp::arg("closed")=(bool)(false) )) );
        bp::implicitly_convertible< QRect const &, QPolygon >();
        QPolygon_exposer.def( bp::init< int, int const * >(( bp::arg("nPoints"), bp::arg("points") )) );
        QPolygon_exposer.def( bp::init< int >(( bp::arg("asize") )) );
        bp::implicitly_convertible< int, QPolygon >();
        { //::QPolygon::boundingRect
        
            typedef ::QRect ( ::QPolygon::*boundingRect_function_type )(  ) const;
            
            QPolygon_exposer.def( 
                "boundingRect"
                , boundingRect_function_type( &::QPolygon::boundingRect ) );
        
        }
        { //::QPolygon::containsPoint
        
            typedef bool ( ::QPolygon::*containsPoint_function_type )( ::QPoint const &,::Qt::FillRule ) const;
            
            QPolygon_exposer.def( 
                "containsPoint"
                , containsPoint_function_type( &::QPolygon::containsPoint )
                , ( bp::arg("pt"), bp::arg("fillRule") ) );
        
        }
        { //::QPolygon::intersected
        
            typedef ::QPolygon ( ::QPolygon::*intersected_function_type )( ::QPolygon const & ) const;
            
            QPolygon_exposer.def( 
                "intersected"
                , intersected_function_type( &::QPolygon::intersected )
                , ( bp::arg("r") ) );
        
        }
        { //::QPolygon::point
        
            typedef void ( ::QPolygon::*point_function_type )( int,int *,int * ) const;
            
            QPolygon_exposer.def( 
                "point"
                , point_function_type( &::QPolygon::point )
                , ( bp::arg("i"), bp::arg("x"), bp::arg("y") ) );
        
        }
        { //::QPolygon::point
        
            typedef ::QPoint ( ::QPolygon::*point_function_type )( int ) const;
            
            QPolygon_exposer.def( 
                "point"
                , point_function_type( &::QPolygon::point )
                , ( bp::arg("index") ) );
        
        }
        { //::QPolygon::putPoints
        
            typedef void ( ::QPolygon::*putPoints_function_type )( int,int,::QPolygon const &,int ) ;
            
            QPolygon_exposer.def( 
                "putPoints"
                , putPoints_function_type( &::QPolygon::putPoints )
                , ( bp::arg("index"), bp::arg("nPoints"), bp::arg("from"), bp::arg("fromIndex")=(int)(0) ) );
        
        }
        { //::QPolygon::setPoint
        
            typedef void ( ::QPolygon::*setPoint_function_type )( int,int,int ) ;
            
            QPolygon_exposer.def( 
                "setPoint"
                , setPoint_function_type( &::QPolygon::setPoint )
                , ( bp::arg("index"), bp::arg("x"), bp::arg("y") ) );
        
        }
        { //::QPolygon::setPoint
        
            typedef void ( ::QPolygon::*setPoint_function_type )( int,::QPoint const & ) ;
            
            QPolygon_exposer.def( 
                "setPoint"
                , setPoint_function_type( &::QPolygon::setPoint )
                , ( bp::arg("index"), bp::arg("pt") )
                , "\n  Misc. QPolygon functions\n" );
        
        }
        { //::QPolygon::subtracted
        
            typedef ::QPolygon ( ::QPolygon::*subtracted_function_type )( ::QPolygon const & ) const;
            
            QPolygon_exposer.def( 
                "subtracted"
                , subtracted_function_type( &::QPolygon::subtracted )
                , ( bp::arg("r") ) );
        
        }
        { //::QPolygon::translate
        
            typedef void ( ::QPolygon::*translate_function_type )( int,int ) ;
            
            QPolygon_exposer.def( 
                "translate"
                , translate_function_type( &::QPolygon::translate )
                , ( bp::arg("dx"), bp::arg("dy") ) );
        
        }
        { //::QPolygon::translate
        
            typedef void ( ::QPolygon::*translate_function_type )( ::QPoint const & ) ;
            
            QPolygon_exposer.def( 
                "translate"
                , translate_function_type( &::QPolygon::translate )
                , ( bp::arg("offset") ) );
        
        }
        { //::QPolygon::translated
        
            typedef ::QPolygon ( ::QPolygon::*translated_function_type )( int,int ) const;
            
            QPolygon_exposer.def( 
                "translated"
                , translated_function_type( &::QPolygon::translated )
                , ( bp::arg("dx"), bp::arg("dy") ) );
        
        }
        { //::QPolygon::translated
        
            typedef ::QPolygon ( ::QPolygon::*translated_function_type )( ::QPoint const & ) const;
            
            QPolygon_exposer.def( 
                "translated"
                , translated_function_type( &::QPolygon::translated )
                , ( bp::arg("offset") ) );
        
        }
        { //::QPolygon::united
        
            typedef ::QPolygon ( ::QPolygon::*united_function_type )( ::QPolygon const & ) const;
            
            QPolygon_exposer.def( 
                "united"
                , united_function_type( &::QPolygon::united )
                , ( bp::arg("r") ) );
        
        }
    }

}

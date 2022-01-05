// This file has been generated by Py++.

#include "boost/python.hpp"
#include "wrappable_v3d.h"
#include "QVector_QPoint.pypp.hpp"

namespace bp = boost::python;

void register_QVector_QPoint_class(){

    { //::QVector< QPoint >
        typedef bp::class_< QVector< QPoint >, boost::noncopyable > QVector_QPoint_exposer_t;
        QVector_QPoint_exposer_t QVector_QPoint_exposer = QVector_QPoint_exposer_t( "QVector_QPoint", bp::no_init );
        bp::scope QVector_QPoint_scope( QVector_QPoint_exposer );
        QVector_QPoint_exposer.def( bp::init< >() );
        QVector_QPoint_exposer.def( bp::init< QVector< QPoint > const & >(( bp::arg("v") )) );
        QVector_QPoint_exposer.def( bp::init< int >(( bp::arg("asize") )) );
        bp::implicitly_convertible< int, QVector< QPoint > >();
        QVector_QPoint_exposer.def( bp::init< int, QPoint const & >(( bp::arg("asize"), bp::arg("t") )) );
    }

}

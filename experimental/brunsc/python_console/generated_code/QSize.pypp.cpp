// This file has been generated by Py++.

#include "boost/python.hpp"
#include "wrappable_v3d.h"
#include "QSize.pypp.hpp"

namespace bp = boost::python;

void register_QSize_class(){

    bp::class_< QSize >( "QSize", "\n Copyright (C) 2011 Nokia Corporation and/or its subsidiary(-ies).\n All rights reserved.\n Contact: Nokia Corporation (qt-info@nokia.com)\n\n This file is part of the QtCore module of the Qt Toolkit.\n\n $QT_BEGIN_LICENSE:LGPL$\n Commercial Usage\n Licensees holding valid Qt Commercial licenses may use this file in\n accordance with the Qt Commercial License Agreement provided with the\n Software or, alternatively, in accordance with the terms contained in\n a written agreement between you and Nokia.\n\n GNU Lesser General Public License Usage\n Alternatively, this file may be used under the terms of the GNU Lesser\n General Public License version 2.1 as published by the Free Software\n Foundation and appearing in the file LICENSE.LGPL included in the\n packaging of this file.  Please review the following information to\n ensure the GNU Lesser General Public License version 2.1 requirements\n will be met: http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.\n\n In addition, as a special exception, Nokia gives you certain additional\n rights.  These rights are described in the Nokia Qt LGPL Exception\n version 1.1, included in the file LGPL_EXCEPTION.txt in this package.\n\n GNU General Public License Usage\n Alternatively, this file may be used under the terms of the GNU\n General Public License version 3.0 as published by the Free Software\n Foundation and appearing in the file LICENSE.GPL included in the\n packaging of this file.  Please review the following information to\n ensure the GNU General Public License version 3.0 requirements will be\n met: http://www.gnu.org/copyleft/gpl.html.\n\n If you have questions regarding the use of this file, please contact\n Nokia at qt-info@nokia.com.\n $QT_END_LICENSE$\n", bp::init< >("\n  QSize inline functions\n") )    
        .def( bp::init< int, int >(( bp::arg("w"), bp::arg("h") )) )    
        .def( 
            "boundedTo"
            , (::QSize ( ::QSize::* )( ::QSize const & ) const)( &::QSize::boundedTo )
            , ( bp::arg("otherSize") ) )    
        .def( 
            "expandedTo"
            , (::QSize ( ::QSize::* )( ::QSize const & ) const)( &::QSize::expandedTo )
            , ( bp::arg("otherSize") ) )    
        .def( 
            "height"
            , (int ( ::QSize::* )(  ) const)( &::QSize::height ) )    
        .def( 
            "isEmpty"
            , (bool ( ::QSize::* )(  ) const)( &::QSize::isEmpty ) )    
        .def( 
            "isNull"
            , (bool ( ::QSize::* )(  ) const)( &::QSize::isNull ) )    
        .def( 
            "isValid"
            , (bool ( ::QSize::* )(  ) const)( &::QSize::isValid ) )    
        .def( bp::self *= bp::other< qreal >() )    
        .def( bp::self += bp::self )    
        .def( bp::self -= bp::self )    
        .def( bp::self /= bp::other< qreal >() )    
        .def( 
            "rheight"
            , (int & ( ::QSize::* )(  ) )( &::QSize::rheight )
            , bp::return_value_policy< bp::copy_non_const_reference >() )    
        .def( 
            "rwidth"
            , (int & ( ::QSize::* )(  ) )( &::QSize::rwidth )
            , bp::return_value_policy< bp::copy_non_const_reference >() )    
        .def( 
            "scale"
            , (void ( ::QSize::* )( int,int,::Qt::AspectRatioMode ) )( &::QSize::scale )
            , ( bp::arg("w"), bp::arg("h"), bp::arg("mode") ) )    
        .def( 
            "scale"
            , (void ( ::QSize::* )( ::QSize const &,::Qt::AspectRatioMode ) )( &::QSize::scale )
            , ( bp::arg("s"), bp::arg("mode") ) )    
        .def( 
            "setHeight"
            , (void ( ::QSize::* )( int ) )( &::QSize::setHeight )
            , ( bp::arg("h") ) )    
        .def( 
            "setWidth"
            , (void ( ::QSize::* )( int ) )( &::QSize::setWidth )
            , ( bp::arg("w") ) )    
        .def( 
            "transpose"
            , (void ( ::QSize::* )(  ) )( &::QSize::transpose ) )    
        .def( 
            "width"
            , (int ( ::QSize::* )(  ) const)( &::QSize::width ) );

}

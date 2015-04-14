# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SwiftWorldFiler
                                 A QGIS plugin
 Updates world files dynamically to have a real time georeferencing
                              -------------------
        begin                : 2015-04-08
        git sha              : $Format:%H$
        copyright            : (C) 2015 by Olivier Dalang
        email                : olivier.dalang@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

from PyQt4.QtCore import *
from PyQt4.QtGui import *

from qgis.core import *

import os.path
import math
import numpy
import random


class SwiftComputer:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        self.iface = iface

        self.pointsPairs = [] # contains correspondance pairs (pixels/world) of QgsPoints

        self.a = 0
        self.b = 0
        self.c = 0
        self.d = 0
        self.e = 0
        self.f = 0

    ########################################################
    # Load / Save pointfiles
    ########################################################

    def loadPointfile(self, path):
        QgsMessageLog.logMessage("Loading point file")
        self.debug()
        try:
            f = open(path, 'r+')
            lines = f.readlines()
            QgsMessageLog.logMessage("lines read : "+str(lines))
            f.close()
        except Exception, e:
            lines = []

        self.pointsPairs = []

        for line in lines[1:]:

            s = line.split(',')

            pixelP = QgsPoint(float(s[0]),float(s[1]))
            worldP = QgsPoint(float(s[2]),float(s[3]))

            self.pointsPairs.append( [pixelP,worldP] )



    def savePointfile(self, path):
        f = open(path, 'w+')
        f.write( 'pixel_x,pixel_y,geo_x,geo_y\n' )
        
        for pixelP,worldP in self.pointsPairs:
            f.write( '{},{},{},{}\n'.format(pixelP.x(),pixelP.y(),worldP.x(),worldP.y()) )
        f.close()


    ########################################################
    # Load / Generate memory ladyer
    ########################################################

    def loadMemoryLayer(self, layer):
        self.pointsPairs = []

        for feature in layer.getFeatures():
            geom = feature.geometry().asPolyline()

            pixelP = self.utmToPixel(geom[0])
            worldP = geom[1]

            self.pointsPairs.append( [pixelP,worldP] )

    def regenerateMemoryLayer(self, memoryLayer):

        provider = memoryLayer.dataProvider()

        # delete all features
        ids = [ f.id() for f in provider.getFeatures()]
        provider.deleteFeatures( ids )

        # add all features

        features = []
        for pixelP,worldP in self.pointsPairs:

            feature = QgsFeature()
            feature.setGeometry(QgsGeometry.fromPolyline([self.pixelToUtm(pixelP),worldP]))
            features.append( feature )
        provider.addFeatures(features)
        


    ########################################################
    # Generate worldfile
    ########################################################

    def generateWorldfile(self, path):

        f = open(path, 'w+')
        f.write( str( self.a )+"\n" )
        f.write( str( self.d )+"\n" )
        f.write( str( self.b )+"\n" )
        f.write( str( self.e )+"\n" )
        f.write( str( self.c )+"\n" )
        f.write( str( self.f )+"\n" )
        f.close()



    def updateTransform(self, extent, rasterWidth, rasterHeight):
        # E : pixel size in the y direction (usually negative)

        useMostCentralThreePoints = True # set to false to use leastSquareTentative

        if len(self.pointsPairs) == 0:
            # We compute the values so that the raster fits in the canvas

            # A : pixel size in X direction
            self.a = extent.width() / rasterWidth

            # B : rotation about X axis
            self.b = 0

            # C : x coord of the upper left pixel
            self.c = extent.xMinimum()

            # D : rotation about Z axis
            self.d = 0

            # E : pixel size in the y direction (usually negative)
            #self.e = -extent.height() / rasterHeight
            self.e = -extent.width() / rasterWidth

            # F : y coord of the upper left pixel
            self.f = extent.yMaximum()

        elif len(self.pointsPairs) == 1:
            # We compute the value so that the size reflect's the canvas size, but we pin the point

            p1 = self.pointsPairs[0]
            p1Pix = p1[0]
            p1Wld = p1[1]

            # A : pixel size in X direction
            self.a = extent.width() / rasterWidth

            # B : rotation about X axis
            self.b = 0

            # C : x coord of the upper left pixel
            self.c = p1Wld.x()-self.a*p1Pix.x()

            # D : rotation about Z axis
            self.d = 0

            # E : pixel size in the y direction (usually negative)
            #self.e = -extent.height() / rasterHeight
            self.e = -self.a

            # F : y coord of the upper left pixel
            self.f = p1Wld.y()-self.e*p1Pix.y()

        elif len(self.pointsPairs) == 2:
            # We compute the values taking the hypothesis pixels are square (A=-E and B=D), so we solve a four variables equation system

            # Input points
            p1 = self.pointsPairs[0]
            p1PixX = p1[0].x()
            p1PixY = p1[0].y()
            p1WldX = p1[1].x()
            p1WldY = p1[1].y()

            p2 = self.pointsPairs[1]
            p2PixX = p2[0].x()
            p2PixY = p2[0].y()
            p2WldX = p2[1].x()
            p2WldY = p2[1].y()

            # We know A=-E and B=D, so we can link both variable by going to the (C;F) point (rotation origin)
            #1) p1WldX - ( A*p1PixX+B*p1PixY ) = p2WldX - ( A*p2PixX+B*p2PixY ) = C
            #2) p1WldY - ( B*p1PixX+(-A)*p1PixY ) = p2WldY - ( B*p2PixX+(-A)*p2PixY ) = F

            # We replace varible names to solve using http://www.numberempire.com/equationsolver.php:
            i,j,k,l = p1[0].x(),p1[0].y(),p1[1].x(),p1[1].y()
            m,n,o,p = p2[0].x(),p2[0].y(),p2[1].x(),p2[1].y()

            #1) (C=)  k - ( A*i+B*j ) = o - ( A*m+B*n )
            #2) (F=)  l - ( B*i+(-A)*j ) = p - ( B*m+(-A)*n )

            # A : pixel size in X direction
            self.a = (j*(p-l)+n*(l-p)+(m-i)*o-k*m+i*k)/(n**2-2*j*n+m**2-2*i*m+j**2+i**2)
            # B : rotation about X axis
            self.b = (m*(p-l)+i*(l-p)+(n-j)*o-k*n+j*k)/(n**2-2*j*n+m**2-2*i*m+j**2+i**2)
            # C : x coord of the upper left pixel
            self.c = k - ( self.a*i+self.b*j )
            # D : rotation about Z axis
            self.d = self.b
            # E : pixel size in the y direction (usually negative)
            self.e = -self.a
            # F : y coord of the upper left pixel
            self.f = l - ( self.b*i+(-self.a)*j )


        elif len(self.pointsPairs) == 3:
            # We solve a six variables equation system

            p1 = self.pointsPairs[0]
            p2 = self.pointsPairs[1]
            p3 = self.pointsPairs[2]

            # Input points
            p1PixX = p1[0].x()
            p1PixY = p1[0].y()
            p1WldX = p1[1].x()
            p1WldY = p1[1].y()

            p2PixX = p2[0].x()
            p2PixY = p2[0].y()
            p2WldX = p2[1].x()
            p2WldY = p2[1].y()

            p3PixX = p3[0].x()
            p3PixY = p3[0].y()
            p3WldX = p3[1].x()
            p3WldY = p3[1].y()


            # We can link those variables by going to the (C;F) point (rotation origin)
            #1) C = p1WldX - ( A*p1PixX + B*p1PixY )
            #2) F = p1WldY - ( D*p1PixX + E*p1PixY )
            #3) C = p2WldX - ( A*p2PixX + B*p2PixY )
            #4) F = p2WldY - ( D*p2PixX + E*p2PixY )
            #5) C = p3WldX - ( A*p3PixX + B*p3PixY )
            #6) F = p3WldY - ( D*p3PixX + E*p3PixY )

            #7) C = p4WldX - ( A*p4PixX + B*p4PixY )
            #8) F = p4WldY - ( D*p4PixX + E*p4PixY )


            # We replace varible names to solve using http://www.numberempire.com/equationsolver.php:
            i,j,k,l = p1[0].x(),p1[0].y(),p1[1].x(),p1[1].y()
            m,n,o,p = p2[0].x(),p2[0].y(),p2[1].x(),p2[1].y()
            q,r,s,t = p3[0].x(),p3[0].y(),p3[1].x(),p3[1].y()

            #1) C = k - ( A*i + B*j )
            #2) F = l - ( D*i + E*j )
            #3) C = o - ( A*m + B*n )
            #4) F = p - ( D*m + E*n )
            #5) C = s - ( A*q + B*r )
            #6) F = t - ( D*q + E*r )

            # A : pixel size in X direction
            self.a = -(j*(s-o)-n*s+o*r+k*(n-r))/(i*(r-n)-m*r+n*q+j*(m-q))
            # B : rotation about X axis
            self.b = (i*(s-o)-m*s+o*q+k*(m-q))/(i*(r-n)-m*r+n*q+j*(m-q))
            # C : x coord of the upper left pixel
            self.c = (i*(o*r-n*s)+j*(m*s-o*q)+k*(n*q-m*r))/(i*(r-n)-m*r+n*q+j*(m-q))
            # D : rotation about Z axis
            self.d = (n*t+j*(p-t)+l*(r-n)-p*r)/(i*(r-n)-m*r+n*q+j*(m-q))
            # E : pixel size in the y direction (usually negative)
            self.e = -(m*t+i*(p-t)+l*(q-m)-p*q)/(i*(r-n)-m*r+n*q+j*(m-q))
            # F : y coord of the upper left pixel
            self.f = (i*(p*r-n*t)+j*(m*t-p*q)+l*(n*q-m*r))/(i*(r-n)-m*r+n*q+j*(m-q)) 


        else:
            # We do a least square fit


            # Definition of one point
            # pWldX1 = A*pPixX1 + B*pPixY1 + C 
            # pWldY1 = D*pPixX1 + E*pPixY1 + F

            # Definition of one point with errors
            # pWldXn = A*pPixXn + B*pPixYn + C + eXn
            # pWldYn = D*pPixXn + E*pPixYn + F + eYn

            # Definition of one point's error
            # eXn = (pWldXn - A*pPixXn - B*pPixYn - C)
            # eYn = (pWldYn - D*pPixXn - E*pPixYn - F)

            # We want the sum of errors squared to converge towards 0

            #1) DER_A( SUM( eXn**2 )) = 0
            #2) DER_B( SUM( eXn**2 )) = 0
            #3) DER_C( SUM( eXn**2 )) = 0

            #4) DER_D( SUM( eYn**2 )) = 0
            #5) DER_E( SUM( eYn**2 )) = 0
            #6) DER_F( SUM( eYn**2 )) = 0


            # Terms :
            # o = pWldXn
            # p = pWldYn
            # x = pPixXn
            # y = pPixYn

            ## eXn =
            # (o - A*x - B*y - C)
            ## eYn =
            # (p - D*x - E*y - F)

            ## eXn^2 =
            # C^2  +  2*y*B*C  +  2*x*A*C  -  2*o*C  +  y^2*B^2  +  2*x*y*A*B  -  2*o*y*B  +  x^2*A^2  -  2*o*x*A  +  o^2
            ## eYn^2 =
            # F^2  +  2*y*E*F  +  2*x*D*F  -  2*p*F  +  y^2*E^2  +  2*x*y*D*E  -  2*p*y*E  +  x^2*D^2  -  2*p*x*D  +  p^2

            # Sums :
            # o = SUM(pWldXn)
            # p = SUM(pWldYn)
            # x = SUM(pPixXn)
            # y = SUM(pPixYn)
            # u = SUM(pPixXn^2) # x^2
            # v = SUM(pPixYn^2) # y^2
            # w = SUM(pPixXn*pPixYn) # x*y
            # j = SUM(pWldXn*pPixXn) # o*x
            # k = SUM(pWldXn*pPixYn) # o*y
            # l = SUM(pWldYn*pPixXn) # p*x
            # m = SUM(pWldyn*pPixYn) # p*y
            # n = SUM(1)


            ## SUM(eXn^2) =
            # n*C^2  +  2*y*B*C  +  2*x*A*C  -  2*o*C  +  v*B^2  +  2*w*A*B  -  2*k*B  +  u*A^2  -  2*j*A  +  n*o^2
            ## SUM(eYn^2) =
            # n*F^2  +  2*y*E*F  +  2*x*D*F  -  2*p*F  +  v*E^2  +  2*w*D*E  -  2*m*E  +  u*D^2  -  2*l*D  +  n*p^2


            ## DERIVATIVE FOR A OF SUM(eXn^2)
            ## DERIVATIVE FOR B OF SUM(eXn^2)
            ## DERIVATIVE FOR C OF SUM(eXn^2)
            ## DERIVATIVE FOR D OF SUM(eYn^2)
            ## DERIVATIVE FOR E OF SUM(eYn^2)
            ## DERIVATIVE FOR F OF SUM(eYn^2)

            #1) 2*x*C+2*w*B+2*u*A-2*j = 0
            #2) 2*y*C+2*v*B+2*w*A-2*k = 0
            #3) 2*n*C+2*y*B+2*x*A-2*o = 0

            #4) 2*x*F+2*w*E+2*u*D-2*l = 0
            #5) 2*y*F+2*v*E+2*w*D-2*m = 0
            #6) 2*n*F+2*y*E+2*x*D-2*p = 0

            #System1 (A;B;C) : 2*x*C+2*w*B+2*u*A-2*j = 0, 2*y*C+2*v*B+2*w*A-2*k = 0, 2*n*C+2*y*B+2*x*A-2*o = 0
            #System2 (D;E;F) : 2*x*F+2*w*E+2*u*D-2*l = 0, 2*y*F+2*v*E+2*w*D-2*m = 0, 2*n*F+2*y*E+2*x*D-2*p = 0

            # Results
            # A = -(j*(n*v-y^2)+w*(o*y-k*n)+x*(k*y-o*v))/(u*(y^2-n*v)-2*w*x*y+v*x^2+n*w^2)
            # B = (u*(o*y-k*n)+x*(-j*y-o*w)+k*x^2+j*n*w)/(u*(y^2-n*v)-2*w*x*y+v*x^2+n*w^2)
            # C = -(u*(o*v-k*y)+j*w*y+(k*w-j*v)*x-o*w^2)/(u*(y^2-n*v)-2*w*x*y+v*x^2+n*w^2)
            # D = -(l*(n*v-y^2)+w*(p*y-m*n)+x*(m*y-p*v))/(u*(y^2-n*v)-2*w*x*y+v*x^2+n*w^2)
            # E = (u*(p*y-m*n)+x*(-l*y-p*w)+m*x^2+l*n*w)/(u*(y^2-n*v)-2*w*x*y+v*x^2+n*w^2)
            # F = -(u*(p*v-m*y)+l*w*y+(m*w-l*v)*x-p*w^2)/(u*(y^2-n*v)-2*w*x*y+v*x^2+n*w^2)

            o,p,x,y,u,v,w,j,k,l,m,n = 0,0,0,0,0,0,0,0,0,0,0,0
            for pair in self.pointsPairs:
                o += pair[1].x()
                p += pair[1].y()
                x += pair[0].x()
                y += pair[0].y()
                u += pair[0].x()**2
                v += pair[0].y()**2
                w += pair[0].x()*pair[0].y()
                j += pair[1].x()*pair[0].x() # o*x
                k += pair[1].x()*pair[0].y() # o*y
                l += pair[1].y()*pair[0].x() # p*x
                m += pair[1].y()*pair[0].y() # p*y
                n += 1

            self.a = -(j*(n*v-y**2)+w*(o*y-k*n)+x*(k*y-o*v))/(u*(y**2-n*v)-2*w*x*y+v*x**2+n*w**2)
            self.b = (u*(o*y-k*n)+x*(-j*y-o*w)+k*x**2+j*n*w)/(u*(y**2-n*v)-2*w*x*y+v*x**2+n*w**2)
            self.c = -(u*(o*v-k*y)+j*w*y+(k*w-j*v)*x-o*w**2)/(u*(y**2-n*v)-2*w*x*y+v*x**2+n*w**2)
            self.d = -(l*(n*v-y**2)+w*(p*y-m*n)+x*(m*y-p*v))/(u*(y**2-n*v)-2*w*x*y+v*x**2+n*w**2)
            self.e = (u*(p*y-m*n)+x*(-l*y-p*w)+m*x**2+l*n*w)/(u*(y**2-n*v)-2*w*x*y+v*x**2+n*w**2)
            self.f = -(u*(p*v-m*y)+l*w*y+(m*w-l*v)*x-p*w**2)/(u*(y**2-n*v)-2*w*x*y+v*x**2+n*w**2)


    def utmToPixel(self, point):
        x = (self.e*point.x()-self.b*point.y()+self.b*self.f-self.e*self.c)/(self.a*self.e-self.d*self.b)
        y = (-self.d*point.x()+self.a*point.y()+self.d*self.c-self.a*self.f)/(self.a*self.e-self.d*self.b)
        return QgsPoint(x,y)


    def pixelToUtm(self, point):
        x = self.a*point.x()+self.b*point.y()+self.c
        y = self.d*point.x()+self.e*point.y()+self.f
        return QgsPoint(x,y)

    def debug(self):
        QgsMessageLog.logMessage(str(self.pointsPairs))

        








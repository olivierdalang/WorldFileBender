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
from qgis.gui import *

from swift_computer import SwiftComputer

import os.path


class SwiftWorldFiler:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor."""
        # Save reference to the QGIS interface
        self.iface = iface

        self.computer = SwiftComputer(self.iface)

        self.memoryLayer = None
        self.rasterLayer = None
        self.pointsFile = None
        self.worldFile = None




    ########################################################
    # Loading / Unloading plugin
    ########################################################

    def unload(self):
        self.unloadProject()
        del self.toolbar

    def initGui(self):
        self.comboBox = QComboBox()
        self.checkBox = QCheckBox('square pixels')

        self.toolbar = self.iface.addToolBar(u'SwiftWorldFiler')
        self.toolbar.addWidget( QLabel('Swift World Filer raster : ') )
        self.toolbar.addWidget( self.comboBox )
        self.toolbar.addWidget( self.checkBox )

        self.iface.newProjectCreated.connect( self.loadProject )
        self.iface.projectRead.connect( self.loadProject )
        self.loadProject()


    ########################################################
    # Unloading / Loading project
    ########################################################

    def unloadProject(self):
        self.unsetupRaster()

        try:
            self.comboBox.currentIndexChanged.disconnect( self.setupRaster )
        except Exception, e:
            pass
        self.comboBox.clear()       
        self.comboBox.addItem('- No project loaded -', None) 

    def loadProject(self):
        QgsMessageLog.logMessage( 'loadingProject' )
        self.unloadProject()

        self.comboBox.clear()
        self.comboBox.addItem('- Disabled -', None)

        for layer in self.iface.legendInterface().layers():
            if layer.type() == QgsMapLayer.RasterLayer:
                self.comboBox.addItem(layer.name(),layer.id())

        self.comboBox.currentIndexChanged.connect( self.setupRaster )


    ########################################################
    # Unsetup / Setup raster
    ########################################################

    def unsetupRaster(self):
        try:
            self.iface.mapCanvas().extentsChanged.disconnect(self.recomputeWorldfile)
        except Exception, e:
            pass

        try:
            # This can throw an exception since it could be delete already by C++ (on project unload) or the layer could not exist (not created yet)
            QgsMapLayerRegistry.instance().removeMapLayer(self.memoryLayer.id())
        except Exception, e:
            pass
            

        self.memoryLayer = None
        self.rasterLayer = None
        self.pointsFile = None
        self.worldFile = None

    def setupRaster(self, comboBoxId):
        self.unsetupRaster()

        layerId = self.comboBox.itemData(comboBoxId)
        if layerId is None:
            return

        # Set rasterlayer
        self.rasterLayer = QgsMapLayerRegistry.instance().mapLayers()[layerId]
        root = os.path.splitext(self.rasterLayer.dataProvider().dataSourceUri())[0]

        # Set pointfile
        self.pointsFile = root+'.swfPoints'

        # Set worldfile
        self.worldFile = root+'.wld'

        # Load the computer
        self.computer.loadPointfile( self.pointsFile )

        # Create the memory layer
        self.memoryLayer = QgsVectorLayer("Linestring", "SWF - pairs", "memory")
        self.memoryLayer.setCrs( self.iface.mapCanvas().mapRenderer().destinationCrs() )
        self.computer.regenerateMemoryLayer( self.memoryLayer )
        self.memoryLayer.editingStopped.connect( self.memoryLayerChanged )
        self.memoryLayer.loadNamedStyle(os.path.join(os.path.dirname(__file__),"PairStyle.qml"))
        QgsMapLayerRegistry.instance().addMapLayer(self.memoryLayer)

        # Connect the extents changed trigger
        self.iface.mapCanvas().extentsChanged.connect(self.recomputeWorldfile)

        self.recomputeWorldfile()


    ########################################################
    # Actions
    ########################################################

    def memoryLayerChanged(self):
        # Reload the computer
        QgsMessageLog.logMessage("memoryLayerChanged")
        self.computer.loadMemoryLayer( self.memoryLayer )
        self.computer.savePointfile( self.pointsFile )
        self.computer.debug()
        self.recomputeWorldfile()

    def recomputeWorldfile(self):

        if not self.memoryLayer.isEditable():
            QgsMessageLog.logMessage("recomputeWorldfile")
            self.computer.updateTransform( self.iface.mapCanvas().extent(), self.rasterLayer.width(), self.rasterLayer.height() )
            self.computer.generateWorldfile( self.worldFile )
            self.computer.regenerateMemoryLayer( self.memoryLayer )        


        # reloading the raster, not sure of the best way... it doesn't reload the extent (while it takes the worldfile into accout !)
        #self.rasterLayer.triggerRepaint()
        #self.rasterLayer.dataProvider().reload()
        #self.iface.mapCanvas().refresh()

        # We remove and readd it... Quite brutal but effective
        f = self.rasterLayer.dataProvider().dataSourceUri()
        n = self.rasterLayer.name()
        QgsMapLayerRegistry.instance().removeMapLayer(self.rasterLayer.id())
        self.rasterLayer = QgsRasterLayer(f, n)
        QgsMapLayerRegistry.instance().addMapLayer(self.rasterLayer, False)
        QgsProject.instance().layerTreeRoot().addLayer(self.rasterLayer)
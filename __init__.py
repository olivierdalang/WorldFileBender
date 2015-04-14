# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SwiftWorldFiler
                                 A QGIS plugin
 Updates world files dynamically to have a real time georeferencing
                             -------------------
        begin                : 2015-04-08
        copyright            : (C) 2015 by Olivier Dalang
        email                : olivier.dalang@gmail.com
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load SwiftWorldFiler class from file SwiftWorldFiler.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .swift_world_filer import SwiftWorldFiler
    return SwiftWorldFiler(iface)

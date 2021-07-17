# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Created on Wed Nov 28 13:49:33 2018

@author: barnafi
"""


class Markers:
    """
    Class containing a list of integers which define markers
    """

    def __init__(self, markers, neumannSolidMarkers, neumannFluidMarkers, noSlipMarkers, none=None):

        self.markers = markers  # Markers dolfin object

        asList = lambda _x: _x if isinstance(_x, list) else [_x]  # use lists only
        self.neumannSolidMarkers = asList(neumannSolidMarkers)
        self.neumannFluidMarkers = asList(neumannFluidMarkers)
        self.noSlipMarkers = asList(noSlipMarkers)

        # Set null flag
        if none:
            self.NONE = none
        else:
            self.NONE = int(markers.array().max() + 42) 

class MarkersSolid:
    """
    Class containing a list of integers which define markers
    """

    def __init__(self, markers, neumannMarkers, robinMarkers, none=None):

        self.markers = markers  # Markers dolfin object

        asList = lambda _x: _x if isinstance(_x, list) else [_x]  # use lists only
        self.neumannMarkers = asList(neumannMarkers)
        self.robinMarkers = asList(robinMarkers)

        # Set null flag
        if none:
            self.NONE = none
        else:
            self.NONE = int(markers.array().max() + 42)
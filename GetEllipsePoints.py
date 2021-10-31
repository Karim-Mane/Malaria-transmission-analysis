#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 10:30:48 2020

@author: karim
"""

import math
import numpy as np
import sys
from shapely.geometry import LineString  #from shapely.geometry import LineString
from shapely import affinity
from haversine import haversine
import pandas as pd


def GetEllipseAxisLength(p1_lat, p1_lng, p2_lat, p2_lng,radius_in_meters):
    c2 = haversine((p1_lat, p1_lng), (p2_lat, p2_lng)) * 1000.0
    if radius_in_meters < c2:
        radius_in_meters=c2
        #raise ValueError("Please specify radius larger than the distance between the two input points.")
    a = radius_in_meters / 2.0
    b = math.sqrt(pow(a, 2) - pow(c2 / 2.0, 2))
    return a, b
    #print(a, b)


def  GetEllipsePointInMeters(a, b, num_points):
    """
    :param a: length of "horizontal" axis in meters
    :param b: length of "vertical" axis in meters
    :param num_points: (half the) number of points to draw
    :return: List of tuples of perimeter points on the ellipse, 
             centered around (0,0), in m.
    """
    x_points = list(np.linspace(-a, a, num_points))#[1:-1]
    y_points_pos = [math.sqrt(pow(a, 2) - pow(x, 2)) * (float(b) / float(a)) for x in x_points]
    y_points_neg = [-y for y in y_points_pos]

    perimeter_points_in_meters = [tuple([-a, 0])] + [tuple([x, y]) for x, y in zip(x_points, y_points_pos)] + [tuple([a, 0])] + list(reversed([tuple([x, y]) for x, y in zip(x_points, y_points_neg)]))
    return perimeter_points_in_meters
    #print(perimeter_points_in_meters)
    


def AddMetersToPoint(center_lng, center_lat, dx, dy):
    """
    :param center_lng, center_lat: GPS coordinates of the center  
           between the two input points.
    :param dx: distance to add to x-axis (lng) in meters
    :param dy: distance to add to y-axis (lat) in meters
    """
    R_EARTH = 6371000
    new_x = (center_lng + (dx / R_EARTH) * (180 / math.pi) / np.cos(center_lat * math.pi/180))
    new_y = center_lat + (dy / R_EARTH) * (180 / math.pi)
    return tuple([new_x, new_y])
    

def main():
    """
    Enter ellipse centers in lat-lng and ellipse perimeter points    
    around the origin (0,0), and get points on the perimeter of the 
    ellipse around the centers in lat-lng.
    :param p1_lat: lat coordinates of center point 1
    :param p1_lng: lng coordinates of center point 1
    :param p2_lat: lat coordinates of center point 2
    :param p2_lng: lng coordinates of center point 2
    :param perimeter_points_in_meters: List of tuples of perimeter 
           points on the ellipse, centered around (0,0), in m.
    :return: List of the points we really want, tuples of (lat,lng)
    """
    
    allArguments = sys.argv
    #print(allArguments)
    argumentList = allArguments[1:]
    p1_lat = float(argumentList[0])
    p1_lng = float(argumentList[1])
    p2_lat = float(argumentList[2])
    p2_lng = float(argumentList[3])
    radius_in_meters = float(argumentList[4])
    number_of_points = int(argumentList[5])
    outDir = argumentList[6]
    outFile = argumentList[7]
    courbature = int(argumentList[8])
    #perimeter_points_in_meters = list(argumentList[4])
    
    val1 = GetEllipseAxisLength(p1_lat, p1_lng, p2_lat, p2_lng,radius_in_meters)
    val = val1[1] + courbature  #increase this to amplify the curve of the line 
    val2 = GetEllipsePointInMeters(val1[0], val, number_of_points) #val2 = GetEllipsePointInMeters(val1[0], val1[1], number_of_points)
    
    center_lng = (p1_lng + p2_lng) / 2.0
    center_lat = (p1_lat + p2_lat) / 2.0
    perimeter_points_in_lng_lat = [AddMetersToPoint(center_lng, center_lat, p[0], p[1]) for p in val2]
    ellipse = LineString(perimeter_points_in_lng_lat)
    
    angle = np.degrees(math.atan2(p2_lat - p1_lat, p2_lng - p1_lng))
    ellipse_rotated = affinity.rotate(ellipse, angle)

    ellipse_points_lng_lat = list(ellipse_rotated.coords)
    ellipse_points = [tuple([p[1], p[0]]) for p in ellipse_points_lng_lat]
    #return ellipse_points
    df = pd.DataFrame(ellipse_points, columns=['lat', 'long'])
    df.to_csv(outDir+'/' + outFile + '.csv', index=False, header=True)
    #print(df)
    
    
    
if __name__ == '__main__':
    main()
    
    

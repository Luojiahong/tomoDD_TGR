# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 20:59:37 2024

@author: luo
"""
import numpy as np
# Convert from local Cartesian coordinates to latitude & longitude
def para(olat,olon,rota):
    # 计算变量
    rad = np.pi / 180.0
    olat *= 60.0  # minutes N
    olon *= 60.0  # minutes W
    rotate = rota * rad
    
    # WGS72 constants
    rearth = 6378.135
    ellip = 298.26

    # Calculate RLATC
    phi = olat * rad / 60.0  # geographical latitude
    beta = phi - np.sin(2 * phi) / ellip  # geocentric latitude
    rlatc = np.tan(beta) / np.tan(phi)
    
    # Calculate aa & bb
    lat1 = np.arctan(rlatc * np.tan(olat * rad / 60.0))  # geocentric lat for OLAT
    lat2 = np.arctan(rlatc * np.tan((olat + 1.0) * rad / 60.0))  # geocentric lat for (OLAT + 1 min)
    dela = lat2 - lat1
    r = rearth * (1.0 - (np.sin(lat1) ** 2) / ellip)  # radius for lat = OLAT
    aa = r * dela  # 1 min geographical lat
    delb = np.arccos(np.sin(lat1) ** 2 + np.cos(rad / 60.0) * np.cos(lat1) ** 2)
    bc = r * delb  # 1 min geographical lon
    bb = r * delb / np.cos(lat1)
    
    sint = np.sin(rotate)
    cost = np.cos(rotate)
    return aa,bb,rlatc,sint,cost
    
def redist(xkm, ykm,olat=31.0,olon=110.5,rota=0):
    # 包含的常量，这些需要从 "geocoord.inc" 文件中获取
    aa,bb,rlatc,sint,cost = para(olat, olon, rota)
    rad = np.pi / 180.0
    olat *= 60.0  # minutes N
    olon *= 60.0  # minutes W
    
    xx = xkm
    yy = ykm

    # 逆时针旋转坐标
    y = yy * cost - xx * sint
    x = yy * sint + xx * cost

    if abs(aa) < 1e-7:
        raise ValueError(f"subr. redist: aa={aa:.5f} bb={bb:.5f} division by zero, run stops here")

    q = y / aa
    lat = (q + olat) / 60.0
    xlat = q + olat - 60.0 * lat
    yp = 60.0 * lat + xlat
    lat1 = np.arctan(rlatc * np.tan(yp * rad / 60.0))
    lat2 = np.arctan(rlatc * np.tan(olat * rad / 60.0))
    lat3 = (lat1 + lat2) / 2.0
    clat1 = np.cos(lat3)
    bcl = bb * clat1

    if abs(bcl) < 1e-6:
        raise ValueError(f"subr. redist: aa={aa:.5f} bb={bb:.5f} cos(lat1)={clat1:.7f} division by zero, run stops here")

    p = x / (bb * clat1)
    lon = (p + olon) / 60.0
    xlon = p + olon - 60.0 * lon
    xlat = lat + xlat / 60.0
    xlon = lon + xlon / 60.0

    return xlat, xlon
# Convert latitude and longitude to kilometers relative
# to center of coordinates by short distance conversion.
def dist(xlat, xlon, olat=31.0, olon=110.5,rota=0):
    
    # Input: xlat, xlon
    # Output: xkm, ykm
    aa,bb,rlatc,sint,cost = para(olat, olon, rota)
    rad = np.pi / 180.0
    olat *= 60.0  # minutes N
    olon *= 60.0  # minutes W

    # Set up short distance conversion
    q = 60 * xlat - olat
    yp = q + olat
    lat1 = np.arctan(rlatc * np.tan(rad * yp / 60.0))
    lat2 = np.arctan(rlatc * np.tan(rad * olat / 60.0))
    lat3 = (lat2 + lat1) / 2.0
    xx = 60 * xlon - olon
    q = q * aa
    xx = xx * bb * np.cos(lat3)
    
    if rota != 0:
        # Rotate coordinate system anticlockwise
        yp = cost * q + sint * xx
        xx = cost * xx - sint * q
        q = yp

    xkm = xx
    ykm = q

    return xkm, ykm

def xy2ll_sdc(ds_vel,clat,clon,rota):
    # convert the XY to lonlat
    lons = []
    lats = []
    for x1 in ds_vel.longitude.values:
        _,lon = redist(x1, 0,clat,clon,rota)
        lons.append(lon)
    for y1 in ds_vel.latitude.values:
        lat,_ = redist(0, y1,clat,clon,rota)
        lats.append(lat)
    yn = lats
    xn = lons
    ds_vel['longitude'] = xn
    ds_vel['latitude'] = yn
    return ds_vel
if __name__ == '__main__':
    lat,lon = redist(0, 30,olat=31.0,olon=110.5,rota=0)
    xkm,ykm = dist(xlat=lat,xlon=lon)
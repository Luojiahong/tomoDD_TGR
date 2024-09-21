import pygmt 
import pandas as pd 
import geopandas as gpd
from shapely.geometry import LineString,box
import numpy as np
from obspy import geodetics
import warnings
warnings.filterwarnings("ignore")

def read_qml(qmlfile):
    lines = open(qmlfile,'r').readlines()
    geo_categories = {}
    for line in lines:
        tmp = line.strip().split()
        if tmp[0] == '<category':
            symbol = [a for a in tmp if a.find('symbol') >= 0][0].split('"')[1]
            name = [a for a in tmp if a.find('value') >= 0][0].split('"')[1]
            label = [a for a in tmp if a.find('label') >= 0][0].split('"')[1]
            label = label.replace('&lt;sup>','@+').replace('&lt;/sup>','@+').replace('&lt;sub>','@-').replace('&lt;/sub>','@-')
            geo_categories.setdefault(symbol,{}).setdefault('name',name)
            geo_categories.setdefault(symbol,{}).setdefault('label',label)
            if symbol == '':continue
        if tmp[0] == '<symbol':
            symbol = [a for a in tmp if a.find('name') >= 0][0].split('"')[1]
            if symbol == '':continue
        if line.find('<Option') >=0 and line.find('name="color"') >= 0 and symbol != '':
            color = [a for a in tmp if a.find('value') >= 0][0].split('"')[1].split(',')
            c = '/'.join([color[0],color[1],color[2]])
            geo_categories[symbol].setdefault('color',c)
    geo_strata = {}
    for i,symbol in enumerate(list(geo_categories)):
        name = geo_categories[symbol]['name']
        label = geo_categories[symbol]['label']
        color = geo_categories[symbol]['color']
        geo_strata.setdefault(name,{}).setdefault('name',name)
        geo_strata.setdefault(name,{}).setdefault('label',label)
        geo_strata.setdefault(name,{}).setdefault('color',color)
    return geo_categories,geo_strata
def get_profile(model,grid_topo,P1,spacing,wave='vp'):
    lon1,lat1,lon2,lat2 = P1
    # 生成测线
    points = pygmt.project(center=[lon1,lat1],endpoint=[lon2,lat2],generate=spacing[0], unit=True)
    ele = pygmt.grdtrack(points=points, grid=grid_topo,newcolname='elevation')

    depths = model.depth.values
    depths = np.arange(0,model.depth.values.max()+spacing[1],spacing[1])
    for i,depth in enumerate(depths):
        grid = model[wave].interp(depth=depth,method='linear')
        track = pygmt.grdtrack(grid=grid,points=points,newcolname='vs')
        track['depth'] = [depth]*len(track)
        if i ==0:
            profile = track.copy()
        else:
            profile = pd.concat([profile,track])
    p_region = pygmt.info(data=profile[['p','depth']],per_column=True,spacing=spacing)
    p_grid = pygmt.xyz2grd(data=profile[['p','depth','vs']],spacing=spacing,region=p_region)
    return points,ele,p_grid,p_region,profile
def get_intersection(gdf,line):
    # 提取交点
    lon1,lat1,lon2,lat2 = line
    sline = LineString([(lon1,lat1),(lon2,lat2)])
    # 计算交点
    intersection_points = []
    names = []
    for i in range(len(gdf)):
        geom = gdf.iloc[i]
        intersection = geom.geometry.intersection(sline)
        #intersection = intersection.buffer(1)
        if not intersection.is_empty:
            try:
                names.append(geom['name'])
            except:
                names.append(i)
            intersection_points.append(intersection)

    # 转换为 GeoDataFrame
    intersection_gdf = gpd.GeoDataFrame(geometry=intersection_points)
    intersection_gdf['name'] = names
    return intersection_gdf
def get_line(topo,gdf,line):
    intersection_gdf = get_intersection(gdf,line)
    lon_begin,lat_begin,lon_end,lat_end = line
    idx = 0
    geo_dict = {}
    data_points = []
    for i in range(len(intersection_gdf)):
        geom = intersection_gdf.iloc[i]
        geom_type = geom.geometry.geom_type
        data = []
        if geom_type == 'Point':
            xy_coords = list(geom.geometry.coords)
            lon1,lat1 = xy_coords[0][0:2]
            data_points.append([lon1,lat1])
            idx = idx+1
        if geom_type == 'LineString':
            xy_coords = list(geom.geometry.coords)
            lon1,lat1 = xy_coords[0][0:2]
            lon2,lat2 = xy_coords[1][0:2]
            data.append([lon1,lat1])
            data.append([lon2,lat2])
            dist = 0
            if lon_begin != lon1:
                dist0 = pygmt.project(center=[lon_begin,lat_begin],endpoint=[lon1,lat1],generate=0.01, unit=True)
                dist = dist0.p.values[-1]
            points = pygmt.project(center=[lon1,lat1],endpoint=[lon2,lat2],generate=0.01, unit=True)
            ele0 = pygmt.grdtrack(points=points, grid=topo,newcolname='elevation')
            idx = idx+1
            geo_dict[idx] = {'coords':data,'name':geom['name'],'track':ele0,'dist':dist}
        if geom_type == 'MultiLineString':
            for line in geom.geometry.geoms:
                xy_coords = list(line.coords)
                lon1,lat1 = xy_coords[0][0:2]
                lon2,lat2 = xy_coords[1][0:2]
                data.append([lon1,lat1])
                data.append([lon2,lat2])
                dist = 0
                if lon_begin != lon1:
                    dist0 = pygmt.project(center=[lon_begin,lat_begin],endpoint=[lon1,lat1],generate=0.01, unit=True)
                    dist = dist0.p.values[-1]
                points = pygmt.project(center=[lon1,lat1],endpoint=[lon2,lat2],generate=0.01, unit=True)
                ele0 = pygmt.grdtrack(points=points, grid=topo,newcolname='elevation')
                idx = idx+1
                geo_dict[idx] = {'coords':data,'name':geom['name'],'track':ele0,'dist':dist}
        if geom_type == 'MultiPolygon':
            print('Luo')
    if len(data_points)>0:
        data_points = pd.DataFrame(data_points,columns=['longitude','latitude'])
        data_points = pygmt.grdtrack(points=data_points, grid=topo,newcolname='elevation')
        track = pygmt.project(data=data_points,center=[lon_begin,lat_begin],endpoint=[lon_end,lat_end],length=None,width=None,unit=True)
        track.columns = ['x','y','z','p','q','r','s']
        geo_dict[idx] = {'coords':data_points,'name':'fault','track':track,'dist':0}
    return geo_dict
def project_data(catalog,cap,P_A,width=15):
    lon1,lat1,lon2,lat2 = P_A
    reloc = pygmt.project(data=catalog[['longitude','latitude','depth']],
                      center=[lon1,lat1],endpoint=[lon2,lat2],
                      length='w',width=[-width,width],convention='pzq',unit=True)
    reloc.columns = ['p','z','q']
    meca = pygmt.project(data=cap[['longitude','latitude']],center=[lon1,lat1],endpoint=[lon2,lat2],width=None,unit=True)
    cap.loc[:,'p'] =  meca[2]
    cap.loc[:,'q'] = meca[3]
    meca = cap[abs(cap.q) <= width]
    return reloc,meca
def plot_meca(fig,cap,P1,width=15,depth=10):
    lon1,lat1,lon2,lat2 = P1
    with pygmt.helpers.GMTTempFile() as temp_file:
        with open(temp_file.name, 'w') as f:
            for i in range(len(cap)):
                line = cap.iloc[i]
                f.write(f'{line.longitude} {line.latitude} {line.depth} {line.strike1} {line.dip1} {line.rake1} {line.magnitude}\n')
        with pygmt.clib.Session() as session:
            cmd = f'{temp_file.name} -Sa0.5 -Aa{lon1}/{lat1}/{lon2}/{lat2}/90/{width}/0/{depth} -Gblack -Q'
            session.call_module('coupe', cmd)
def plot_profile_location(fig,P,label=['A','B'],offset=['0.0c/-0.1c',None]):
    lon1,lat1,lon2,lat2 = P
    dist,az,baz = geodetics.gps2dist_azimuth(lat1=lat1,lon1=lon1,lat2=lat2,lon2=lon2)
    track = pygmt.project(data=None,center=[lon1,lat1],endpoint=[lon2,lat2],generate=0.01,unit=True)
    length = track.p.max()
    fig.plot(x=track.r,y=track.s,pen='0.5p,black,--')
    fig.text(x=track.r.values[0],y=track.s.values[0],text=label[0],justify='BL',fill='white',offset=offset[0],font='10p,blue',no_clip=True)
    fig.text(x=track.r.values[-1],y=track.s.values[-1],text=label[1],justify='BL',fill='white',offset=offset[1],font='10p,blue',no_clip=True)
    #data = [[(lon1+lon2)/2,(lat1+lat2)/2,az,length,2*width]]
    #fig.plot(data=data,style='J',pen='0.2p,--')
    return length
def plot_profile_value(fig,model,topo,catalog,meca,P,scale,label=['N','S'],colorbar=False,interval=100,width=3):
    lon1,lat1,lon2,lat2 = P
    _,ele,vp,_,_ = get_profile(model,topo,P,[0.1,0.1],wave='vp')
    _,_,vs,_,_ = get_profile(model,topo,P,[0.1,0.1],wave='vs')
    _,_,dwsp,_,_ = get_profile(model,topo,P,[0.1,0.1],wave='dwsp')
    _,_,dwss,_,_ = get_profile(model,topo,P,[0.1,0.1],wave='dwss')

    fig.plot(x=ele.p, y=-1*ele.elevation, pen='0.1p',close='+y0', fill='gray')
    track = pygmt.project(data=catalog[['longitude','latitude','depth']],center=[lon1,lat1],endpoint=[lon2,lat2],length=None,width=None,unit=True)
    catalog['p'] = track[3]
    catalog['q'] = track[4]
    reloc = catalog[(catalog.q>=-width)&(catalog.q<=width)]
    fig.grdimage(grid=vp,cmap='velp.cpt')
    fig.grdcontour(grid=dwsp,interval=interval,limit=[interval,interval],pen='0.7p,white')
    if colorbar:fig.colorbar(position='JMR+w0c/0.2c+o0.2c/0c+e',frame=['xaf+lP-wave velocity (km/s)'],cmap='velp.cpt')
    fig.plot(x=reloc.p,y=reloc.depth,style='cc',size=0.06*(reloc.magnitude+1),fill=reloc.reltime,cmap='time.cpt',pen='0.2p')
    if len(meca) >=1:fig.plot(x=meca.p,y=meca.depth,style='h0.1i',fill='red',pen='0.1p')
    fig.text(position='cBL',offset='0.1c/0.1c',text=label[0],font='10p',no_clip=True)
    fig.text(position='cBR',offset='-0.1c/0.1c',text=label[1],font='10p',no_clip=True)

    fig.shift_origin(yshift='-{}c'.format(8*scale*1.0+0.2))
    fig.grdimage(grid=vs,cmap='vels.cpt',frame=['xaf','yaf+lDepth (km)','Wsen'])
    fig.grdcontour(grid=dwss,interval=interval,limit=[interval,interval],pen='0.7p,white')
    if colorbar:fig.colorbar(position='JMR+w0c/0.2c+o0.2c/0c+e',frame=['xaf+lS-wave velocity (km/s)'],cmap='vels.cpt')
    fig.plot(x=reloc.p,y=reloc.depth,style='cc',size=0.06*(reloc.magnitude+1),fill=reloc.reltime,cmap='time.cpt',pen='0.2p')
    if len(meca) >=1:fig.plot(x=meca.p,y=meca.depth,style='h0.1i',fill='red',pen='0.1p')

    fig.shift_origin(yshift='-{}c'.format(8*scale*1.0+0.2))
    fig.grdimage(grid=vp/vs,cmap='vpvs.cpt',frame=['xaf+lDistance along the profile (km)','yaf+lDepth (km)','WSen'])
    fig.grdcontour(grid=dwsp,interval=interval,limit=[interval,interval],pen='0.7p,white')
    if colorbar:fig.colorbar(position='JMR+w0c/0.2c+o0.2c/0c+e',frame=['xaf+lVp/Vs ratio'],cmap='vpvs.cpt')
    fig.plot(x=reloc.p,y=reloc.depth,style='cc',size=0.06*(reloc.magnitude+1),fill=reloc.reltime,cmap='time.cpt',pen='0.2p')
    if len(meca) >=1:fig.plot(x=meca.p,y=meca.depth,style='a0.1i',fill='red',pen='0.1p')
def plot_elevation(fig,topo,TGR,geo_strata,fault,geo,line):
    lon1,lat1,lon2,lat2 = line
    points = pygmt.project(center=[lon1,lat1],endpoint=[lon2,lat2],generate=0.01, unit=True)
    elevation = pygmt.grdtrack(points=points, grid=topo,newcolname='elevation')

    gdf_fault = get_line(topo,fault,line)
    gdf_strata = get_line(topo,geo,line)
    gdf_TGR = get_line(topo,TGR,line)
    fig.plot(x=elevation.p,y=elevation.elevation/1000.,pen='0.5p',close='+y0', fill='gray')
    for i in list(gdf_strata):
        data = gdf_strata[i]
        ele,name,dist = data['track'],data['name'],data['dist']
        color = geo_strata[name]['color']
        fig.plot(x=ele.p+dist,y=ele.elevation/1000,pen=f'1p,{color}',close='+y0', fill=color)
    for i in list(gdf_TGR):
        data = gdf_TGR[i]
        ele,dist = data['track'],data['dist']
        #fig.plot(x=ele.p+dist,y=ele.elevation/1000,pen='1p,0/112/255')
    try:
        gdf_fault = gdf_fault[list(gdf_fault)[0]]['track']
        fig.plot(x=gdf_fault.p,y=gdf_fault.z/1000.,style='h0.1c',pen='0.3p',fill='red')
    except:
        pass
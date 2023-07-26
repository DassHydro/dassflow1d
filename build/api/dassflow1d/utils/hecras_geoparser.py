# coding: utf-8
from __future__ import print_function

import math
import numpy as np
#from osgeo import ogr
import matplotlib.pyplot as plt
from shapely.geometry import LineString, mapping


class CenterLine:
  
  def __init__(self, x, y):
    
    self.x = x
    self.y = y
    self.__compute_length__()
    
  def __compute_length__(self):
    
    self.length = 0.0
    for ix in range(1, self.x.size):
      dx = self.x[ix] - self.x[ix-1]
      dy = self.y[ix] - self.y[ix-1]
      self.length += np.sqrt(dx**2 + dy**2)
    
    
class Reach:
  
  def __init__(self, name, centerline=None):
    
    self.name = name
    self.centerline = None
    if isinstance(centerline, tuple):
      self.centerline = CenterLine(centerline[0], centerline[1])
    elif isinstance(centerline, CenterLine):
      self.centerline = centerline
    elif centerline is not None:
      raise ValueError("'centerline' must be a tuple (x,y) or a Centerline")
    
    self.length = None
    if self.centerline:
      self.__compute_length__()
      
    self.sections = []
      
      
  def ogr_geometry(self, dxnodes=None, with_xnodes=False):
    
    if with_xnodes and dxnodes:
      xnodes = []
  
    geometry = ogr.Geometry(ogr.wkbLineString)
    if dxnodes is None:
      for ix in range(0, self.centerline.x.size):
        geometry.AddPoint(self.centerline.x[ix], self.centerline.y[ix])
        
    else:
      
      length = 0.0
      length_prev = 0.0
      npoints = 0
      for ix in range(1, self.centerline.x.size):
        dx = self.centerline.x[ix] - self.centerline.x[ix-1]
        dy = self.centerline.y[ix] - self.centerline.y[ix-1]
        dl = np.sqrt(dx**2 + dy**2)
        if length + dl > (0.5 + float(npoints)) * dxnodes:
          alpha = (length + dl - (0.5 + float(npoints)) * dxnodes) / dl
          #print("xnode=", (0.5 + float(npoints)) * dx, npoints, dxnodes)
          #print("alpha=", alpha, (length, (0.5 + float(npoints)) * dxnodes, length+dl), length_prev)
          #choice = raw_input()
          x = self.centerline.x[ix-1] + alpha * dx
          y = self.centerline.y[ix-1] + alpha * dy
          geometry.AddPoint(x,y)
          if with_xnodes:
            xnodes.append(length + alpha * dl)
          length_prev = (0.5 + float(npoints)) * dxnodes
          npoints += 1
        length += dl
      
    if with_xnodes and dxnodes:
      return geometry, xnodes
    else:
      return geometry

      
  def shapely_geometry(self):
  
    x = self.centerline.x
    y = self.centerline.y
    geometry = LineString([(x[i], y[i]) for i in range(0, x.size)])
      
    return geometry
      
    
    
  def __compute_length__(self):
    
    self.length = 0.0
    #print "__compute_length__:nx=%i" % (self.centerline.x.size)
    #print "                   x=%f,y=%f" % (self.centerline.x[0], self.centerline.y[0])
    for ix in range(1, self.centerline.x.size):
      dx = self.centerline.x[ix] - self.centerline.x[ix-1]
      dy = self.centerline.y[ix] - self.centerline.y[ix-1]
      self.length += np.sqrt(dx**2 + dy**2)
      #print "                   x=%f,y=%f (length=%f)" % (self.centerline.x[ix], 
                                                          #self.centerline.y[ix], 
                                                          #self.length)
    
    
class CrossSection:
  
  def __init__(self, name, dist, t, z, banks=None, cut_line=None, 
               description=None):
    
    self.name = name
    self.dist = dist
    self.banks = banks
    self.t = t
    self.z = z
    self.cut_line = cut_line
    self.description = description
    
    if cut_line:
      self.__compute_t_cut_line__()
    if banks is not None:
      ileft = np.searchsorted(t, banks[0])
      iright = np.searchsorted(t, banks[1])
      self.ibanks = [ileft, iright]
      self.tbanks = [t[ileft], t[iright]]
      self.zbanks = [z[ileft], z[iright]]
      
      
  def get_top_width_segment(self, zs, epsilon=0.01, valid_plot=False):
    
    if self.banks is not None:
      if zs > self.zbanks[0] + epsilon:
        i = 0
      else:
        i = self.ibanks[0]
      if zs > self.zbanks[1] + epsilon:
        ilast = len(self.z)
      else:
        ilast = self.ibanks[1]+1
    else:
      i = 0
      imax = len(self.z)
    while self.z[i] > zs - epsilon:
      i+=1
    imin = i
    
    # Find point under water surface
    while i < ilast:
      if self.z[i] < zs - epsilon:
        imax = i
      i = i + 1
    imax = min(imax+1, ilast)
      
    #while self.z[i] < zs:
      #i+=1
    #imax = i
    
    x = []
    y = []
    if valid_plot: t = []
    if imin == 0:
      x.append(self.xt[0])
      y.append(self.xt[0])
      if valid_plot: t.append(self.t[0])
      alpha = 0.0
      ifirst = imin+1
    else:
      alpha = (zs - self.z[imin]) / (self.z[imin-1] - self.z[imin])
      x.append(self.xt[imin] + alpha * (self.xt[imin-1] - self.xt[imin]))
      y.append(self.yt[imin] + alpha * (self.yt[imin-1] - self.yt[imin]))
      if valid_plot: t.append(self.t[imin] + alpha * (self.t[imin-1] - self.t[imin]))
      if alpha > 1e-12:
        ifirst = imin+1
      else:
        ifirst = imin
    for i in range(ifirst, imax):
      x.append(self.xt[i])
      y.append(self.yt[i])
      if valid_plot: t.append(self.t[i])
    if imax < ilast-1:
      alpha = (zs - self.z[imax-1]) / (self.z[imax] - self.z[imax-1])
      if alpha > 1e-12:
        x.append(self.xt[imax-1] + alpha * (self.xt[imax] - self.xt[imax-1]))
        y.append(self.yt[imax-1] + alpha * (self.yt[imax] - self.yt[imax-1]))
        if valid_plot: t.append(self.t[imax-1] + alpha * (self.t[imax] - self.t[imax-1]))
        
    if valid_plot:
      plt.plot(self.t, self.z)
      plt.plot(t, zs * np.ones(len(t)))
      if self.banks is not None:
        plt.plot(self.tbanks, self.zbanks, 'ro')
      plt.show()
        
    return x, y
      
      
  def ogr_geometry(self):
  
    x = self.cut_line[0]
    y = self.cut_line[1]
    geometry = ogr.Geometry(ogr.wkbLineString)
    for ix in range(0, x.size):
      geometry.AddPoint(x[ix], y[ix])
      
    return geometry
      
      
  def shapely_geometry(self, z=False):
  
    if z is True:
        x = self.xt
        y = self.yt
        z = self.z
        geometry = LineString([(x[i], y[i], z[i]) for i in range(0, x.size)])
    else:
        x = self.cut_line[0]
        y = self.cut_line[1]
        geometry = LineString([(x[i], y[i]) for i in range(0, x.size)])
      
    return geometry
      
      
  def __compute_t_cut_line__(self):
  
    x = self.cut_line[0]
    y = self.cut_line[1]
    alpha1 = np.zeros(x.size)
    alpha2 = np.zeros(self.t.size)
    dist = 0.0
    for ix in range(1, x.size):
      dist += np.sqrt((x[ix] - x[ix-1])**2 + (y[ix] - y[ix-1])**2)
      alpha1[ix] = dist
    dist = 0.0
    for it in range(1, self.t.size):
      dist += self.t[it] - self.t[it-1]
      alpha2[it] = dist
    #print alpha1
    #print alpha2
      
    self.xt = np.interp(alpha2, alpha1, x)
    self.yt = np.interp(alpha2, alpha1, y)
    
    
class HecRasGeoParser:

  def __init__(self, geomfile, parse_only_real_sections=True, debug=False):
    
    self.parse_only_real_sections = parse_only_real_sections
    self.reaches = []
    self.sections = []
    
    geomIO = open(geomfile, "r")
    
    self.__parse_header__(geomIO)
    
    line = geomIO.readline().strip()
    words = line.split("=")
    while words[0] == "Junct Name":
        
      if debug:
        print("words=", words)
        
        # TODO parse junctions !
      for i in range(0, 8):
        line = geomIO.readline().strip()
      line = geomIO.readline().strip()
      words = line.split("=")
      
    while words[0] == "River Reach":
        
      if debug:
        print("Parse reach: %s ..." % words[1])
    
      reach = self.__parse_reach__(geomIO, line)
      self.reaches.append(reach)
      line = geomIO.readline().strip()
      words = line.split("=")
        
      while words[0] == "Type RM Length L Ch R ":
            
        #if debug:
            #print("Parse section...")
        
        section = self.__parse_section__(geomIO, line)
        
        if section is not None:
            self.sections.append(section)
            self.reaches[-1].sections.append(section)
        
        line = geomIO.readline().strip()
        words = line.split("=")
      
    #line = geomIO.readline().strip()
    #words = line.split("=")
    #print "words:", words
    geomIO.close()
    
    
  def ogr_channel_geometry(self):
    
    left_geometry = ogr.Geometry(ogr.wkbLineString)
    right_geometry = ogr.Geometry(ogr.wkbLineString)
    for section in self.sections:
      
      banks = section.banks
      ileft = np.searchsorted(section.t, banks[0])
      left_geometry.AddPoint(section.xt[ileft], section.yt[ileft])
      iright = np.searchsorted(section.t, banks[1])
      right_geometry.AddPoint(section.xt[iright], section.yt[iright])
      #print "ileft=%i, iright=%i (banks=%f,%f)" % (ileft, iright, banks[0], banks[1])
      
    return left_geometry, right_geometry
    
    
  def __parse_header__(self, geomIO):
    
    # Parse title
    line = geomIO.readline().strip()
    words = line.split("=")
    if words[0] != "Geom Title":
      raise RuntimeError("Expected keyword 'Geom Title' in geometry file"
                         "(current line is %s)" % line)
    self.title = words[1]
    
    # Parse program version
    line = geomIO.readline().strip()
    words = line.split("=")
    if words[0] != "Program Version":
      raise RuntimeError("Expected keyword 'Program Version' in geometry file"
                         "(current line is %s)" % line)
    self.program_version = float(words[1])
    
    # Parse viewing rectangle
    line = geomIO.readline().strip()
    words = line.split("=")
    if words[0] != "Viewing Rectangle":
      raise RuntimeError("Expected keyword 'Viewing Rectangle' in geometry file"
                         "(current line is %s)" % line)
    self.viewing_rectangle = map(float, words[1].split(','))
    
    # Skip blank line
    line = geomIO.readline().strip()
    
    
  def __parse_reach__(self, geomIO, line):
    
    # Parse reach name
    words = line.split("=")
    name = words[1]
    
    # Parse number of xy
    line = geomIO.readline().strip()
    words = line.split("=")
    n = int(words[1])
    
    # Parse number x,y tuples
    x, y = self.__parse_xy_tuples__(geomIO, n)
    #print "x.size=", x.size, n
    #x = []
    #y = []
    #while len(x) < nxy:
      #line = geomIO.readline().strip()
      ##print line
      
      #words0 = line.split()
      #words = []
      #for word in words0:
        #if len(word) > 16:
          #for i in range(0, int(math.ceil(len(word)/16.0))):
            #words.append(word[i*16:i*16+16])
        #else:
          #words.append(word)
          
      ##print words
          
      #for i in range(0, len(words)/2):
        #x.append(float(words[i*2]))
        #y.append(float(words[i*2+1]))
        
      ##print x
      ##print "len(x)=", len(x), nxy
      ##choice = raw_input()
      
    reach = Reach(name, centerline=(x, y))
        
    # Parse remaining lines
    line = geomIO.readline().strip()
    while len(line) > 0:
      line = geomIO.readline().strip()
      
    return reach
    
    
  def __parse_section__(self, geomIO, line):
    
    # Parse header
    words = line.split("=")
    dist = words[1].split(',')[1]
    if dist.strip()[-1] != "*":
      real_section = True
      dist = float(words[1].split(',')[3])
      #print "Real CS : ", dist
    else:
      real_section = False
      dist = float(words[1].split(',')[3])
      
    description = None
    line = geomIO.readline().strip()
    if line.strip() == "BEGIN DESCRIPTION:":
      line = geomIO.readline().strip()
      description = line
      while line.strip() != "END DESCRIPTION:":
        line = geomIO.readline().strip()
        description = description + "\n%s" % line
      line = geomIO.readline().strip()
      
    # Parse node name
    words = line.split("=")
    if words[0] == "Node Name":
      name = words[1]
      line = geomIO.readline().strip()
    else:
      name = None
      
    # Parse cut line
    words = line.split("=")
    if words[0] == "XS GIS Cut Line": # Optional
      #raise RuntimeError("Expected 'XS GIS Cut Line' keyword in line %s" % line)
      n = int(words[1])
      x, y = self.__parse_xy_tuples__(geomIO, n)
      cutline = (x, y)
      line = geomIO.readline().strip()
      words = line.split("=")
    else:
      # TODO orthogonal profile ?
      cutline = None
    
    if words[0] != "Node Last Edited Time":
      raise RuntimeError("Expected 'Node Last Edited Time' keyword in line %s" % line)
    
    line = geomIO.readline().strip()
    words = line.split("=")
    if words[0] != "#Sta/Elev":
      raise RuntimeError("Expected '#Sta/Elev' keyword in line %s" % line)
    
    n = int(words[1])
    t, z = self.__parse_xy_tuples__(geomIO, n, maxsize=8)
    
    # Skip #Mann lines
    line = geomIO.readline().strip()
    words = line.split("=")
    if words[0] != "#Mann":
      raise RuntimeError("Expected '#Mann' keyword in line %s" % line)
    line = geomIO.readline().strip()
    
    # Parse  bank stations
    line = geomIO.readline().strip()
    words = line.split("=")
    if words[0] != "Bank Sta":
      raise RuntimeError("Expected 'Bank Sta' keyword in line %s" % line)
    banks = list(map(float, words[1].split(',')))
    
    # Skip XS Rating Curve line
    line = geomIO.readline().strip()
    words = line.split("=")
    if words[0] != "XS Rating Curve":
      raise RuntimeError("Expected 'XS Rating Curve' keyword in line %s" % line)
    
    # Skip Exp/Cntr line
    line = geomIO.readline().strip()
    words = line.split("=")
    if words[0] != "Exp/Cntr":
      raise RuntimeError("Expected 'Exp/Cntr' keyword in line %s" % line)
    
    section = CrossSection(name, dist, t, z, cut_line=cutline, banks=banks,
                           description=description)
    
    # Skip blank line
    line = geomIO.readline().strip()
    #print "line=[%s]" % line
    
    if not real_section and self.parse_only_real_sections:
      section = None
      
    return section
    
    
  def __parse_xy_tuples__(self, geomIO, n, maxsize=16):
    
    # Parse number x,y tuples
    x = []
    y = []
    while len(x) < n:
      line = geomIO.readline().rstrip()
      
      count = int(math.ceil(len(line) / float(maxsize)))
      words = []
      for i in range(0, count):
        words.append(line[i*maxsize:i*maxsize+maxsize])
        
      ##print "test:", len(x), n, len(x) < n
      
      #words0 = line.split()
      #words = []
      #for word in words0:
        #if len(word) > maxsize:
          #word = " "+word
          #for i in range(0, int(math.ceil(len(word)/float(maxsize)))):
            #words.append(word[i*maxsize:i*maxsize+maxsize])
        #else:
          #words.append(word)
          
      for i in range(0, len(words)//2):
        x.append(float(words[i*2]))
        y.append(float(words[i*2+1]))
        
      #if maxsize == 8:
        #print line, len(x), n
        
    return np.array(x), np.array(y)

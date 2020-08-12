#!/usr/bin/env python3
#===---------------------------------------------------------------------------------------------===
# reentry.py - Simple euler solver for spacecraft reentries.
#
# Created on 2019-03-07 by Amy Parent <amy@amyparent.com>
# Copyright (c) 2018 Amy Parent
# Licensed under the MIT License
# =^•.•^=
#===---------------------------------------------------------------------------------------------===
import sys
import json
import numpy as np
from matplotlib import pyplot as plt

class Planet(object):
  """ Base class for planetary bodies.
  
  allows us to define more complex atmospheric models down the road, and to do coordinate system
  conversions.
  """
  radius = 0
  grav_k = 0
  name = ""
  
  def __init__(self, radius, grav_k, name):
    self.radius = radius
    self.grav_k = grav_k
    self.name = name
  
  def gravity(self, r):
    return  -(self.grav_k / (r**2.0))
  
  def altitude(self, pos):
    return np.linalg.norm(pos) - self.radius
  
  def cartesian(self, lat, alt):
    r = alt + self.radius
    return (r * np.cos(lat), r * np.sin(lat))

  def polar(self, x, y):
    alt = np.sqrt(x**2 + y**2) - self.radius
    return (np.arctan2(x, y), alt)
    
  def density(self, alt):
    return 0

# Scale Height is defined for T=15C (288.15K)
class Earth(Planet):
  def __init__(self):
    Planet.__init__(self, 6371e3, 3.986004418e14, "Earth")
  
  def density(self, alt):
    return 1.221 * np.exp(- alt / 8.43e3)

class Mars(Planet):
  def __init__(self):
    Planet.__init__(self, 3389e3, 4.282837e13, "Mars")
    
  def density(self, alt):
    return 0.015 * np.exp(- alt / 12e3)

def sim_run(sim, planet, craft):
  """ Calculates the trajectory & loads of a craft reentering atmosphere on a planetary body
  
  Uses a **very** basic Euler integrator and would probably benefit from a better solver.
  However, the degree of simplification of the equation system makes it unlikely the work would
  be worth it, given we probably have larger-magnitude errors induced by the model
  
  The model equations are simple:
  
   - g = (mu * mass) / r^2      -> towards body centre
   - Drag = 0.5*rho*v^2 / Beta  -> opposite velocity vector
   - Lift = L/D * Drag          -> normal to velocity vector
  """
  # Initialise our position, velocity, and acceleration components
  max_it = sim['max_it']
  dt = sim['delta_t']
  fpa = np.radians(sim['fpa'])
  
  p = np.array([0, sim['entry_interface'] + planet.radius])
  x = np.zeros(max_it)
  y = np.zeros(max_it)
  
  v = sim['velocity'] * np.array([np.cos(fpa), np.sin(fpa)])
  vx = np.zeros(max_it)
  vy = np.zeros(max_it)
  
  a = np.array([0, 0])
  ax = np.zeros(max_it)
  ay = np.zeros(max_it)
  
  t = np.arange(0, max_it * dt, dt)
  
  beta = craft['ballistic_coef']
  ld = craft['lift_drag']
  
  # Very basic Euler integrator, running until we reach ~10km
  # (doesn't take parachute into acount -- decent would slow down then)
  k = 0
  for _ in range(0, max_it):
    p = p + v * dt
    x[k], y[k] = p
    
    r = np.linalg.norm(p)
    rho = planet.density(planet.altitude(p))
    v_mag = np.linalg.norm(v)
    normal = np.array([v[1], v[0]])
    
    # Calculate drag and lift, this is pretty basic continuum equations
    # Because we use the ballistic coefficient, this isn't really a force but the acceleration
    # cause by aerodynamic forces
    aero_accel = 0.5*rho*v_mag * (ld*normal/beta - v/beta)
    gravity_accel = planet.gravity(r) * (p/r)
    
    a =  aero_accel + gravity_accel
    ax[k], ay[k] = a
    
    v = v + a * dt
    vx[k], vy[k] = v
    
    k += 1
    if planet.altitude(p) <= sim['stop_alt']:
      print('done in %d iterations' % k)
      break
      
  return (
    np.resize(x, k),
    np.resize(y, k),
    np.resize(vx, k),
    np.resize(vy, k),
    np.resize(ax, k),
    np.resize(ay, k),
    np.resize(t, k)
  )

def do_plot(xlabel, x, ylabel, y, label, title, fname):
  """ Basic utility function to simplify plotting
  """
  plt.plot(x, y, label=label)
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.legend()
  plt.tight_layout()
  plt.show()
  #plt.savefig(fname, dpi=300)

if __name__ == '__main__':
  bodies = {
    "Earth": Earth(),
    "Mars": Mars()
  }

  if len(sys.argv) != 2:
    print('usage: %s parameter_file' % sys.argv[0])
    exit(1)

  with open(sys.argv[1], 'r') as f:
    params = json.load(f)

  craft = params['craft']
  planet = bodies[params['planet']]
  sim = params['sim']

  # Run the simulation
  x, y, vx, vy, ax, ay, t = sim_run(sim, planet, craft)

  title = r'%s ($\beta=%.0f\mathrm{kg}/\mathrm{m}^2$) @ %s' % (craft['name'], craft['ballistic_coef'], planet.name)
  label = r'$fpa=%.2f, v_\mathrm{EI}=%.0f\mathrm{ms}^{-1}$' % (sim['fpa'], sim['velocity'])

  # convert to useful cooridnates
  lat, alt = planet.polar(x, y)
  downrange = (lat) * planet.radius
  vn = np.linalg.norm([vx, vy])

  # Compute the velocity magnitude
  v = np.sqrt(vx**2.0 + vy**2.0)

  # Get the axial load ((ax,ay) projected onto (vx,vy))
  aa = np.abs((ax*vx + ay*vy)/v)

  # Time and distance to go
  tti = np.max(t) - t
  dtg = np.max(downrange) - downrange


  f1 = plt.figure(figsize=(6, 2))
  do_plot(
    'downrange (km)', downrange/1e3,
    'altitude (km)', alt/1e3,
    label, title, '%s-traj.png' % craft['name']
  
  )

  f2 = plt.figure(figsize=(4, 3))
  do_plot(
    'axial loads (g)', aa/9.81,
    'altitude (km)', alt/1e3,
    label, title, '%s-load_alt.png' % craft['name']
  
  )

  f3 = plt.figure(figsize=(4, 3))
  do_plot(
    'time since EI (s)', t,
    'axial loads (g)', aa/9.81,
    label, title, '%s-load_time.png' % craft['name']
  )

  f4 = plt.figure(figsize=(4, 3))
  do_plot(
    'distance to splashdown (km)', dtg/1e3,
    'time to parachute deploy (s)', tti/60.0,
    label, title, '%s-dtg.png' % craft['name']
  )

  f5 = plt.figure(figsize=(4, 3))
  do_plot(
    'velocity (km/s)', v/1e3,
    'altitude (km)', alt/1e3,
    label, title, '%s-vel.png' % craft['name']
  )
  
  plt.close()
  

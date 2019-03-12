#!/usr/bin/env python3
from pylab import *
import py_param_phi

sphere_radius=1
dz = 10
Cx, Cy = .5, .8660254

#figure(1)
#for Cx, Cy, color in [(.5, .8660254, 'red'), (.85, .9660254, 'blue')]:
#    f = lambda pphi, ptheta: py_param_phi.m_py_param_phi.py_iterative_phi_theta_from_param_phi_and_param_theta(sphere_radius, dz, Cx, Cy, pphi, ptheta)
#
#    for ptheta in [-1,0,1]:
#        phi,theta=np.array([ f(pphi, ptheta) for pphi in  linspace(-2,2,1e5)]).T
#        plot(rad2deg(theta), rad2deg(phi), color=color)
#
#    for pphi in [-2,-1,1,2]:
#        phi,theta=np.array([ f(pphi, ptheta) for ptheta in  linspace(-1,1,1e5)]).T
#        plot(rad2deg(theta), rad2deg(phi), color=color)

f = lambda pphi, ptheta: py_param_phi.m_py_param_phi.py_iterative_phi_theta_from_param_phi_and_param_theta(sphere_radius, dz, Cx, Cy, pphi, ptheta)

figure(num=2, figsize=(8,4))
clf()
Cx, Cy = .8, 0.9660254
for ptheta in [0,]:
    phi,theta=np.array([ f(pphi, ptheta) for pphi in  linspace(-2,2,1e5)]).T
    plot(rad2deg(theta), rad2deg(phi))

for pphi in [-1,-2,1,2]:
    phi,theta=np.array([ f(pphi, ptheta) for ptheta in  linspace(-1,1,1e5)]).T
    plot(rad2deg(theta), rad2deg(phi))


gca().set_title(r'$\rm \phi^*,\ \theta^*\ Phase\ Diagram\ for\ C=({:.2f}, {:.2f})$'.format(Cx, Cy))
gca().set_xlabel(r'${\rm zenith\ \theta\ [deg]}$', fontsize='large')
gca().set_ylabel(r'${\rm azimuth\ \phi\ [deg]}$', fontsize='large')

gca().set_xlim(0,90)
gca().set_ylim(-60,75)

plt.locator_params(axis='y', nbins=10)
plt.locator_params(axis='x', nbins=6)

props = dict(boxstyle='round', facecolor='gray', alpha=0.1)
txtprops = dict(fontsize=10, verticalalignment='top', bbox=props)

textstr = '\n'.join((
    r'Base:  D',
    r'Left:  D',
    r'Right: D',
    ))
gca().text(2, 30, textstr, **txtprops)

textstr = '\n'.join((
    r'Base:  S',
    r'Left:  D',
    r'Right: D',
    ))
gca().text(30, 30, textstr, **txtprops)

textstr = '\n'.join((
    r'Base:  S',
    r'Left:  D',
    r'Right: S',
    ))
gca().text(70, -30, textstr, **txtprops)

textstr = '\n'.join((
    r'Base:  S',
    r'Left:  S',
    r'Right: D',
    ))
gca().text(75, 61, textstr, **txtprops)

phi, theta = np.rad2deg(f(.5, 0))
gca().annotate(r'$\theta^* = 0$', xy=(theta, phi), xytext=(2,-40),
            bbox=props,
            arrowprops=dict(facecolor='black', arrowstyle = "->",))

phi, theta = np.rad2deg(f(1, .5))
gca().annotate(r'$\phi^* = 1$', xy=(theta, phi), xytext=(60,-0),
            bbox=props,
            arrowprops=dict(facecolor='black', arrowstyle = "->",))
phi, theta = np.rad2deg(f(-1, .5))
gca().annotate(r'$\phi^* = 2$', xy=(theta, phi), xytext=(60,15),
            bbox=props,
            arrowprops=dict(facecolor='black', arrowstyle = "->",))

savefig('param_phi_param_theta_phase_diagram.pdf', bbox_inches='tight')

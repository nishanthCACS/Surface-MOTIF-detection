# -*- coding: utf-8 -*-
"""
Created on %Dec-2017 or before

@author: %A.Nishanth C00294860
"""

'''  plot the resukts and save it
'''
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')

## The triangles in parameter space determine which x, y, z points are
## connected by an edge
##ax.plot_trisurf(x, y, z, triangles=tri.triangles, cmap=plt.cm.Spectral)
ax.plot_trisurf(coordinates[:,0],coordinates[:,1], coordinates[:,2], triangles=tri.simplices, cmap=plt.cm.Spectral)
plt.show()
fig.savefig('Only_convex_x_y_z.png')
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')

## The triangles in parameter space determine which x, y, z points are
## connected by an edge
##ax.plot_trisurf(x, y, z, triangles=tri.triangles, cmap=plt.cm.Spectral)
ax.plot_trisurf(coordinates[:,1],coordinates[:,2], coordinates[:,0], triangles=tri.simplices, cmap=plt.cm.Spectral)
plt.show()
fig.savefig('Only_convex_y_z_x.png')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')

## The triangles in parameter space determine which x, y, z points are
## connected by an edge
##ax.plot_trisurf(x, y, z, triangles=tri.triangles, cmap=plt.cm.Spectral)
ax.plot_trisurf(coordinates[:,1],coordinates[:,0], coordinates[:,2], triangles=tri.simplices, cmap=plt.cm.Spectral)
plt.show()
fig.savefig('Only_convex_y_x_z.png')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')

## The triangles in parameter space determine which x, y, z points are
## connected by an edge
##ax.plot_trisurf(x, y, z, triangles=tri.triangles, cmap=plt.cm.Spectral)
ax.plot_trisurf(coordinates[:,2],coordinates[:,1], coordinates[:,0], triangles=tri.simplices, cmap=plt.cm.Spectral)
plt.show()
fig.savefig('Only_convex_z_y_x.png')
#%%
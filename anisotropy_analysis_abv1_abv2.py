# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 14:05:34 2025

@author: cmilano
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
np.set_printoptions(precision=4, suppress=False)

# Description: This script demonstrates the implementation of anisotropic analysis
###########################################
#############INPUT SECTION#################
###########################################

#If you want to save a figure, add the command: plt.savefig('directory to figure/figure.png', dpi=300) under the desired figure

# Read data, skip the first line (header)
tau_data1 = np.loadtxt("./Example_data/Reynolds_stress_tensor_ABV1.txt", skiprows=1)

# Assign values
tau111, tau221, tau331, tau121, tau131, tau231 = tau_data1


# Read data, skip the first line (header)
tau_data2 = np.loadtxt("./Example_data/Reynolds_stress_tensor_ABV2.txt", skiprows=1)

# Assign values
tau112, tau222, tau332, tau122, tau132, tau232 = tau_data2

# Read data, skip the first line (header)
uvw_data1 = np.loadtxt("./Example_data/mean_velocity_ABV1.txt", skiprows=1)

# Assign values
u_mean1, v_mean1, w_mean1 = uvw_data1


# Read data, skip the first line (header)
uvw_data2 = np.loadtxt("./Example_data/mean_velocity_ABV2.txt", skiprows=1)

# Assign values
u_mean2, v_mean2, w_mean2 = uvw_data2


##########################################
##########################################
##########################################
##########################################
###########################################
#############OUTPUT SECTION#################
###########################################


# Construct symmetric 3x3 matrix
tau_matrix1 = np.array([
    [tau111, tau121, tau131],
    [tau121, tau221, tau231],
    [tau131, tau231, tau331]
])

print("Reynolds Stress tensor at ABV1:")
print(tau_matrix1)

# Construct symmetric 3x3 matrix
tau_matrix2 = np.array([
    [tau112, tau122, tau132],
    [tau122, tau222, tau232],
    [tau132, tau232, tau332]
])

print("Reynolds Stress tensor at ABV2:")
print(tau_matrix2)


vel_mean1= np.array([u_mean1, v_mean1, w_mean1])

vel_mean2= np.array([u_mean2, v_mean2, w_mean2])

print("Mean velocity at ABV1:")
print(", ".join([f"{v:.6e}" for v in vel_mean1]))

print("Mean velocity at ABV2:")
print(", ".join([f"{v:.6e}" for v in vel_mean2]))



#Reorders RS in a tensor
rij1 = np.array([tau111, tau121, tau131, 
                      tau121, tau221, tau231,
                      tau131, tau231, tau331])

rij1 = rij1.reshape(3, 3)

rij2 = np.array([tau112, tau122, tau132, 
                      tau122, tau222, tau232,
                      tau132, tau232, tau332])

rij2 = rij2.reshape(3, 3)


#Calculates k
k1 = 0.5*(tau111 + tau221 + tau331)
k2 = 0.5*(tau112 + tau222 + tau332)

print("TKE at ABV1:")
print(k1)

print("TKE at ABV2:")
print(k2)

#Calculates anisotropic RS tensor
bij1 = (rij1 - 2 * k1 * np.eye(3, 3) / 3) / (2 * k1)
bij2 = (rij2 - 2 * k2 * np.eye(3, 3) / 3) / (2 * k2)

# Extract upper triangle (including diagonal)
upper_bij1 = bij1[np.triu_indices(3)]
upper_bij2 = bij2[np.triu_indices(3)]

# Write to text file with header
with open("Anisotropic_reynolds_stress_tensor_ABV1.txt", "w") as f:
    f.write("# Elements of the Anisotropic_reynolds_stress_tensor_ABV1\n")
    f.write("# Order: b11, b12, b13, b22, b23, b33\n")
    for value in upper_bij1:
        f.write(f"{value:.6e}\t")


# Write to text file with header
with open("Anisotropic_reynolds_stress_tensor_ABV2.txt", "w") as f:
    f.write("# Elements of the Anisotropic_reynolds_stress_tensor_ABV2\n")
    f.write("# Order: b11, b12, b13, b22, b23, b33\n")
    for value in upper_bij1:
        f.write(f"{value:.6e}\t")
        
#Calculates eigenvalues and eigenvectors of the RS and anisotropic RS
lambda_rij1, eigvecs_rij1 = np.linalg.eig(rij1)
lambda_bij1, eigvecs_bij1 = np.linalg.eig(bij1)

lambda_rij2, eigvecs_rij2 = np.linalg.eig(rij2)
lambda_bij2, eigvecs_bij2 = np.linalg.eig(bij2)

#Orders eigenvalues from largest to smallest and prints their values
lambda_rij_sorted1 = np.sort(lambda_rij1)[::-1]
lambda_bij_sorted1 = -1/3 + lambda_rij1/(2*k1)

lambda_rij_sorted2 = np.sort(lambda_rij2)[::-1]
lambda_bij_sorted2 = -1/3 + lambda_rij2/(2*k2)


# Write to text file with header
with open("Reynolds_stress_Eigenvalues_ABV1.txt", "w") as f:
    f.write("# Eigevalues of the Reynolds_stress tensor at ABV1\n")
    f.write("# Order: lambdau1, lambdau2, lambdau3\n")
    for value in lambda_rij_sorted1:
        f.write(f"{value:.6e}\t")


# Write to text file with header
with open("Reynolds_stress_Eigenvalues_ABV2.txt", "w") as f:
    f.write("# Eigevalues of the Reynolds_stress tensor at ABV2\n")
    f.write("# Order: lambdau1, lambdau2, lambdau3\n")
    for value in lambda_rij_sorted2:
        f.write(f"{value:.6e}\t")


# Write to text file with header
with open("Anisotropic_Reynolds_stress_Eigenvalues_ABV1.txt", "w") as f:
    f.write("# Eigevalues of the Anisotropic Reynolds_stress tensor at ABV1\n")
    f.write("# Order: lambdab1, lambdab2, lambdab3\n")
    for value in lambda_bij_sorted1:
        f.write(f"{value:.6e}\t")


# Write to text file with header
with open("Anisotropic_Reynolds_stress_Eigenvalues_ABV2.txt", "w") as f:
    f.write("# Eigevalues of the Anisotropic Reynolds_stress tensor at ABV2\n")
    f.write("# Order: lambdab1, lambdab2, lambdab3\n")
    for value in lambda_bij_sorted2:
        f.write(f"{value:.6e}\t")
        
#Calculates invariants for Lumley triangle
eta1 = np.sqrt((lambda_bij_sorted1[0]**2 + lambda_bij_sorted1[1]**2 +lambda_bij_sorted1[0]*lambda_bij_sorted1[1])/3)
eta2 = np.sqrt((lambda_bij_sorted2[0]**2 + lambda_bij_sorted2[1]**2 +lambda_bij_sorted2[0]*lambda_bij_sorted2[1])/3)


csi1 = np.cbrt(-0.5*lambda_bij_sorted1[0]*lambda_bij_sorted1[1]*(lambda_bij_sorted1[0]+lambda_bij_sorted1[1]))
csi2 = np.cbrt(-0.5*lambda_bij_sorted2[0]*lambda_bij_sorted2[1]*(lambda_bij_sorted2[0]+lambda_bij_sorted2[1]))


pi_b1 = -3*eta1**2
tri_b1= 2*csi1**3

pi_b2 = -3*eta2**2
tri_b2= 2*csi2**3

inv1 = np.array([eta1, csi1, pi_b1, tri_b1])
inv2 = np.array([eta2, csi2, pi_b2, tri_b2])

# Write to text file with header
with open("Invariants_ABV1.txt", "w") as f:
    f.write("# Invariants of the Reynols stress tensor at ABV1\n")
    f.write("# Order: eta, csi, II_b, III_b\n")
    for value in inv1:
        f.write(f"{value:.6e}\t")


# Write to text file with header
with open("Invariants_ABV2.txt", "w") as f:
    f.write("# Invariants of the Reynols stress tensor at ABV2\n")
    f.write("# Order: eta, csi, II_b, III_b\n")
    for value in inv2:
        f.write(f"{value:.6e}\t")

###########################################
###########################################
# Lumley triangle
# Description: This script demonstrates the plots of the Lumley triangle 

#Define axes limits for Lumley triangle
xneg = np.linspace(-1/6,0,100)
xpos = np.linspace(1/3,0,100)
xtot = np.linspace(-1/6, 1/3, 100)

yneg=-xneg
ypos = xpos
ytot = np.sqrt(1/27 +2*xtot**3)

#Plot Lumley triangle
plt.figure(figsize=(8,7))
plt.plot(xneg,yneg, 'k')
plt.plot(xpos,ypos, 'k')
plt.plot(xtot, ytot, 'k')
plt.plot(csi1, eta1, 'ro', label='ABV1')
plt.plot(csi2, eta2, 'go', label='ABV2')
plt.xlabel(r'$\xi$', fontsize=14)
plt.ylabel(r'$\eta$', fontsize=14)
plt.rcParams['figure.dpi'] = 300
plt.title('Lumley triangle')
plt.legend(fontsize=14)
plt.savefig(r'./lumley.png', dpi=300)


print("Lumley triangle correctly generated.")

#Define axes limits for anisotropic invariant map (AIM)

xneg = np.linspace(-1/108,0,100)
xpos = np.linspace(0,2/27,100)
xtot = np.linspace(-1/108,2/27,100)



yneg = 3*(np.abs(xneg)/2)**(2/3)
ypos= 3*(np.abs(xpos)/2)**(2/3)
ytot = 3*xtot + 1/9


plt.figure(figsize=(6,5))

plt.plot(xneg, yneg, 'k')
plt.plot(xpos, ypos, 'k')
plt.plot(xtot, ytot, 'k')
plt.plot(tri_b1, -pi_b1, 'ro', label='ABV1')
plt.plot(tri_b2, -pi_b2, 'go', label='ABV2')
plt.xlim(-0.02,0.08)
plt.ylim(-0.1,0.4)
plt.xlabel('$III_b$')
plt.ylabel('$-II_b$')
plt.title('Anisotropic Invariant Map (AIM)')
plt.legend(fontsize=14)
plt.savefig(r'./aim.png', dpi=300)

print("Anisotropic invariant map correctly generated.")


#Define eigenvectors of anisotropic RS tensor
xr11 = eigvecs_bij1[:,1]
xr21 = eigvecs_bij1[:,2]
xr31 = eigvecs_bij1[:,0]


#Normalize the velocity vector
vel_mean1 = vel_mean1/np.linalg.norm(vel_mean1)

xr12 = eigvecs_bij2[:,1]
xr22 = eigvecs_bij2[:,2]
xr32 = eigvecs_bij2[:,0]


#Normalize the velocity vector
vel_mean2 = vel_mean2/np.linalg.norm(vel_mean2)
# # Calculate the dot product between the anisotropic RS tensor eigenvectors and the velocity vector
# dot_product1 = np.dot(xr1, vel_mean)


# angle1 = (180/np.pi)*math.acos(dot_product1)
# if angle1 > 90:
#     angle1 = 180 - angle1

# print('Angle between eigenvector 1 and velocity vector: ', angle1)

# dot_product2 = np.dot(xr2, vel_mean)

# angle2 = (180/np.pi)*math.acos(dot_product2)
# if angle2 > 90:
#     angle2 = 180 - angle2

# print('Angle between eigenvector 1 and velocity vector: ', angle2)

# dot_product3 = np.dot(xr3, vel_mean)

# angle3 = (180/np.pi)*math.acos(dot_product3)
# if angle3 > 90:
#     angle3 = 180 - angle3

# print('Angle between eigenvector 1 and velocity vector: ', angle3)




# Data for full RS matrix not normalized
categories = [r'$\overline{uu}$', r'$\overline{vv}$', r'$\overline{ww}$', r'$\overline{uv}$', r'$\overline{uw}$', r'$\overline{vw}$']

# FDPTV_values = 1.531*1.531*np.array([0.001493, 0.003164, 0.003423, -0.001198, 0.0005082, -0.0006702])
RS_values1 = np.array([tau111, tau221, tau331, tau121, tau131, tau231])
RS_values2 = np.array([tau112, tau222, tau332, tau122, tau132, tau232])

# DES_values = 1.531*1.531*np.array([0.001182, 0.00386, 0.006548, -0.0009414, 0.001784, 0.001103]) old values


# Set the width of the bars and additional offset
bar_width = 0.35
bar_offset = 0.1

plt.figure()
# Create the bar chart
index = 1.35*np.arange(len(categories))
plt.bar(index, RS_values1, bar_width)


# Add spaces after every 3 bars
plt.xticks(index + bar_width , categories)

# Customize the chart
plt.title('Reynolds stess tensor components: ABV1')
plt.ylabel(r'Magnitude $[m^2/s^2]$')
plt.grid()
plt.ylim(-0.005,0.025)
plt.savefig(r'./Reynolds_stress_abv1.png', dpi=300)

print("Reynolds stresses plot for ABV1 correctly generated.")


plt.figure()
# Create the bar chart
index = 1.35*np.arange(len(categories))
plt.bar(index, RS_values2, bar_width)


# Add spaces after every 3 bars
plt.xticks(index + bar_width , categories)

# Customize the chart
plt.title('Reynolds stess tensor components: ABV2')
plt.ylabel(r'Magnitude $[m^2/s^2]$')
plt.grid()
plt.ylim(-0.005,0.025)

plt.savefig(r'./Reynolds_stress_abv2.png', dpi=300)

print("Reynolds stresses plot for ABV2 correctly generated.")


####################################################################
# Data for full bij matrix
categories = [r'${b_{11}}$', r'${b_{22}}$',r'${b_{33}}$',r'${b_{12}}$',r'${b_{13}}$',r'${b_{23}}$',]

Anisotropic_RS_values1 = np.array([bij1[0,0], bij1[1,1], bij1[2,2], bij1[0,1], bij1[0,2], bij1[1,2]])
Anisotropic_RS_values2 = np.array([bij2[0,0], bij2[1,1], bij2[2,2], bij2[0,1], bij2[0,2], bij2[1,2]])


# Set the width of the bars
bar_width = 0.35
plt.figure()

# Create the bar chart
index = 1.35*np.arange(len(categories))
plt.figure()
plt.bar(index, Anisotropic_RS_values1, bar_width)


plt.ylabel('Magnitude $[m^2/s^2]$', fontsize=12)

plt.xticks(index + bar_width , categories)

plt.grid()
plt.title('Anisotropic Reynolds stress tensor components: ABV1')
plt.ylim(-0.14,0.14)

plt.savefig(r'./Anisotropic_Reynolds_Stress_ABV1.png', dpi=300)


print("Anisotropic Reynolds stresses plot for ABV1 correctly generated.")


plt.figure()
plt.bar(index, Anisotropic_RS_values2, bar_width)


plt.ylabel('Magnitude $[m^2/s^2]$', fontsize=12)

plt.xticks(index + bar_width , categories)

plt.grid()
plt.ylim(-0.14,0.14)

plt.title('Anisotropic Reynolds stress tensor components: ABV2')
plt.savefig(r'./Anisotropic_Reynolds_Stress_ABV2.png', dpi=300)

print("Anisotropic Reynolds stresses plot for ABV2 correctly generated.")



####################################################################
# Data for diagonal RS not normalized
categories = [r'$\lambda_{u_1}$', r'$\lambda_{u_2}$',r'$\lambda_{u_3}$']

Diagonal_RS_values1 = lambda_rij_sorted1
Diagonal_RS_values2 = lambda_rij_sorted2


# # Set the width of the bars

bar_width = 0.35

# Create the bar chart
index = 1.35*np.arange(len(categories))
plt.figure()
plt.bar(index, Diagonal_RS_values1, bar_width, label='RANS EARSM ABV1 S4')


# Customize the chart
plt.title('Reynolds stress tensor component in principal axes: ABV1')
plt.ylabel('Magnitude $[m^2/s^2]$')
plt.xticks(index , categories, fontsize=14)
plt.ylim(0.0,0.03)

plt.grid()
plt.tight_layout()  # Adjust the layout
plt.savefig(r'./Reynolds_stress_Eigenvalues_abv1.png', dpi=300)
# Display the chart

print("Reynolds stresses eigenvalues plot for ABV1 correctly generated.")

plt.figure()
plt.bar(index, Diagonal_RS_values2, bar_width, label='RANS EARSM ABV1 S4')


# Customize the chart
plt.title('Reynolds stress tensor component in principal axes: ABV2')
plt.ylabel('Magnitude $[m^2/s^2]$')
plt.xticks(index , categories, fontsize=14)
plt.ylim(0.0,0.03)

plt.grid()
plt.tight_layout()  # Adjust the layout
plt.savefig(r'./Reynolds_stress_Eigenvalues_abv2.png', dpi=300)
# Display the chart
print("Reynolds stresses eigenvalues plot for ABV2 correctly generated.")

#######################################################################################
# Data for the bar chart principal axes bij
categories = [r'$\lambda_{b_1}$', r'$\lambda_{b_2}$',r'$\lambda_{b_3}$']

Diagonal_anisotropic_RS_values1 = lambda_bij_sorted1
Diagonal_anisotropic_RS_values2 = lambda_bij_sorted2


# Set the width of the bars
bar_width = 0.35
plt.figure()

# Create the bar chart
index = 1.35*np.arange(len(categories))
plt.bar(index, Diagonal_anisotropic_RS_values1, bar_width)

# Customize the chart
plt.title(r'${b_{ij}}$ in principal axes: ABV1')
# plt.xlabel('Categories')
plt.ylabel('Magnitude')
plt.xticks(index , categories, fontsize=14)

plt.ylim(-0.2,0.2)

plt.grid()
plt.tight_layout()  # Adjust the layout
plt.savefig(r'./Anisotropic_Reynolds_stress_Eigenvalues_abv1.png', dpi=300)

print("Anisotropic Reynolds stresses eigenvalues plot for ABV1 correctly generated.")


plt.figure()

# Create the bar chart
index = 1.35*np.arange(len(categories))
plt.bar(index, Diagonal_anisotropic_RS_values2, bar_width)

# Customize the chart
plt.ylim(-0.2,0.2)

plt.title(r'${b_{ij}}$ in principal axes: ABV2')
# plt.xlabel('Categories')
plt.ylabel('Magnitude')
plt.xticks(index , categories, fontsize=14)


plt.grid()
plt.tight_layout()  # Adjust the layout
plt.savefig(r'./Anisotropic_Reynolds_stress_Eigenvalues_abv2.png', dpi=300)


print("Anisotropic Reynolds stresses eigenvalues plot for ABV2 correctly generated.")




#RS ellipsoid

# Define ellipsoid parameters
center = [u_mean1,v_mean1,w_mean1]  # coordinates of the center
center = center
radii = np.sqrt(lambda_rij_sorted1)  # semi-axes lengths (a, b, c)
end = np.array([2*u_mean1,2*v_mean1,2*w_mean1]) #end of the velocity vector
	

# Define principal axes vectors
v1 = eigvecs_rij1[:,1]  # first principal axis
v2 = eigvecs_rij1[:,2]   # second principal axis
v3 = eigvecs_rij1[:,0]   # third principal axis

# Create ellipsoid meshgrid
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 50)
x = radii[0] * np.outer(np.cos(u), np.sin(v))
y = radii[1] * np.outer(np.sin(u), np.sin(v))
z = radii[2] * np.outer(np.ones(np.size(u)), np.cos(v))

# Define transformation matrix
T = np.column_stack((v1, v2, v3))

# Transform the ellipsoid
ellipsoid_points = np.stack([x.flatten(), y.flatten(), z.flatten()], axis=0)
transformed_points = np.dot(T, ellipsoid_points)

# Extract the transformed coordinates
x_transformed = transformed_points[0, :].reshape(x.shape)
y_transformed = transformed_points[1, :].reshape(y.shape)
z_transformed = transformed_points[2, :].reshape(z.shape)

# Translate the transformed ellipsoid
x_transformed += center[0]
y_transformed += center[1]
z_transformed += center[2]

# Normalize z_transformed to range from 0 to 1
z_norm = (z_transformed - np.min(z_transformed)) / (np.max(z_transformed) - np.min(z_transformed))

# Plot the transformed ellipsoid
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x_transformed, y_transformed, z_transformed, rstride=1, cstride=1,
                facecolors=plt.cm.jet(z_norm), alpha=0.6, linewidth=0.2, edgecolors='k', shade=False)


# Plot the principal axes vectors
#You might need to adjust k1,k2,k3 to fit the vectors in the plot area
k1 = 5
k2 = 6
k3 = 3
origin = np.array(center)
ax.quiver(*origin, *v1/k1, arrow_length_ratio=0.13,color='r', label='Principal Axis 1')
ax.quiver(*origin, *v2/k2,arrow_length_ratio=0.13, color='g', label='Principal Axis 2')
ax.quiver(*origin, *v3/k2,arrow_length_ratio=0.13, color='b', label='Principal Axis 3')
ax.quiver(*origin, *end/k3,arrow_length_ratio=0.13, color='purple', label='Velocity vector')
ax.set_xlim(u_mean1-0.1,u_mean1 +0.1)
ax.set_ylim(v_mean1-0.1,v_mean1+ 0.1)
ax.set_zlim(w_mean1-0.1, w_mean1+0.1)

# Set plot labels and title
ax.set_xlabel('U')
ax.set_ylabel('V')
ax.set_zlabel('W')
ax.set_title('Reynolds stress ellipsoid: ABV1')



# Set the view angles, might need to be adjustes to have optimal view.
ax.view_init(elev=20, azim=-20)

plt.legend()
plt.savefig(r'./Reynolds_stress_ellipsoid_abv1.png', dpi=300)

print("Reynolds stresses ellipsoid for ABV1 correctly generated.")




#RS ellipsoid

# Define ellipsoid parameters
center = [u_mean2,v_mean2,w_mean2]  # coordinates of the center
center = center
radii = np.sqrt(lambda_rij_sorted2)  # semi-axes lengths (a, b, c)
end = np.array([2*u_mean2,2*v_mean2,2*w_mean2]) #end of the velocity vector
	

# Define principal axes vectors
v1 = eigvecs_rij2[:,1]  # first principal axis
v2 = eigvecs_rij2[:,2]   # second principal axis
v3 = eigvecs_rij2[:,0]   # third principal axis

# Create ellipsoid meshgrid
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 50)
x = radii[0] * np.outer(np.cos(u), np.sin(v))
y = radii[1] * np.outer(np.sin(u), np.sin(v))
z = radii[2] * np.outer(np.ones(np.size(u)), np.cos(v))

# Define transformation matrix
T = np.column_stack((v1, v2, v3))

# Transform the ellipsoid
ellipsoid_points = np.stack([x.flatten(), y.flatten(), z.flatten()], axis=0)
transformed_points = np.dot(T, ellipsoid_points)

# Extract the transformed coordinates
x_transformed = transformed_points[0, :].reshape(x.shape)
y_transformed = transformed_points[1, :].reshape(y.shape)
z_transformed = transformed_points[2, :].reshape(z.shape)

# Translate the transformed ellipsoid
x_transformed += center[0]
y_transformed += center[1]
z_transformed += center[2]

# Normalize z_transformed to range from 0 to 1
z_norm = (z_transformed - np.min(z_transformed)) / (np.max(z_transformed) - np.min(z_transformed))

# Plot the transformed ellipsoid
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x_transformed, y_transformed, z_transformed, rstride=1, cstride=1,
                facecolors=plt.cm.jet(z_norm), alpha=0.6, linewidth=0.2, edgecolors='k', shade=False)


# Plot the principal axes vectors
#You might need to adjust k1,k2,k3 to fit the vectors in the plot area
k1 = 4
k2 = 4
k3 = 1
origin = np.array(center)
ax.quiver(*origin, *v1/k1, arrow_length_ratio=0.13,color='r', label='Principal Axis 1')
ax.quiver(*origin, *v2/k2,arrow_length_ratio=0.13, color='g', label='Principal Axis 2')
ax.quiver(*origin, *v3/k2,arrow_length_ratio=0.13, color='b', label='Principal Axis 3')
ax.quiver(*origin, *end/k3,arrow_length_ratio=0.13, color='purple', label='Velocity vector')
ax.set_xlim(u_mean2-0.15,u_mean2 +0.15)
ax.set_ylim(v_mean2-0.15,v_mean2+ 0.15)
ax.set_zlim(w_mean2-0.15, w_mean2+0.15)

# Set plot labels and title
ax.set_xlabel('U')
ax.set_ylabel('V')
ax.set_zlabel('W')
ax.set_title('Reynolds stress ellipsoid: ABV2')



# Set the view angles, might need to be adjustes to have optimal view.
ax.view_init(elev=20, azim=-20)

plt.legend()
plt.savefig(r'./Reynolds_stress_ellipsoid_abv2.png', dpi=300)

print("Reynolds stresses ellipsoid for ABV2 correctly generated.")


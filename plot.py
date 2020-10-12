import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math
import argparse

HALF_PI = 1.5707963267948966192313216916397514420985846996875529104874722961
yticks = [-0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
yticks_extended = [-1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
colors = ['r', 'm', 'k', 'b', 'w', 'c']

def parse_space_sep_file(filename):
	# res = []
	return np.genfromtxt(filename, delimiter=' ',
		names=['t', 'E_r', 'E_i', 'Var_r', 'Var_i'])

def plot_gates(num_gates):
	data = parse_space_sep_file('sigma_z.out')
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.set_title("Time evolution")    
	ax1.set_xlabel('time')
	ax1.set_ylabel('P_z')
	ax1.set_ylim(bottom=-0.1, top=1.1)
	ax1.set_xlim(left=0, right=num_gates*HALF_PI)
	ax1.set_yticks(yticks)

	# plot data
	ax1.plot(data['t'], data['E_r'], color='r', label='the data')
	# draw 0, 0.5, 1.0 probability lines
	ax1.plot(data['t'], len(data['t'])*[0.5], color='g', linestyle='dashed')
	ax1.plot(data['t'], len(data['t'])*[0.0], color='g', linestyle='dashed')
	ax1.plot(data['t'], len(data['t'])*[1.0], color='g', linestyle='dashed')
	# draw gate boundaries
	for i in range(1, num_gates+1):
		ax1.plot(len(yticks)*[i*HALF_PI], yticks, color='b', linestyle='dashed')
	plt.show()

def plot_bloch_sphere(num_gates, solid, plot_open):
	if plot_open:
		x = parse_space_sep_file('x_open')
		y = parse_space_sep_file('y_open')
		z = parse_space_sep_file('z_open')
	else:
		x = parse_space_sep_file('x')
		y = parse_space_sep_file('y')
		z = parse_space_sep_file('z')
	fig = plt.figure(figsize = (10,10))
	ax = fig.add_subplot(111, projection='3d')

	u = np.linspace(0, 2 * np.pi, 100)
	v = np.linspace(0, np.pi, 100)

	sphere_x = 1 * np.outer(np.cos(u), np.sin(v))
	sphere_y = 1 * np.outer(np.sin(u), np.sin(v))
	sphere_z = 1 * np.outer(np.ones(np.size(u)), np.cos(v))
	marker_size=100
	if solid:
		marker_size=500
		ax.plot_surface(sphere_x, sphere_y, sphere_z, rstride=4, cstride=4, color='w')
	else:
		ax.plot_wireframe(sphere_x, sphere_y, sphere_z, rstride=7, cstride=7, color='gray', linestyle='dashed', linewidth=0.5)
	ax.plot3D(x['E_r'][1:-1], y['E_r'][1:-1], z['E_r'][1:-1], 'g')
	# really messy
	handles = []
	sz = (len(x) - 1) / float(num_gates)
	for i in range(num_gates+1):
		idx = int(round(sz * i))
		handles.append(ax.scatter3D(x['E_r'][idx], y['E_r'][idx], z['E_r'][idx], s=marker_size, marker='D', color=colors[i], edgecolor='k'))
	# mark the center of the sphere
	ax.scatter3D(0, 0, 0, s=marker_size, marker='X', color='c', edgecolor='k')

	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	# plt.legend(handles=handles)
	plt.show()

# extract r, theta and phi from x, y, z since:
# x = sin(theta)*cos(phi),
# y = sin(theta)*sin(phi),
# z = cos(theta)
def plot_r_theta_phi(num_gates, superimpose_open, skip_unitary_plot):
	# parse and prepare data
	if not skip_unitary_plot:
		x = parse_space_sep_file('x')
		y = parse_space_sep_file('y')
		z = parse_space_sep_file('z')
		phi, theta, r = [], [], []
		for i in range(len(x['E_r'])):
			# r = sqrt(x^2+y^2+z^2)
			x_sq_y_sq = y['E_r'][i]**2 + x['E_r'][i]**2
			r.append(math.sqrt(z['E_r'][i]**2 + x_sq_y_sq))
			theta.append(HALF_PI - math.atan2(z['E_r'][i], math.sqrt(x_sq_y_sq)))
			phi.append(math.atan2(y['E_r'][i], x['E_r'][i]))

		theta = map(lambda x: x / (2*HALF_PI), theta)
		phi = map(lambda x: x / (2*HALF_PI), phi)

	# parse and prepare data for open system
	if superimpose_open:
		x = parse_space_sep_file('x_open')
		y = parse_space_sep_file('y_open')
		z = parse_space_sep_file('z_open')
		phi_open, theta_open, r_open = [], [], []
		for i in range(len(x['E_r'])):
			# r_open = sqrt(x^2+y^2+z^2)
			x_sq_y_sq = y['E_r'][i]**2 + x['E_r'][i]**2
			r_open.append(math.sqrt(z['E_r'][i]**2 + x_sq_y_sq))
			theta_open.append(HALF_PI - math.atan2(z['E_r'][i], math.sqrt(x_sq_y_sq)))
			phi_open.append(math.atan2(y['E_r'][i], x['E_r'][i]))

		theta_open = map(lambda x: x / (2*HALF_PI), theta_open)
		phi_open = map(lambda x: x / (2*HALF_PI), phi_open)

	# # plot data
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.set_title("Time evolution")
	ax1.set_xlabel('time')
	ax1.set_ylabel('normalized values')
	ax1.set_ylim(bottom=-0.1, top=1.1)
	ax1.set_xlim(left=0, right=num_gates*HALF_PI)
	ax1.set_yticks(yticks_extended)

	# plot data
	handles = []
	if not skip_unitary_plot:
		handles.append(ax1.plot(x['t'], theta, color='r', label='theta')[0])
		handles.append(ax1.plot(x['t'], phi, color='m', label='phi')[0])
		handles.append(ax1.plot(x['t'], r, color='k', label='r')[0])
	if superimpose_open:
		handles.append(ax1.plot(x['t'], theta_open, color='g', label='theta_open')[0])
		handles.append(ax1.plot(x['t'], phi_open, color='deepskyblue', label='phi_open')[0])
		handles.append(ax1.plot(x['t'], r_open, color='darkorange', label='r_open')[0])
	# draw -1.0, -0.5, 0, 0.5, 1.0
	ax1.plot(x['t'], len(x['t'])*[-1.0], color='g', linestyle='dashed', linewidth=0.3)
	ax1.plot(x['t'], len(x['t'])*[-0.5], color='g', linestyle='dashed', linewidth=0.3)
	ax1.plot(x['t'], len(x['t'])*[0.5], color='g', linestyle='dashed', linewidth=0.3)
	ax1.plot(x['t'], len(x['t'])*[0.0], color='g', linestyle='dashed', linewidth=0.3)
	ax1.plot(x['t'], len(x['t'])*[1.0], color='g', linestyle='dashed', linewidth=0.3)
	# draw gate boundaries
	for i in range(1, num_gates+1):
		ax1.plot(len(yticks_extended)*[i*HALF_PI], yticks_extended, color='b', linestyle='dashed', linewidth=0.3)
	plt.legend(handles=handles)
	plt.show()

def main(args):
	args = parser.parse_args()
	if args.mode == 'gate':
		plot_gates(int(args.num_gates))
	elif args.mode == 'tmp':
		plot_r_theta_phi(int(args.num_gates), args.plot_open=='1', args.skip_unitary_plot=='1')
	else:
		plot_bloch_sphere(int(args.num_gates), args.solid_sphere=='1', args.plot_open=='1')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--mode', default='bloch', help='bloch or gate')
	parser.add_argument('--solid_sphere', default='0', help='bloch sphere solid or wire, 1/0 (meaning True/False)')
	parser.add_argument('--num_gates', default='3', help='how many gates were applied?')
	parser.add_argument('--plot_open', default='0', help='should plot open system? 1/0 (meaning True/False)')
	parser.add_argument('--skip_unitary_plot', default='0', help='should plot unitary system? 1/0 (meaning True/False)')
	args = parser.parse_args()
	main(args)

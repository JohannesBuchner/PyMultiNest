
import numpy
from StringIO import StringIO
import re

class Analyzer(object):
	"""
		Class for accessing the output of MultiNest.
		
		After the run, this class provides a mean of reading the result
		in a structured way.
		
	"""
	def __init__(self, n_params, outputfiles_basename = "chains/1-"):
		self.outputfiles_basename = outputfiles_basename
		self.n_params = n_params
		"""
		[root].txt
		Compatable with getdist with 2+nPar columns. Columns have 
		sample probability, -2*loglikehood, samples. Sample
		probability is the sample prior mass multiplied by its 
		likelihood & normalized by the evidence.
		"""
		self.data_file = "%s.txt" % self.outputfiles_basename
		print '  analysing data from %s' % self.data_file

		"""[root]post_separate.dat
		This file is only created if mmodal is set to T. Posterior 
		samples for modes with local log-evidence value
		greater than nullZ, separated by 2 blank lines. Format is the 
		same as [root].txt file."""
		self.post_file = "%spost_seperate.dat" % self.outputfiles_basename

		"""[root]stats.dat
		Contains the global log-evidence, its error & local log-evidence 
		with error & parameter means & standard
		deviations as well as the  best fit & MAP parameters of each of 
		the mode found with local log-evidence > nullZ.
		"""
		self.stats_file = "%sstats.dat" % self.outputfiles_basename

		"""
		[root]post_equal_weights.dat
		Contains the equally weighted posterior samples
		"""
		self.equal_weighted_file = "%spost_equal_weights.dat" % self.outputfiles_basename
	
	def get_data(self):
		"""
			fetches self.data_file
		"""
		return numpy.loadtxt(self.data_file, ndmin=2)
	
	def _read_error_line(self, l):
		#print '_read_error_line', l
		name, values = l.split('  ', 1)
		name = name.strip(': ').strip()
		v, error = values.split(" +/- ")
		return name, float(v), float(error)
	def _read_error_into_dict(self, l, d):
		name, v, error = self._read_error_line(l)
		d[name.lower()] = v
		d['%s error' % name.lower()] = error
	def _read_table(self, txt, d = None, title = None):
		if title is None:
			title, table = txt.split("\n", 1)
		else:
			table = txt		
		header, table = txt.split("\n", 1)
		data = numpy.loadtxt(StringIO(table), ndmin=2)
		if d is not None:
			d[title.strip().lower()] = data
		return data

	def get_stats(self):
		"""
			information about the modes found:
			mean, sigma, maximum a posterior in each dimension
		"""
		lines = file(self.stats_file).readlines()
		text = "".join(lines)
		parts = text.split("\n\n\n")
		del parts[0]
		stats = {
			'modes':[]
		}
		# Global Evidence
		self._read_error_into_dict(lines[0], stats)
		i = 0
		for p in parts:
			modelines = p.split("\n\n")
			mode = {
				'index':i
			}
			i = i + 1
			modelines1 = modelines[0].split("\n")
			# Strictly local evidence
			self._read_error_into_dict(modelines1[1], mode)
			self._read_error_into_dict(modelines1[2], mode)
			t = self._read_table(modelines[1], title = "Parameters")
			mode['mean'] = t[:,1].tolist()
			mode['sigma'] = t[:,2].tolist()
			mode['maximum'] = self._read_table(modelines[1])[:,1].tolist()
			mode['maximum a posterior'] = self._read_table(modelines[1])[:,1].tolist()
			stats['modes'].append(mode)
		return stats

import scipy
import scipy.interpolate
import scipy.special
import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib.cm as cm
import itertools

class PlotMarginal(object):
	"""
		This class can be used to 
		plot marginal and conditional likelihoods.
		
		@param analyser: A Analyzer instance
	"""
	def __init__(self, analyser):
		self.analyser = analyser

	def plot_conditional(self, dim1, dim2 = None, 
		with_ellipses = True, with_points = True,
		only_interpolate = False, use_log_values = False,
		grid_points = 40, marginalization_type='sum'
	):
		"""
			Generate a conditional/marginal probability plot.
			 (marginalize all but two/one dimensions).
		
			@param dim1: first dimension to use
			@param dim2: second dimension to use (set to None for marginal plot)
		
			@param with_ellipses: Show ellipses for the resulting modes
		
			@param with_points: Show the sampled points in the plot
		
			@param only_interpolate: Use a interpolation of the points 
				instead of the points.
		
			@param use_log_values: Log-plot
		
			@param grid_points: how many bins the plot shall have
			
			@param marginalization_type: how should the "marginal" or "conditional" be calculated:
			       can be one of: 
			       
			       - **sum**  ... real marginalization, and the default
			       - **max**
			       - **mean**
			       
		"""
		data = self.analyser.get_data()
		stats = self.analyser.get_stats()
		n_params = self.analyser.n_params
		modes = stats['modes']
		# determining min/max
		min1 = min([mode['mean'][dim1] - 3*mode['sigma'][dim1] for mode in modes])
		max1 = max([mode['mean'][dim1] + 3*mode['sigma'][dim1] for mode in modes])
		if dim2 is not None:
			min2 = min([mode['mean'][dim2] - 3*mode['sigma'][dim2] for mode in modes])
			max2 = max([mode['mean'][dim2] + 3*mode['sigma'][dim2] for mode in modes])
		#print dim1, min1, max1, dim2, min2, max2
		# create grid
		n = grid_points
		if dim2 is not None:
			m = n
			grid_x, grid_y = numpy.mgrid[min1:max1:n*1j, min2:max2:n*1j]
			binsize2 = (max2 - min2) / m
		else:
			m = 1
			grid_x = numpy.mgrid[min1:max1:n*1j]
			grid_y = [0]
		
		binsize1 = (max1 - min1) / n
		
		dim1_column = data[:,2 + dim1]
		if dim2 is not None:
			dim2_column = data[:,2 + dim2]
			coords = numpy.array([dim1_column, dim2_column]).transpose()
		else:
			coords = dim1_column.transpose()
		values = data[:,0]
		if use_log_values:
			values = numpy.log(values)
		grid_z = numpy.zeros((n,m))
		minvalue = values.min()
		maxvalue = values.max()
		
		# for each grid item, find the matching points and put them in.
		for row, col in itertools.product(range(len(grid_x)), range(len(grid_y))):
			if dim2 is not None:
				xc = grid_x[row,col]
				here_x = numpy.abs(dim1_column - xc) < binsize1 / 2.
				yc = grid_y[row,col]
				here_y = numpy.abs(dim2_column - yc) < binsize2 / 2.
			else:
				xc = grid_x[row]
				here_x = numpy.abs(dim1_column - xc) < binsize1 / 2.
				here_y = True
			
			bin = values[numpy.logical_and(here_x, here_y)]
			if bin.size != 0:
				if marginalization_type == 'max':
					grid_z[row,col] = bin.max()
				elif marginalization_type == 'sum':
					grid_z[row,col] = bin.sum()
				elif marginalization_type == 'mean':
					grid_z[row,col] = bin.mean()
				elif marginalization_type == 'count':
					grid_z[row,col] = bin.size
				else:
					assert False, "marginalization_type should be mean, sum or max"
			else:
				grid_z[row,col] = minvalue
		
		#print 'maxima', values.max(), grid_z.max()
		# plot gridded data
		if only_interpolate:
		#   version A: interpolated -- may look weird because of the 
		#              loss of dimensions
			assert dim2 is not None
			grid_z = scipy.interpolate.griddata(coords, values, (grid_x, grid_y), method='cubic')
		
		if dim2 is not None:
			plt.xlim((min1,max1))
			plt.ylim((min2,max2))
			plt.imshow(grid_z.transpose(), origin='lower', aspect='auto',
				cmap=cm.gray_r, alpha = 0.8, extent=(min1,max1,min2,max2))
			plt.colorbar()
		else:
			#plt.xlim(min1, max1)
			plt.plot(grid_x, grid_z[:,0], '-', color='grey', drawstyle='steps')
		# add contours
		if use_log_values:
			levels = [maxvalue, maxvalue - .5, maxvalue - 1.0, maxvalue - 2.0]
		else:
			levels = [maxvalue, maxvalue / 3, maxvalue / 10, maxvalue / 100]
		leveltitles = ['max', 'max/3', 'max/10', 'max/100']
		
		if dim2 is not None:
			plt.contour(grid_x, grid_y, grid_z, levels, linewidths=0.5, colors='k')
		elif marginalization_type == 'max':
			# for the other types, these levels are pretty useless.
			for i in range(len(levels)):
				plt.plot([grid_x.min(), grid_x.max()], [levels[i]] * 2, '--', color='grey', label=leveltitles[i])
		
		# add points
		if with_points and dim2 is not None:
			plt.scatter(dim1_column, dim2_column, marker='+', color='black', s=1, alpha=0.4)
		# add ellipses
		i = 0
		for mode in modes:
			if not with_ellipses or dim2 is None:
				break
			el_xcenter = mode['mean'][dim1]
			el_xsize = mode['sigma'][dim1] * 2
			el_ycenter = mode['mean'][dim2]
			el_ysize = mode['sigma'][dim2] * 2
			
			ell = matplotlib.patches.Ellipse((el_xcenter, el_ycenter), 
				el_xsize, el_ysize, facecolor='#AAAAFF', 
				alpha=0.7, edgecolor='#3333AA')
			plt.axes().add_artist(ell)
			# xtext = xcenter + xsize / 2.
			# ytext = ycenter + ysize / 2.
			# put the text at 45 degrees or so
			ellipse_x = el_xcenter + el_xsize / 2. / 2. ** 0.5
			ellipse_y = el_ycenter + el_ysize / 2. / 2. ** 0.5
	
			textx = 20 * numpy.sin(numpy.pi / 8. * i)
			texty = 20 * numpy.cos(numpy.pi / 8. * i)
			text = 'Mode %d' % i
			arrow = dict(facecolor='blue', width=0.1, headwidth=0.9, 	
				linewidth=0.1)
			plt.annotate(text, xy=(ellipse_x, ellipse_y),
				xycoords='data', textcoords='offset points', size=4,
				xytext=(textx, texty),
				arrowprops=arrow, horizontalalignment='left', verticalalignment='bottom',
			)
			i = i + 1
		
	def plot_modes_marginal(self, dim1, cumulative = False, grid_points = 200, 
		with_ellipses = True, with_points = True):
		"""
			generate a marginal probability plot, visualizing the calculated
			modes.
		"""
		
		data = self.analyser.get_data()
		stats = self.analyser.get_stats()
		n_params = self.analyser.n_params
		modes = stats['modes']
		dim1_column = data[:,2 + dim1]
		values = data[:,0]
		min1 = min([mode['mean'][dim1] - 3*mode['sigma'][dim1] for mode in modes])
		max1 = max([mode['mean'][dim1] + 3*mode['sigma'][dim1] for mode in modes])
		x = numpy.linspace(min1, max1, grid_points)
		# 1-d. just need to plot
		# Z_mode / sqrt(2*pi) / sigma * exp(-0.5 * (x-u)^2 / sigma^2 )
		# sum{ Z_mode * 0.5 * [1 + erf( (x-u)/sigma/sqrt(2))]}
		f = numpy.sqrt(2 * numpy.pi)
		def norm(x, u, sigma):
			z = (x - u) / sigma
			return 1. / f / sigma * numpy.exp(-0.5 * z**2)

		def normcum(x, u, sigma):
			z = (x - u) / sigma
			return 0.5 * (1 + scipy.special.erf(z / numpy.sqrt(2)))

		pdff = lambda x: sum([m['strictly local evidence'] * norm(x,m['mean'][dim1],m['sigma'][dim1]) for m in modes])
		pdf = numpy.vectorize(pdff)(x)
		cdff = lambda x: sum([m['strictly local evidence'] * normcum(x,m['mean'][dim1],m['sigma'][dim1]) for m in modes])
		cdf = numpy.vectorize(cdff)(x)
		
		if cumulative:
			y = cdf
		else:
			y = pdf	
		plt.plot(x, y, '-')
		if with_points:
			plt.plot(dim1_column, values, 'x')
		
		# can add ellipses too:
		# x-mean, x-sigma is from mode parameters
		# y-height = Z_mode, y is cumulative (sum of previous Z_modes)
		cumy = 0
		i = 0
		for m in sorted(modes, key=lambda x:x['mean'][dim1]):
			if not with_ellipses:
				break
			el_xcenter = m['mean'][dim1]
			el_xsize = 2 * m['sigma'][dim1]
			el_ysize = m['strictly local evidence']
			if cumulative:
				# set to cdf step size
				y = cdff(el_xcenter)
				y = cumy
				cumy += el_ysize
			else:
				# set to pdf peak
				el_ysize /= numpy.sqrt(2 * numpy.pi) * m['sigma'][dim1]
				y = 0
			el_ycenter = el_ysize / 2. + y
			
			ell = matplotlib.patches.Ellipse((el_xcenter, el_ycenter), 
				el_xsize, el_ysize, facecolor='#AAAAFF', 
				alpha=0.7, edgecolor='#3333AA')
			plt.axes().add_artist(ell)
			# xtext = xcenter + xsize / 2.
			# ytext = ycenter + ysize / 2.
			# put the text at 45 degrees or so
			ellipse_x = el_xcenter + el_xsize / 2. / 2. ** 0.5
			ellipse_y = el_ycenter + el_ysize / 2. / 2. ** 0.5
	
			textx = 20 * numpy.sin(numpy.pi / 8. * i)
			texty = 20 * numpy.cos(numpy.pi / 8. * i)
			text = 'Mode %d' % m['index']

			arrow = dict(facecolor='blue', width=0.1, headwidth=0.9, 	
				linewidth=0.1)
			plt.annotate(text, xy=(ellipse_x, ellipse_y),
				xycoords='data', textcoords='offset points', size=4,
				xytext=(textx, texty),
				arrowprops=arrow, horizontalalignment='left',
				verticalalignment='bottom',
			)
			i = i + 1

		

	"""
	create n-dimensional interpolation using scipy.interpolate.griddata
	
	create 2D marginalization plots
	add modes as 2D-ellipses; colors should add up cumulatively
	add colourmap
	add contours -- sort cubes by value; add while evidence is < Jeffrey's Factors
	
	create 1D marginalization plots
	add mode parameter estimations; as (cumulative) bars
	
	"""

	"""
	marginal is the same as conditional, call it with dim2 = None for the 
	marginal behaviour.
	"""
	plot_marginal = plot_conditional

	


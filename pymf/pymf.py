# -*- coding: utf-8 -*-
# 
#  This file is part of PyMF.
# 
#  PyMF is free software; you can redistribute it and/or modify
#  it under the terms of the MIT License.
# 
#  PyMF is distributed in the hope that it will be useful,but 
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
#  the provided copy of the MIT License for more details.

import numpy as np


def corr2cov(corr, var):
	'''Converts a given correlation matrix into a covariance matrix 

	Parameters
	----------
	corr: 2D float array
		Correlation matrix of n components
	var: 1D float array
		Variance of the n components
		
	Returns
	-------
	corr: 2D float array
		Covariance matrix
	'''
    
	x_var, y_var = np.meshgrid(np.sqrt(var), np.sqrt(var))

	cov = corr * x_var * y_var
		    
	return(cov)


def rad_profile(image, pixel_size_arcmin, return_k=False):
	'''Computes the azimuthally-averaged radial profile of a given map.

	Parameters
	----------
	image: 2D float or complex array
		Input image.
	pixel_size_arcmin: float
		Pixel size in arcmin.
	return_k: bool, optional
		If set to True, the provided image is assumed to be a power spectrum
		and the x-coordinate of the returned data is converted to the 
		two-dimensional spatial frequency k. Default: None	

	Returns
	-------
	kradialprofile: float array
		x-coordinate of the radial profile. Either the radial separation 
		from the image center or spatial frequency k, depending on the 
		value of the variable return_k.
	radialprofile: float or complex array
		azimuthally averaged radial profile.

	''' 
	
	nxpix, nypix = image.shape[0], image.shape[1]

	YY, XX = np.indices((image.shape))
	r = np.sqrt((XX - nxpix//2)**2 + (YY - nypix//2.)**2)

	if return_k is True:
		k = 360/(pixel_size_arcmin/60.)*np.sqrt(((XX-nxpix//2)/nxpix)**2 + ((YY-nypix//2)/nypix)**2)
	else:
		k = np.copy(r)
		
	#r = np.round(r).astype(np.int)
	r_int = r.astype(np.int).ravel()
	k = k.ravel()
	image = image.ravel()

	uniq = np.unique(r_int.ravel())

	kradialprofile = np.zeros(len(uniq))
	radialprofile = np.zeros(len(uniq), dtype=image.dtype)

	for i in np.arange(len(uniq)):
		index = np.where(r_int == uniq[i])[0]
		kradialprofile[i] = np.mean(k[index])
		radialprofile[i] = np.mean(image[index])
				
	return(kradialprofile, radialprofile)


def power_spec(image, pixel_size_arcmin, return_k=False):
	'''Computes the azimuthally-averaged power spectrum of a given map.

	Parameters
	----------
	image: 2D float array
		Input image
	pixel_size_arcmin: float
		Pixel size in arcmin.
	return_k: bool, optional
		If set to True, the provided image is assumed to be a power spectrum
		and the x-coordinate of the returned data is converted to the 
		two-dimensional spatial frequency k. Default: None	

	Returns
	-------
	k: float array
		x-coordinate of the radial profile. Either the radial separation 
		from the image center or spatial frequency k, depending on the 
		value of the variable return_k.
	Pk: float or complex array
		azimuthally-averaged power spectrum
	''' 

	npix = image.shape[0] * image.shape[1]
	
	Fk = np.fft.fftshift((np.fft.fft2(image, norm=None))) / npix
	ps=(np.absolute((Fk))**2) 
	
	k, Pk = rad_profile(ps, pixel_size_arcmin, return_k=return_k)
	
	return(k, Pk)


def cross_spec(images, pixel_size_arcmin, return_k=False, global_cross=False):
	'''Computes the azimuthally-averaged cross-spectrum of a given set of maps.

	Parameters
	----------
	images: 3D float array
		Data cube containing the input images
	pixel_size_arcmin: float
		Pixel size in arcmin.
	return_k: bool, optional
		If set to True, the provided image is assumed to be a powerspectrum
		and the x-coordinate of the returned data is converted to the 
		two-dimensional spatial frequency k. Default: None
	global_cross: bool, optional
		If set to True, the cross power of the maps is computed globally
		(i.e. averaged over the full maps / all scales). If False, the 
		cross power is computed at each spatial frequency. Default: False

	Returns
	-------
	cc_maps: float array
		Data cube of dimensions n_nu * n_nu containing the azimuthally-averaged 
		cross-spectra maps.
	''' 
	
	nf = images.shape[0]
	npix = images.shape[1] * images.shape[2]

	cc_maps = np.array(np.zeros((nf,nf,npix)), dtype=np.complex128)

	if global_cross is False:

		for i in np.arange(nf):
			FK1 = np.fft.fftshift((np.fft.fft2(images[i,:,:], norm=None))) / npix
			for j in np.arange(nf):
				FK2 = np.fft.fftshift((np.fft.fft2(images[j,:,:], norm=None))) / npix
				cc = FK1*np.conj(FK2)
				k, Ck = rad_profile(cc, pixel_size_arcmin, return_k=return_k)
				cc_map = make_filter_map(cc, k, Ck)
				cc_maps[i,j,:] = cc_map.reshape(npix)

	else:

		corr = np.corrcoef(images.reshape(nf,npix))

		power = np.zeros((nf,npix))
		for i in np.arange(nf):
			k, Pk = power_spec(images[i,:,:],1)
			power[i,:] = make_filter_map(images[0,:,:], k, Pk).reshape(npix)
      
		for i in np.arange(npix):
			cc_maps[:,:,i] = corr2cov(corr, power[:,i])
		
	return(cc_maps)


def make_filter_map(reference_image, x, profile):
	'''Creates a 2D map from an azimuthally-averaged radial profile. 
       This is useful for the computing 2D window functions from profiles.

	Parameters
	----------
	reference_image: 2D float or complex array
		Input reference image. Only used to determine the desired image 
		dimensions.
	x: float array
		x-coordinate of the azimuthally-averaged radial profile.
		Radial separation from the image center in usints of pixels.
	window: float array
		Azimuthally averaged radial profile to be used to construct the 
		2D map	

	Returns
	-------
	image: 2D float or complex array
		2D map created from the provided azimuthally-averaged radial 
		profile.

	'''
	
	nxpix, nypix = reference_image.shape[0], reference_image.shape[1]

	YY, XX = np.indices((reference_image.shape))

	r = np.sqrt((XX - nxpix//2.)**2 + (YY - nypix//2.)**2)
	r=r.ravel()

	image = np.interp(r, x, profile)
	image = image.reshape(nxpix,nypix)

	return(image)


def filter_map_mf(image, source, noise_map=None, sigma_noise=None):
	'''Computes and applies a matched filter to a given map. The filter is build from a 
       provided source template and the power spectrum of the map.

	Parameters
	----------
	image: 2D float array
		Input image.
	source: 2D float array
		Source template. Needs to have identical dimensions as the provided image.
	noise_map: 2D float array, optional
		The provided noise map will be used to compute noise power spectrum. Otherwise
		the noise power spectrum is directly computed from the image. Default: None
	sigma_noise: float, optional
		Standard deviation of the map noise. Does allow analytic computation of the
		noise power spectrum. Useful for tests with white noise. Default: None

	Returns
	-------
	filtered_image: 2D float array
		Image convolved with the computed matched filter.
	filter_ft: 2D float array
		2D window function of the matched filter.
	noise: float
		Noise level of the filtered map.

	''' 

	#determine number of pixels
	npix = image.shape[0] * image.shape[1]

	#compute FFTs
	source_ft = np.fft.fftshift(np.fft.fft2(source)) / npix
	image_ft = np.fft.fftshift(np.fft.fft2(image)) / npix

	#compute noise power spectrum
	if noise_map is not None:
		k, noise_ps = power_spec(noise_map, 1)
	elif sigma_noise is not None:
		noise_ft = np.ones_like(image) * sigma_noise**2 / npix
		k, noise_ps = rad_profile(noise_ft, 1)
	else:
		k, noise_ps = power_spec(image, 1)

	ps_map = make_filter_map(image, k, noise_ps)

	#compute filter
	filter_ft = (source_ft/ps_map) / np.sum(source_ft**2/ps_map)

	#compute results
	filtered_image = np.fft.ifftshift(np.real(np.fft.ifft2(np.fft.ifftshift(filter_ft*image_ft)))) * npix
	noise = np.sqrt(np.real(np.sum(filter_ft**2*ps_map)))

	return(filtered_image, filter_ft, noise)


def filter_map_cmf(image, templates, response, noise_map=None, sigma_noise=None):
	'''Computes and applies a constrained matched filter to a given map. The filter is 
       build from a provided source template and the power spectrum of the map.

	Parameters
	----------
	image: 2D float array
		Input image.
	templates: 3D float array
		Data cube containing the individual source templates that will be used to
		construct the filter. The dimensions are supposed to be n * n_x * n_y, 
		where n is the number of templates and n_x and n_y are the number of pixels
		along the * and y axis of the images. The latter have to be identical to
		the dimensions of image.
	response: float array
		Array containing the desired response of the filter to each provided source 
        template.
	noise_map: 2D float array, optional
		The provided noise map will be used to compute noise power spectrum. Otherwise
		the noise power spectrum is directly computed from the image. Default: None
	sigma_noise: float, optional
		Standard deviation of the map noise. Does allow analytic computation of the
		noise power spectrum. Useful for tests with white noise. Default: None

	Returns
	-------
	filtered_image: 2D float array
		Image convolved with the computed constrained matched filter.
	filter_ft: 2D float array
		2D window function of the constrained matched filter.
	noise: float
		Noise level of the filtered map.

	'''

	#determine dimensions
	n_temp, nxpix, nypix= templates.shape[0], templates.shape[1], templates.shape[2]
	npix = nxpix * nypix

	#compute FFTs
	template_ft = np.empty((n_temp, nxpix, nypix), dtype=np.complex128)
	for i in np.arange(n_temp):
		template_ft[i,:,:] = np.fft.fftshift(np.fft.fft2(templates[i,:,:])) / npix
	
	image_ft = np.fft.fftshift(np.fft.fft2(image)) / npix

	#compute noise power spectrum
	if noise_map is not None:
		k, noise_ps = power_spec(noise_map, 1)
	elif sigma_noise is not None:
		noise_ft = np.ones_like(image) * sigma_noise**2 / npix
		k, noise_ps = rad_profile(noise_ft, 1)
	else:
		k, noise_ps = power_spec(image, 1)

	ps_map = make_filter_map(image, k, noise_ps)

	#compute filter
	norm_matrix = np.empty((n_temp,n_temp), dtype=np.complex128)
	for i in np.arange(n_temp):
		for j in np.arange(n_temp):
			norm_matrix[i,j] = np.sum(template_ft[i,:,:] / ps_map * template_ft[j,:,:])

	filter_matrix = abs(template_ft.reshape(n_temp, npix))
	for i in np.arange(n_temp):
		filter_matrix[i,:] /= ps_map.reshape(npix)

	filter_ft = (response@np.linalg.inv(norm_matrix)@filter_matrix)
	filter_ft = filter_ft.reshape(nxpix, nypix)

	#compute results
	filtered_image = np.real(np.fft.ifft2(np.fft.ifftshift(filter_ft*image_ft))) * npix
	noise = np.sqrt(np.real(np.sum(filter_ft**2*ps_map)))

	return(filtered_image, filter_ft, noise)


def filter_map_mmf(images, source, spec, noise_maps=None, sigma_noise=None, no_cross = False, global_cross=False):
	'''Computes and applies a matched multifilter to a set of multi-frequency map. The 
       filter is build from provided spatial and spectral source templates and the cross 
       spectrum of a set of maps at different frequencies.

	Parameters
	----------
	images: 3D float array
		Data cube containing the input images that are provided at n_f frequencies.
		The dimensions are supposed to be n_f * n_x * n_y, where n_x and n_y are 
		the number of pixels along the * and y axis of the images. 
	source: 2D or 3D float array
		Source template. Can either be a single image or alternatively a data cube
		if each image comes at a different spatial resolution. Needs to have identical 
		dimensions as images if a cube is provided. Otherwise the single source image 
		needs to match the dimensions of a single image.
	spec: float array
		SED of the source. Needs to be provided at the same n_f frequencies as the 
		input images.
	noise_maps: 3D float array, optional
		Data cube containing provided noise maps that will be used to compute noise 
		power spectrum matrix. Otherwise the noise power spectrum matrix is directly 
		computed from the images. noise_maps need to have the same dimensions as images. 
		Default: None
	sigma_noise: float array, optional
		Array containing the standard deviation of the noise for the n_f maps. Does 
		allow the analytic computation of the noise power spectrum. Useful for tests 
		with white noise. Default: None
	no_cross: bool, optional
		If True the frequency-to-frequency cross terms of the n_f maps will be ignored,
		i.e. non-diagonal elements of the noise power spectrum matrix are set to zero. 
		Default: False
	global_cross: bool, optional
		If set to True, the cross power of the maps is computed globally
		(i.e. averaged over the full maps / all scales). If False, the 
		cross power is computed at each spatial frequency. Default: False

	Returns
	-------
	filtered_image: 2D float array
		Linear combination of the images convolved with the computed matched multifilter.
	filter_ft: 3D float array
		Data cube containing the 2D window functions of the matched multifilter.
	noise: float
		Noise level of the filtered map.

	'''

	#determine dimensions
	nf, nxpix, nypix = images.shape[0], images.shape[1], images.shape[2]
	npix = nxpix*nypix

	#compute FFTs
	tau = np.zeros((nf, npix),dtype=np.complex128)
	for i in np.arange(nf):
		if source.shape[0] == nf:
			source_ft = np.fft.fftshift(np.fft.fft2(source[i,:,:])) / npix
		else: 
			source_ft = np.fft.fftshift(np.fft.fft2(source)) / npix
		source_ft = source_ft.reshape(npix)
		tau[i,:] = spec[i]*source_ft

	#compute noise power spectrum
	if noise_maps is not None:
		cc_maps = np.real(cross_spec(noise_maps, 1, global_cross=global_cross))
	elif sigma_noise is not None:
		cc_maps = np.zeros((nf,nf,npix))
		for i in np.arange(nf):
			noise_ps = np.ones_like(images[0,:,:]) * sigma_noise[i]**2 / npix
			cc_maps[i,i,:] = noise_ps.reshape(npix)
	else:
		cc_maps = np.real(cross_spec(images, 1, global_cross=global_cross))

	#compute filters
	filters = np.zeros((nf, npix),dtype=np.complex128)
	norm = np.zeros(npix,dtype=np.complex128)
	for i in np.arange(npix):
		try:
			if no_cross is True:
				C_inverse = np.linalg.inv(cc_maps[:,:,i]*np.identity(6))
			else:
				C_inverse = np.linalg.inv(cc_maps[:,:,i])

			F = tau[:,i]
			filters[:,i] = (C_inverse@F)
			norm[i] = (np.transpose(F)@C_inverse@F)

		except:
			#print('Matrix inversion failed at pixel: ', i)
			filters[:,i] = filters[:,i-1]
			norm[i] = norm[i-1]

		if i == nxpix*nypix/2 + nxpix/2:
			filters[:,i] = filters[:,i-1]
			norm[i] = norm[i-1]

	#compute results
	noise = np.real(np.sqrt(1/np.sum(norm)))
	filter_ft = filters.reshape((nf,nxpix,nypix)) / np.sum(norm)
	filter_ft[:,0,0] = filter_ft[:,0,nypix-1]

	result = np.zeros((nxpix,nypix),dtype=np.complex128)
	for i in np.arange(nf):
		image_ft = np.fft.fftshift(np.fft.fft2(images[i,:,:])) / npix
		result += filter_ft[i,:,:]*image_ft
	
	filtered_image = np.fft.fftshift(np.real(np.fft.ifft2(np.fft.ifftshift(result)))) * npix

	if (no_cross is True) or (global_cross is True):
		noise = np.std(filtered_image)

	return(filtered_image, filter_ft, noise)


def filter_map_cmmf(images, sources, spec, response, noise_maps=None, sigma_noise=None, no_cross = False, global_cross=False):
	'''Computes and applies a constrained matched multifilter to a set of multi-frequency 
       map. The filter is build from provided spatial and spectral source templates and 
       the cross spectrum of a set of maps at different frequencies.

	Parameters
	----------
	images: 3D float array
		Data cube containing the input images that are provided at n_f frequencies.
		The dimensions are supposed to be n_f * n_x * n_y, where n_x and n_y are 
		the number of pixels along the * and y axis of the images. 
	sources: 3D or 4D float array
		Source templates of the n_c components. Can either be a n_c * n_x * n_y 
		data cube or alternatively a hypercube of dimensions n_c * n_f * n_x * n_y
		if each image comes at a different spatial resolution. 
	spec: 2D float array
		SEDs of the sources. Needs to be provided at the same n_f frequencies as the 
		input images. Dimensions are n_c * n_f.
	response: float array
		Array containing the desired response of the filter to each provided source 
		template.
	noise_maps: 3D float array, optional
		Data cube containing provided noise maps that will be used to compute noise 
		power spectrum matrix. Otherwise the noise power spectrum matrix is directly 
		computed from the images. noise_maps need to have the same dimensions as images. 
		Default: None
	sigma_noise: float array, optional
		Array containing the standard deviation of the noise for the n_f maps. Does 
		allow the analytic computation of the noise power spectrum. Useful for tests 
		with white noise. Default: None
	no_cross: bool, optional
		If True the frequency-to-frequency covariance of the n_f maps will be ignored,
		i.e. non-diagonal elements of the noise power spectrum matrix are set to zero. 
		Default: False
	global_cross: bool, optional
		If set to True, the cross power of the maps is computed globally
		(i.e. averaged over the full maps / all scales). If False, the 
		cross power is computed at each spatial frequency. Default: False

	Returns
	-------
	filtered_image: 2D float array
		Linear combination of the images convolved with the computed constrained matched 
		multifilter.
	filter_ft: 3D float array
		Data cube containing the 2D window functions of the constrained matched multifilter.
	noise: float
		Noise level of the filtered map.

	'''

	#determine dimensions
	nf, nxpix, nypix = images.shape[0], images.shape[1], images.shape[2]
	n_components = len(response)
	npix = nxpix*nypix

	#compute FFTs
	tau = np.zeros((n_components, nf, npix),dtype=np.complex128)
	for i in np.arange(n_components):
		for j in np.arange(nf):
			if sources.shape[1] == nf:
				source_ft = np.fft.fftshift(np.fft.fft2(sources[i,j,:,:])) / npix
			else: 
				source_ft = np.fft.fftshift(np.fft.fft2(sources[i,:,:])) / npix
			source_ft = source_ft.reshape(npix)
			tau[i,j,:] = spec[i,j]*source_ft

	#compute noise power spectrum
	if noise_maps is not None:
		cc_maps =  np.real(cross_spec(noise_maps, 1, global_cross=global_cross))
	elif sigma_noise is not None:
		cc_maps = np.zeros((nf,nf,npix))
		for i in np.arange(nf):
			noise_ps = np.ones_like(images[0,:,:]) * sigma_noise[i]**2 / npix
			cc_maps[i,i,:] = noise_ps.reshape(npix)
	else:
		cc_maps = np.real(cross_spec(images, 1, global_cross=global_cross))

	#compute filters
	norm_matrix = np.zeros((n_components,n_components,npix),dtype=np.complex128)
	for i in np.arange(npix):
		try:
			if no_cross is True:
				C_inverse = np.linalg.inv(cc_maps[:,:,i]*np.identity(6))
			else:
				C_inverse = np.linalg.inv(cc_maps[:,:,i])
		
			F = tau[:,:,i]
			norm_matrix[:,:,i] = F@C_inverse@np.transpose(F)

		except:
			#print('Matrix inversion failed at pixel: ', i)
			norm_matrix[:,:,i] = norm_matrix[:,:,i-1]

		if i == nxpix*nypix/2 + nxpix/2:
			norm_matrix[:,:,i] = norm_matrix[:,:,i-1]

	norm_matrix = np.sum(norm_matrix, axis=2)

	filters = np.zeros((nf, npix),dtype=np.complex128)
	variance = 0
	for i in np.arange(npix):
		try:
			if no_cross is True:
				C_inverse = np.linalg.inv(cc_maps[:,:,i]*np.identity(6))
			else:
				C_inverse = np.linalg.inv(cc_maps[:,:,i])
			F = tau[:,:,i]
			filters[:,i] = response@np.linalg.inv(norm_matrix)@F@C_inverse

		except:
			#print('Matrix inversion failed at pixel: ', i)
			filters[:,i] = filters[:,i-1]

		if i == nxpix*nypix/2 + nxpix/2:
			filters[:,i] = filters[:,i-1]

		variance += np.transpose(filters[:,i])@cc_maps[:,:,i]@filters[:,i]

	#compute results
	noise = np.real(np.sqrt(variance))
	filter_ft = filters.reshape((nf,nxpix,nypix))
	filter_ft[:,0,0] = filter_ft[:,0,nypix-1]

	result = np.zeros((nxpix,nypix),dtype=np.complex128)
	for i in np.arange(nf):
		image_ft = np.fft.fftshift(np.fft.fft2(images[i,:,:])) / npix
		result += image_ft*filter_ft[i,:,:]
	
	filtered_image = np.fft.fftshift(np.real(np.fft.ifft2(np.fft.ifftshift(result)))) * npix

	if (no_cross is True) or (global_cross is True):
		noise = np.std(filtered_image)

	return(filtered_image, filter_ft, noise)
